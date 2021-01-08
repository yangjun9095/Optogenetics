function SingleTrace_Analysis_Optogenetics

% Author : Yang Joon Kim : Jan, 2018

% The goal is to plot the MS2 spot / Nuclear fluorescence trace of single
% particle / nucleus with info (ROI, AP bin, particle #, etc.)

% Define the Prefix (dataset), and directory
Prefix='2018-05-09-3A3-MCP-mCherry-2x-LEXY-mCi-Bcd2';
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(Prefix);
% Note that info from the MovieDatabase, I need to include the code from pipeline to
% recall that info.

%Load all the information ( about CompiledParticles, CompiledNuclei)
% Note that CompiledParticles, and CompiledNuclei are structure, and
% Particles, Nuclei are .CompiledParticles/CompiledNuclei field in that structure (cell)

CompiledParticles = load([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'])
Particles = CompiledParticles.CompiledParticles;
CompiledNuclei = load([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat'])
Nuclei = CompiledNuclei.CompiledNuclei;

load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])
load([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'])
load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'])

%Create folder to store plots of single traces
singleTracesFolder = ([DropboxFolder,filesep,Prefix,filesep,'SingleTraceAnalysis']);
mkdir(singleTracesFolder)

%% Extract the information (Spot fluo, Nuclear fluo, and corresponding time frames)

% First, we need to find the index of CompiledParticles of interest(ex.schnitzcells)
all_Nuclei = extractfield(Particles,'Nucleus');

% Define the index of schnitzcells, for indexing of nucleus containing
% CompiledParticles later
for j=1:length(Nuclei)
    schnitzIndex(j) = Nuclei(j).schnitz;
end
    
for i=1:length(all_Nuclei)
    clear TotalFrames
    clear TotalIndex
    clear TotalFluo
    clear NucFluo
    clear NucFrames
    
    NucIndex(i) = all_Nuclei(i); 
    
    % For each nucleus, find the MS2 particle that is inside
    TotalFrames = Particles(i).Frame;
    TotalIndex = Particles(i).Index;
    
    % Extract the fluo intensity for each particle at all time frames
    % Use the z-stack where the Spots.brighestZ is the same as Spots.z
    % Spots(Frame).Fits(Index).FixedAreaIntensity(.z == .brightestZ)
%     for j=1:length(TotalFrames)
%         clear TempFits
%         TempFits=Spots(TotalFrames(j)).Fits(TotalIndex(j));
%         TotalFluo(j)=TempFits.FixedAreaIntensity (TempFits.z == TempFits.brightestZ);
%     
%     end

    % Extract the MS2 Spot Fluorescence, identical but easier than the
    % script above, (.z == .brightestZ)
    % : just grabbing CompiledCompiledParticles.Fluo field.
    TotalFluo = Particles(i).Fluo;
    
    % ROI and AP bin info
    % First, define ROI based on the y position (yPos field), this should be fixed more
    % systematically later. (Using xmlread to extract the ROI positional
    % info)
    if nanmean(Particles(i).yPos)<200
        ROI = 0; % non-ROI
    elseif nanmean(Particles(i).yPos)>300
        ROI = 1; % ROI
    else 
        ROI = 2; % gray area between ROI and non-ROI
    end
    
    % Second, define the AP bin based on the MeanAP field
    % referred APfilter (CompileCompiledParticles)
    APbinID = 0:0.025:1;
    APbin = max(find(APbinID<=Particles(i).MeanAP));
    
    % Extract the Nuclear fluorescence, time frames
    % Basically, we use the CompiledCompiledNuclei(.schnitz == NucIndex(i)).FluoMax
    % I need to think about the case when schnitzcell is not in the
    % CompiledCompiledNuclei since it's on the boundary or some other reasons.
    if ismember(NucIndex(i),schnitzIndex)
        NucFluo = Nuclei(find(schnitzIndex == NucIndex(i))).FluoMax;
        NucFrames = Nuclei(schnitzIndex == NucIndex(i)).Frames;
    else 
        NucFluo = nan;
        NucFrames = nan;
    end
    
    AllTraces(i).Nucleus = NucIndex(i);
    AllTraces(i).NucFluo = NucFluo;
    AllTraces(i).NucFrames = NucFrames;
    
    AllTraces(i).ParticleFrames = TotalFrames;
    AllTraces(i).Fluo = TotalFluo;
    AllTraces(i).nc = Particles(i).nc;
    AllTraces(i).ROI = ROI;
    AllTraces(i).AP = APbin;

end
        
%% Plot single traces (at specific AP bin, nuclear cycle)
% The idea here is that I am plotting the single traces on top of the mean
% value (MeanVectorAP, with error bars), to see how individual CompiledNuclei
% behaves differently from the averaged. In the end, I want to plot both
% ROI and non-ROI region in the same plot.

% setup figure parameters
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultAxesFontWeight','bold')

%setup scale: find the maximum and the minimum particle fluorescence
MaxParticleFluo = nanmax([Particles.Fluo]);
MinParticleFluo = nanmin([Particles.Fluo]);
%MaxNuclearFluo = nanmax([Nuclei.FluoMax]);
%MinNuclearFluo = nanmin([Nuclei.FluoMax]);
MaxMovieTime = (FrameInfo(length(FrameInfo)).Time)/60; %in minutes
set(0, 'defaultaxesylim', [MinParticleFluo MaxParticleFluo])

% Info about the Averaged Particles/ Nuclei
nc12 = CompiledParticles.nc12;
nc13 = CompiledParticles.nc13;
nc14 = CompiledParticles.nc14;

% Plot the single traces with the Mean
for i=1:length(AllTraces)
    % I am only getting particles which have more than 5 frames, and also
    % apply some conditions, such as nuclear cycle, ROI,AP  and the nuclear
    % fluo is not Nans.
    if length(AllTraces(i).ParticleFrames) > 3 && (AllTraces(i).nc == 12 | AllTraces(i).nc == 13)...   
                &&AllTraces(i).ROI==1&& ~isnan(mean(AllTraces(i).NucFluo))%&&AllTraces(i).AP == 13
        % Setting for ROI or not, I need to change the MeanVectorAP_ROI as
        % we
%         ROI string
        if AllTraces(i).ROI==1
            ROI = 'ROI';
        elseif AllTraces(i).ROI==0
            ROI = 'nonROI';
        else
            ROI = 'gray area';
        end
        %First, deal with ROI, and save them separately

         Folder = [singleTracesFolder,'/',ROI];
         
         % Define the AP bin
         APbin = AllTraces(i).AP;
         
        hold on
        % First, plot MS2 spot fluorescence
        % Averaged trace (over one AP bin)
        shadedErrorBar(nc12:nc13,CompiledParticles.MeanVectorAP_nonROI(nc12:nc13,APbin),...
            CompiledParticles.SDVectorAP_nonROI(nc12:nc13,APbin)./CompiledParticles.NParticlesAP_nonROI(nc12:nc13,APbin),...
            'lineProps',{'Color',[0.9 0 0],'LineWidth',2})
        yyaxis left
        % single trace
        plot(AllTraces(i).ParticleFrames,AllTraces(i).Fluo,'-or','MarkerFaceColor','r','LineWidth',1,'MarkerSize',3)

        xlabel('Frame')
        ylabel('transcript intensity (a.u.)','color','r') 
        ylim([MinParticleFluo MaxParticleFluo])
        
        % Second, plot nuclear fluorescence
        % Averaged trace (over one AP bin)
        shadedErrorBar(nc12:nc13,CompiledNuclei.MeanVectorAP_nonROI(nc12:nc13,APbin),...
            CompiledNuclei.SDVectorAP_nonROI(nc12:nc13,APbin)./CompiledNuclei.NParticlesAP_nonROI(nc12:nc13,APbin),...
            'lineProps',{'Color',[0 0.9 0],'LineWidth',2})
        yyaxis right
        ylabel('input protein intensity (a.u.)','color','g');
        % single trace
        plot(AllTraces(i).NucFrames,AllTraces(i).NucFluo,'-og','MarkerFaceColor','g','LineWidth',1,'MarkerSize',3)
        % old code to get the nuclear fluo        
%         plot(schnitzcells(NucIndex(i)).frames,...
%                 max(schnitzcells(NucIndex(i)).Fluo,[],2),'-og','MarkerFaceColor','g','LineWidth',1,'MarkerSize',3)
        %ylim([MinNuclearFluo MaxNuclearFluo])
        ylim([0 200])
        hold off
        
        title(['Particle ',num2str(i),' AP =',num2str(AllTraces(i).AP),' nc=',num2str(AllTraces(i).nc),' ',ROI])
        %legend('MS2 Spot Fluo','Nuclear Fluo')
        %pause
        save figure
        saveas(gcf, [Folder, '\', 'Particle',num2str(i),'AP',num2str(APbin),'nc=',num2str(AllTraces(i).nc),ROI, '.png'])%[nc14Folder '\' num2str(nucleus) '.png'])
        close all
    end
end

%% Compare the ROI and non-ROI (spot fluo and nuc fluo) at the same AP bin
% Assuem that nc13,nc14, and CompiledParticles, and CompiledNuclei are all
% defined above.

% Define the AP bin
APbin = 17;
hold on
% plot ROI first
    shadedErrorBar(nc13:nc14,CompiledParticles.MeanVectorAP_ROI(nc13:nc14,APbin),...
        CompiledParticles.SDVectorAP_ROI(nc13:nc14,APbin)./CompiledParticles.NParticlesAP_ROI(nc13:nc14,APbin),...
        'lineProps',{'Color',[0.9 0 0],'LineWidth',2})
    % Nuclear fluo
    shadedErrorBar(nc13:nc14,CompiledNuclei.MeanVectorAP_ROI(nc13:nc14,APbin),...
        CompiledNuclei.SDVectorAP_ROI(nc13:nc14,APbin)./CompiledNuclei.NParticlesAP_ROI(nc13:nc14,APbin),...
        'lineProps',{'Color',[0 0.9 0],'LineWidth',2})
% plot non-ROI
    shadedErrorBar(nc13:nc14,CompiledParticles.MeanVectorAP(nc13:nc14,APbin),...
        CompiledParticles.SDVectorAP(nc13:nc14,APbin)./CompiledParticles.NParticlesAP(nc13:nc14,APbin),...
        'lineProps',{'Color',[0.5 0 0],'LineWidth',2})
    % Nuclear fluo
    shadedErrorBar(nc13:nc14,CompiledNuclei.MeanVectorAP(nc13:nc14,APbin),...
        CompiledNuclei.SDVectorAP(nc13:nc14,APbin)./CompiledNuclei.NParticlesAP(nc13:nc14,APbin),...
        'lineProps',{'Color',[0 0.5 0],'LineWidth',2})
    
    title('')
    xlabel('')
    ylabel('')
    legend('')
hold off
%% Nuclear florescence over AP (ROI and non-ROI), or over time


    


end



