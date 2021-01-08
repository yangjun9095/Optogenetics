function SingleTracePlots_Optogenetics(varargin)

% Author : Simon Alamos, edited by YJK : 1/8/2018

% The goal is to plot the MS2 spot trace and nuclear fluo of single
% particle / nucleus with info (ROI, AP bin, particle #, etc.)

%Information about about folders
[~,~,DefaultDropboxFolder,~,~]=...
    DetermineLocalFolders;

% read arguments
if isempty(varargin)%looks for the folder to analyze
    FolderTemp=uigetdir(DefaultDropboxFolder,'Select folder with data to analyze');
    Dashes=strfind(FolderTemp,filesep);
    Prefix=FolderTemp((Dashes(end)+1):end);
else
    Prefix=varargin{1};
    for i=2:length(varargin)
%         if strcmp(varargin{i},'ForceAP')
%             ForceAP=1;
%         elseif strcmp(varargin{i},'SkipTraces')
%             SkipTraces=1;
%         elseif strcmp(varargin{i},'SkipFluctuations')
%             SkipFluctuations=1;
%         elseif strcmp(varargin{i},'SkipFits')    
%             SkipFits=1;
%         elseif strcmp(varargin{i},'SkipMovie')    
%             SkipMovie=1;
%         elseif strcmp(varargin{i},'SkipAll')        
%             SkipTraces=1;
%             SkipFluctuations=1;
%             SkipFits=1;
%             SkipMovie=1;
%         elseif strcmp(varargin{i},'ApproveAll')    
%             ApproveAll=1;
%         elseif strcmp(varargin{i},'MinParticles')
%             if ~isnumeric(varargin{i+1})
%                 error('Wrong input parameters. After ''MinParticles'' you should input the desired minimum number of particles per approved AP bin')
%             else
%                 MinParticles=varargin{i+1};
%             end
%         elseif strcmp(varargin{i},'MinTime')
%             if ~isnumeric(varargin{i+1})
%                 error('Wrong input parameters. After ''MinTime'' you should input the desired minimum number of frames per particle.')
%             else
%                 minTime=varargin{i+1};
%             end
%         end
    end
end

FilePrefix=[Prefix,'_'];

%What type of experiment are we dealing with? Get this out of
%MovieDatabase.xlsx
[~,~,DropboxFolder,~, PreProcPath,...
    ~, ~, ~, ~, ~,~] = readMovieDatabase(Prefix);

%Note that some of this information is redundant given what we get out of
%readMovieDatabase above. We'll have to integrate this better.
[~,XLSTxt,XLSRaw]=xlsread([DefaultDropboxFolder,filesep,'MovieDatabase.xlsx']);
ExperimentTypeColumn=find(strcmp(XLSRaw(1,:),'ExperimentType'));
ExperimentAxisColumn=find(strcmp(XLSRaw(1,:),'ExperimentAxis'));
APResolutionColumn = find(strcmp(XLSRaw(1,:),'APResolution'));

DataFolderColumn=find(strcmp(XLSRaw(1,:),'DataFolder'));
Dashes=findstr(Prefix,'-');
PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));
    if isempty(PrefixRow)
        PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'/',Prefix(Dashes(3)+1:end)]));
        if isempty(PrefixRow)
            error('Could not find data set in MovieDatabase.XLSX. Check if it is defined there.')
        end
    end
        
if isempty(PrefixRow)
    error('Entry not found in MovieDatabase.xlsx')
end

ExperimentType=XLSRaw{PrefixRow,ExperimentTypeColumn};
ExperimentAxis=XLSRaw{PrefixRow,ExperimentAxisColumn};
APResolution = XLSRaw{PrefixRow,APResolutionColumn};


%Load all the information
MS2Particles = load([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'])
Particles = MS2Particles.CompiledParticles;
Nuclei = load([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat'])
Nuclei = Nuclei.CompiledNuclei;

load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])
load([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'])
load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'])



%Create folder to store plots of single traces
singleTracesFolder = ([DropboxFolder,filesep,Prefix,filesep,'SingleTraces']);
mkdir(singleTracesFolder)

%% Extract the information (Spot fluo, Nuclear fluo, and corresponding time frames)

% First, we need to find the index of particles of interest(ex.schnitzcells)
all_nuclei = extractfield(Particles,'Nucleus');

for i=1:length(all_nuclei)
    clear TotalFrames
    clear TotalIndex
    clear TotalFluo
    clear NucFluo
    clear NucFrames
    
    NucIndex(i) = all_nuclei(i); 
    
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

    % Extract the MS2 Spot Fluorescence, same, but easier than above,
    % just grabbing CompiledParticles.Fluo field.
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
    % referred APfilter (CompileParticles)
    APbinID = 0:0.025:1;
    APbin = max(find(APbinID<=Particles(i).MeanAP));
    
%     % Extract the Nuclear fluorescence, time frames
%     % Basically, we use the CompiledNuclei(.schnitz == NucIndex(i)).FluoMax
%     NucFluo = Nuclei(find(Nuclei.schnitz == NucIndex(i))).FluoMax;
%     NucFrames = Nuclei(Nuclei.schnitz == NucIndex(i)).
    
     AllTraces(i).Nucleus = NucIndex(i);
%     AllTraces(i).NucFluo = NucFluo;
%     AllTraces(i).NucFrames = 
    
    AllTraces(i).ParticleFrames = TotalFrames;
    AllTraces(i).Fluo = TotalFluo;
    AllTraces(i).nc = Particles(i).nc;
    AllTraces(i).ROI = ROI;
    AllTraces(i).AP = APbin;

end


        

%% Plot single traces
% setup figure params
set(0,'DefaultAxesFontSize',16)
set(0,'DefaultAxesFontWeight','bold')

%setup scale: find the maximum and the minimum particle fluorescence
MaxParticleFluo = nanmax([Particles.Fluo]);
MinParticleFluo = nanmin([Particles.Fluo]);
MaxMovieTime = (FrameInfo(length(FrameInfo)).Time)/60; %in minutes
set(0, 'defaultaxesylim', [MinParticleFluo MaxParticleFluo])

% ROI string
%         if AllTraces(i).ROI==1
%             ROI = 'ROI';
%         elseif AllTraces(i).ROI==0
%             ROI = 'nonROI';
%         else
%             ROI = 'gray area';
%         end
% First, deal with ROI, and save them separately
Folder = [singleTracesFolder,'\','GrayArea'];
for i=1:length(AllTraces)
    if length(AllTraces(i).Frames) > 5 && AllTraces(i).nc == 13 && AllTraces(i).ROI==2
        ROI = 'Gray Area';
        plot(AllTraces(i).Frames,AllTraces(i).Fluo,'-or','MarkerFaceColor','r','LineWidth',1,'MarkerSize',3)
        hold on
        plot(schnitzcells(NucIndex(i)).frames,...
                max(schnitzcells(NucIndex(i)).Fluo,[],2),'-og','MarkerFaceColor','g','LineWidth',1,'MarkerSize',3)
        hold off
        
        title(['Particle ',num2str(i),' AP =',num2str(AllTraces(i).AP),' nc=',num2str(AllTraces(i).nc),' ',ROI])
        xlabel('Frame')
        ylabel('Fluo')
        ylim([MinParticleFluo MaxParticleFluo])
        legend('MS2 Spot Fluo','Nuclear Fluo')
        %pause
        save figure
        saveas(gcf, [Folder, '\', 'Particle',num2str(i),'AP',num2str(AllTraces(i).AP),'nc=',num2str(AllTraces(i).nc),ROI, '.png'])%[nc14Folder '\' num2str(nucleus) '.png'])
        close all
    end
end

% %loop over particles to plot each of them
% 
% for i = 1:length(CompiledParticles)
%     
%     %get relevant data
%     particleFluo = CompiledParticles(i).Fluo;
%     particleError = CompiledParticles(i).FluoError;
%     particleOffset = CompiledParticles(i).Off;
%     particleFrames = CompiledParticles(i).Frame;
%     particleTimes = MovieTimes(particleFrames);
%     particleTimes = particleTimes;
%     particleNC = CompiledParticles(i).nc;
%convert frames to time in min
% MovieTimes = [FrameInfo.Time]/60;
%     if length(particleFrames) > 5 && particleNC == 13 
%         % Choose the particles only if the particle lasts more than 5 frames, and also is from nc 13
% 
%         %plot
%         figure
%         shadedErrorBar(MovieTimes,MeanVectorAll,SDVectorAll,'lineProps',{'Color',[0.7 0.7 0.7],'LineWidth',2})
%         hold on
%         plot(particleTimes,particleFluo,'ro-','MarkerFaceColor','r','LineWidth',1,'MarkerSize',3)
%         hold off
%         %legend('Standard deviation','mean field of view activity','','something')
%         ylim([MinParticleFluo MaxParticleFluo])
%         xlim([0 MaxMovieTime])
%         ylabel('Particle Fluorescence (a.u)')
%         xlabel('Time (min)')
%         title(['Particle #' num2str(i)])
%         pause
%         %save figure
%         saveas(gcf, [singleTracesFolder '\' num2str(i) '.png'])%[nc14Folder '\' num2str(nucleus) '.png'])
%         close all
%     end
% %end
% 
% %plot total transcriptional activity per frame


    


end



