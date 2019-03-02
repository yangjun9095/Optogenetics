% Here, I want to compare two optogenetics datasets with or without
% illumination with the blue light.
% What I ultimately want to have is Figure 4, Garcia 2013 analysis for
% these optogenetic Bcd constructs driving the hbP2P.V1-MS2-lacZ
% expression.

%% Load the datasets (using LoadMS2sets)
NoBlueData = LoadDataSets('noBlue')
BlueData = LoadDataSets('Blue')

%% Plot the accumulated mRNA over AP, in nc13 and nc 14.
% I will use the IntegratemRNA.m to calculate the accumulated mRNA

[TotalProd_NoBlue,TotalProdError_NoBlue,TotalProdN_NoBlue,...
    MeanTotalProd_NoBlue,SDTotalProd_NoBlue,SETotalProd_NoBlue] = IntegratemRNA([NoBlueData.Particles],1,1);

[TotalProd_Blue,TotalProdError_Blue,TotalProdN_Blue,...
    MeanTotalProd_Blue,SDTotalProd_Blue,SETotalProd_Blue] = IntegratemRNA([BlueData.Particles],1,1);

%% Plot to compare different quantities

% 1) Total mRNA per nucleus at nc 13
figure(1)
hold on
errorbar(0:0.025:1, nanmean(TotalProd_NoBlue(1:3,:,13),1),nanmean(TotalProdError_NoBlue(1:3,:,13),1),'r')%,TotalProdError_NoBlue(1,:,13)./sqrt(TotalProdN_NoBlue(1,:,13)),'r')
errorbar(0:0.025:1, nanmean(TotalProd_Blue(1:3,:,13),1),nanmean(TotalProdError_Blue(1:3,:,13),1),'b')%,TotalProdError_Blue(1,:,13)./sqrt(TotalProdN_Blue(1,:,13)),'b')

hold off
xlabel('AP bin')
ylabel('Integrated mRNA (AU)')
title('Integrated mRNA over AP @ nc 13')
%legend('Control','Export')
% plot settings, like a title, x,y labels, etc.

% 2) Total mRNA per nucleus at nc 14
figure(2)
hold on
errorbar(0:0.025:1, nanmean(TotalProd_NoBlue(1:3,:,14),1),nanmean(TotalProdError_NoBlue(1:3,:,14),1),'r')%,TotalProdError_NoBlue(1,:,14)./sqrt(TotalProdN_NoBlue(1,:,14)),'r')
errorbar(0:0.025:1, nanmean(TotalProd_Blue(1:3,:,14),1),nanmean(TotalProdError_Blue(1:3,:,14),1),'b')%,TotalProdError_Blue(1,:,14)./sqrt(TotalProdN_Blue(1,:,14)),'b')
hold off
xlabel('AP bin')
ylabel('Integrated mRNA (AU)')
title('Integrated mRNA over AP @ nc 14')
legend('Control','Export')

% 3) 
%% Plot single embryos
figure(1)
hold on
%errorbar(0:0.025:1, nanmean(TotalProd_NoBlue(1:3,:,13),1),nanmean(TotalProdError_NoBlue(1:3,:,13),1),'r')%,TotalProdError_NoBlue(1,:,13)./sqrt(TotalProdN_NoBlue(1,:,13)),'r')
errorbar(0:0.025:1, TotalProd_NoBlue(1,:,13),TotalProdError_NoBlue(1,:,13),'r')%,TotalProdError_NoBlue(1,:,13)./sqrt(TotalProdN_NoBlue(1,:,13)),'r')
errorbar(0:0.025:1, TotalProd_NoBlue(2,:,13),TotalProdError_NoBlue(2,:,13),'r')
errorbar(0:0.025:1, TotalProd_NoBlue(3,:,13),TotalProdError_NoBlue(3,:,13),'r')

%errorbar(0:0.025:1, nanmean(TotalProd_Blue(1:3,:,13),1),nanmean(TotalProdError_Blue(1:3,:,13),1),'b')%,TotalProdError_Blue(1,:,13)./sqrt(TotalProdN_Blue(1,:,13)),'b')
errorbar(0:0.025:1, TotalProd_Blue(1,:,13),TotalProdError_Blue(1,:,13),'b')%,TotalProdError_NoBlue(1,:,13)./sqrt(TotalProdN_NoBlue(1,:,13)),'r')
errorbar(0:0.025:1, TotalProd_Blue(2,:,13),TotalProdError_Blue(2,:,13),'b')
errorbar(0:0.025:1, TotalProd_Blue(3,:,13),TotalProdError_Blue(3,:,13),'b')

xlabel('AP bin')
ylabel('Integrated mRNA (AU)')
title('Integrated mRNA over AP @ nc 14')

figure(2)
hold on
%errorbar(0:0.025:1, nanmean(TotalProd_NoBlue(1:3,:,13),1),nanmean(TotalProdError_NoBlue(1:3,:,13),1),'r')%,TotalProdError_NoBlue(1,:,13)./sqrt(TotalProdN_NoBlue(1,:,13)),'r')
errorbar(0:0.025:1, TotalProd_NoBlue(1,:,14),TotalProdError_NoBlue(1,:,14),'r')%,TotalProdError_NoBlue(1,:,13)./sqrt(TotalProdN_NoBlue(1,:,13)),'r')
errorbar(0:0.025:1, TotalProd_NoBlue(2,:,14),TotalProdError_NoBlue(2,:,14),'r')
errorbar(0:0.025:1, TotalProd_NoBlue(3,:,14),TotalProdError_NoBlue(3,:,14),'r')

%errorbar(0:0.025:1, nanmean(TotalProd_Blue(1:3,:,13),1),nanmean(TotalProdError_Blue(1:3,:,13),1),'b')%,TotalProdError_Blue(1,:,13)./sqrt(TotalProdN_Blue(1,:,13)),'b')
errorbar(0:0.025:1, TotalProd_Blue(1,:,14),TotalProdError_Blue(1,:,14),'b')%,TotalProdError_NoBlue(1,:,13)./sqrt(TotalProdN_NoBlue(1,:,13)),'r')
errorbar(0:0.025:1, TotalProd_Blue(2,:,14),TotalProdError_Blue(2,:,14),'b')
errorbar(0:0.025:1, TotalProd_Blue(3,:,14),TotalProdError_Blue(3,:,14),'b')

xlabel('AP bin')
ylabel('Integrated mRNA (AU)')
title('Integrated mRNA over AP @ nc 14')

%% RNAP II Loading rate (by FitMeanAPSymmetric)
% Also average this rate, and calculate the total SD for all embryos.
noBlueFit = load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2018-02-05-P2PV1-MCP-mCherry-LEXY-mCi-Bcd_noBlue\MeanFits.mat')
BlueFit = load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2018-02-05-P2PV1-MCP-mCherry-LEXY-mCi-Bcd_Blue2\MeanFits.mat')

noBlueFit_nc13 = nan(1,41);
noBlueFit_nc14 = nan(1,41);
BlueFit_nc13 = nan(1,41);
BlueFit_nc14 = nan(1,41);

noBlueFit_nc13_SD = nan(1,41);
noBlueFit_nc14_SD = nan(1,41);
BlueFit_nc13_SD = nan(1,41);
BlueFit_nc14_SD = nan(1,41);


for i=8:22
    noBlueFit_nc13(i) = noBlueFit.FitResults(i,2).RateFit;
    noBlueFit_nc13_SD(i)=noBlueFit.FitResults(i,2).SDRateFit;
    
    noBlueFit_nc14(i) = noBlueFit.FitResults(i,3).RateFit;
    noBlueFit_nc14_SD(i) = noBlueFit.FitResults(i,3).SDRateFit;
    
    BlueFit_nc13(i) = BlueFit.FitResults(i,2).RateFit;
    BlueFit_nc13_SD(i) = BlueFit.FitResults(i,2).SDRateFit;
    
    BlueFit_nc14(i) = BlueFit.FitResults(i,3).RateFit;
    BlueFit_nc14_SD(i) = BlueFit.FitResults(i,3).SDRateFit;
end

%% Plot RNAP II loading rate
hold on
errorbar(0:0.025:1,noBlueFit_nc13,noBlueFit_nc13_SD,'r')
errorbar(0:0.025:1,BlueFit_nc13,BlueFit_nc13_SD,'b')
hold off
xlabel('AP bin')
ylabel('RNAP Loading Rate (AU)')
title('RNAP Loading Rate @ nc 13')
legend('Control','Export')

figure(2)
hold on
errorbar(0:0.025:1,noBlueFit_nc14,noBlueFit_nc14_SD,'r')
errorbar(0:0.025:1,BlueFit_nc14,BlueFit_nc14_SD,'b')

xlabel('AP bin')
ylabel('RNAP Loading Rate (AU)')
title('RNAP Loading Rate @ nc 14')
legend('Control','Export')

%% Things to check

% 1) offset of MCP-mCherry
% noBlue = load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2018-01-11-P2PV1-MCP-mCherry-LEXY-mCi-Bcd_noBlue\CompiledParticles.mat')
% Blue = load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2018-01-15-P2PV1-MCP-mCherry-LEXY-mCi-Bcd_Blue\CompiledParticles.mat')
% 
% hold on
% plot(noBlue.ElapsedTime(noBlue.nc13:noBlue.nc14),noBlue.MeanOffsetVector(noBlue.nc13:noBlue.nc14))
% plot(Blue.ElapsedTime(Blue.nc13-21:Blue.nc14-21),Blue.MeanOffsetVector(Blue.nc13:Blue.nc14))

% 2) Bicoid concentration

% 3) Averaging the Rate fit for multiple embryos with proper error bar

% 4) Find the AP bins where the Bicoid concentration profile is similar
% (for Blue and no Blue), then see how the output looks like.



%% Nuclear fluorescence     % This also needs averaging over multiple embryos.
% Bcd concentration (to check whether it was homozygous)
% noBlueNuc = load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2018-02-05-P2PV1-MCP-mCherry-LEXY-mCi-Bcd_noBlue\CompiledNuclei.mat')
% BlueNuc = load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2018-02-05-P2PV1-MCP-mCherry-LEXY-mCi-Bcd_Blue2\CompiledNuclei.mat')

% Define the time points (for synchronization)
for i=1:length(NoBlueData)
    noBlue_nc13(i) = NoBlueData(i).Particles.nc13;
    noBlue_nc14(i) = NoBlueData(i).Particles.nc14;
    NucFluo_noBlue{i} = NoBlueData(i).Nuclei.MeanVectorAP(noBlue_nc13(i):end,:);
    NucFluo_STD_noBlue{i} = NoBlueData(i).Nuclei.SDVectorAP(noBlue_nc13(i):end,:);
    NucFluo_N_noBlue{i} = NoBlueData(i).Nuclei.NParticlesAP(noBlue_nc13(i):end,:);
end

for i=1:length(BlueData)
    Blue_nc13(i) = BlueData(i).Particles.nc13;
    Blue_nc14(i) = BlueData(i).Particles.nc14;
    NucFluo_Blue{i} = BlueData(i).Nuclei.MeanVectorAP(Blue_nc13(i):end,:);
    NucFluo_STD_Blue{i} = BlueData(i).Nuclei.SDVectorAP(Blue_nc13(i):end,:);
    NucFluo_N_Blue{i} = BlueData(i).Nuclei.NParticlesAP(Blue_nc13(i):end,:);
end

% Averaging the Nuc fluo over multiple embryos
% For all timepoints, for all embryos

% No Blue
NucFluo_NoBlueTotal=zeros(71,41);
NucFluo_Var_NoBlue=zeros(71,41);
for i=1:length(NoBlueData)
    NucFluo_NoBlueTotal = NucFluo_NoBlueTotal + NucFluo_noBlue{i}(1:71,:);
    NucFluo_Var_NoBlue = NucFluo_Var_NoBlue + NucFluo_STD_noBlue{i}(1:71,:).^2;
end
    NucFluo_NoBlueAvg = NucFluo_NoBlueTotal./length(NoBlueData);
    NucFluo_NoBlueSTD = sqrt(NucFluo_Var_NoBlue);

% Blue
NucFluo_BlueTotal=zeros(69,41);
NucFluo_Var_Blue=zeros(69,41);
for i=1:length(BlueData)
    NucFluo_BlueTotal = NucFluo_BlueTotal + NucFluo_Blue{i}(1:69,:);
    NucFluo_Var_Blue = NucFluo_Var_Blue + NucFluo_STD_Blue{i}(1:69,:).^2;
end
    NucFluo_BlueAvg = NucFluo_BlueTotal./length(BlueData);
    NucFluo_BlueSTD = sqrt(NucFluo_Var_Blue);

% Define the AP bin
AP = 10;
hold on
errorbar((1:71)*0.66,NucFluo_NoBlueAvg(:,AP),NucFluo_NoBlueSTD(:,AP),'r')
errorbar((1:69)*0.66,NucFluo_BlueAvg(:,AP),NucFluo_BlueSTD(:,AP),'b')

xlabel('Time (min)')
ylabel('Nuclear Fluorescence (AU)')
title(['Nuclear Fluorescence over Time',' @ ',num2str((AP-1)*2.5),'% AP'])
legend('Control','Export')

% 
% hold on
% errorbar(noBlueNuc.ElapsedTime(noBlue_nc13:end)-noBlueNuc.ElapsedTime(noBlue_nc13),...
%         noBlueNuc.MeanVectorAP(noBlue_nc13:end,AP),...
%         (noBlueNuc.SDVectorAP(noBlue_nc13:end,AP)./sqrt(noBlueNuc.NParticlesAP(noBlue_nc13:end,AP))),'r')
% errorbar(BlueNuc.ElapsedTime(Blue_nc13:end)-BlueNuc.ElapsedTime(Blue_nc13),...
%         BlueNuc.MeanVectorAP(Blue_nc13:end,AP),...
%         (BlueNuc.SDVectorAP(Blue_nc13:end,AP)./sqrt(BlueNuc.NParticlesAP(Blue_nc13:end,AP))),'b')
% 


%% Averaging multiple embryos
% As of 2018-02-12, we have three embryos for each conditions (either the
% blue light was on or off for the whole time)

% BlueData
% First, averaging the MeanVectorAP with appropriate error bar of SEM.
% To do this, we have to be careful about synchronization.
% Plot the Mean Spot fluo of all datasets individually.

AP=15;
figure(1) % BlueData
ylim([0 500])
hold on
for i=1:length(BlueData)
%     SpotFluo_Blue(i,:,:) = BlueData(i).Particles.MeanVectorAP;
%     SDSpotFluo_Blue(i,:,:) = BlueData(i).Particles.SDVectorAP;
%     NSpots_Blue(i,:,:) = BlueData(i).Particles.NParticlesAP;
    plot(BlueData(i).Particles.ElapsedTime(BlueData(i).Particles.nc13:end)-...
            BlueData(i).Particles.ElapsedTime(BlueData(i).Particles.nc13),...
            BlueData(i).Particles.MeanVectorAP(BlueData(i).Particles.nc13:end,AP),'b')
    pause
end

%figure(2) % NoBlueData

hold on
for i=1:length(NoBlueData)

    plot(NoBlueData(i).Particles.ElapsedTime(NoBlueData(i).Particles.nc13:end)-...
            NoBlueData(i).Particles.ElapsedTime(NoBlueData(i).Particles.nc13),...
            NoBlueData(i).Particles.MeanVectorAP(NoBlueData(i).Particles.nc13:end,AP),'r')
    pause
end

title(['Mean Spot Fluorescence',' @ ','AP = ',num2str((AP-1)*2.5),'%'])
xlabel('Time (min)')
ylabel('Mean Spot Fluorescence (AU)')

%% Averaging the mean spot fluo
% First, define Nan matrix to plug in the MeanVectorAP values from the different
% embryos. We need this because different embryos have different length of
% the time series.

for i=1:length(BlueData)
    [dummy,BlueLength(i)] = size(BlueData(i).Particles.ElapsedTime(BlueData(i).Particles.nc13:end));
end

MaxBlueLength = max(BlueLength);

VarSum = zeros(MaxBlueLength,41);
Sum = zeros(MaxBlueLength,41);
TotalNumber = zeros(MaxBlueLength,41);
Denominator = zeros(MaxBlueLength,41);

for i=1:length(BlueData)
    SpotFluo_Blue = zeros(BlueLength(i),41);
    NParticles_Blue = zeros(BlueLength(i),41);
    SDSpotFluo_Blue = zeros(BlueLength(i),41);
    
    % nc
    nc13 = BlueData(i).Particles.nc13;
    nc14 = BlueData(i).Particles.nc14;
    
    % Remove Nans
    BlueData(i).Particles.MeanVectorAP(isnan(BlueData(i).Particles.MeanVectorAP)) = 0;
    BlueData(i).Particles.NParticles(isnan( BlueData(i).Particles.NParticlesAP)) = 0;
    BlueData(i).Particles.SDVectorAP(isnan(BlueData(i).Particles.SDVectorAP)) = 0;
    
    MeanSpotFluo = BlueData(i).Particles.MeanVectorAP(nc13:end,:);
    NParticlesAP = BlueData(i).Particles.NParticlesAP(nc13:end,:);
    SDSpotFluo = BlueData(i).Particles.SDVectorAP(nc13:end,:);
    
    if BlueLength(i)==MaxBlueLength
        SpotFluo_Blue= MeanSpotFluo;
        NParticles_Blue = NParticlesAP;
        SDSpotFluo_Blue = SDSpotFluo;
    else 
        SpotFluo_Blue= MeanSpotFluo;
        SpotFluo_Blue(BlueLength(i)+1:MaxBlueLength,1:41)=0;
        NParticles_Blue = NParticlesAP;
        NParticles_Blue(BlueLength(i)+1:MaxBlueLength,1:41)=0;
        SDSpotFluo_Blue= SDSpotFluo;
        SDSpotFluo_Blue(BlueLength(i)+1:MaxBlueLength,1:41)=0;
    end
    
    
    Sum = Sum + (SpotFluo_Blue.*NParticles_Blue);
    
    VarSum = VarSum + SDSpotFluo_Blue.^2 .* (NParticles_Blue);
    TotalNumber = TotalNumber + NParticles_Blue;
    Denominator = Denominator + NParticles_Blue;
end


BlueMeanSpotFluo = Sum./TotalNumber;
BlueSEMSpotFluo = sqrt(VarSum./Denominator);

%%
% No Blue Data
for i=1:length(NoBlueData)
    [dummy,NoBlueLength(i)] = size(NoBlueData(i).Particles.ElapsedTime(NoBlueData(i).Particles.nc13:end));
end

MaxNoBlueLength = max(NoBlueLength);

VarSum = zeros(MaxNoBlueLength,41);
Sum = zeros(MaxNoBlueLength,41);
TotalNumber = zeros(MaxNoBlueLength,41);
Denominator = zeros(MaxNoBlueLength,41);

for i=1:length(NoBlueData)
    SpotFluo_NoBlue = zeros(NoBlueLength(i),41);
    NParticles_NoBlue = zeros(NoBlueLength(i),41);
    SDSpotFluo_NoBlue = zeros(NoBlueLength(i),41);
    
    % nc
    nc13 = NoBlueData(i).Particles.nc13;
    nc14 = NoBlueData(i).Particles.nc14;
    
    % Remove Nans
    NoBlueData(i).Particles.MeanVectorAP(isnan(NoBlueData(i).Particles.MeanVectorAP)) = 0;
    NoBlueData(i).Particles.NParticles(isnan(NoBlueData(i).Particles.NParticlesAP)) = 0;
    NoBlueData(i).Particles.SDVectorAP(isnan(NoBlueData(i).Particles.SDVectorAP)) = 0;
    
    MeanSpotFluo = NoBlueData(i).Particles.MeanVectorAP(nc13:end,:);
    NParticlesAP = NoBlueData(i).Particles.NParticlesAP(nc13:end,:);
    SDSpotFluo = NoBlueData(i).Particles.SDVectorAP(nc13:end,:);
    
    if NoBlueLength(i)==MaxNoBlueLength
        SpotFluo_NoBlue= MeanSpotFluo;
        NParticles_NoBlue = NParticlesAP;
        SDSpotFluo_NoBlue = SDSpotFluo;
    else 
        SpotFluo_NoBlue= MeanSpotFluo;
        SpotFluo_NoBlue(NoBlueLength(i)+1:MaxNoBlueLength,1:41)=0;
        NParticles_NoBlue = NParticlesAP;
        NParticles_NoBlue(NoBlueLength(i)+1:MaxNoBlueLength,1:41)=0;
        SDSpotFluo_NoBlue= SDSpotFluo;
        SDSpotFluo_NoBlue(NoBlueLength(i)+1:MaxNoBlueLength,1:41)=0;
    end
    
    
    Sum = Sum + (SpotFluo_NoBlue.*NParticles_NoBlue);
    
    VarSum = VarSum + SDSpotFluo_NoBlue.^2 .* (NParticles_NoBlue-1);
    TotalNumber = TotalNumber + NParticles_NoBlue;
    Denominator = Denominator + NParticles_NoBlue-1;
end


NoBlueMeanSpotFluo = Sum./TotalNumber;
NoBlueSEMSpotFluo = sqrt(VarSum./Denominator);

hold on
AP = 10;
errorbar((1:length(BlueMeanSpotFluo))*0.66,BlueMeanSpotFluo(:,AP),BlueSEMSpotFluo(:,AP),'b')
errorbar((1:length(NoBlueMeanSpotFluo))*0.66,NoBlueMeanSpotFluo(:,AP),NoBlueSEMSpotFluo(:,AP),'r')

title(['Mean Spot Fluorescence',' @ ','AP = ',num2str((AP-1)*2.5),'%'])
xlabel('Time (min)')
ylabel('Spot Fluorescence (AU)')
legend('Export','Control')

%% Fraction ON
