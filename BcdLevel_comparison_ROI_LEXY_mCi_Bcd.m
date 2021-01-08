% script for generating plots of comparing LEXY-mCitrine-Bcd Bicoid level
% for ROI vs non-ROI.
function BcdLevel_comparison_ROI_LEXY_mCi_Bcd
%% Description
% There can be multiple ways to generate the plots of comparing the Bcd
% level for the export/non-export. But, here, we will take a single embryo,
% then compare the ROI (Region-Of-Interest) illuminated with the blue light
% (458 nm laser), and non-ROI. 

%% Load the datasets
% First, define the file paths
DataPath = 'S:\YangJoon\Dropbox\Optogenetics'
% DataName = '2017-10-23-P2PV1-MCP-mCherry-LEXY-mCi-Bcd_ROI';
DataName = '2018-05-27-P2PV1-MCP-mCherry-LEXY-mCi-Bcd-BcdE1-ROI2';
FigPath = 'S:\YangJoon\Dropbox\Optogenetics\Figures';

Data = load([DataPath, filesep, DataName, filesep, 'CompiledNuclei.mat']);

%% Extract useful fields
Time = Data.ElapsedTime;
nc13 = Data.nc13;
nc14 = Data.nc14;

Fluo_ROI = Data.MeanVectorAP_ROI;
Fluo_SE_ROI = Data.SDVectorAP_ROI./Data.NParticlesAP_ROI;

Fluo_nonROI = Data.MeanVectorAP_nonROI;
Fluo_SE_nonROI = Data.SDVectorAP_nonROI./Data.NParticlesAP_nonROI;

% 2017-10-23 dataset, as at that time, we did define the MeanVectorAP for
% MeanVectorAP_nonROI.
% Fluo_nonROI = Data.MeanVectorAP;
% Fluo_SE_nonROI = Data.SDVectorAP./Data.NParticlesAP;
%% generate plots
% pick an AP bin
APbin = 15; % X %

hold on
errorbar(Time(nc13:end) - Time(nc13), Fluo_ROI(nc13:end, APbin), Fluo_SE_ROI(nc13:end, APbin))
errorbar(Time(nc13:end) - Time(nc13), Fluo_nonROI(nc13:end, APbin), Fluo_SE_nonROI(nc13:end, APbin))
xline(Time(nc14) - Time(nc13), '--')

xlabel('time (min)')
ylabel('Bcd fluorescence (AU)')
legend('ROI','non-ROI')

box on
StandardFigure(gcf,gca)

saveas(gcf,[FigPath,filesep,'LEXY-mCi-Bcd_level_comparison_ROI-nonROI_',num2str((APbin-1)*2.5),'%.tif']); 
saveas(gcf,[FigPath,filesep,'LEXY-mCi-Bcd_level_comparison_ROI-nonROI_',num2str((APbin-1)*2.5),'%.pdf']); 
saveas(gcf,[FigPath,filesep,'LEXY-mCi-Bcd_level_comparison_ROI-nonROI_',num2str((APbin-1)*2.5),'%.fig']); 
end
