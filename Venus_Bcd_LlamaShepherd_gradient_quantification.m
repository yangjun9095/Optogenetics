% Compare the Venus-Bcd gradient in the absence/presence of Llama Shepherd.

% I took the intensity profile from ImageJ
% Let's smoothen the intensity profile


Venus_Bcd_BcdE1_smooth = movmean(Venus_Bcd_BcdE1,50);
Venus_Bcd_BcdE1_smooth_STD = movstd(Venus_Bcd_BcdE1,50);

Venus_Bcd_BcdE1_LlamaShepherd_smooth = movmean(Venus_Bcd_BcdE1_LlamaShepherd,50);
Venus_Bcd_BcdE1_LlamaShepherd_smooth_STD = movstd(Venus_Bcd_BcdE1_LlamaShepherd,50);

Venus_Bcd_BcdE1_LlamaShepherd2_smooth = movmean(Venus_Bcd_BcdE1_LlamaShepherd2,50);
Venus_Bcd_BcdE1_LlamaShepherd2_smooth_STD = movstd(Venus_Bcd_BcdE1_LlamaShepherd2,50);

hold on
errorbar(Venus_Bcd_BcdE1(1:50:end,1)/max(Venus_Bcd_BcdE1(1:50:end,1)),...
       Venus_Bcd_BcdE1_smooth(1:50:end,2)/max(Venus_Bcd_BcdE1_smooth(1:50:end,2)),...
       Venus_Bcd_BcdE1_smooth_STD(1:50:end,2)/max(Venus_Bcd_BcdE1_smooth(1:50:end,2)))

errorbar(Venus_Bcd_BcdE1_LlamaShepherd(1:50:end,1)/max(Venus_Bcd_BcdE1_LlamaShepherd(1:50:end,1)),...
        Venus_Bcd_BcdE1_LlamaShepherd_smooth(1:50:end,2)/max(Venus_Bcd_BcdE1_LlamaShepherd_smooth(1:50:end,2)),...
        Venus_Bcd_BcdE1_LlamaShepherd_smooth_STD(1:50:end,2)/max(Venus_Bcd_BcdE1_LlamaShepherd_smooth(1:50:end,2)))

errorbar(Venus_Bcd_BcdE1_LlamaShepherd2(1:50:end,1)/max(Venus_Bcd_BcdE1_LlamaShepherd2(1:50:end,1)),...
        Venus_Bcd_BcdE1_LlamaShepherd2_smooth(1:50:end,2)/max(Venus_Bcd_BcdE1_LlamaShepherd2_smooth(1:50:end,2)),...
        Venus_Bcd_BcdE1_LlamaShepherd2_smooth_STD(1:50:end,2)/max(Venus_Bcd_BcdE1_LlamaShepherd2_smooth(1:50:end,2)))

title('Venus-Bcd gradient')
xlabel('AP axis')
ylabel('Normalized Venus-Bcd intensity')
legend('Control','Llama-Shepherd','Llama-Shepherd')
    
StandardFigure(gcf,gca)
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Optogenetics\LlamaShephereds\Bcd_gradient_Quantification';
saveas(gcf,[FigPath,'Venus-Bcd_BcdE1_with_without_LlamaShepherd_comparison.tif'])
saveas(gcf,[FigPath,'Venus-Bcd_BcdE1_with_without_LlamaShepherd_comparison.pdf']) 