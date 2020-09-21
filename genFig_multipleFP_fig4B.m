% genFig_multipleFP_fig4B.m
%
% 2020/09/11 SHP
% Multiple face patch partially overlapping population
%   making figure 4B: bar graph of average correlation for each ROI of average map for each recording site

clear all;

dirFig = '/projects/parksh/NeuroMRI/_labNote/_figs';

load('/procdata/parksh/_macaque/Art/Clustering_CorrMap_4FPs_Movie123_ArtRHROI_set01_probability.mat', 'Clustering_meanROI')

fid = fopen('/procdata/parksh/_macaque/Art/ROIs/i64_colorscale.pal');
A = fscanf(fid, '%s');
fclose(fid);
matRGB = sscanf(A(8:end), '#%2x%2x%2x', [3 inf])';
matRGB = flipud(unique(matRGB, 'rows', 'stable'));
matRGB = matRGB./255;

orderROI = [1 2 22 3 4 35 34 12 13 14 29 30 6 7 8 36 9 10 11 32 15 5 23 26 37 27 28 16 17 33 18 19 20 21 31 24 25]; % 1:37;


for iArea = 1:length(Clustering_meanROI.setArea)
    matR_ROI_area{iArea} = mean(Clustering_meanROI.matR(Clustering_meanROI.catAreaID==iArea, orderROI));
end

setArea_order = [4 1 2 3]; % ML-AF-AM-AAM

fig4B = figure;
set(fig4B, 'Color', 'w', 'PaperPositionMode', 'auto');
for iArea = 1:4
    idArea = setArea_order(iArea);
    sp(iArea) = subplot(4, 1, iArea);
    b = bar(matR_ROI_area{idArea});
    b.FaceColor = 'flat';
    b.CData = matRGB(orderROI, :);
    b.BarWidth = 1;
    b.EdgeColor = 'none';
end
set(sp(:), 'Box', 'off', 'TickDir', 'out')
set(sp(:), 'XColor', 'none')
axis(sp(:), 'tight')