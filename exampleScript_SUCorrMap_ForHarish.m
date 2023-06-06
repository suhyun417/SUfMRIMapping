% exampleScript_SUCorrMap.m
%
% 2023/06/05 SHP
%  - load correlation matrix between cells and functional ROIs
%  - lines to reproduce Park et al. (2022) Fig. S8 (D)
%
% corrMap_multipleFP.matR_SUfROI_org : a matrix containing Spearman correlations
% between responses of each single unit (n=389) and fMRI responses of each
% functional ROIs (n=37)
% 
% infoSU.catSubjID : monkey IDs for each single unit
% infoSU.catAreaID : recorded area IDs for each single unit (use it with infoSU.setArea)
% infoSU.setArea : set of names of recorded areas (use it with infoSU.catAreaID)
% infoSU.catChanID : channel IDs for each single unit% 
% infoSU.catCorrMapClusterID : cluster IDs for each single unit, based on
%                           their whole-brain correlation maps
% infoSU.orderSU_figS8 : order of single units as shown in Fig. S8 (D) of
%                       Park et al. (2022) Science Advances
% 
% infofROIs.nameROI : name of each funtional ROI
% infofROIs.orderROI_figS8 : order of funtional ROIs as shown in Fig. S8 (D) of 
%                       Park et al. (2022) Science Advances


clear all;

load('/nifvault/procdata/parksh/_macaque/_Harish/multipleFP_corrMapSUfROI.mat')

% replicate Fig.S8 panel D
fig3b2 = figure;
set(fig3b2, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 350 1000]);

sp1 = subplot('Position', [0.08 0.05 0.62 0.9]);
imagesc(corrMap_multipleFP.matR_SUfROI_org(infoSU.orderSU_figS8, infofROIs.orderROI_figS8))
set(gca, 'CLim', [-1 1].*0.5)
locDiff = find(abs(diff(infoSU.catCorrMapClusterID(infoSU.orderSU_figS8)))>0);
set(sp1, 'YTick', [], 'YTickLabel', []) 
set(sp1, 'XTick', []) 
% colormap(sp1, cMap_corrSUMA)
set(sp1,  'YColor', 'none', 'Box', 'off', 'XColor', 'none')
L = line(repmat([-0.8 37.5]', 1, length(locDiff)), [locDiff+0.5 locDiff+0.5]', 'Color', 'k',...
    'LineWidth', 2, 'LineStyle', ':');
set(sp1, 'XLim', [-0.8 37.5]);

sp2 = subplot('Position', [0.7 0.05 0.1 0.9]);
imagesc(infoSU.catAreaID(infoSU.orderSU_figS8)) 
set(sp2, 'XColor', 'none') 
set(sp2, 'YTick', locDiff+0.5, 'YTickLabel', [])
cMap_Area = [179 226 205; 141 160 203; 252 141 98; 231 41 138]./255; % 
colormap(sp2, cMap_Area)
set(sp2,  'YColor', 'none', 'Box', 'off', 'XColor', 'none', 'XTick', [])
L2 = line(repmat([0.5;1.5], 1, length(locDiff)), [locDiff+0.5 locDiff+0.5]', 'Color', 'k',...
    'LineWidth', 2, 'LineStyle', ':');
set(sp2, 'XLim', [0.5 1.5]);




