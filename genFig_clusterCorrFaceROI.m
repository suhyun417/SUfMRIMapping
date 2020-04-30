% genFig_fig3.m


%% Settings
ss = pwd;
if ~isempty(strfind(ss, 'Volume')) % if it's local
    dirProjects = '/Volumes/PROJECTS';
    dirProcdata = '/Volumes/PROCDATA';
    dirLibrary = '/Volumes/LIBRARY';
else % on virtual machine
    dirProjects = '/projects';
    dirProcdata = '/procdata';
    dirLibrary = '/library';
end
    
% Add necessary toolbox % Should be 
addpath(fullfile(dirLibrary, 'matlab_utils')) % for convolution
addpath(fullfile(dirProjects, 'parksh/_toolbox/'))
addpath(fullfile(dirProjects, 'parksh/_toolbox/hslcolormap'))

% Set directories 
nameSubjNeural = 'Tor';
nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataHome = fullfile(dirProcdata, 'parksh');
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% Directory for saving figures as graphic files
dirFig = fullfile(dirProjects, 'parksh/NeuralBOLD/_labNote/_figs/');

% Load data
load(fullfile(dirDataNeural, 'ClusterCorrMap_FaceROI.mat'))

%
cMap = [0 0 0; 230 159 0; 86 180 233; 0 158 115; 240 228 66; 0 114 178; 213 94 0; 204 121 167]./255;
marker = {'o', '*', 'x', 's', 'd', '+', '^'};


%%
setCluster = ClusterCorr_faceROI.indNewClusterOrder;
meanCorr_ROI = ClusterCorr_faceROI.meanCorr_ROI;

matMeanCorr_selROI = mean(cat(3, meanCorr_ROI(:,ClusterCorr_faceROI.setROIs_right), meanCorr_ROI(:,ClusterCorr_faceROI.setROIs_left)), 3);


setROI = {'AM', 'AD', 'AL+AF', 'MF', 'ML', 'PL'};

fig_faceROICorr=figure;
set(fig_faceROICorr, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [600 600 410 500])
set(gca, 'ColorOrder', cMap); hold on;
plot(matMeanCorr_selROI', 'o-', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2)
hold on
set(gca, 'TickDir', 'out', 'Box', 'off')
line([0.5 6.5], [0 0], 'LineStyle', ':', 'Color', 'k')
set(gca, 'YLim', [-0.4 0.4], 'XLim', [0.5 6.5])
set(gca, 'FontSize', 15)
set(gca, 'XTick', 1:6, 'XTickLabel', setROI)
ylabel('Correlation (\rho)')

% fig_faceROICorr=figure;
% set(fig_faceROICorr, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [600 600 410 500])
% set(gca, 'ColorOrder', cMap); hold on;
% plot(meanCorr_ROI(:,ClusterCorr_faceROI.setROIs_right)', 'o-', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 2)
% hold on
% set(gca, 'TickDir', 'out', 'Box', 'off')
% line([0.5 6.5], [0 0], 'LineStyle', ':', 'Color', 'k')
% set(gca, 'YLim', [-0.4 0.4], 'XLim', [0.5 6.5])
% set(gca, 'FontSize', 15)
% set(gca, 'XTick', 1:6, 'XTickLabel', setROI)



