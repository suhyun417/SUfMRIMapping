% saveSDF.m
% SHP 11/26/17

dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';

nameSubjNeural = 'Mat'; % 'Dex';% 'Mat';
% Load data files
dirDataHome = '/procdata/parksh/_macaque';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);

filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
load(fullfile(dirDataNeural, filenameNeural), 'paramSDF') % to get the cell ID

setCellID = paramSDF.setCellIDs;
setMovie = [1 2 3];
% sampleCell = {'065a', '075a', '118a', '082a'}; %'/procdata/parksh/Tor/
% win = [0 300]; % seconds


%% Sample cell SDF to show diversity
sigma = 1000; % sigma (for SDF computation) in msec
S = createCellRegressor_indMov(dirDataNeural, setCellID, setMovie, sigma); %{'065a', '082a', '089a'}, 1);

tempC = bone(15); % jet(size(S(1).matsdf{1}, 2))*0.75; % temporary ColorOrder for line plot

figSDF = figure;
set(figSDF, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 940 640]); %
% case of five example neurons
startY = 0.09;
height = 0.24;
gap = 0.07;
for iCell = 1:length(setCellID)
    figure(figSDF); clf;
    
    for iMovie = 1:length(setMovie)
        SP(iMovie) = subplot('position', [0.1 startY+height*(iMovie-1)+gap*(iMovie-1) 0.8 height]);
    end
    set(SP, 'ColorOrder', tempC, 'NextPlot', 'replacechildren')
    for iMovie = 1:length(setMovie)
        plot(SP(iMovie), S(iCell, setMovie(iMovie)).matsdf{1}, 'LineWidth', 2)
        title(SP(iMovie), sprintf('Movie #%d, Cell %s', setMovie(iMovie), char(setCellID{iCell})))
    end
    
    set(SP, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 15)
    set(SP(2:3), 'XTickLabel', [])
    set(SP(1), 'XTickLabel', 0:50:300)
    xlabel(SP(1), 'Time (sec)')
    ylabel(SP(2), 'Spike / sec')
    
    input('')

%     % save as eps
%     print(figSDF, sprintf(fullfile(dirFig, '%s_Cell%sMovie%s_SDF'), nameSubjNeural, char(setCellID{iCell}), num2str(setMovie, '%d%d%d')), '-depsc')
    
end


