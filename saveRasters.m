% saveRasters.m
%
% save raster plots for each neuron
% 2017/10/30 SHP

dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';

nameSubjNeural = 'Mat';
% Load data files
dirDataHome = '/procdata/parksh/';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);

filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
load(fullfile(dirDataNeural, filenameNeural), 'paramSDF') % to get the cell ID

setCellID = paramSDF.setCellIDs;
% sampleCell = {'065a', '075a', '118a', '082a'}; %'/procdata/parksh/Tor/
win = [0 300]; % seconds
setMovie = [4 5 6]; %

figRaster = figure;
set(figRaster, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 940 640]) %[100 100 1200 200])
for iCell = 1:length(setCellID)
    figure(figRaster); clf;
    for iMov = 1:3
        
        idMov = setMovie(iMov);
        fName = fullfile(dirDataNeural, sprintf('%smov%dsig%s.mat', lower(nameSubjNeural(1)), idMov, char(setCellID{iCell})));
        load(fName)
        
        switch lower(dat.h.units)
            case 'sec'
                win = [0 300]; % seconds or milliseconds
            case 'ms'
                win = [0 300*1000];
        end
        spikes = {};
        for t=1:length( dat.s)
            ts = dat.s{t};
            spikes{t} = ts(find((ts>=win(1)) & ts<=win(2)));
        end
        
        figure(figRaster); %clf;
        subplot(3, 1, iMov) %subplot('Position', [0 0 1 1])   
        iS=1;
        set(gca,'YLim',[iS-1 iS]);
        line([spikes{iS} spikes{iS}]', repmat([iS-1 iS]', 1, length(spikes{iS})), 'Color', 'k')
        
        for iS=1:length(spikes)
            set(gca,'YLim',[iS-1 iS]);
            line([spikes{iS} spikes{iS}]', repmat([iS-1 iS]', 1, length(spikes{iS})), 'Color', 'k')
            hold on;
            %     yline(spikes{iS},prop,val);
        end
        
        set(gca,'YLim',[0 length(dat.s)]);
        set(gca,'XLim', win);
        set(gca,'YTick', [], 'XTick', 0:50000:300000, 'XTickLabel', 0:50:300)
        xlabel('Time (sec)')
        title(sprintf('Movie #%d, Cell %s', idMov, char(setCellID{iCell})))
        
    end
    % save as eps
    print(figRaster, sprintf(fullfile(dirFig, '%s_cell%s_movie%s'), nameSubjNeural, char(setCellID{iCell}), num2str(setMovie, '%d%d%d')), '-depsc')
    
end