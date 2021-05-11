% genFig_multipleFP_matTS.m
%
% 2021/04/13 SHP
% Added parts to show time series of neurons from same electrode
%
% 2021/03/03 SHP
% Added parts to apply face-selectivity-index to draw plot using only
% face-selective cells
%
% 2020/10/20 SHP
% Create a plot to show neuronal time series to movie and show example
% cells with similar time courses from different recording sites

load('/procdata/parksh/_macaque/matSDF_Movie123_allCells.mat', 'infoTS_subj', 'matTS_FP')
dirFig = '/projects/parksh/NeuroMRI/_labNote/_figs';

% % get the Spice cells with everything 
% locCell = find(cellfun(@numel, infoTS_subj(4).validChanID)<3);
% get only face selective cells
load('/procdata/parksh/_macaque/multipleFP_fsi.mat')
locFaceCell =  1:389; %find(fsi.matFSI(:,1)>0.33); % find(abs(fsi.matFSI(:,1))>0.33);

matFR_TR = matTS_FP.matFR_TR(:, locFaceCell);
matTS_norm = zscore(matFR_TR); 
catAreaID = matTS_FP.catAreaID(locFaceCell);
catChanID = matTS_FP.catChanID(locFaceCell);

% Colormap for recording sites (in the order of "AF-pAM-aAM-ML")
% cMap_Area = [91 148 203; 237 28 35; 248 148 29; 6 177 102]./255; % from Kenji's schematic
% cMap_Area(4, :) = cMap_Area(4, :).*0.7; % make the green a bit darker

% cMap_Area = [194 165 207; 166 219 160; 0 136 55; 123 50 148]./255; %from colorbrewer2.org, diverging 4 classes
% cMap_Area =  [194 165 207; 123 50 148; 166 219 160; 0 136 55]./255; % swap pAM(2nd) and ML(4th)
% cMap_Area = [179 226 205; 252 141 98; 141 160 203; 231 41 138]./255; % from colorbrewer2.org: mixture of qualitive 3
cMap_Area = [179 226 205; 141 160 203; 252 141 98; 231 41 138]./255; % 

%% PCA
% matTS_norm = zscore(matTS_FP.matFR_TR);
[coeff, score, latent, tsq, explained] = pca(matTS_norm');
[sortedScore, indSortChan] = sort(score(:,1));

%% plot
% Fig S1a
fig_matTS = figure;
set(fig_matTS, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [500 30 550 1030]);
imagesc(matTS_norm(:, indSortChan)')
set(gca, 'XTick', 0:25:375, 'XTickLabel', 0:1:15, 'YTick', [])
colormap(hot); %jet)
set(gca, 'CLim', [-1 1].*5) %[0 5])
xlabel('time (m)')
ylabel('cells')
set(gca, 'YColor', 'none', 'TickDir', 'out', 'Box', 'off')

fig_matTS_sites = figure;
set(fig_matTS_sites, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [500 30 50 1030]);
imagesc(catAreaID(indSortChan)) %magesc(matTS_FP.catAreaID(indSortChan))
colormap(cMap_Area)
set(gca, 'Box', 'off')
axis off
% colormap(cMap_Area_temp)

% figure
% plot(zscore(matTS_FP.matFR_TR(:, indSortChan(389-10:389))))
% axis tight
% figure
% plot(zscore(matTS_FP.matFR_TR(:, indSortChan(1:10))))
% axis tight
% figure
% plot(coeff(:,1))
% axis tight

% Example cells with similar ts
[matR] = corr(matFR_TR, 'type', 'Spearman'); %corr(matTS_FP.matFR_TR);
matR_uni = tril(matR, -1);
vectR_uni = matR_uni(:);
[sortedR, indPair] = sort(vectR_uni, 'descend');
[i, j] = ind2sub(size(matR_uni), indPair(1:10));
catChanID([i j]) %matTS_FP.catChanID([i j])

[iset, jset] = ind2sub(size(matR_uni), indPair);
pairArea = catAreaID([iset jset]); %matTS_FP.catAreaID([iset jset]);
ttt = diff(pairArea, 1, 2);
setPairHighRDiffArea = find(abs(ttt)>0, 20); % pairs from different areas with high (top 10) correlation
[i, j] = ind2sub(size(matR_uni), indPair(setPairHighRDiffArea));  
setPairHighRDiffArea_indChan = [i j];
sortedR(setPairHighRDiffArea)
catChanID([i j]) %

for iPP = 11:20
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 665 145]);
P = plot(2.4:2.4:900, matTS_norm(:, setPairHighRDiffArea_indChan(iPP,:)), 'LineWidth', 2);
P(1).Color = cMap_Area(catAreaID(setPairHighRDiffArea_indChan(iPP,1)), :); %cMap_Area(matTS_FP.catAreaID(setPairHighRDiffArea_indChan(iPP,1)), :);
P(2).Color = cMap_Area(catAreaID(setPairHighRDiffArea_indChan(iPP,2)), :); %cMap_Area(matTS_FP.catAreaID(setPairHighRDiffArea_indChan(iPP,2)), :);
set(gca, 'TickDir', 'out', 'box', 'off')
% legend(catChanID(setPairHighRDiffArea_indChan(iPP,:))) %legend(matTS_FP.catChanID(setPairHighRDiffArea_indChan(iPP,:)))
print(gcf, fullfile(dirFig, sprintf('matTS_FR_TR_exampleTSPair_%d_%s_%s', ...
    iPP, catChanID{setPairHighRDiffArea_indChan(iPP,1)}, catChanID{setPairHighRDiffArea_indChan(iPP,2)})), '-depsc')
end

% Example cell with dissimilar time course from same area
iArea = 4; 
clear sortedR indPair
locArea = find(catAreaID==iArea);
[matR] = corr(matTS_norm(:,catAreaID==iArea), 'type', 'Spearman'); %corr(matTS_FP.matFR_TR);
matR_uni = tril(matR, -1);
vectR_uni = matR_uni(:);
[sortedR, indPair] = sort(abs(vectR_uni), 'ascend');
find(sortedR>0.03,1)
[i, j] = ind2sub(size(matR_uni), indPair(1:10));
catChanID(locArea([i j])) %matTS_FP.catChanID([i j])

setPairs = {'157AFMoc','138AFMoc'; '145AFMoc', '138AFMoc'; '232AFMoc', '137AFMoc'; ...
    '21AMWas', '14AMWas'; '112AMMoc', '111AMMoc'; '27Dan', '18Dan'; '17Dav', '11Dav'; '37Dav', '33Dav'};

for iPP = 1:size(setPairs, 1)
curCellID = setPairs(iPP,:);
tLoc = contains(catChanID, curCellID); 

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 665 145]);
P = plot(2.4:2.4:900, matTS_norm(:, tLoc), 'LineWidth', 1);
P(1).Color = cMap_Area(catAreaID(find(tLoc, 1)), :); %
P(2).Color = cMap_Area(catAreaID(find(tLoc, 1)), :).*0.5;
% set(P(:), 'Color', cMap_Area(catAreaID(find(tLoc, 1)), :)); %cMap_Area(matTS_FP.catAreaID(setPairHighRDiffArea_indChan(iPP,1)), :);
% P(2).LineStyle = '-'; %cMap_Area(catAreaID(setPairHighRDiffArea_indChan(iPP,2)), :); %cMap_Area(matTS_FP.catAreaID(setPairHighRDiffArea_indChan(iPP,2)), :);
set(gca, 'TickDir', 'out', 'box', 'off')
% legend(catChanID(setPairHighRDiffArea_indChan(iPP,:))) %legend(matTS_FP.catChanID(setPairHighRDiffArea_indChan(iPP,:)))
print(gcf, fullfile(dirFig, sprintf('matTS_FR_TR_exampleTSPair_sameFPdiffCell_%s_%s', ...
    curCellID{1}, curCellID{2})), '-depsc')
end



%% example 3 cells with similar time course 
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 665 145]);
P = plot(2.4:2.4:900, matTS_norm(:, [153 74 141]), 'LineWidth', 1); %2);
P(1).Color = cMap_Area(catAreaID(153), :); %cMap_Area(matTS_FP.catAreaID(setPairHighRDiffArea_indChan(iPP,1)), :);
P(2).Color = cMap_Area(catAreaID(74), :); %cMap_Area(matTS_FP.catAreaID(setPairHighRDiffArea_indChan(iPP,2)), :);
P(3).Color = cMap_Area(catAreaID(141), :); 
set(gca, 'TickDir', 'out', 'box', 'off')
% % legend(catChanID(setPairHighRDiffArea_indChan(iPP,:))) %legend(matTS_FP.catChanID(setPairHighRDiffArea_indChan(iPP,:)))
% print(gcf, fullfile(dirFig, sprintf('matTS_FR_TR_exampleTSPair_width1_%s_%s_%s', ...
%     catChanID{153}, catChanID{74}, catChanID{141})), '-depsc')


set(gca, 'XColor', 'none')
ylim([-4 4])
L = line([5 65], [-1.8 -1.8], 'Color', 'k', 'LineWidth', 3);
xlim([-5 900])
print(gcf, fullfile(dirFig, sprintf('matTS_FR_TR_exampleTSPair_width1_1minScaleBar_newColor_%s_%s_%s', ...
    catChanID{153}, catChanID{74}, catChanID{141})), '-depsc')

%% Raster
setIndCell = [153 74 141 200 26];
for iCell = 1:length(setIndCell)
    curCellID = catChanID{setIndCell(iCell)};
    nameSubjNeural = char(curCellID(end-2:end));
    chanID = char(curCellID(1:end-3));
    iMov = 1;
    
    load(sprintf('/procdata/parksh/_macaque/%s/%smov%dsig%s.mat', ...
        nameSubjNeural, lower(nameSubjNeural), iMov, chanID))
    % danmov1sig10.mat
    % mocmov1sig232AF.mat
    % wasmov1sig47AM.mat
    
    % Plotting raster
    fig_raster = figure;
    set(fig_raster, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 665 145]);
    ax_raster = subplot('Position', [0 0 1 1]);
    hold on;
    
    locValidTrial = find(cellfun(@isempty, dat.s)<1);
    for i=1:length(locValidTrial)
        indTrial = locValidTrial(i);
        if length(dat.s{indTrial})>1000
            indx = 1:round(length(dat.s{indTrial})/1000):length(dat.s{indTrial});
            plotdat_x = repmat(dat.s{indTrial}(indx)',[3 1]);
        else
            plotdat_x = repmat(dat.s{indTrial}',[3 1]);
        end
        plotdat_y = repmat([i-1; i; NaN],[1 size(plotdat_x,2)]);
        plot(ax_raster, plotdat_x(:),plotdat_y(:),'-','LineWidth',0.01, 'Color', cMap_Area(catAreaID(setIndCell(iCell)), :));
    end
    xlim([0 300000]);
    ylim([0 length(locValidTrial)]);
    axis off
    box off
    print(fig_raster, fullfile(dirFig, sprintf('multipleFP_FigS_raster_exampleCell_%s_movie%d_newColor', curCellID, iMov)), '-depsc');
end


%% Cells from the same electrode
% Which units?
setExampleCellIDs = { '613aSig', '613cSig'; '030aSpi','030cSpi'; '037aSpi', '037cSpi'};
% setExampleCellIDs = {'013aSig', '013bSig'; '613aSig', '613bSig'; '003aSpi', '003bSpi'; '008aSpi', '008bSpi'; ...
%     '011aSpi', '011bSpi'; '015aSpi', '015bSpi'; '020aSpi', '020bSpi'; '022aSpi', '022bSpi'; '025aSpi', '025bSpi'; ...
%     '030aSpi', '030bSpi'; '031aSpi', '031bSpi'; '033aSpi', '033bSpi'; '034aSpi', '034bSpi'; '036aSpi', '036bSpi'; ...
%     '037aSpi', '037bSpi'; '042aSpi', '042bSpi'; '044aSpi', '044bSpi'; '054aSpi', '054bSpi'};
% setExampleCellIDs = { '613aSig', '613bSig', '613cSig'; '030aSpi', '030bSpi', '030cSpi'; '037aSpi', '037bSpi', '037cSpi'};
% setExampleCellIDs = {'082aTor', '082bTor'; '085aTor','085bTor'; '086aTor', '086bTor'; '088aTor', '088bTor'; '089aTor', '089bTor';...
%       '091aTor', '091bTor'; '096aTor', '096bTor'; '122aTor', '122bTor'};
% setExampleCellIDs = {'17AMWas', '18AMWas'; '19AMWas', '20AMWas'; '22AMWas', '23AMWas'; '28AMWas', '29AMWas'; ...
%     '31AMWas', '32AMWas'; '33AMWas', '34AMWas'; '38AMWas', '39AMWas'; '41AMWas', '42AMWas'; '48AMWas', '49AMWas'};
% setExampleCellIDs = {'108AMMoc', '109AMMoc'; '130AFMoc', '131AFMoc'; '217AFMoc', '218AFMoc';...
%     '242AMMoc', '243AMMoc'; '244AMMoc', '245AMMoc'; '249AMMoc', '250AMMoc'; '257AMMoc', '258AMMoc'; ...
%     '264AMMoc', '265AMMoc'};
% setExampleCellIDs = {'14Spi', '15Spi'; '23Spi', '24Spi'; '42Spi', '43Spi'; '46Spi', '47Spi'}; %
% setExampleCellIDs = {'29Dav', '30Dav'; '130AFMoc', '131AFMoc'; '28AMWas', '29AMWas'; '108AMMoc', '109AMMoc'};
matTS_norm = zscore(matTS_FP.matFR_TR); 

for iPair = 1:size(setExampleCellIDs, 1)
    
    curCellID = setExampleCellIDs(iPair, :); % two cells from one electrode
    tLoc = contains(matTS_FP.catChanID, curCellID); 
    
    figure;
    set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 665 145]);
%     P = plot(2.4:2.4:900, matTS_FP.matFR_TR(:, tLoc), 'LineWidth', 1);
    P = plot(2.4:2.4:900, matTS_norm(:, tLoc), 'LineWidth', 1);
    P(1).Color = ones(1, 3).*0; %cMap_Area(catAreaID(153), :); %cMap_Area(matTS_FP.catAreaID(setPairHighRDiffArea_indChan(iPP,1)), :);
    P(2).Color = ones(1, 3).*0.6; %cMap_Area(catAreaID(74), :); %cMap_Area(matTS_FP.catAreaID(setPairHighRDiffArea_indChan(iPP,2)), :);
    set(gca, 'TickDir', 'out', 'box', 'off')
    % legend(catChanID(setPairHighRDiffArea_indChan(iPP,:))) %legend(matTS_FP.catChanID(setPairHighRDiffArea_indChan(iPP,:)))
    print(gcf, fullfile(dirFig, sprintf('matTS_FR_TR_exampleTSPair_width1_%s_%s', ...
        curCellID{1}, curCellID{2})), '-depsc')
    
%     %% Raster
%     setIndCell = [153 74 141 200 26];
%     for iCell = 1:length(setIndCell)
%         curCellID = catChanID{setIndCell(iCell)};
%         nameSubjNeural = char(curCellID(end-2:end));
%         chanID = char(curCellID(1:end-3));
%         iMov = 1;
%         
%         load(sprintf('/procdata/parksh/_macaque/%s/%smov%dsig%s.mat', ...
%             nameSubjNeural, lower(nameSubjNeural), iMov, chanID))
        % danmov1sig10.mat
        % mocmov1sig232AF.mat
        % wasmov1sig47AM.mat
        
%         % Plotting raster
%         fig_raster = figure;
%         set(fig_raster, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 665 145]);
%         ax_raster = subplot('Position', [0 0 1 1]);
%         hold on;
%         
%         locValidTrial = find(cellfun(@isempty, dat.s)<1);
%         for i=1:length(locValidTrial)
%             indTrial = locValidTrial(i);
%             if length(dat.s{indTrial})>1000
%                 indx = 1:round(length(dat.s{indTrial})/1000):length(dat.s{indTrial});
%                 plotdat_x = repmat(dat.s{indTrial}(indx)',[3 1]);
%             else
%                 plotdat_x = repmat(dat.s{indTrial}',[3 1]);
%             end
%             plotdat_y = repmat([i-1; i; NaN],[1 size(plotdat_x,2)]);
%             plot(ax_raster, plotdat_x(:),plotdat_y(:),'-','LineWidth',0.01, 'Color', cMap_Area(catAreaID(setIndCell(iCell)), :));
%         end
%         xlim([0 300000]);
%         ylim([0 length(locValidTrial)]);
%         axis off
%         box off
%         print(fig_raster, fullfile(dirFig, sprintf('multipleFP_FigS_raster_exampleCell_%s_movie%d', curCellID, iMov)), '-depsc');
%     end
    
end



%%
% plot(zscore(matTS_FP.matFR_SU_1hz(:, setPairHighRDiffArea_indChan(1,:))))
% 
% 
% ttt=load('/procdata/parksh/MovieRegressors/dbtmMriReg.mat'); % Face scale regressor (in TR unit)
% scaleRGR = ttt.reg.xx(7,:)';
% 
% scaleRGR_norm = (scaleRGR-nanmean(scaleRGR))./nanstd(scaleRGR);
% 
% figure;
% plot(1:2.4:900, zscore(coeff(:,1)), 'b')
% hold on
% plot(1:2.4:900, scaleRGR_norm, 'm')
% xlabel('Time (s)')
% ylabel('Normalized amplitude')
% legend('1st PC', 'Face Area')
% set(gca, 'TickDir', 'out', 'Box', 'off')
% set(gcf, 'Color', 'w')
% 
% 
% % Example cells with similar ts
% [matR] = corr(matTS_FP.matFR_TR);
% matR_uni = tril(matR, -1);
% vectR_uni = matR_uni(:);
% [sortedR, indPair] = sort(vectR_uni, 'descend');
% [i, j] = ind2sub(size(matR_uni), indPair(1:10));
% matTS_FP.catChanID([i j])
% 
% 
% 
% flagSM = 1; % flag for compression and smoothing
% setMovie = [1 2 3];
% fullRGR4fps = createMovieRGR_4fps_indMov(setMovie, flagSM); %createFullMovieRegressors_4fps_indMov(setMovID); %
% 
% % full regressors
% catRGRfull=[];matRGRfull=[];
% for iMov=1:length(setMovie)
%     m = setMovie(iMov);
%     matCurRGR = fullRGR4fps(m).smoRegressors; %fullRGR4fps(iMov).regressors(:,indValidRGR); %fullRGR4fps(iMov).regressors;
%     catRGRfull = cat(1, catRGRfull, matCurRGR); % concatenation across movies
% end
% catRGRfull_TRresolution = resample(catRGRfull, 0.25*100, 2.4*100);
% 
% 
% % scaleRGR_resampled = resample(scaleRGR, 2.4*100, 0.25*100); %matRGR = resample(catRGR, 0.25*100, 2.4*100);
% matRGRfull = catRGRfull; %cat(2, catRGRfull, scaleRGR_resampled); %cat(2, matRGR, scaleRGR);
% varnamesfull = fullRGR4fps(1).features; % cat(1, fullRGR4fps(1).features, {'Face size'});
% 
% % subset of regressors
% indValidRGR = [1 2 6 7 3 21 20 32 22 31 25]; %[1, 3, 9, 20, 21, 22, 25]; 
% % 1: 'Luminance', 2: 'Contrast', 6: Low spatial Frequency 7: High spatial frequencty 3: 'Motion (Speed)', 
% % 21: 'One face', 20: 'Number of faces', 32: 'Face size', 22: 'Body parts', 31: 'Hands', 25: 'Any animal'
% % matRGRvalid = matRGRfull(:,indValidRGR);
% varnamesvalid = varnamesfull(indValidRGR);
% 
% % Compute correlation between neural TS and feature TS
% R_ClusterMovieRGRfull=NaN(size(matRGRfull,2), size(meanFRCluster4fps,2));
% R_SUmovieRGRfull=NaN(size(matRGRfull,2), size(matFR4fps,2));
% for iRGR = 1:size(matRGRfull,2)
%     r_c=[]; r_su=[];
%     
% %     % averaged TS in each cluster
% %     r_c = corr(matRGRfull(:,iRGR), meanFRCluster4fps, 'rows', 'complete', 'type', 'Spearman');
%     % single unit TS
%     r_su = corr(matRGRfull(:,iRGR), matFR4fps, 'rows', 'complete', 'type', 'Spearman');
%     
% %     R_ClusterMovieRGRfull(iRGR, :) = r_c;
%     R_SUmovieRGRfull(iRGR, :) = r_su;    
%     
% end
% 
% R_SUmovieRGRvalid = R_SUmovieRGRfull(indValidRGR,:);
% R_ClusterMovieRGRvalid = R_ClusterMovieRGRfull(indValidRGR,:);
% 
% %%%%
% % tempDMRGR = resample(matSizeRGR, 4, 30);  % down sample in 4 hz
% 
% % Version 1. upsampling from the original 71 features in 4fps (high-level) and 10fps (low-level)
% [fullRGR30fps] = createFullMovieRegressors_30fps_indMov(setMovie);
% catFullRGR = cat(1, fullRGR30fps.regressors);
% 
% % Version 2. upsampling from the 11 features selected features in 4fps
% flagSM = 1; % flag for compression and smoothing
% fullRGR4fps = createMovieRGR_4fps_indMov(setMovie, flagSM); %createFullMovieRegressors_4fps_indMov(setMovID); %
% 
% catRGRfull_30fps=[];matRGRfull_30fps=[];
% for iMov=1:length(setMovie)
%     m = setMovie(iMov);
%     matCurRGR = fullRGR4fps(m).smoRegressors; %fullRGR4fps(iMov).regressors(:,indValidRGR); %fullRGR4fps(iMov).regressors;
%     for iRGR = 1:size(matCurRGR, 2)
%         tRGR = resample(matCurRGR(:,iRGR), 30, 4); 
%         matCurRGR_30fps(:,iRGR) = tRGR;
%     end
%     catRGRfull_30fps = cat(1, catRGRfull_30fps, matCurRGR_30fps); %matCurRGR); % concatenation across movies
% end
% % scaleRGR_resampled = resample(scaleRGR, 2.4*100, 0.25*100); %matRGR = resample(catRGR, 0.25*100, 2.4*100);
% matRGRfull_30fps = catRGRfull_30fps; %cat(2, catRGRfull, scaleRGR_resampled); %cat(2, matRGR, scaleRGR);
% varnamesfull = fullRGR4fps(1).features; % cat(1, fullRGR4fps(1).features, {'Face size'});
% 
% % subset of regressors
% indValidRGR = [1 2 6 7 3 21 20 32 22 31 25]; %[1, 3, 9, 20, 21, 22, 25]; 
% % 1: 'Luminance', 2: 'Contrast', 6: Low spatial Frequency 7: High spatial frequencty 3: 'Motion (Speed)', 
% % 21: 'One face', 20: 'Number of faces', 32: 'Face size', 22: 'Body parts', 31: 'Hands', 25: 'Any animal'
% matRGRvalid_30fps = matRGRfull_30fps(:,indValidRGR);
% varnamesvalid = varnamesfull(indValidRGR);
% 
% BRfeatureParams.version1.featureNames = fullRGR30fps(1).features;
% BRfeatureParams.version1.notes = 'upsampled 71 regressors (including DM''s 5 size regressors) in 30fps from 4fps (high-level) and 10fps (low-level) using "createFullMovieRegressors_30fps_indMov.m"';
% BRfeatureParams.version2.featureNames = varnamesvalid;
% BRfeatureParams.version2.notes = 'upsampled 11 regressors in 30fps from 4fps';