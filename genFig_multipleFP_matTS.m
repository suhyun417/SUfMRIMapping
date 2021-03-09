% genFig_multipleFP_matTS.m
%
% 2020/10/20 SHP
% Create a plot to show neuronal time series to movie and show example
% cells with similar time courses from different recording sites

load('/procdata/parksh/_macaque/matSDF_Movie123_allCells.mat', 'infoTS_subj', 'matTS_FP')

% % get the Spice cells with everything 
% locCell = find(cellfun(@numel, infoTS_subj(4).validChanID)<3);
% get only face selective cells
load('/procdata/parksh/_macaque/multipleFP_fsi.mat')
locFaceCell = find(abs(fsi.matFSI(:,1))>0.33);

matFR_TR = matTS_FP.matFR_TR(:, locFaceCell);
catAreaID = matTS_FP.catAreaID(locFaceCell);
catChanID = matTS_FP.catChanID(locFaceCell);

% cMap_Area = [91 148 203; 237 28 35; 248 148 29; 6 177 102]./255; % from Kenji's schematic
% cMap_Area(4, :) = cMap_Area(4, :).*0.7; % make the green a bit darker

cMap_Area = [194 165 207; 166 219 160; 0 136 55; 123 50 148]./255; %from colorbrewer2.org, diverging 4 classes
cMap_Area_temp = cat(1,cMap_Area, zeros(6, 3));

% matTS_norm = zscore(matTS_FP.matFR_TR);
matTS_norm = zscore(matFR_TR); % only face-selective neurons
[coeff, score, latent, tsq, explained] = pca(matTS_norm');
[sortedScore, indSortChan] = sort(score(:,1));
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

plot(zscore(matTS_FP.matFR_SU_1hz(:, setPairHighRDiffArea_indChan(1,:))))


ttt=load('/procdata/parksh/MovieRegressors/dbtmMriReg.mat'); % Face scale regressor (in TR unit)
scaleRGR = ttt.reg.xx(7,:)';

scaleRGR_norm = (scaleRGR-nanmean(scaleRGR))./nanstd(scaleRGR);

figure;
plot(1:2.4:900, zscore(coeff(:,1)), 'b')
hold on
plot(1:2.4:900, scaleRGR_norm, 'm')
xlabel('Time (s)')
ylabel('Normalized amplitude')
legend('1st PC', 'Face Area')
set(gca, 'TickDir', 'out', 'Box', 'off')
set(gcf, 'Color', 'w')


% Example cells with similar ts
[matR] = corr(matTS_FP.matFR_TR);
matR_uni = tril(matR, -1);
vectR_uni = matR_uni(:);
[sortedR, indPair] = sort(vectR_uni, 'descend');
[i, j] = ind2sub(size(matR_uni), indPair(1:10));
matTS_FP.catChanID([i j])



flagSM = 1; % flag for compression and smoothing
setMovie = [1 2 3];
fullRGR4fps = createMovieRGR_4fps_indMov(setMovie, flagSM); %createFullMovieRegressors_4fps_indMov(setMovID); %

% full regressors
catRGRfull=[];matRGRfull=[];
for iMov=1:length(setMovie)
    m = setMovie(iMov);
    matCurRGR = fullRGR4fps(m).smoRegressors; %fullRGR4fps(iMov).regressors(:,indValidRGR); %fullRGR4fps(iMov).regressors;
    catRGRfull = cat(1, catRGRfull, matCurRGR); % concatenation across movies
end
catRGRfull_TRresolution = resample(catRGRfull, 0.25*100, 2.4*100);


% scaleRGR_resampled = resample(scaleRGR, 2.4*100, 0.25*100); %matRGR = resample(catRGR, 0.25*100, 2.4*100);
matRGRfull = catRGRfull; %cat(2, catRGRfull, scaleRGR_resampled); %cat(2, matRGR, scaleRGR);
varnamesfull = fullRGR4fps(1).features; % cat(1, fullRGR4fps(1).features, {'Face size'});

% subset of regressors
indValidRGR = [1 2 6 7 3 21 20 32 22 31 25]; %[1, 3, 9, 20, 21, 22, 25]; 
% 1: 'Luminance', 2: 'Contrast', 6: Low spatial Frequency 7: High spatial frequencty 3: 'Motion (Speed)', 
% 21: 'One face', 20: 'Number of faces', 32: 'Face size', 22: 'Body parts', 31: 'Hands', 25: 'Any animal'
% matRGRvalid = matRGRfull(:,indValidRGR);
varnamesvalid = varnamesfull(indValidRGR);

% Compute correlation between neural TS and feature TS
R_ClusterMovieRGRfull=NaN(size(matRGRfull,2), size(meanFRCluster4fps,2));
R_SUmovieRGRfull=NaN(size(matRGRfull,2), size(matFR4fps,2));
for iRGR = 1:size(matRGRfull,2)
    r_c=[]; r_su=[];
    
%     % averaged TS in each cluster
%     r_c = corr(matRGRfull(:,iRGR), meanFRCluster4fps, 'rows', 'complete', 'type', 'Spearman');
    % single unit TS
    r_su = corr(matRGRfull(:,iRGR), matFR4fps, 'rows', 'complete', 'type', 'Spearman');
    
%     R_ClusterMovieRGRfull(iRGR, :) = r_c;
    R_SUmovieRGRfull(iRGR, :) = r_su;    
    
end

R_SUmovieRGRvalid = R_SUmovieRGRfull(indValidRGR,:);
R_ClusterMovieRGRvalid = R_ClusterMovieRGRfull(indValidRGR,:);

%%%%
% tempDMRGR = resample(matSizeRGR, 4, 30);  % down sample in 4 hz

% Version 1. upsampling from the original 71 features in 4fps (high-level) and 10fps (low-level)
[fullRGR30fps] = createFullMovieRegressors_30fps_indMov(setMovie);
catFullRGR = cat(1, fullRGR30fps.regressors);

% Version 2. upsampling from the 11 features selected features in 4fps
flagSM = 1; % flag for compression and smoothing
fullRGR4fps = createMovieRGR_4fps_indMov(setMovie, flagSM); %createFullMovieRegressors_4fps_indMov(setMovID); %

catRGRfull_30fps=[];matRGRfull_30fps=[];
for iMov=1:length(setMovie)
    m = setMovie(iMov);
    matCurRGR = fullRGR4fps(m).smoRegressors; %fullRGR4fps(iMov).regressors(:,indValidRGR); %fullRGR4fps(iMov).regressors;
    for iRGR = 1:size(matCurRGR, 2)
        tRGR = resample(matCurRGR(:,iRGR), 30, 4); 
        matCurRGR_30fps(:,iRGR) = tRGR;
    end
    catRGRfull_30fps = cat(1, catRGRfull_30fps, matCurRGR_30fps); %matCurRGR); % concatenation across movies
end
% scaleRGR_resampled = resample(scaleRGR, 2.4*100, 0.25*100); %matRGR = resample(catRGR, 0.25*100, 2.4*100);
matRGRfull_30fps = catRGRfull_30fps; %cat(2, catRGRfull, scaleRGR_resampled); %cat(2, matRGR, scaleRGR);
varnamesfull = fullRGR4fps(1).features; % cat(1, fullRGR4fps(1).features, {'Face size'});

% subset of regressors
indValidRGR = [1 2 6 7 3 21 20 32 22 31 25]; %[1, 3, 9, 20, 21, 22, 25]; 
% 1: 'Luminance', 2: 'Contrast', 6: Low spatial Frequency 7: High spatial frequencty 3: 'Motion (Speed)', 
% 21: 'One face', 20: 'Number of faces', 32: 'Face size', 22: 'Body parts', 31: 'Hands', 25: 'Any animal'
matRGRvalid_30fps = matRGRfull_30fps(:,indValidRGR);
varnamesvalid = varnamesfull(indValidRGR);

BRfeatureParams.version1.featureNames = fullRGR30fps(1).features;
BRfeatureParams.version1.notes = 'upsampled 71 regressors (including DM''s 5 size regressors) in 30fps from 4fps (high-level) and 10fps (low-level) using "createFullMovieRegressors_30fps_indMov.m"';
BRfeatureParams.version2.featureNames = varnamesvalid;
BRfeatureParams.version2.notes = 'upsampled 11 regressors in 30fps from 4fps';