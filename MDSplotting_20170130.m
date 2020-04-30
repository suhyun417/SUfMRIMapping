
setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi'};

for iSubj = 1:length(setNameSubjNeural)
    nameSubjNeural = setNameSubjNeural{iSubj};
    dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
    filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
    load(fullfile(dirDataNeural, filenameNeural))
    setMovie = [1 2 3];
    switch lower(nameSubjNeural)
        case 'spi'
            excChanIndex = [10 13 22 27 30 49]; % cells were not same acrossd two days
            validC = setdiff(1:length(paramSDF.setCellIDs), excChanIndex)';
        otherwise
            [indDataMat, CellID, movieID] = genDataMatrix_SU(nameSubjNeural, 0); % data matrix
            validC = find(indDataMat*ismember(movieID, setMovie)>0); % valid channel with movie [1 2 3]
    end
    indMovieNeuron = find(ismember(paramSDF.setMovIDs, setMovie)>0);
    nt = 375;
    matNeuralRGR = NaN(nt, length(validC));
    for iChan = 1:length(validC) % compute correlation channel-by-channel
        neuralrgrs=[];
        neuralrgrs = cat(1, S(validC(iChan), indMovieNeuron).mnFR);
        matNeuralRGR(:,iChan) = neuralrgrs';
    end
    tempMat(iSubj).matMnFR = matNeuralRGR;
end
tempMat
tempMat(1)
tempMat(2)
catMatSDF = cat(2, tempMat.matMnFR);
size(catMatSDF)
D = pdist(catMatSDF, 'euclidean');
[Y2,stress,disparities] = mdscale(D,2);
figure
plot(Y2(:,1), Y2(:,2), 'ro', 'LineWidth', 2, 'MarkerSize', 5)
text(Y2(:,1), Y2(:,2), cat(1, paramClustering.validChanID))
cat(1, paramClustering.validChanID)
size(Y2)
D = pdist(catMatSDF', 'euclidean');
[Y2,stress,disparities] = mdscale(D,2);
plot(Y2(:,1), Y2(:,2), 'ro', 'LineWidth', 2, 'MarkerSize', 5)
size(Y2)
text(Y2(:,1), Y2(:,2), cat(1, paramClustering.validChanID))
size(matR_SU_all_moviemask)
% title('2D MDS of correlation based on K=6 voxel clustering')
% figure(11)
% title('2D MDS of correlation based on K=13 voxel clustering')
% xlim([-1.3 1.3])
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% dirFig
% print(gcf, fullfile(dirFig, '2DMDS_CorrMap_13ROIs'), '-depsc')
% print(gcf, fullfile(dirFig, '2DMDS_CorrMap_6ROIs'), '-depsc')
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
title('2D MDS of time series in TR resolution (no MION convolution)')
size(catMatSDF)
help zscore
catMatSDF_norm = zscore(catMatSDF);
D = pdist(catMatSDF_norm', 'euclidean');
[Y2,stress,disparities] = mdscale(D,2);
figure
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
plot(Y2(:,1), Y2(:,2), 'ro', 'LineWidth', 2, 'MarkerSize', 5)
text(Y2(:,1), Y2(:,2), cat(1, paramClustering.validChanID))
title('2D MDS of time series in TR resolution (no MION convolution): centered')
catMatSDF_centered = catMatSDF - repmat(mean(catMatSDF), 375, 1);
size(catMatSDF_centered)
% D = pdist(catMatSDF_norm', 'euclidean');
D = pdist(catMatSDF_centered', 'euclidean');
[Y2,stress,disparities] = mdscale(D,2);
figure
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
plot(Y2(:,1), Y2(:,2), 'ro', 'LineWidth', 2, 'MarkerSize', 5)
text(Y2(:,1), Y2(:,2), cat(1, paramClustering.validChanID))
title('2D MDS of time series in TR resolution (no MION convolution): z-scored')
title('2D MDS of time series in TR resolution (no MION convolution): centered')
print(gcf, fullfile(dirFig, '2DMDS_SDF_TR_centered'), '-depsc')
print(gcf, fullfile(dirFig, '2DMDS_SDF_TR_zscored'), '-depsc')
D = pdist(matR_SU_all_moviemask', 'euclidean');
[Y2,stress,disparities] = mdscale(D,2);
figure
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
plot(Y2(:,1), Y2(:,2), 'ro', 'LineWidth', 2, 'MarkerSize', 5)
title('2D MDS of time series in TR resolution (no MION convolution): centered')
title('2D MDS of correlation based on movie-driven voxels')
text(Y2(:,1), Y2(:,2), cat(1, paramClustering.validChanID))