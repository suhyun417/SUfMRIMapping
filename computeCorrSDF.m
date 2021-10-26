% computeCorrSDF.m
%
% 2017/02/28 SHP
% compute correlation bewteen neuronal time courses
% across neurons within each cluster

setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi'};
nameSubjNeural = 'Spi';
dirDataHome = '/procdata/parksh/';
dirDataNeural = fullfile('/procdata/parksh/', nameSubjNeural);

load(fullfile(dirDataNeural, sprintf('matSDF_%sMovie123.mat', cell2mat(setNameSubjNeural))), 'matSDF')
matFR_SU = cat(2, matSDF.matFR_SU);
matFR_SU_norm = cat(2, matSDF.matFR_SU_norm);
matFR_TR = cat(2, matSDF.matFR_TR);
matNeuralRGR = cat(2, matSDF.matNeuralRGR);


load(fullfile(dirDataNeural, 'Clustering_SU_allCells_validVoxels_critCorr1_7Means_prob.mat'))

orderIndCluster = [4 2 3 6 7 5 1];%
for iK = 1:length(orderIndCluster)
    [r, p] = corr(matFR_SU(:, resultProbClustering(orderIndCluster(iK)).validIndCells), 'rows', 'complete', 'type', 'Spearman');
    temp = tril(r, -1);
    setR = temp(abs(temp)>0);
    FR_SU.setRSDF(iK).matR_org = r;
    FR_SU.setRSDF(iK).setR = setR;
    FR_SU.setRSDF(iK).meanR = mean(setR);
    FR_SU.setRSDF(iK).medianR = median(setR);
    
    [r, p] = corr(matFR_SU_norm(:, resultProbClustering(orderIndCluster(iK)).validIndCells), 'rows', 'complete', 'type', 'Spearman');
    temp = tril(r, -1);
    setR = temp(abs(temp)>0);
    FR_SU_norm.setRSDF(iK).matR_org = r;
    FR_SU_norm.setRSDF(iK).setR = setR;
    FR_SU_norm.setRSDF(iK).meanR = mean(setR);
    FR_SU_norm.setRSDF(iK).medianR = median(setR);
    
    [r, p] = corr(matFR_TR(:, resultProbClustering(orderIndCluster(iK)).validIndCells), 'rows', 'complete', 'type', 'Spearman');
    temp = tril(r, -1);
    setR = temp(abs(temp)>0);
    FR_TR.setRSDF(iK).matR_org = r;
    FR_TR.setRSDF(iK).setR = setR;
    FR_TR.setRSDF(iK).meanR = mean(setR);
    FR_TR.setRSDF(iK).medianR = median(setR);
    
    [r, p] = corr(matNeuralRGR(:, resultProbClustering(orderIndCluster(iK)).validIndCells), 'rows', 'complete', 'type', 'Spearman');
    temp = tril(r, -1);
    setR = temp(abs(temp)>0);
    FR_TR_MION.setRSDF(iK).matR_org = r;
    FR_TR_MION.setRSDF(iK).setR = setR;
    FR_TR_MION.setRSDF(iK).meanR = mean(setR);
    FR_TR_MION.setRSDF(iK).medianR = median(setR);
    
    
%     figure(100);
%     subplot(7,1,iK)
%     plot(zscore(matFR_SU(:, resultProbClustering(orderIndCluster(iK)).validIndCells)))
end

validIndCells = cat(1, resultProbClustering(orderIndCluster).validIndCells);
matFR_SU_reorder = matFR_SU(:, validIndCells);
matFR_SU_norm_reorder = matFR_SU_norm(:, validIndCells);
matFR_TR_reorder = matFR_TR(:, validIndCells);
matNeuralRGR_reorder = matNeuralRGR(:, validIndCells);


[r_5_7, p] = corr(matFR_SU(:, resultProbClustering(orderIndCluster(5)).validIndCells), matFR_SU(:, resultProbClustering(orderIndCluster(7)).validIndCells),...
    'rows', 'complete', 'type', 'Spearman');

[r_mean_5_7, p] = corr(mean(matFR_SU(:, resultProbClustering(orderIndCluster(5)).validIndCells), 2), mean(matFR_SU(:, resultProbClustering(orderIndCluster(7)).validIndCells), 2),...
    'rows', 'complete', 'type', 'Spearman');

%% raster plots of cell group 2
setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi'};
dirDataHome = '/procdata/parksh/';

indCell_monkey{1} = 1:48;
indCell_monkey{2} = 49:53;
indCell_monkey{3} = 54:69;
indCell_monkey{4} = 70:135;
% indCell_monkey(1,:) = [1 49 55 70];
% indCell_monkey(2,:) = [48 54 69 135];

load('/procdata/parksh/Spi/Clustering_SU_allCells_validVoxels_critCorr1_7Means_prob.mat')

orderIndCluster = [4 2 3 6 7 5 1];%
analTrialCluster = struct([]);
for iCG = 1:length(orderIndCluster) %2;
 
    curIndCells = resultProbClustering(orderIndCluster(iCG)).validIndCells;
    
    cellRaster_cluster = struct([]);
    for iCell = 1: length(curIndCells)
        
        curCell = curIndCells(iCell);
        
        for iMonk = 1:4
            if sum(ismember(curCell, indCell_monkey{iMonk}))<1
                continue;
            else
                indMonkey = iMonk;
                [a, indCell] = ismember(curCell, indCell_monkey{iMonk});
            end
        end
%         [i, j] = find(curCell>indCell_monkey);
%         if isempty(j)
%             j = 1;
%         end
%         indMonkey = max(j);
        
        % Load the data
        nameSubjNeural = setNameSubjNeural{indMonkey};
        dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
        filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
        load(fullfile(dirDataNeural, filenameNeural))
        fprintf(1, 'Loading single unit data of %s: %s .... \n', nameSubjNeural, filenameNeural)
        
        setMovie=[1 2 3];
        switch lower(nameSubjNeural)
            case 'spi'
                excChanIndex = [10 13 22 27 30 49]; % cells were not same acrossd two days
                validC = setdiff(1:length(paramSDF.setCellIDs), excChanIndex)';
            otherwise
                [indDataMat, CellID, movieID] = genDataMatrix_SU(nameSubjNeural, 0); % data matrix
                validC = find(indDataMat*ismember(movieID, setMovie)>0); % valid channel with movie [1 2 3]
        end
        
        tempCat = cat(1, S(validC(indCell), setMovie).matFR);
        ntrial = [];
        for iM = setMovie
            ntrial(iM,1) = size(tempCat{iM}, 2);
        end
        
        minNumTrial = min(ntrial);
        
        curCellRaster = [];
        for iM = setMovie
            curCellRaster = cat(1, curCellRaster, tempCat{iM}(:, 1:minNumTrial));
        end
        
        [r, p] = corr(curCellRaster, 'rows', 'complete', 'type', 'Spearman');
        temp = tril(r, -1);
        setR = temp(abs(temp)>0);
        
        cellRaster_cluster(iCell).curCellRaster = curCellRaster;
        cellRaster_cluster(iCell).meanR_acrossTrials = mean(setR);
    end
    
    analTrialCluster(iCG).cellRaster_cluster = cellRaster_cluster;
    analTrialCluster(iCG).setMeanR_acrossTrials = cat(1, cellRaster_cluster.meanR_acrossTrials);
    analTrialCluster(iCG).meanAvgR_acrossTrials = mean(cat(1, cellRaster_cluster.meanR_acrossTrials));
end

[catSetR{1:7}] = deal(analTrialCluster.setMeanR_acrossTrials);



