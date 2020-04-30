% quickCompareMaps.m

nameSubjBOLD = 'Art';
setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi'}; 
% nameSubjNeural = 'Sig'; % 'Rho'; % 'Spi'; %'Tor';

numSubject = size(setNameSubjNeural, 2);
matR_SU_all_MION = cell([1 3]);
% valid voxels
load(sprintf('/procdata/parksh/%s/%s_MaskArrays.mat', nameSubjBOLD, nameSubjBOLD), 'movieDrivenAmp');
[nx ny nz] = size(movieDrivenAmp.mask_amp1);
nVox = nx*ny*nz;
moviemask_vec = reshape(movieDrivenAmp.mask_amp1, nVox, 1);

catChanID=[];
for typeMION = 1:3
    matR_SU_all=[];
    for iSubj = 1:numSubject
        nameSubjNeural = setNameSubjNeural{iSubj}; %'Spi'; %'Tor';
%         dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
        
        switch typeMION
            case 1
                load(sprintf('/procdata/parksh/%s/CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjNeural, nameSubjBOLD), ...
                    'matR_SU', 'paramCorr');
            otherwise
                load(sprintf('/procdata/parksh/%s/CorrMap_SU_%s%sMovie123_new_kernel%d.mat', nameSubjNeural, nameSubjNeural, nameSubjBOLD, typeMION), ...
                    'matR_SU', 'paramCorr');
        end
        
        % 2017/01/19: Excluded some units, in case of Spice's data which hasn't been went through this procedure yet
        switch lower(nameSubjNeural)
            case 'spi'
                excChanIndex = [10 13 22 27 30 49]; % cells were not same acrossd two days
                validChanIndex_clustering = setdiff(paramCorr.validChanIndex, excChanIndex);
                validChanID_clustering = cat(2, paramCorr.validChanID(validChanIndex_clustering,:), num2str(ones(size(validChanIndex_clustering)).*iSubj)); %paramCorr.validChanID(validChanIndex_clustering,:);
            otherwise
                validChanIndex_clustering = (1:length(paramCorr.validChanIndex))';
                validChanID_clustering = cat(2, paramCorr.validChanID, num2str(ones(size(paramCorr.validChanIndex)).*iSubj)); %paramCorr.validChanID;
        end
        
        matR_SU_valid = matR_SU(:, validChanIndex_clustering);
        matR_SU_all = cat(2, matR_SU_all, matR_SU_valid);
        
        clear matR_SU matR_SU_valid
        
        catChanID = cat(1, catChanID, validChanID_clustering);
    end
    matR_SU_all_MION{typeMION} = matR_SU_all(moviemask_vec,:); 
end

% Comparison

% 1. same cell, different kernels, map comparison
combMat = nchoosek([1 2 3], 2); % all possible pairs

rVal_BtwnKernels=NaN(size(matR_SU_all_MION{1}, 2), size(combMat, 1));
for iPair = 1:size(combMat, 1)
    for iChan = 1:size(matR_SU_all_MION{1}, 2)
        [r, p] = corr(matR_SU_all_MION{combMat(iPair, 1)}(:,iChan), matR_SU_all_MION{combMat(iPair, 2)}(:,iChan), 'rows', 'complete', 'type', 'Spearman');
        rVal_BtwnKernels(iChan, iPair) = r;
    end
end
figure;
plot(rVal_BtwnKernels)
mean(rVal_BtwnKernels)

% 2. across cells, for each kernel, map comparison 
for iK = 1:3
    [r, p] = corr(matR_SU_all_MION{iK}, 'rows', 'complete', 'type', 'Spearman');
    matR_acrossSUMaps{iK} = r;
end

% 2-1. matrix of correlation between pairs of single unit maps for each kernel
dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';

for iK=1:3
    figure;
    set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 575 545])
    imagesc(matR_acrossSUMaps{iK})
    set(gca, 'CLim', [-1 1])
    title(sprintf('Correlation between single unit maps: kernel %d', iK))
    print(gcf, fullfile(dirFig, sprintf('HRFComparison_CorrMatrix_kernel%d', iK)), '-dtiff', '-r120');
end



% 2-2. MDS plot for each kernel result
load(sprintf('/procdata/parksh/%s/Clustering_%s%sMovie123_new_masked_probability_critCorr1.mat', nameSubjBOLD, cell2mat(setNameSubjNeural), nameSubjBOLD),...
    'paramClustering*')  
matR_SU_all_MION_valid = cell([1 3]);
for iK = 1:3
    matR_SU_all_MION_valid{iK}  = matR_SU_all_MION{iK}(paramClustering_global.locValidVox, :);
    D = pdist(matR_SU_all_MION_valid{iK}', 'euclidean');
    [Y2,stress,disparities] = mdscale(D,2);
    % [Y3,stress,disparities] = mdscale(D,3);
    
    setD_kernel{iK} = D;
    setY2_kernel{iK} = Y2;
end

catChanID = cat(1, paramClustering.validChanID);
indMonkey = str2num(catChanID(:,5));
marker = {'o', '^', 'square', 'diamond'}; 


cMap = [        0         0    1.0000
         0    0.5000         0
    1.0000         0         0];

figMDS=figure;
set(figMDS, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1300 600 700 700])

% for iCell = 1:length(setY2_kernel{iK})
%     plot(setY2_kernel{iK}(iCell, 1), setY2_kernel{iK}(iCell, 2), 'o', 'Marker', marker{indMonkey(iCell)},...
%         'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 2)
%     hold on
% end
for iK = 1:3
    figure(figMDS);
    text(setY2_kernel{iK}(:, 1), setY2_kernel{iK}(:, 2), catChanID, 'Color', cMap(iK, :)) %paramCorr.validChanID(curChan,:))
    hold on
end
axis square
xlim([-38 38]) %xlim([-35 35]) %xlim([-25 30])
ylim([-15 15]) %ylim([-12 12])
% set(gca, 'XTick', [], 'YTick', [], 'LineWidth', 2,  'Box', 'off')


figure
plot(matRBetweenKernels)
mean(matRBetweenKernels)
[rOrgBetweenChan, p] = corr(matR_SU_org_moviemask, 'rows', 'complete', 'type', 'Spearman');
figure
imagesc(rOrgBetweenChan)
set(gca, 'CLim', [0 1])
a = tril(rOrgBetweenChan, -1);
mean(a(abs(a)>0))


