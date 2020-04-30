% test tSNE

addpath('/projects/parksh/_toolbox/bhtsne/')
addpath('/projects/parksh/_toolbox/tSNE/')

load('/procdata/parksh/Tor/CorrMap_SU_TorArtMovie123.mat', 'matR_SU')
load('/procdata/parksh/Tor/Clustering_TorArtMovie123_new_masked.mat') % get cell IDs from another file


% tSNE parameters
no_dims = 2;
initial_dims = 48; %30; %1000;
perplexity = 7; %30;
% theta = 0.5;

matIndClust_SU = cat(2, Clustering_moviemask.resultKMeans.SU_indCluster);

targetK = 7;
[indClust, indSortChan_tsne]=sort(matIndClust_SU(:,targetK-1));
[aa, bb] = sort(indSortChan_tsne);
labels = indClust(bb);

mappedX = tsne(matR_SU', labels, no_dims, 48, perplexity);

cMap = jet(targetK);
figure;
set(gcf, 'Color', 'w')
scatter(mappedX(:,1), mappedX(:,2), 50, cMap(labels, :), 'filled')

% diary 'tSNE_result.txt'
for iT = 1:20
    [tX, cost] = tsne(matR_SU', [], no_dims, initial_dims, perplexity);
    setX{iT} = tX;
end
% diary off
%     

% 
% [indClustSDFRGR, indSortChanSDFRGR_tsne]=sort(matIndClustSDFRGR_SU(:,targetK-1));
% [aaa, bbb] = sort(indSortChanSDFRGR_tsne);
% labelsSDFRGR = indClustSDFRGR(bbb);
% 
% mappedXSDFRGR = tsne(R_SUmovieRGR', labelsSDFRGR, no_dims, 30, perplexity);
% figure;
set(gcf, 'Color', 'w')
scatter(mappedXSDFRGR(:,1), mappedXSDFRGR(:,2), 50, cMap(labelsSDFRGR, :), 'filled')