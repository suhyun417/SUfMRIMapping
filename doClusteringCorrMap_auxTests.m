% doClusteringCorrMap_auxTests.m

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
% addpath(fullfile(dirProjects, 'parksh/_toolbox/Boot_Time_Series'))

% Set directories 
nameSubjNeural = 'Spi'; %'Tor';
nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataHome = fullfile(dirProcdata, 'parksh');
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% Directory for saving figures as graphic files
dirFig = fullfile(dirProjects, 'parksh/NeuralBOLD/_labNote/_figs/');

%% Load the data
load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'paramCorr', 'matR_SU') % get cell IDs from another file
load(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD))) % get cell IDs from another file
% load(fullfile(dirDataNeural, sprintf('ClusteringSDF_%sMovie123.mat', nameSubjNeural)))

%%
setK = Clustering.setK;

matWSS=[];
matExpVar=[];
for iK = 1:length(setK)
    curK = setK(iK);
    indClust = Clustering.resultKMeans(iK).SU_indCluster;
    [sortedClust, indSortedChan]=sort(indClust);
    
%     tExpVar=[];
%     for ii = 1:curK
%         tExpVar(ii,1) = Clustering.resultKMeans(iK).SU_sumD(ii)/(2*sum(sortedClust==ii));
%     end
    
    matWSS(iK,1) = sum(Clustering.resultKMeans(iK).SU_sumD);
%     matExpVar(iK,1) = sum(tExpVar);

end

[a, c, totalSS] = kmeans(matR_SU', 1);
betweenSS = totalSS-matWSS;
% totalVar = totalSS/(2*size(matR_SU,2)); %totalD/(2*size(matR_SU,2));

propExplained = (totalSS-matWSS)./totalSS; %matExpVar./totalSS;
figure;
plot(setK, propExplained, 'o-')

%%

% Silhouette plot
for i=1:length(Clustering_moviemask.resultKMeans)
figure(100);
[sil, h] = silhouette(matR_SU_moviemask', Clustering_moviemask.resultKMeans(i).SU_indCluster, 'sqEuclidean');
drawnow;
title(sprintf('K = %d: avg silh val = %2.4f', i+1, mean(sil)))
fprintf(1, '\nK = %d: avg silh val = %2.4f', i+1, mean(sil))
setMeanSil(i,1) = mean(sil);
input('')
end


opts = statset('Display','iter');
for iK = 1:length(Clustering.setK)
K = Clustering.setK(iK);
[IDX_SU, C, SUMD_SU, totalD] = kmeans(matR_SU', K, 'Replicates', Clustering.numReplicates, 'Options', opts);
cc(iK).ind = IDX_SU;
cc(iK).c = C;
cc(iK).sumd = SUMD_SU;
cc(iK).totald = totalD;
end

