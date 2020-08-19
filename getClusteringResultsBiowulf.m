% getClusteringRestulsBiowulf.m
%
% 2020/08/14 SHP
% Load and merge clustering results across multiple *.mat files from
% parallel processing using Swarm in Biowulf

clear all;

%% Setting
nameSubjBOLD = 'Art';
dirDataOrg = sprintf('/procdata/parksh/_macaque/%s/clustering_multipleFP', nameSubjBOLD);
dirSaveData = sprintf('/procdata/parksh/_macaque/%s/', nameSubjBOLD);
saveFileName = 'Clustering_CorrMap_4FPs_Movie123_probability.mat';

setK = 2:40;

%% Retrieve necessary results into new struct and rename the struct
for iK = 1:length(setK)
    
    load(fullfile(dirDataOrg, sprintf('Clustering_CorrMap_4FPs_Movie123_probability_K%d.mat', setK(iK))));
    
    C_b.resultKMeans(iK) = Clustering_brainmask.resultKMeans;
    if isfield(Clustering_moviemask, 'resultKMeans') % for some (larger K) cases, the movie mask KMeans didn't finish
        C_m.resultKMeans(iK) = Clustering_moviemask.resultKMeans;
    else
        C_m.resultKMeans(iK).SU_indCluster = [];
        C_m.resultKMeans(iK).SU_sumD = [];
        C_m.resultKMeans(iK).Vox_indCluster = [];
        C_m.resultKMeans(iK).Vox_sumD = [];
    end
    
end

C_b.totalSS_SU = Clustering_brainmask.totalSS_SU;
C_b.totalSS_Vox = Clustering_brainmask.totalSS_Vox;
C_m.totalSS_SU = Clustering_moviemask.totalSS_SU;
C_m.totalSS_Vox = Clustering_moviemask.totalSS_Vox;

clear Clustering*

Clustering_brainmask = C_b;
Clustering_moviemask = C_m;

clear C_*

paramClustering_global.setK = setK;


load(sprintf('/procdata/parksh/_macaque/%s/matR4clusteringmultiplepatches.mat', nameSubjBOLD));

Clustering_brainmask.matR = matR_brainmask;
Clustering_moviemask.matR = matR_moviemask;

Clustering_brainmask.infoCells = infoCells;
Clustering_moviemask.infoCells = infoCells;

%% save
save(fullfile(dirSaveData, saveFileName), 'Clustering*', 'param*')