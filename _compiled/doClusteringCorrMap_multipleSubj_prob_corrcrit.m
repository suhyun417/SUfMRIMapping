function [] = doClusteringCorrMap_multipleSubj_prob_corrcrit(flagParallel, flagSave)
% 2017/01/25 SHP
% Written for Biowulf environment
% Apply K-means clustering to correlation map of single units from multiple
% subjects
% Modified from 'doClusteringCorrMap_voxel_multipleSubjects.m'


% %% Settings
% ss = pwd;
% if ~isempty(strfind(ss, 'Volume')) % if it's local
%     dirProjects = '/Volumes/PROJECTS';
%     dirProcdata = '/Volumes/PROCDATA';
%     dirLibrary = '/Volumes/LIBRARY';
% else % on virtual machine
%     dirProjects = '/projects';
%     dirProcdata = '/procdata';
%     dirLibrary = '/library';
% end
    
% % Add necessary toolbox 
% addpath(fullfile(dirLibrary, 'matlab_utils')) % for convolution

% Set directories 
setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi'};
nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataHome = '/data/parks20/procdata/NeuroMRI/'; %fullfile(dirProcdata, 'parksh');
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% load the masks (movie-driven & brain-only mask)
load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp'); %, 'brainMask_BlockAna3D');


%% Load the corr map and select valid channel: subject by subject
numSubject = size(setNameSubjNeural, 2);
matR_SU_all = [];

for iSubj = 1:numSubject
    nameSubjNeural = setNameSubjNeural{iSubj}; %'Spi'; %'Tor';
    dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
    load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'matR_SU', 'paramCorr');
    
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
    
    paramClustering(iSubj).nameSubj = nameSubjNeural;
    paramClustering(iSubj).validChanIndex = validChanIndex_clustering;
    paramClustering(iSubj).validChanID = validChanID_clustering;
    
    matR_SU_valid = matR_SU(:, validChanIndex_clustering);
    matR_SU_all = cat(2, matR_SU_all, matR_SU_valid);
    
     clear matR_SU matR_SU_valid
end


%% Clustering on masked map 
numRepeat = 100; % number of repetition for entire clustering

setK = 2:20; %15;
opts = statset('Display','final');
numReplicates = 5; %


paramClustering_global.methods = 'KMeans';
paramClustering_global.setK = setK;
paramClustering_global.numReplicates = numReplicates;
paramClustering_global.descriptions = 'Clustering of neurons based on the masked whole-brain correlation maps: multiple repetitions (numRepeat) of an execution of kmeans function, which had its own multiple "replicates" (numReplicates)';
paramClustering_global.numRepeat = numRepeat;

if flagParallel
    pool = parpool;                      % Invokes workers
    stream = RandStream('mlfg6331_64');  % Random number stream
    opts = statset('UseParallel',1,'UseSubstreams',1, 'Streams',stream,...
        'MaxIter', 1000, 'Display','final');
    paramClustering_global.parallel = opts;
end


[nx ny nz] = size(movieDrivenAmp.mask_amp1);
nVox = nx*ny*nz;

% Apply movie-driven mask to correlation matrix 
moviemask_vec = reshape(movieDrivenAmp.mask_amp1, nVox, 1); % change the 3D mask to 1D
matR_SU_all_moviemask = matR_SU_all(moviemask_vec,:); %matR_SU(moviemask_vec,:); % 15495 voxels

% Select voxels: voxels that showed correlation higher than 0.3 with any one
% of the neurons
critCorr = 0.3; 
critNumNeuron = 6; %13;
matValidVox = abs(matR_SU_all_moviemask)>critCorr;
locValidVox = sum(matValidVox, 2)>critNumNeuron;
matR_SU_all_moviemask_valid  = matR_SU_all_moviemask(locValidVox, :);


paramClustering_global.critCorr = critCorr;
paramClustering_global.critNumNeuron = critNumNeuron;
paramClustering_global.locValidVox = locValidVox;

% brainmask_vec = reshape(movieDrivenAmp.map_sm_brain>0, nVox, 1); % change the 3D mask to 1D
% matR_SU_all_brainmask = matR_SU_all(brainmask_vec,:); %matR_SU(brainmask_vec,:); % 27113 voxels

% % totalSS
[a, c, totalSS] = kmeans(matR_SU_all_moviemask_valid', 1); % single unit
Clustering_moviemask_valid.totalSS_SU = totalSS;
[a, c, totalSS] = kmeans(matR_SU_all_moviemask_valid, 1); % voxel
Clustering_moviemask_valid.totalSS_Vox = totalSS;

% 2017/01/19 Record the results of each replicate (number of iterations,
% sum of the distance)
% diary(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new_masked_ReplicatesRecords', nameSubjNeural, nameSubjBOLD)))
for iK = 1:length(setK)
    
    K = setK(iK);                
    
    SU_indCluster_moviemask = NaN(size(matR_SU_all_moviemask_valid, 2), numRepeat);
    SU_sumD_moviemask = NaN(K, numRepeat);
    
    Vox_indCluster_moviemask = NaN(sum(locValidVox), numRepeat);
    Vox_sumD_moviemask = NaN(K, numRepeat);

    for iRep = 1:numRepeat
        
        fprintf(1, ':: K = %d; Movie-driven mask ::\n', K);
        
        % Cluster single units based on whole brain correlation
        [IDX_SU, C_SU, SUMD_SU] = kmeans(matR_SU_all_moviemask_valid', K,...
            'Replicates', numReplicates, 'Options', opts);
        
        SU_indCluster_moviemask(:, iRep) = IDX_SU;
        SU_sumD_moviemask(:, iRep) = SUMD_SU;
        
        % Cluster voxels based on 50 singel unit correlation
        [IDX_Vox, C_Vox, SUMD_Vox] = kmeans(matR_SU_all_moviemask_valid, K,...
            'Replicates', numReplicates, 'Options', opts);
        
        Vox_indCluster_moviemask(:, iRep) = IDX_Vox;
        Vox_sumD_moviemask(:, iRep) = SUMD_Vox;
        
    end

    Clustering_moviemask_valid.resultKMeans(iK).SU_indCluster = SU_indCluster_moviemask;
    Clustering_moviemask_valid.resultKMeans(iK).SU_sumD = SU_sumD_moviemask;
    Clustering_moviemask_valid.resultKMeans(iK).Vox_indCluster = Vox_indCluster_moviemask;
    Clustering_moviemask_valid.resultKMeans(iK).Vox_sumD = Vox_sumD_moviemask;
    
    if flagSave
        save(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked_probability_critCorr1.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)), ... %nameSubjNeural, nameSubjBOLD)),...
            'Clustering*', 'paramClustering*');
        fprintf(1, ':: K = %d; Movie-driven mask :: Results saved \n', K);
    end
end

%% Perform gap statistics to get the optimal number of clusters
% Gap statistics to get the optimal number
%fprintf(1, 'Performing gap statistics on SU clustering... \n');
%eva_clusteringSU = evalclusters(matR_SU_all_moviemask_valid', 'kmeans', 'gap', 'KList', 
%paramClustering_global.setK, 'Distance', 'sqEuclidean');
%save(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked_probability_critCorr1_eval.mat', 
%cell2mat(setNameSubjNeural), nameSubjBOLD)), ... %nameSubjNeural, nameSubjBOLD)),...
%        'eva_clusteringSU');
%fprintf(1, '...done! Evaluation results saved. \n')
%
%fprintf(1, 'Performing gap statistics on Voxel clustering... \n');
%eva_clusteringVox = evalclusters(matR_SU_all_moviemask_valid, 'kmeans', 'gap', 'KList', 
%paramClustering_global.setK, 'Distance', 'sqEuclidean'); 
%fprintf(1, '...done! Evaluation results saved. \n')
%save(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked_probability_critCorr1_eval.mat', 
%cell2mat(setNameSubjNeural), nameSubjBOLD)), ... %nameSubjNeural, nameSubjBOLD)),...
%        'eva_clusteringSU', 'eva_clusteringVox');

if flagParallel
    delete(pool);
end


%%%   
%% plot Clustering results
%figure;
%set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
%for iK=1:length(setK)
%    plot(iK+1, mean(Clustering_moviemask.resultKMeans(iK).Vox_sumD), 'ko', 'MarkerSize', 10, 'LineWidth', 2)
%    hold on
%end
%xlim([1 16])
%set(gca, 'LineWidth', 2, 'FontSize', 12)
%box off
%title('Clustering of voxels')
%xlabel('Number of cluster')
%ylabel('Within-cluster distance')
%% print(gcf, fullfile(dirFig, 'kmeans_distanceElbowPlot_Vox'), '-depsc');


% %% Gap statistics
% tempRand = -1 + (1+1).*rand(size(matR_SU)); % uniformly distributed random samples in interval [-1, 1]
% for iK = 1:length(setK)
%     
%     K = setK(iK);
%     % Cluster single units based on whole brain correlation
%     [IDX_SU_Rand, C, SUMD_SU_Rand] = kmeans(tempRand', K, 'Replicates', numReplicates, 'Options', opts); 
% %     % Cluster voxels based on 50 singel unit correlation
% %     [IDX_Vox, C, SUMD_Vox] = kmeans(matR_SU, K, 'Replicates', numReplicates, 'Options', opts);
%     
%     
%     Clustering_Rand.resultKMeans(iK).SU_indCluster = IDX_SU_Rand;
%     Clustering_Rand.resultKMeans(iK).SU_sumD = SUMD_SU_Rand;
% %     Clustering.resultKMeans(iK).Vox_indCluster = IDX_Vox;
% %     Clustering.resultKMeans(iK).Vox_sumD = SUMD_Vox;
%     
% end



