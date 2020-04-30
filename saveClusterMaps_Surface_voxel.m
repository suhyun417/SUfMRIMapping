% saveClusterMaps_Surface.m
%
% 1) load fMRI correlation maps computed based on Toroid's cells and Artemis
% or Avalanche's fMRI data
% 2) load clustering results based on various values (e.g. whole brain corr maps, corr between cells and movie
% regressors, cell spike density function etc.)
% 3) for each cluster, average the corr maps
% 4) save it as AFNI brik/head format so that it can be mapped onto the surface anatomy

clear all;


nameSubjNeural = 'Tor'; %'Spi'; % 'Tor'; % 'Sig'; %'Tor';
nameSubjBOLD = 'Art'; %'Ava'; %'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';

%
addpath('/library/matlab_utils/')

% Load files
dirDataHome = '/procdata/parksh/';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% 1) fMRI correlation maps
load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'matR_SU', 'paramCorr')

% 2) Clustering results
% load(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new_masked_significant.mat', nameSubjNeural, nameSubjBOLD)))  % cluste
% % load(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD))) %
load(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new_masked.mat', nameSubjNeural, nameSubjBOLD)))  % clustering based on the maps with movie-driven mask applied
% % load(fullfile(dirDataNeural, 'Clustering_TorArtAvaMovie123.mat')) % based on two sets of whole brain maps from two different monkeys
% % Clustering = ClusteringAll;

% 3) fMRI movie-driven activity mask
load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp', 'brainMask_BlockAna3D');

% % 4) Single-unit time courses
% load(fullfile(dirDataNeural, sprintf('%s_movieTS_SU_indMov.mat', nameSubjNeural)))

% filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
% fileNameNeural_BLP = [nameSubjNeural, '_movieTS_BLPLFP_indMov.mat'];
% filenameBOLD = [nameSubjBOLD, '_movieTS_fMRI_indMov.mat'];
%
% fprintf(1, '\nLoading single unit data of %s: %s ....', nameSubjNeural, filenameNeural)
% load(fullfile(dirDataNeural, filenameNeural))
% load(fullfile(dirDataNeural, fileNameNeural_BLP))
% fprintf(1, '\nLoading fMRI data of %s: %s ....\n', nameSubjBOLD, filenameBOLD)
% load(fullfile(dirDataBOLD, filenameBOLD))

% % Get movie IDs common in two dataset
% commonSetMovie = intersect(paramBOLD.unimov, paramSDF.setMovIDs);
% dataBOLD.mvoltc = voltcIndMov(commonSetMovie);
% dataBOLD.unimov = commonSetMovie;
%
% setMovie = [1 2 3];

% 2017/01/20 add Spice's case where valid channels were selected at
% the clustering stage
switch lower(nameSubjNeural)
    case 'spi'
        matR_SU_org = matR_SU;
        clear matR_SU
        matR_SU = matR_SU_org(:, paramClustering.validChanIndex);
end

%% Get indices of voxel clusters
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% for iK = 1:length(paramClustering.setK)
%     plot(paramClustering.setK(iK), Clustering_moviemask.resultKMeans(iK).Vox_sumD, 'ko', 'LineWidth', 2, 'MarkerSize', 10)
%     hold on
% end
% xlim([1 16])
% set(gca, 'LineWidth', 2, 'FontSize', 12)
% box off
% title('Clustering of voxels')
% xlabel('Number of cluster')
% ylabel('Within-cluster distance')

% matIndClust_Vox = cat(2, Clustering_moviemask.resultKMeans.Vox_indCluster);
for iK = 1:length(paramClustering.setK)
    sortTargetK = paramClustering.setK(iK);
    % [sortedClust, indSortVoxel]=sort(matIndClust_Vox(:,sortTargetK-1));
    
    [nx ny nz] = size(movieDrivenAmp.mask_amp1);
    nVox = nx*ny*nz;
    moviemask_vec = reshape(movieDrivenAmp.mask_amp1, nVox, 1);
    
    voxelClusterMap = zeros(size(moviemask_vec));
    voxelClusterMap(moviemask_vec) = Clustering_moviemask.resultKMeans(sortTargetK-1).Vox_indCluster;
    
    % voxelClusterMap_vol = reshape(voxelClusterMap./sortTargetK, [nx, ny, nz]); % interim work-around for labeling of different clusters in AFNI
    voxelClusterMap_vol = reshape(voxelClusterMap, [nx, ny, nz]); % each voxel now has cluster identity (integers)
    
    
    
    pname = [dirDataBOLD, '/tempSURF/'];
    dirSPEC = [dirDataBOLD, '/Anatomy/_suma/'];
    
    cd /projects/parksh/NeuralBOLD/analysis/BlockAna/
    blockana;
    S_neuralRegressor(nameSubjNeural, nameSubjBOLD);
    global STDPATH DSP DATA GH
    
    % convert the map to the surface
    % for iMask = 1:2
    % %     iMask = 2;
    %
    %     for iK = 1:sortTargetK
    
    cd /projects/parksh/NeuralBOLD/analysis/BlockAna/
    
    %         switch iMask
    %             case 1 % unmasked (everything within brain)
    %                 fileHead = sprintf('new_%s', nameSubjNeural);
    %                 DSP.proc.fncvol_3d = mapR_Cluster(:,:,:,iK).*brainMask_BlockAna3D;
    %             case 2 % masked
    fileHead = sprintf('new_masked_%s%s', nameSubjNeural, nameSubjBOLD); %sprintf('new_maskedSignificant_%s', nameSubjNeural);% sprintf('new_masked_%s', nameSubjNeural);
    DSP.proc.fncvol_3d = voxelClusterMap_vol; %mapR_Cluster(:,:,:,iK).*movieDrivenAmp.mask_amp1; %reshape(mapR, [nx, ny, nz]).*movieDrivenAmp.mask_amp1;
    %         end
    
    %         fname = sprintf('%s_Cluster%d_%dMeans_Art_AvgCorrMapMovie123_noFiltering+orig.BRIK', fileHead, iK, sortTargetK);%
    fname = sprintf('%s_VoxClusters_%d_Movie123_noFiltering+orig.BRIK', fileHead,sortTargetK);%
    
    vol = single(DSP.proc.fncvol_3d);
    
    %     fname = sprintf('%s_%s_Movie123_noFiltering_%s+orig.BRIK', fileHead, nameSubjNeural, fileTail);%
    %
    %     vol = single(DSP.proc.fncvol_3d);  %single(mapR_Cluster(:,:,:,iK)); % vol = single(DSP.proc.fncvol_3d);
    % %     cellID = validChanID(iUnit,:);
    % %     fprintf(1, 'Unit # %d, Cell ID: %s, %s \n', iUnit, cellID);
    
    
    % Do dumpFunctionalBrik
    Info = DATA.Info;
    Info.DATASET_RANK = DATA.Info.DATASET_RANK;
    Info.DATASET_DIMENSIONS = DATA.Info.DATASET_DIMENSIONS;
    Info.DELTA = DATA.Info.DELTA;
    Info.BYTEORDER_STRING = DATA.Info.BYTEORDER_STRING;
    Info.TYPESTRING = DATA.Info.TYPESTRING;
    Info.SCENE_DATA = DATA.Info.SCENE_DATA;
    Info.ORIGIN = DATA.Info.ORIGIN;
    Info.TAXIS_FLOATS = DATA.Info.TAXIS_FLOATS;
    Info.TAXIS_OFFSETS = DATA.Info.TAXIS_OFFSETS;
    Info.IDCODE_STRING = DATA.Info.IDCODE_STRING;
    Info.IDCODE_DATE = DATA.Info.IDCODE_DATE;
    Info.HISTORY_NOTE = DATA.Info.HISTORY_NOTE;
    
    Info.DATASET_DIMENSIONS = [size(vol) 0 0];
    Info.DATASET_RANK(2) = 1;
    Info.BRICK_TYPES = 3;  %3 = float
    Info.BRICK_FLOAT_FACS = 0;
    Info.BRICK_STATS = [min(vol(:)) max(vol(:))];
    Info.BRICK_LABS = DATA.Info.BRICK_LABS(1:2);
    Info.BRICK_KEYWORDS = DATA.Info.BRICK_KEYWORDS;
    Info.TAXIS_NUMS(1) = 1;
    Info.TypeName = 'float';
    Info.TypeBytes = 4;
    Info.ByteOrder = 'ieee-le';
    Info.Orientation = DATA.Info.Orientation;
    Info.FileFormat = DATA.Info.FileFormat;
    Info.RootName = DATA.Info.RootName;
    Info.Extension_1D = DATA.Info.Extension_1D;
    
    code0 = getAFNIOrientationCode('LR');
    code1 = getAFNIOrientationCode('AP');
    code2 = getAFNIOrientationCode('DV');
    
    Info.ORIENT_SPECIFIC = [code0 code1 code2];
    
    doti = strfind(fname,'+');
    prefix = fname(1:doti-1);
    Opt.Prefix = prefix;
    Opt.NoCheck = 0;
    
    fprintf('Writing %s as BRIK and HEAD files to %s\n',prefix,pname);
    curdir = cd(pname);
    
    %
    % Compensating for an AFNI inconsistency
    %
    if isfield(Info,'WARP_TYPE'),
        Info = rmfield(Info,'WARP_TYPE');
    end
    
    WriteBrik(vol,Info,Opt);  % needs to have a modified version to overwrite
    % an existing file with "w+"
    
    %%%% ADDED BY BER 4/5/13 for alignment purposes.
    % copy data to new file so original dimensions are saved.
    copyDataCMD = ['3dcopy  ' prefix ' ' prefix '_refit'];
    disp(copyDataCMD);
    system(copyDataCMD);
    
    % refit the x y z dimensioned based on strange AFNI conventions
    refitDataCMD = ['3drefit -xdel 1.5 -ydel 1.5 -zdel 1.5 ' prefix '_refit+orig'];
    disp(refitDataCMD);
    system(refitDataCMD);
    
    monk=DSP.afni_prefix(1:3);
    % now that data is refit, align the center with the surface anatomy center
    alignDataCMD = ['@Align_Centers -base ' STDPATH.afni_uw '/Anatomy/' monk '_t1_refit_shft+orig. -dset ' prefix '_refit+orig.'];
    disp(alignDataCMD);
    system(alignDataCMD);
    
    % finally shift the origin slightly because of strange ANFI stuff
    shiftDataCMD = ['3drefit -dyorigin -0.5 -dzorigin 0.5 ' prefix '_refit_shft+orig.'];
    disp(shiftDataCMD);
    system(shiftDataCMD);
    %%%%
    
    % Copy relevant files to where the surface files are
    [s,mess,messID]=movefile([pname, '*_refit_shft+orig*'],  dirSPEC);
    if ~s
        fprintf(1,' %s \n', mess);
    end
end
% end
% end


