% saveClusterMaps_Surface_multipleFP_voxelbased.m
%
% 2020/08/15 SHP
% Modified from "saveClusterMaps_Surface_multipleFP_ROIbased.m" to load voxel-based clustering
% results and save the maps for each cluster
% - Load clustering results based on whole brain corr maps
% - For each cluster, average the corr maps
% - Save it as AFNI brik/head format so that it can be mapped onto the surface anatomy

clear all;

nameSubjBOLD = 'Art'; %'Ava';

%
% addpath('/library/matlab_utils/')

% Load files
dirDataHome = '/procdata/parksh/_macaque/';
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% 1) fMRI correlation maps
% load(sprintf('/procdata/parksh/_macaque/CorrMap_SU_AllCells%s_corticalFPMerged.mat', nameSubjBOLD), 'info*', 'corrMap_merged_FP'); %, 'info*', 'corrMap_Area', 'corrMap_merged');

% 2) Clustering results (voxel-based clustering)
load(sprintf('/procdata/parksh/_macaque/%s/Clustering_CorrMap_4FPs_Movie123_probability.mat', nameSubjBOLD), 'Clustering_brainmask')

% 3) fMRI movie-driven activity mask
load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp', 'brainMask_BlockAna3D');


%% 1. Average fMRI maps for each cluster
setK = 2:20; %paramClustering_global.setK; %Clustering.setK;

matWSS=[];
matExpVar=[];
for iK = 1:length(setK)
    curK = setK(iK);
    matWSS(:,iK) = sum(Clustering_brainmask.resultKMeans(iK).SU_sumD); %sum(Clustering.resultKMeans(iK).SU_sumD);
end

totalSS = Clustering_brainmask.totalSS_SU;
propExplained = (totalSS-matWSS)./totalSS; %matExpVar./totalSS;

for sortTargetK = 10 %9:10 %5; %6:8 %4; %:5 %sortTargetK = 4; %8; %6; %7; %6; %7;
    
    
    curK = sortTargetK; %9; %6; %7;
    locMode = find(propExplained(:,curK-1)==mode(propExplained(:,curK-1)));
    locMin = find(propExplained(:,curK-1)==min(propExplained(:,curK-1)));
    [sortedClust, indSortChan] = sort(Clustering_brainmask.resultKMeans(curK-1).SU_indCluster(:, locMode(1)));
    
    % matIndClust_SU = cat(2, Clustering_moviemask.resultKMeans.SU_indCluster); % cat(2, Clustering.resultKMeans.SU_indCluster);
    % matIndClust_SU =
    % [sortedClust, indSortChan]=sort(matIndClust_SU(:,sortTargetK-1));
    
    % average maps for each cluster
    nx = 40; ny = 64; nz = 32;
    nVox = nx*ny*nz;
    brainmask_vec = reshape(movieDrivenAmp.map_sm_brain>0, nVox, 1); % change the 3D mask to 1D
    mapR_Cluster=[];
    for iK = 1:sortTargetK
        R_avgCluster = zeros(size(brainmask_vec));
        R_avgCluster(brainmask_vec) = mean(Clustering_brainmask.matR(:,indSortChan(sortedClust==iK)), 2);
        tempClustMapR = reshape(R_avgCluster, [nx, ny, nz]);
        mapR_Cluster = cat(4, mapR_Cluster, tempClustMapR);
    end
    
    % % compute residual
    % GRAND average map
    % mapR_grandAvg = reshape(mean(matR_SU,2), [nx, ny, nz]);
    % R_grandAvg = mean(matR_SU,2);
    % R_avgCluster_residual = R_avgCluster - repmat(R_grandAvg, 1, sortTargetK);
    % mapR_Cluster_residual = [];
    % for iK = 1:sortTargetK
    %     tempClustMapR_residual = reshape(R_avgCluster_residual(:,iK), [nx, ny, nz]);
    %     mapR_Cluster_residual = cat(4, mapR_Cluster_residual, tempClustMapR_residual);
    % end
    
    
    pname = [dirDataBOLD, '/tempSURF/'];
    dirSPEC = [dirDataBOLD, '/Anatomy/_suma/'];
    
    cd /projects/parksh/_toolbox/BlockAna/
    blockana;
    S_neuralRegressor('Tor', nameSubjBOLD);
    global STDPATH DSP DATA GH
    
    % convert the map to the surface
    
    for iK = 1:sortTargetK
        
        cd /projects/parksh/_toolbox/BlockAna/

        fname = sprintf('Clustering4FP_ArtBrainVoxelBasedMovie123_%dMeans_Cluster%d_AvgCorrMap+orig.BRIK', sortTargetK, iK);%
        
        DSP.proc.fncvol_3d = mapR_Cluster(:,:,:,iK); %.*brainMask_BlockAna3D;        
        vol = single(DSP.proc.fncvol_3d);
               
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
end


