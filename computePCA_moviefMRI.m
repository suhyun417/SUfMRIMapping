% computePCA_moviefMRI.m
%
% 2019/1/8 SHP
%         - Modify the part where the residual is computed: before residual was only computed for "masked" maps but now it is done for unmasked maps
% 2018/10 SHP
%         - Perform PCA on fMRI voxel time series for each movie 
%         - Save the PCs and residuals
    

clear all;

%% PCA on fMRI data
% load the masks (movie-driven & brain-only mask)
nameSubjBOLD = 'Art'; %'Ava'; %'Art'; 
load(sprintf('/procdata/parksh/_macaque/%s/%s_MaskArrays.mat', nameSubjBOLD, nameSubjBOLD), 'movieDrivenAmp', 'brainMask_BlockAna3D');
[nx, ny, nz] = size(movieDrivenAmp.mask_amp1);
nVox = nx*ny*nz;
moviemask_vec = reshape(movieDrivenAmp.mask_amp1, nVox, 1); % change the 3D mask to 1D

load(sprintf('/procdata/parksh/_macaque/%s/%s_movieTS_fMRI_indMov.mat', nameSubjBOLD, nameSubjBOLD))

indMovieBOLD = [1 2 3];

for iM = 1:3;
    fmritc=[];
    curvoltc = voltcIndMov{iM};
    avgvoltc = repmat(nanmean(curvoltc,4),[1 1 1 size(curvoltc,4)]);
    if ~isempty(find(avgvoltc==0, 1))
        avgvoltc(avgvoltc==0) = realmin; % get rid of zeros because it causes NaNs in percent signals
    end
    pcvoltc = ((curvoltc - avgvoltc)./avgvoltc)*100;
    fmritc = pcvoltc(:,:,:,8:125); %pcvoltc;
    
    [nx, ny, nz, nt] = size(fmritc);
    nVox = nx*ny*nz;
    
    matTS = reshape(fmritc, nVox, nt);
    
    %% PCA on the whole brain
    [coeff, score, latent, tsquared, explained] = pca(zscore(matTS), 'Economy', 'off', 'Centered', 'on');
    [residuals] = pcares(zscore(matTS), 1);
    %%% x = zscore(matTS);
    %%% ndim = 1;
    %%% reconstructed = repmat(mean(x,1),n,1) + score(:,1:ndim)*coeff(:,1:ndim)';
    %%% residuals = x - reconstructed;
    
    resultsPCA(iM).coeff = coeff(:,1:10);
    resultsPCA(iM).score = score(:,1:10);
    resultsPCA(iM).explained = explained;
    paramPCA.option = {'Economy', 'off', 'Centered', 'on'};
    paramPCA.flag_zscore = 1;
    
    resultsPCAres(iM).residuals = residuals;
    paramPCAres.ndim = 1;
    paramPCAres.flag_zscore = 1;
    
    
    %% PCA on movie-masked one
    matTS_moviemask = matTS(moviemask_vec,:); %matR_SU(moviemask_vec,:); % 15495 voxels
    
    [coeff, score, latent, tsquared, explained] = pca(zscore(matTS_moviemask), 'Economy', 'off', 'Centered', 'on');
    [residuals] = pcares(zscore(matTS_moviemask), 1);
    %%% x = zscore(matTS);
    %%% ndim = 1;
    %%% reconstructed = repmat(mean(x,1),n,1) + score(:,1:ndim)*coeff(:,1:ndim)';
    %%% residuals = x - reconstructed;
    
    resultsPCA_moviemask(iM).coeff = coeff(:,1:10);
    resultsPCA_moviemask(iM).score = score(:,1:10);
    resultsPCA_moviemask(iM).explained = explained;
    paramPCA_moviemask.option = {'Economy', 'off', 'Centered', 'on'};
    paramPCA_moviemask.flag_zscore = 1;
    
    resultsPCAres_moviemask(iM).residuals = residuals;
    paramPCAres_moviemask.ndim = 1;
    paramPCAres_moviemask.flag_zscore = 1;
    
    
end

% concatenate the first principal component across movies
catCoeff = cat(1, resultsPCA.coeff);
pc1_fMRI = [];
for iM = 1:3
    pc1_fMRI = cat(1, pc1_fMRI, NaN(7,1), catCoeff(118*(iM-1)+1:118*iM, 1));
end

catCoeff_moviemask = cat(1, resultsPCA_moviemask.coeff);
pc1_fMRI_moviemask = [];
for iM = 1:3
    pc1_fMRI_moviemask = cat(1, pc1_fMRI_moviemask, NaN(7,1), catCoeff_moviemask(118*(iM-1)+1:118*iM, 1));
end
% pc1_fMRI = catCoeff(:,1);

save(sprintf('/procdata/parksh/_macaque/%s/%s_movieTS_fMRI_Movie123_PCA.mat', nameSubjBOLD, nameSubjBOLD), 'paramPCA*', 'resultsPCA*', 'pc1_fMRI*')

clear fmritc pcvoltc avgvoltc curvoltc voltcIndMov

% some evaluation on PCs
catVar = cat(2, resultsPCA_moviemask.explained);
cumsum(catVar(1:10, :))


%% save the surface map for score of PCs
% to visualize scores for each voxels related to each principal component
load('/procdata/parksh/Art/Art_movieTS_fMRI_Movie123_PCA.mat', 'resultsPCA');

nameSubjBOLD = 'Art'; 
load(sprintf('/procdata/parksh/%s/%s_MaskArrays.mat', nameSubjBOLD, nameSubjBOLD), 'movieDrivenAmp', 'brainMask_BlockAna3D');

for iPC = 2:3 %1;
for iMovie = [1 2 3];
% for iPC = 5; %1:nPC
    
    [nx ny nz] = size(movieDrivenAmp.mask_amp1);
    nVox = nx*ny*nz;
    moviemask_vec = reshape(movieDrivenAmp.mask_amp1, nVox, 1);
    
%     tempPC_valid = zeros(size(paramClustering_global.locValidVox));
%     tempPC_valid(paramClustering_global.locValidVox) = coeff_e(:, iPC); %Clustering_moviemask_valid.resultKMeans(sortTargetK-1).Vox_indCluster;

    PCMap = NaN(size(moviemask_vec));
    PCMap(moviemask_vec) = resultsPCA_moviemask(iMovie).score(:,iPC);
    
%     voxelClusterMap = zeros(size(moviemask_vec));
%     voxelClusterMap(moviemask_vec) = Clustering_moviemask.resultKMeans(sortTargetK-1).Vox_indCluster;
    
    % voxelClusterMap_vol = reshape(voxelClusterMap./sortTargetK, [nx, ny, nz]); % interim work-around for labeling of different clusters in AFNI
    PCMap_vol = reshape(PCMap, [nx, ny, nz]); % each voxel now has cluster identity (integers)
    
    
    dirDataBOLD = sprintf('/procdata/parksh/%s', nameSubjBOLD);
    pname = [dirDataBOLD, '/tempSURF/'];
    dirSPEC = [dirDataBOLD, '/Anatomy/_suma/'];
    
    cd /projects/parksh/_toolbox/BlockAna/ %/projects/parksh/NeuralBOLD/analysis/BlockAna/
    blockana;
    S_neuralRegressor('Tor', nameSubjBOLD); % random neural subject to get the global field going
    global STDPATH DSP DATA GH
    
    % convert the map to the surface
    % for iMask = 1:2
    % %     iMask = 2;
    %
    %     for iK = 1:sortTargetK
    
%     cd /projects/parksh/_toolbox/BlockAna/ %/projects/parksh/NeuralBOLD/analysis/BlockAna/
    
    %         switch iMask
    %             case 1 % unmasked (everything within brain)
    %                 fileHead = sprintf('new_%s', nameSubjNeural);
    %                 DSP.proc.fncvol_3d = mapR_Cluster(:,:,:,iK).*brainMask_BlockAna3D;
    %             case 2 % masked
%     fileHead = sprintf('%s', nameSubjBOLD); %sprintf('new_maskedSignificant_%s', nameSubjNeural);% sprintf('new_masked_%s', nameSubjNeural);
    DSP.proc.fncvol_3d = PCMap_vol; %mapR_Cluster(:,:,:,iK).*movieDrivenAmp.mask_amp1; %reshape(mapR, [nx, ny, nz]).*movieDrivenAmp.mask_amp1;
    %         end
    
    %         fname = sprintf('%s_Cluster%d_%dMeans_Art_AvgCorrMapMovie123_noFiltering+orig.BRIK', fileHead, iK, sortTargetK);%
    fname = sprintf('%s_fMRI_PC%dScore_Movie%d_noFiltering+orig.BRIK', nameSubjBOLD, iPC, iMovie);%
    
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
end

