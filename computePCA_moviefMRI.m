% computePCA_moviefMRI.m
%
% 2020/7/29 SHP
%         - Add PCA on brain-only voxels
%         - Add computing r and r-squared for each voxel with 1st PC
% 2020/7/18 SHP
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
catTS = []; catTS_nan = [];
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
    
    catTS = cat(2, catTS, matTS);
    catTS_nan = cat(2, catTS, NaN(size(matTS, 1), 7), matTS);
    
%     %% PCA on the whole slice package (including voxels outside of the brain)
%     [coeff, score, latent, tsquared, explained] = pca(zscore(matTS), 'Economy', 'off', 'Centered', 'on');
%     [residuals] = pcares(zscore(matTS), 1);
%     %%% x = zscore(matTS);
%     %%% ndim = 1;
%     %%% reconstructed = repmat(mean(x,1),n,1) + score(:,1:ndim)*coeff(:,1:ndim)';
%     %%% residuals = x - reconstructed;
%     
%     resultsPCA(iM).coeff = coeff(:,1:10);
%     resultsPCA(iM).score = score(:,1:10);
%     resultsPCA(iM).explained = explained;
%     paramPCA.option = {'Economy', 'off', 'Centered', 'on'};
%     paramPCA.flag_zscore = 1;
%     
%     resultsPCAres(iM).residuals = residuals;
%     paramPCAres.ndim = 1;
%     paramPCAres.flag_zscore = 1;
%     
%     %% PCA on brain-masked one
%     brainmask_vec = reshape(movieDrivenAmp.map_sm_brain>0, nVox, 1); % change the 3D mask to 1D
%     matTS_brainmask = matTS(brainmask_vec, :); %matR_SU(brainmask_vec,:); % 27113 voxels
%     
%     [coeff, score, latent, tsquared, explained] = pca(zscore(matTS_brainmask), 'Economy', 'off', 'Centered', 'on');
%     [residuals] = pcares(zscore(matTS_brainmask), 1);
%     %%% x = zscore(matTS);
%     %%% ndim = 1;
%     %%% reconstructed = repmat(mean(x,1),n,1) + score(:,1:ndim)*coeff(:,1:ndim)';
%     %%% residuals = x - reconstructed;
%     
%     resultsPCA_brainmask(iM).coeff = coeff(:,1:10);
%     resultsPCA_brainmask(iM).score = score(:,1:10);
%     resultsPCA_brainmask(iM).explained = explained;
%     resultsPCA_brainmask(iM).matTS = matTS_brainmask;
%     paramPCA_brainmask.brainmask_1d = brainmask_vec;
%     paramPCA_brainmask.option = {'Economy', 'off', 'Centered', 'on'};
%     paramPCA_brainmask.flag_zscore = 1;
%     
%     resultsPCAres_brainmask(iM).residuals = residuals;
%     paramPCAres_brainmask.ndim = 1;
%     paramPCAres_brainmask.flag_zscore = 1;
    
    
%     %% PCA on movie-masked one
%     matTS_moviemask = matTS(moviemask_vec,:); %matR_SU(moviemask_vec,:); % 15495 voxels
%     
%     [coeff, score, latent, tsquared, explained] = pca(zscore(matTS_moviemask), 'Economy', 'off', 'Centered', 'on');
%     [residuals] = pcares(zscore(matTS_moviemask), 1);
%     %%% x = zscore(matTS);
%     %%% ndim = 1;
%     %%% reconstructed = repmat(mean(x,1),n,1) + score(:,1:ndim)*coeff(:,1:ndim)';
%     %%% residuals = x - reconstructed;
%     
%     resultsPCA_moviemask(iM).coeff = coeff(:,1:10);
%     resultsPCA_moviemask(iM).score = score(:,1:10);
%     resultsPCA_moviemask(iM).explained = explained;
%     paramPCA_moviemask.option = {'Economy', 'off', 'Centered', 'on'};
%     paramPCA_moviemask.flag_zscore = 1;
%     
%     resultsPCAres_moviemask(iM).residuals = residuals;
%     paramPCAres_moviemask.ndim = 1;
%     paramPCAres_moviemask.flag_zscore = 1;
    
    
end


%% concatenate the first principal component across movies
catCoeff = cat(1, resultsPCA.coeff);
pc1_fMRI = [];
for iM = 1:3
    pc1_fMRI = cat(1, pc1_fMRI, NaN(7,1), catCoeff(118*(iM-1)+1:118*iM, 1));
end

catCoeff_brainmask = cat(1, resultsPCA_brainmask.coeff);
pc1_fMRI_brainmask = [];
for iM = 1:3
    pc1_fMRI_brainmask = cat(1, pc1_fMRI_brainmask, NaN(7,1), catCoeff_brainmask(118*(iM-1)+1:118*iM, 1));
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
catVar = cat(2, resultsPCA_brainmask.explained);
cumsum(catVar(1:10, :))


%% PCA on concatenated TS
brainmask_vec = reshape(movieDrivenAmp.map_sm_brain>0, nVox, 1); % change the 3D mask to 1D
catTS_brainmask = catTS(brainmask_vec, :); %matR_SU(brainmask_vec,:); % 27113 voxels
% catTS_nan_brainmask = catTS_nan(brainmask_vec, :);

[coeff, score, latent, tsquared, explained] = pca(zscore(catTS_brainmask), 'Economy', 'off', 'Centered', 'on');
[residuals] = pcares(zscore(catTS_brainmask), 1);

resultsPCA_concat_brainmask.coeff = coeff(:,1:10);
resultsPCA_concat_brainmask.score = score(:,1:10);
resultsPCA_concat_brainmask.explained = explained;
paramPCA_concat_brainmask.option = {'Economy', 'off', 'Centered', 'on'};
paramPCA_concat_brainmask.flag_zscore = 1;
paramPCA_concat_brainmask.setMovie = [1 2 3];

resultsPCAres_concat_brainmask.residuals = residuals;
resultsPCAres_concat_brainmask.ndim = 1;
resultsPCAres_concat_brainmask.flag_zscore = 1;


%% compute r-squared for concated PC
for iPC = 1:3
    catR=[];
    catR = corr(catTS', resultsPCA_concat_brainmask.coeff(:,iPC), 'type', 'Spearman');
    resultsCorrPC_concat_brainmask.indPC(iPC).matR = catR;
    resultsCorrPC_concat_brainmask.indPC(iPC).matRsq = catR.^2; % r squared
end


%% Compute R2 for each voxel with the 1st PC
clear all;

nameSubjBOLD = 'Art'; %'Art'; %'Ava'; %'Art';
load(sprintf('/procdata/parksh/_macaque/%s/%s_movieTS_fMRI_indMov.mat', nameSubjBOLD, nameSubjBOLD))
load(sprintf('/procdata/parksh/_macaque/%s/%s_movieTS_fMRI_Movie123_PCA.mat', nameSubjBOLD, nameSubjBOLD), 'resultsPCA*')

indMovieBOLD = [1 2 3];
% first, for each movie separately
catTS = [];
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
    
    % correlation with 1st PC from the whole slice package
    matR = corr(matTS', resultsPCA(iM).coeff(:,1), 'type', 'Spearman');
    resultsCorrPC_wholebrain.indMovie(iM).matR = matR;
    resultsCorrPC_wholebrain.indMovie(iM).matRsq = matR.^2; % r squared
    
    % correlation with 1st PC from the brain
    matR = corr(matTS', resultsPCA_brainmask(iM).coeff(:,1), 'type', 'Spearman');
    resultsCorrPC_brainmask.indMovie(iM).matR = matR;
    resultsCorrPC_brainmask.indMovie(iM).matRsq = matR.^2; % r squared
    
    % correlation with 1st PC from the movie-masked brain
    matR = corr(matTS', resultsPCA_moviemask(iM).coeff(:,1), 'type', 'Spearman');
    resultsCorrPC_moviemask.indMovie(iM).matR = matR;
    resultsCorrPC_moviemask.indMovie(iM).matRsq = matR.^2; % r squared
    
    catTS = cat(2, catTS, matTS);
end

% using concatenated TS
catCoeff = []; catR = [];
catCoeff = cat(1, resultsPCA.coeff);
catR = corr(catTS', catCoeff(:,1), 'type', 'Spearman');
resultsCorrPC_wholebrain.matR_catTS = catR;
resultsCorrPC_wholebrain.matRsq_catTS = catR.^2; % r squared

catCoeff = []; catR = [];
catCoeff = cat(1, resultsPCA_brainmask.coeff);
catR = corr(catTS', catCoeff(:,1), 'type', 'Spearman');
resultsCorrPC_brainmask.matR_catTS = catR;
resultsCorrPC_brainmask.matRsq_catTS = catR.^2; % r squared

catCoeff = []; catR = [];
catCoeff = cat(1, resultsPCA_moviemask.coeff);
catR = corr(catTS', catCoeff(:,1), 'type', 'Spearman');
resultsCorrPC_moviemask.matR_catTS = catR;
resultsCorrPC_moviemask.matRsq_catTS = catR.^2; % r squared


save(sprintf('/procdata/parksh/_macaque/%s/%s_movieTS_fMRI_Movie123_PCA.mat', nameSubjBOLD, nameSubjBOLD), 'resultsCorr*', '-append')

%% save the surface map for some results of PCA analysis
% to visualize scores for each voxels related to each principal component
nameSubjBOLD = 'Art';
load(sprintf('/procdata/parksh/_macaque/%s/%s_movieTS_fMRI_Movie123_PCA.mat', nameSubjBOLD, nameSubjBOLD), 'paramPCA_brainmask', 'resultsPCA*', 'resultsCorr*')
% load(sprintf('/procdata/parksh/%s/%s_MaskArrays.mat', nameSubjBOLD, nameSubjBOLD), 'movieDrivenAmp', 'brainMask_BlockAna3D');

dirDataBOLD = sprintf('/procdata/parksh/_macaque/%s', nameSubjBOLD);
pname = [dirDataBOLD, '/tempSURF/'];
dirSPEC = [dirDataBOLD, '/Anatomy/_suma/'];

cd /projects/parksh/_toolbox/BlockAna/ %/projects/parksh/NeuralBOLD/analysis/BlockAna/
blockana;
S_neuralRegressor('Tor', nameSubjBOLD); % random neural subject to get the global field going
global STDPATH DSP DATA GH

% save the score maps
for iPC = 1:3 %1;
    for iMovie = [1 2 3];
        
        cd /projects/parksh/_toolbox/BlockAna/ %/projects/parksh/NeuralBOLD/analysis/BlockAna/
                     
        nx = 40; ny = 64; nz = 32;
        nVox = nx*ny*nz;
        
        PCMap = NaN(size(paramPCA_brainmask.brainmask_1d));
        PCMap(paramPCA_brainmask.brainmask_1d) = resultsPCA_brainmask(iMovie).score(:,iPC);
        PCMap_vol = reshape(PCMap, [nx, ny, nz]);
        
%         PCMap_vol = reshape(resultsPCA(iMovie).score(:,iPC), [nx, ny, nz]);
        
        % convert the map to the surface
        %     fileHead = sprintf('%s', nameSubjBOLD); %sprintf('new_maskedSignificant_%s', nameSubjNeural);% sprintf('new_masked_%s', nameSubjNeural);
        DSP.proc.fncvol_3d = PCMap_vol; %mapR_Cluster(:,:,:,iK).*movieDrivenAmp.mask_amp1; %reshape(mapR, [nx, ny, nz]).*movieDrivenAmp.mask_amp1;
        fname = sprintf('%s_fMRI_PC%dScore_Movie%d_brainmask+orig.BRIK', nameSubjBOLD, iPC, iMovie);%
        
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
        
        % Compensating for an AFNI inconsistency
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

%% Save the r squared map
clear all;
nameSubjBOLD = 'Art';
load(sprintf('/procdata/parksh/_macaque/%s/%s_movieTS_fMRI_Movie123_PCA.mat', nameSubjBOLD, nameSubjBOLD), 'paramPCA_brainmask', 'resultsPCA*', 'resultsCorr*')
% load(sprintf('/procdata/parksh/%s/%s_MaskArrays.mat', nameSubjBOLD, nameSubjBOLD), 'movieDrivenAmp', 'brainMask_BlockAna3D');

dirDataBOLD = sprintf('/procdata/parksh/_macaque/%s', nameSubjBOLD);
pname = [dirDataBOLD, '/tempSURF/'];
dirSPEC = [dirDataBOLD, '/Anatomy/_suma/'];

cd /projects/parksh/_toolbox/BlockAna/ %/projects/parksh/NeuralBOLD/analysis/BlockAna/
blockana;
S_neuralRegressor('Tor', nameSubjBOLD); % random neural subject to get the global field going
global STDPATH DSP DATA GH

for iType = 4 %1:3
       
    switch iType
        case 1 % no mask (whole slice package)
            filetail = 'wholeSlice';
            resultsS = resultsCorrPC_wholebrain;            
        case 2 % brain mask (voxels within brain)
            filetail = 'brainmask';
            resultsS = resultsCorrPC_brainmask;
        case 3 % movie mask (only movie-driven voxels)
            filetail = 'moviemask';
            resultsS = resultsCorrPC_moviemask;
        case 4 % brain mask, with concatenated time series
            filetail = 'concat_brainmask';
            resultsS = resultsCorrPC_concat_brainmask;
    end

    for iMovie = [0 1 2 3];
        
        cd /projects/parksh/_toolbox/BlockAna/ %/projects/parksh/NeuralBOLD/analysis/BlockAna/
                     
        nx = 40; ny = 64; nz = 32;
        nVox = nx*ny*nz;
        
        if iType < 4
            if iMovie > 0
                PCMap_vol = reshape(resultsS.indMovie(iMovie).matRsq, [nx, ny, nz]);
                fname = sprintf('%s_fMRI_PC1Rsquared_Movie%d_%s+orig.BRIK', nameSubjBOLD, iMovie, filetail);%
            else
                PCMap_vol = reshape(resultsS.matRsq_catTS, [nx, ny, nz]);
                fname = sprintf('%s_fMRI_PC1Rsquared_catMovie123_%s+orig.BRIK', nameSubjBOLD, filetail);%
            end
        else %  concatenated (15 min) time sereis iType == 4 
            if iMovie > 0
                iPC = iMovie;
                PCMap_vol = reshape(resultsS.indPC(iPC).matRsq, [nx, ny, nz]); 
                fname = sprintf('%s_fMRI_PC%d_catTSMovie123_Rsquared_%s+orig.BRIK', nameSubjBOLD, iPC, filetail);%
            else
                continue;
            end
        end            
        
        % convert the map to the surface
        DSP.proc.fncvol_3d = PCMap_vol;         
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
        
        % Compensating for an AFNI inconsistency
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


