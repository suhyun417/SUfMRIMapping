% saveClusterMaps_Surface.m
%
% 1) load fMRI correlation maps computed based on Toroid's cells and Artemis
% or Avalanche's fMRI data
% 2) load clustering results based on various values (e.g. whole brain corr maps, corr between cells and movie
% regressors, cell spike density function etc.)
% 3) for each cluster, average the corr maps
% 4) save it as AFNI brik/head format so that it can be mapped onto the surface anatomy

clear all;

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
addpath(fullfile(dirProjects, 'parksh/_toolbox/afni_matlab'))
addpath(fullfile(dirProjects, 'parksh/_toolbox/hslcolormap'))

% Set directories
setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi'};
nameSubjNeural = 'Spi'; % 'Tor';
nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataHome = fullfile(dirProcdata, 'parksh');
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% Directory for saving figures as graphic files
dirFig = fullfile(dirProjects, 'parksh/NeuralBOLD/_labNote/_figs/');

%% Plotting parameters
% cMap = [0 0 0; 230 159 0; 86 180 233; 0 158 115; 240 228 66; 0 114 178; 213 94 0; 204 121 167]./255;
cMap = [228	26	28;
    55	126	184;
    77	175	74;
    152	78	163;
    255	127	0;
    255 217 47; % dark yellow %255	255	255; %white %255	255	51; % yellow was too similar to another yellow in mat2
    166	86	40;
    247	129	191;
    ]./255;
marker = {'o', '^', 'square', 'diamond'}; % {'o','*', 'x', 's', 'd', '+', '^'};
% orderIndNewClust = [3 2 7 6 5 1 4];
% cMap_newclust = cMap(orderIndNewClust,:);

%% Load data
% 1) Clustering results
load(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked_probability_critCorr1.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)))
% load(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new_masked.mat', nameSubjNeural, nameSubjBOLD)));
% 2) Movie-driven mask
load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp');
% 3) Valid cells considering likelihood being clustered together
load(fullfile(dirDataNeural, 'Clustering_SU_allCells_validVoxels_critCorr1_7Means_prob.mat'))

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
    matR_SU_valid = matR_SU(:, validChanIndex_clustering);
    matR_SU_all = cat(2, matR_SU_all, matR_SU_valid);
    clear matR_SU matR_SU_valid
end

[nx ny nz] = size(movieDrivenAmp.mask_amp1);
nVox = nx*ny*nz;

% Apply movie-driven mask to correlation matrix
% 2017/01/19 & valid channels
moviemask_vec = reshape(movieDrivenAmp.mask_amp1, nVox, 1); % change the 3D mask to 1D
matR_SU_all_moviemask = matR_SU_all(moviemask_vec,:); %matR_SU(moviemask_vec,:); % 15495 voxels
matR_SU_all_moviemask_valid  = matR_SU_all_moviemask(paramClustering_global.locValidVox, :);

%% Run PCA
[coeff_e, score_e, latent_e, tsquared_e, explained_e] = pca(matR_SU_all_moviemask_valid', 'Economy', false);
figure; 
plot(explained_e(1:20), 'o-');

%% Get indices of voxel clusters
nPC = 4;

% matIndClust_Vox = cat(2, Clustering_moviemask.resultKMeans.Vox_indCluster);
for iPC = 5; %1:nPC
%     sortTargetK = paramClustering_global.setK(iPC);
%     % [sortedClust, indSortVoxel]=sort(matIndClust_Vox(:,sortTargetK-1));
    
    [nx ny nz] = size(movieDrivenAmp.mask_amp1);
    nVox = nx*ny*nz;
    moviemask_vec = reshape(movieDrivenAmp.mask_amp1, nVox, 1);
    
    tempPC_valid = zeros(size(paramClustering_global.locValidVox));
    tempPC_valid(paramClustering_global.locValidVox) = coeff_e(:, iPC); %Clustering_moviemask_valid.resultKMeans(sortTargetK-1).Vox_indCluster;

    PCMap = zeros(size(moviemask_vec));
    PCMap(moviemask_vec) = tempPC_valid;
    
%     voxelClusterMap = zeros(size(moviemask_vec));
%     voxelClusterMap(moviemask_vec) = Clustering_moviemask.resultKMeans(sortTargetK-1).Vox_indCluster;
    
    % voxelClusterMap_vol = reshape(voxelClusterMap./sortTargetK, [nx, ny, nz]); % interim work-around for labeling of different clusters in AFNI
    PCMap_vol = reshape(PCMap, [nx, ny, nz]); % each voxel now has cluster identity (integers)
    
    
    
    pname = [dirDataBOLD, '/tempSURF/'];
    dirSPEC = [dirDataBOLD, '/Anatomy/_suma/'];
    
    cd /projects/parksh/_toolbox/BlockAna/ %/projects/parksh/NeuralBOLD/analysis/BlockAna/
    blockana;
    S_neuralRegressor(nameSubjNeural, nameSubjBOLD);
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
    fileHead = sprintf('new_masked_%s%s', cell2mat(setNameSubjNeural), nameSubjBOLD); %sprintf('new_maskedSignificant_%s', nameSubjNeural);% sprintf('new_masked_%s', nameSubjNeural);
    DSP.proc.fncvol_3d = PCMap_vol; %mapR_Cluster(:,:,:,iK).*movieDrivenAmp.mask_amp1; %reshape(mapR, [nx, ny, nz]).*movieDrivenAmp.mask_amp1;
    %         end
    
    %         fname = sprintf('%s_Cluster%d_%dMeans_Art_AvgCorrMapMovie123_noFiltering+orig.BRIK', fileHead, iK, sortTargetK);%
    fname = sprintf('%s_PrinComp%d_CritCorr1_Movie123_noFiltering+orig.BRIK', fileHead, iPC);%
    
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


