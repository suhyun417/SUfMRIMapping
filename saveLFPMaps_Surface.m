function [] = saveLFPMaps_Surface(nameSubjNeural, nameSubjBOLD)
% saveLFPMaps_Surface.m
%
% 1) load fMRI correlation maps computed based on Toroid's cells and Artemis
% or Avalanche's fMRI data  
% 2) load clustering results based on various values (e.g. whole brain corr maps, corr between cells and movie
% regressors, cell spike density function etc.)
% 3) for each cluster, average the corr maps 
% 4) save it as AFNI brik/head format so that it can be mapped onto the surface anatomy

% nameSubjNeural = 'Tor'; % 'Sig'; %'Rho'; % 'Sig'; %'Tor';
% nameSubjBOLD = 'Art'; %'Ava'; %'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';

% 
addpath('/library/matlab_utils/')

% Load files
dirDataHome = '/procdata/parksh/';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% 1) fMRI correlation maps
load(fullfile(dirDataNeural, sprintf('CorrMap_LFP_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'matR_LFP', 'paramCorr')
 
%  % 2) Clustering results
%  load(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD))) %
% % load(fullfile(dirDataNeural, 'Clustering_TorArtAvaMovie123.mat')) % based on two sets of whole brain maps from two different monkeys
% % Clustering = ClusteringAll;

% 3) fMRI movie-driven activity mask
load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp');
 


%% fMRI map for each SU
% average maps for each cluster
nx = 40; ny = 64; nz = 32;


pname = [dirDataBOLD, '/tempSURF/'];
dirSPEC = [dirDataBOLD, '/Anatomy/_suma/'];

cd /projects/parksh/NeuralBOLD/analysis/BlockAna/
blockana;
S_neuralRegressor
global STDPATH DSP DATA GH 

% convert the map to the surface 
for iFreq = 1:size(matR_LFP,2)
    
    cd /projects/parksh/NeuralBOLD/analysis/BlockAna/
%     cellID = paramCorr.validChanID(iFreq,:);
    rangeFq = paramCorr.blpFreq(iFreq, :);
    fprintf(1, 'iFreq # %d, Frequency range: %d - %d \n', iFreq, rangeFq(1), rangeFq(2));
    
    % case1: unmasked version
    fname = sprintf('LFP_%d%d_%s_Movie123_noFiltering+orig.BRIK', rangeFq(1), rangeFq(2), nameSubjNeural);%sprintf('new_masked_%s_%s_Movie123_noFiltering+orig.BRIK', cellID, nameSubjNeural);%
    DSP.proc.fncvol_3d = reshape(matR_LFP(:,iFreq), [nx, ny, nz]); %.*movieDrivenAmp.mask_amp1;
    
%     % case2: masked version
%     fname = sprintf('new_masked_LFP_%d%d_%s_Movie123_noFiltering+orig.BRIK', rangeFq(1), rangeFq(2), nameSubjNeural);%
%     DSP.proc.fncvol_3d = reshape(matR_LFP(:,iFreq), [nx, ny, nz]).*movieDrivenAmp.mask_amp1;
    
    vol = single(DSP.proc.fncvol_3d);  %single(mapR_Cluster(:,:,:,iK)); % vol = single(DSP.proc.fncvol_3d);  
    
    
    %% Do dumpFunctionalBrik
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


