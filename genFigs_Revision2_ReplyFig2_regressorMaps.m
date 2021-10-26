% genFigs_Revision2_ReplyFig2_regressorMaps.m
%
% 2017/06/02

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

% Set directories
nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataHome = fullfile(dirProcdata, 'parksh');
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);


% % Clustering based on corr maps
% load(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new_masked.mat', nameSubjNeural, nameSubjBOLD)));
% % % Clustering based on time series
% % load(fullfile(dirDataNeural, sprintf('ClusteringSDF_%sMovie123.mat', nameSubjNeural)))

% Directory for saving figures as graphic files
dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';

%% fMRI data
filenameBOLD = [nameSubjBOLD, '_movieTS_fMRI_indMov.mat'];
fprintf(1, '\nLoading fMRI data of %s: %s ....\n', nameSubjBOLD, filenameBOLD)
load(fullfile(dirDataBOLD, filenameBOLD))

% Get movie IDs common in two dataset
setMovie = [1 2 3];
dataBOLD.mvoltc = voltcIndMov(setMovie);
dataBOLD.unimov = setMovie;

% 1. fMRI tc
fmritc=[];
indMovieBOLD = find(ismember(dataBOLD.unimov, setMovie)>0);
for iM = indMovieBOLD %1:length(indMovieBOLD)
    curvoltc = dataBOLD.mvoltc{iM};
    avgvoltc = repmat(nanmean(curvoltc,4),[1 1 1 size(curvoltc,4)]);
    if ~isempty(find(avgvoltc==0, 1))
        avgvoltc(avgvoltc==0) = realmin; % get rid of zeros because it causes NaNs in percent signals
    end
    pcvoltc = ((curvoltc - avgvoltc)./avgvoltc)*100;
    fmritc = cat(4,fmritc,pcvoltc);
end

% Reshape BOLD 4-d data
[nx, ny, nz, nt] = size(fmritc);
nVox = nx*ny*nz;

%%
% Movie regressors
flagSM = 1; % flag for compression and smoothing

fullRGR4fps = createMovieRGR_4fps_indMov(setMovie, flagSM); %createFullMovieRegressors_4fps_indMov(setMovID); %
% ttt=load('/procdata/parksh/MovieRegressors/dbtmMriReg.mat'); % Face scale regressor (in TR unit)
% scaleRGR = ttt.reg.xx(7,:)';

% full regressors
catRGRfull=[];matRGRfull=[];
for iMov=1:length(setMovie)
    m = setMovie(iMov);
    matCurRGR = fullRGR4fps(m).smoRegressors; %fullRGR4fps(iMov).regressors(:,indValidRGR); %fullRGR4fps(iMov).regressors;
    catRGRfull = cat(1, catRGRfull, matCurRGR); % concatenation across movies
end
% scaleRGR_resampled = resample(scaleRGR, 2.4*100, 0.25*100); %matRGR = resample(catRGR, 0.25*100, 2.4*100);
matRGRfull = resample(catRGRfull, 0.25*100, 2.4*100); %catRGRfull; %cat(2, catRGRfull, scaleRGR_resampled); %cat(2, matRGR, scaleRGR);
varnamesfull = fullRGR4fps(1).features; % cat(1, fullRGR4fps(1).features, {'Face size'});

% subset of regressors
indValidRGR = [1 2 6 7 3 21 20 32 22 31 25]; %[1, 3, 9, 20, 21, 22, 25];
% 1: 'Luminance', 2: 'Contrast', 6: Low spatial Frequency 7: High spatial frequencty 3: 'Motion (Speed)',
% 21: 'One face', 20: 'Number of faces', 32: 'Face size', 22: 'Body parts', 31: 'Hands', 25: 'Any animal'
% matRGRvalid = matRGRfull(:,indValidRGR);
varnamesvalid = varnamesfull(indValidRGR);

% Compute correlation between fMRI TS and feature TS
for iRGR = 1:size(matRGRfull,2)
    
    Rvals = [];
    [Rvals, Pvals] = corr(reshape(fmritc, nVox, nt)', matRGRfull(:,iRGR),...
        'rows','complete', 'type', 'Spearman');
    
    corrMap_movieRGRfull(:, iRGR) = Rvals*(-1); % because of MION
    
    
    %     r_c=[]; r_su=[];
    %
    %     % averaged TS in each cluster
    %     r_c = corr(matRGRfull(:,iRGR), meanFRCluster4fps, 'rows', 'complete', 'type', 'Spearman');
    %     % single unit TS
    %     r_su = corr(matRGRfull(:,iRGR), matFR4fps, 'rows', 'complete', 'type', 'Spearman');
    %
    %     R_ClusterMovieRGRfull(iRGR, :) = r_c;
    %     R_SUmovieRGRfull(iRGR, :) = r_su;
    
end

corrMap_movieRGRvalid = corrMap_movieRGRfull(:,indValidRGR);

paramCorr.varnamesfull = varnamesfull;
paramCorr.varnamesvalid = varnamesvalid;
paramCorr.indValidRGR = indValidRGR;

save(fullfile(dirDataBOLD, sprintf('CorrMap_movieRGR_%sMovie123.mat', nameSubjBOLD)),...
    'corrMap_movieRGRfull', 'corrMap_movieRGRvalid', 'paramCorr');


%% save maps to the surface
addpath('/library/matlab_utils/')

nameSubjNeural = 'Tor';
nameSubjBOLD = 'Art';

% Load files
dirDataHome = '/procdata/parksh/';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);


% 1) fMRI correlatio map of movie regressors
load(fullfile(dirDataBOLD, sprintf('CorrMap_movieRGR_%sMovie123.mat', nameSubjBOLD)),...
    'corrMap_movieRGRfull', 'corrMap_movieRGRvalid', 'paramCorr');

% 2) fMRI movie-driven activity mask
load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp');


%% fMRI map for each SU
% average maps for each cluster
nx = 40; ny = 64; nz = 32;


pname = [dirDataBOLD, '/tempSURF/'];
dirSPEC = [dirDataBOLD, '/Anatomy/_suma/'];

cd /projects/parksh/_toolbox/BlockAna/
blockana;
S_neuralRegressor(nameSubjNeural, nameSubjBOLD)
global STDPATH DSP DATA GH

% convert the map to the surface
for iRGR = 1:size(paramCorr.varnamesvalid,1)
    
    for iMask = 1:2
        cd /projects/parksh/_toolbox/BlockAna/
        %     cellID = paramCorr.validChanID(iFreq,:);
        
        % case1: unmasked version
        switch iMask
            case 1
                fname = sprintf('CorrMap_movieRGR_validRGR%02d_%sMovie123+orig.BRIK', iRGR, nameSubjBOLD);
                DSP.proc.fncvol_3d = reshape(corrMap_movieRGRvalid(:,iRGR), [nx, ny, nz]); %.*movieDrivenAmp.mask_amp1;
            case 2
                % case2: masked version
                fname = sprintf('new_masked_CorrMap_movieRGR_validRGR%02d_%sMovie123+orig.BRIK', iRGR, nameSubjBOLD);%
                DSP.proc.fncvol_3d = reshape(corrMap_movieRGRvalid(:,iRGR), [nx, ny, nz]).*movieDrivenAmp.mask_amp1;
        end
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
end

