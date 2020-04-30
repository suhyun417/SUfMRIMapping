% genMovieDrivenBOLDMap_cv.m
%
% compute cofficient of variance (std/mean)
% and its inverse
% and convert it into suma space


nameSubjNeural = 'Tor';
nameSubjBOLD = 'Art';

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
addpath(fullfile(dirLibrary, 'matlab_utils'))
addpath(fullfile(dirProjects, 'parksh/_toolbox/Boot_Time_Series'))

% Set directories 
% nameSubjNeural = 'Tor';
% nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataHome = fullfile(dirProcdata, 'parksh');
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% Directory for saving figures as graphic files
dirFig = fullfile(dirProjects, 'parksh/NeuralBOLD/_labNote/_figs/');


%% Load the original time series
filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
% fileNameNeural_BLP = [nameSubjNeural, '_movieTS_BLPLFP_indMov.mat'];
filenameBOLD = [nameSubjBOLD, '_movieTS_fMRI_indMov.mat'];

% fprintf(1, '\nLoading single unit data of %s: %s ....', nameSubjNeural, filenameNeural)
% load(fullfile(dirDataNeural, filenameNeural))
% load(fullfile(dirDataNeural, fileNameNeural_BLP))
fprintf(1, '\nLoading fMRI data of %s: %s ....\n', nameSubjBOLD, filenameBOLD)
load(fullfile(dirDataBOLD, filenameBOLD))


%% Compute CV and CV inverse (tSNR)
% Get movie IDs common in two dataset
setMovie = [1 2 3];
indMovieBOLD = find(ismember(paramBOLD.unimov, setMovie)>0); 

% 1. fMRI tc in percent signal
[nx, ny, nz, nt] = size(voltcIndMov{1});
nVox = nx*ny*nz;

voltc = []; 
voltc = cat(4, voltcIndMov{indMovieBOLD});
voltc = reshape(voltc, nVox, nt*3)';

avgvoltc = repmat(nanmean(voltc),[size(voltc,1), 1]); %
if ~isempty(find(avgvoltc==0, 1))
    avgvoltc(avgvoltc==0) = realmin; % get rid of zeros because it causes NaNs in percent signals
end

vol_cv = nanstd(voltc)./nanmean(voltc);
vol_cvinv = nanmean(voltc)./nanstd(voltc);

%% Make map
map_cv = reshape(vol_cv, [nx, ny, nz]);
map_cvinv = reshape(1./vol_cvinv, [nx, ny, nz]);

voltc_long = []; 
voltc_long = cat(4, voltcIndMov{:});
voltc_long = reshape(voltc_long, nVox, nt*15)';

% avgvoltc = repmat(nanmean(voltc_long),[size(voltc_long,1), 1]); %
% if ~isempty(find(avgvoltc==0, 1))
%     avgvoltc(avgvoltc==0) = realmin; % get rid of zeros because it causes NaNs in percent signals
% end

vol_cvinv_long = nanmean(voltc_long)./nanstd(voltc_long);

%% Make map
map_cvinv_long = reshape(1./vol_cvinv_long, [nx, ny, nz]);
% Then go do blockana -> run "S_neuralRegressor"
% Or do manual check of voxel positions
% addpath('./BlockAna')
% [avol,Info] = loadBrik(fullfile(dirDataBOLD, 'Anatomy'),'e66_t1_avg_reg2+orig.BRIK');
% [avol_3d,params3d] = canonicalizeFor3DDisplay(avol,Info);


%% Save a flat map
pname = [dirDataBOLD, '/tempSURF/'];
dirSPEC = [dirDataBOLD, '/Anatomy/_suma/'];

cd /projects/parksh/NeuralBOLD/analysis/BlockAna/
blockana;
S_neuralRegressor
global STDPATH DSP DATA GH 

% convert the map to the surface 
% for iUnit = 1:size(matR_SU,2)
    
    cd /projects/parksh/NeuralBOLD/analysis/BlockAna/
    
    fname = sprintf('%s_Movie123_tSNR+orig.BRIK', nameSubjBOLD);%
    

    nx = 40; ny = 64; nz = 32;
        
    DSP.proc.fncvol_3d = map_cvinv; %map_cv; %map_avgAmp; %reshape(matR_SU(:,iUnit), [nx, ny, nz]);
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
% end



















