function [] = saveSUMaps_Surface_pcares_masked(nameSubjNeural, setMovie)
%
% 2018/10/24


% nameSubjNeural = 'Tor'; %'Dav'; %'Spi'; %'Mat'; %'Ava'; %'Mat'; %'Spi'; %'Sig'; %'Rho'; % 'Sig'; %'Tor';
nameSubjBOLD = 'Art'; %'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';

%
addpath('/library/matlab_utils/')

% Load files
dirDataHome = '/procdata/parksh/_macaque';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% setMovie = 1; % [1 2 3]; %[4 5 6]; %[1 2 3]; %[4 5 6]; %[1 2 3 4 5 6]; %
tempS = num2str(setMovie);
MovieStr = tempS(~isspace(tempS));

% 1) fMRI correlation maps
% if sum(strcmpi(nameSubjNeural, {'spice', 'spi'}))
%     dirDataNeural = '/procdata/parksh/Spi/2016Nov_movie/';
% end
load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie%s_pcares_masked.mat', nameSubjNeural, nameSubjBOLD, MovieStr)), 'matR_SU', 'paramCorr')

% 2) fMRI movie-driven activity mask
load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp');



%% fMRI map for each SU
nx = 40; ny = 64; nz = 32;
nVox = nx*ny*nz;

% Apply movie-driven mask to correlation matrix 
% 2017/01/19 & valid channels 
moviemask_vec = reshape(movieDrivenAmp.mask_amp1, nVox, 1); % change the 3D mask to 1D

pname = [dirDataBOLD, '/tempSURF/'];
dirSPEC = [dirDataBOLD, '/Anatomy/_suma/', nameSubjNeural, '/']; %[dirDataBOLD, '/Anatomy/_suma/'];
% if sum(strcmpi(nameSubjNeural, {'spice', 'spi'}))
%     dirSPEC = [dirDataBOLD, '/Anatomy/_suma/', nameSubjNeural,  '/_2016Nov/']; %'/_2018Jan/'];
% end
if ~exist(dirSPEC, 'dir')
    mkdir(dirSPEC)
end
    

cd /projects/parksh/_toolbox/BlockAna/ %/projects/parksh/NeuralBOLD/analysis/BlockAna/
blockana;
S_neuralRegressor(nameSubjNeural, nameSubjBOLD); %S_neuralRegressor
global STDPATH DSP DATA GH


    
% convert the map to the surface
for iUnit = 1:size(matR_SU,2)
    
    cd /projects/parksh/_toolbox/BlockAna/ %/projects/parksh/NeuralBOLD/analysis/BlockAna/
    cellID = paramCorr.validChanID(iUnit,:);
    fileHead = 'new_masked_';
    fname = sprintf('%s%s_%s_Movie%s_1PCres_noFiltering+orig.BRIK', fileHead, nameSubjNeural, cellID, MovieStr);%
    
    fprintf(1, 'Unit # %d, Cell ID: %s, %s \n', iUnit, cellID);
    
    pcaresCorrMap = zeros(size(moviemask_vec));
    pcaresCorrMap(moviemask_vec) = matR_SU(:, iUnit);
    
    DSP.proc.fncvol_3d = reshape(pcaresCorrMap, [nx, ny, nz]); %mapR_Cluster(:,:,:,iK).*movieDrivenAmp.mask_amp1; %reshape(mapR, [nx, ny, nz]).*movieDrivenAmp.mask_amp1;
    vol = single(DSP.proc.fncvol_3d);
    
    
    
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

cd /projects/parksh/NeuralBOLD/analysis



