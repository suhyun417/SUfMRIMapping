% genCritCorr1ValidVoxMap_pcaresAllCells.m
%
% 2019/07/01 SHP
% To visualize the location of the voxels that are used for clustering
% based on the pca-residual maps with all cells, make & save the whole
% brain map of indices of valid 579 voxels

clear all;

%% Load the data
load('/procdata/parksh/_macaque/Art/Clustering_TorRhoSigSpiMatDanMocWasArtMovie123_pcares_masked_critCorr1.mat', 'paramClustering_global')
load('/procdata/parksh/_macaque/Art/Art_MaskArrays.mat', 'movieDrivenAmp');


%% Make a map with indices of 579 valid voxels
[nx, ny, nz] = size(movieDrivenAmp.mask_amp1);
nVox = nx*ny*nz;

moviemask_vec = reshape(movieDrivenAmp.mask_amp1, nVox, 1);

validVoxelMap = zeros(size(moviemask_vec));
validVoxelMap(moviemask_vec) = paramClustering_global.locValidVox; %voxIndCluster_valid;

voxelClusterMap_vol = reshape(validVoxelMap, [nx, ny, nz]); % each voxel now has cluster identity (integers)


%% Save maps to AFNI/SUMA
dirDataBOLD = '/procdata/parksh/_macaque/Art';
pname = [dirDataBOLD, '/tempSURF/'];
dirSPEC = [dirDataBOLD, '/Anatomy/_suma/'];

cd /projects/parksh/_toolbox/BlockAna/ %/projects/parksh/NeuralBOLD/analysis/BlockAna/
blockana;
S_neuralRegressor('Tor', 'Art') %arbitrary 
global STDPATH DSP DATA GH


fname = 'ValidVoxMap_Clustering_TorRhoSigSpiMatDanMocWasArtMovie123_pcares_masked_critCorr1+orig.BRIK';

DSP.proc.fncvol_3d = voxelClusterMap_vol; 
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

cd /projects/parksh/NeuralBOLD/analysis

