% computeCorrMap_sigCI.m

nameSubjNeural = 'Tor';
nameSubjBOLD = 'Art';

dirDataHome = '/procdata/parksh/';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);
dirDataCI = fullfile(dirDataNeural, 'corrMap_resultsCI/resultsCI_eachCell');


% Original correlation map
load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123.mat', nameSubjNeural, nameSubjBOLD)), 'paramCorr', 'matR_SU') % get cell IDs from another file

%% 
% 3) fMRI movie-driven activity mask
load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp', 'brainMask_BlockAna3D');

[a, b, c] = ind2sub(size(movieDrivenAmp.mask_amp1), find(movieDrivenAmp.mask_amp1==1)); % Get indices of AF voxels in EPI 3D coords
mask_amp1 = [a b c];
mask_amp1_sub = sub2ind(size(movieDrivenAmp.mask_amp1), a, b, c);

% convert the map to the surface 
for iUnit = 1:size(matR_SU,2)
    
%     cd /projects/parksh/NeuralBOLD/analysis/BlockAna/
    cellID = paramCorr.validChanID(iUnit,:);
    fprintf(1, 'Unit # %d, Cell ID: %s \n', iUnit, cellID);
    
%     fname = sprintf('%s_%s_Movie123_CI95_noFiltering+orig.BRIK', cellID, nameSubjNeural);%
    
    
    fileNameCI = sprintf('corrcoeffCI_%s%s_cell%d.mat',...
        nameSubjNeural, nameSubjBOLD, iUnit);
    load(fullfile(dirDataCI, fileNameCI))
    
    catCI95 = cat(1, resultBS(1).VoxCI.matCI95);
    catCI99 = cat(1, resultBS(1).VoxCI.matCI99);
    
    CI95_in = catCI95(mask_amp1_sub, :);
    setCI95_in(iUnit, :) = mean(CI95_in);
    % catRho_org = cat(1, resultBS(1).VoxCI.rho_org).*(-1);
    corrMap_org = matR_SU(:,iUnit);
end



    
    aa=cat(2, corrMap_org<catCI95(:,1), corrMap_org>catCI95(:,2));
    sigMask_CI95 = sum(aa,2);
    
    aa2=cat(2, corrMap_org<catCI99(:,1), corrMap_org>catCI99(:,2));
    sigMask_CI99 = sum(aa2,2);
    
    
    
    %% fMRI map for each SU
    % average maps for each cluster
    nx = 40; ny = 64; nz = 32;
    
    corrMap_CI95 = corrMap_org.*sigMask_CI95;
    corrMap_CI99 = corrMap_org.*sigMask_CI99;
    
    corrMap_CI95_3d = reshape(corrMap_CI95, [nx, ny, nz]);
    corrMap_CI99_3d = reshape(corrMap_CI99, [nx, ny, nz]);

    
    
    DSP.proc.fncvol_3d = corrMap_CI95_3d; %corrMap_CI99_3d; %reshape(matR_SU(:,iUnit), [nx, ny, nz]);
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
