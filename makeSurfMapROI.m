% makeSurfMapROI.m

cd /projects/parksh/NeuralBOLD/analysis/BlockAna
blockana;

global DATA DSP 

nameSubjNeural = 'Tor';
nameSubjBOLD = 'Art';

% load fmri ROIs
dirDataHome = '/procdata/parksh/';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

dirROI = fullfile(dirDataBOLD, 'ROIs');
if ~exist(dirROI)
    visROIs=[]; faceROIs=[];
else
    d_vis = dir(fullfile(dirROI, '*VisROIs.mat'));
    d_face = dir(fullfile(dirROI, '*faceROIs2.mat'));
    load(fullfile(dirROI, d_vis.name));
    load(fullfile(dirROI, d_face.name));
    
    % get ROI information into structures with different names
    tempName_visROIs = char(fieldnames(load(fullfile(dirROI, d_vis.name))));
    tempName_faceROIs = char(fieldnames(load(fullfile(dirROI, d_face.name))));
    eval(['visROIs=', tempName_visROIs, ';']) % get visual ROI data into "visROIs"
    eval(['faceROIs=', tempName_faceROIs, ';']) % get face ROI data into "faceROIs"
    
    % clear up
    eval(['clear ' tempName_visROIs])
    eval(['clear ' tempName_faceROIs])
end

% Get ROI names & coordinates in EPI coords 
resizeFactor = DSP.proc.params3d.res./faceROIs(1).params.res; % calculate the scaling factor

% convert the map to the surface 
% First, get all the map together
[dvolx, dvoly, dvolz]  = size(visROIs(1).vol3D); % dimension of volume

catAllVisROIs = cat(4, visROIs.vol3D);
matAllVisROIs = reshape(catAllVisROIs, dvolx*dvoly*dvolz, length(visROIs));
catAllFaceROIs = cat(4, faceROIs.vol3D);
matAllFaceROIs = reshape(catAllFaceROIs, dvolx*dvoly*dvolz, length(faceROIs));

allVisROIsInd= sum(matAllVisROIs,2); %   matAllVisROIs*[1:8]';
% tempAllVisROIs = NaN(size(allVisROIsInd));
% tempAllVisROIs(allVisROIsInd>0) = 1;
allVisROIsInd_vol =  reshape(allVisROIsInd, [dvolx, dvoly, dvolz] );
allVisROIsInd_fun = decimate3D(allVisROIsInd_vol, resizeFactor, .25); % turn anat_res ROIs into func_res ROIs

allFaceROIsInd = matAllFaceROIs*[1:8]'; %sum(matAllFaceROIs,2); %matAllFaceROIs*[1:8]';
allFaceROIsInd_vol =  reshape(allFaceROIsInd, [dvolx, dvoly, dvolz] );
allFaceROIsInd_fun = decimate3D(allFaceROIsInd_vol, resizeFactor, .25); % turn anat_res ROIs into func_res ROIs

targetMapVol{1}.vol = allVisROIsInd_fun;
targetMapVol{1}.fname = sprintf('%s_ROI_%s+orig.BRIK', nameSubjBOLD, 'AllVisROI');%

targetMapVol{2}.vol = allFaceROIsInd_fun;
targetMapVol{2}.fname = sprintf('%s_ROI_%s+orig.BRIK', nameSubjBOLD, 'AllFaceROI');%

% % colormap for 8 different ROIs (from FreeSurferColorLUT)
% cMap = summer(8);
% cMap = cat(2, cMap, [0:1:7]');
% ID = [1:8]';
% nameVisROI = {visROIs(1).name(3:end), visROIs(2).name(3:end), visROIs(3).name(3:end), visROIs(4).name(3:end), visROIs(5).name(3:end),...
%     visROIs(6).name(3:end), visROIs(7).name(3:end), visROIs(8).name(3:end)}';
% nameFaceROI = {faceROIs(1).name(3:end), faceROIs(2).name(3:end), faceROIs(3).name(3:end), faceROIs(4).name(3:end), faceROIs(5).name(3:end),...
%     faceROIs(6).name(3:end), faceROIs(7).name(3:end), faceROIs(8).name(3:end)}';
% 
% T = table(ID, nameVisROI, cMap(:,1), cMap(:,2), cMap(:,3), ones(8,1) , 'VariableNames', {'ID', 'Label', 'R', 'G', 'B', 'A'});
% writetable(T, 'visROIsLUT.txt')    
% T = table(ID, nameFaceROI, cMap(:,1), cMap(:,2), cMap(:,3), ones(8,1), 'VariableNames', {'ID', 'Label', 'R', 'G', 'B', 'A'});
% writetable(T, 'faceROIsLUT.txt')    
 

pname = [dirDataBOLD, '/tempSURF/'];
dirSPEC = [dirDataBOLD, '/Anatomy/_suma/'];


for iMap = 1:length(targetMapVol) % iFR = 1:length(faceROIs)
    
    cd /projects/parksh/NeuralBOLD/analysis/BlockAna/
%     cellID = validChanID(iUnit,:);
%     fprintf(1, 'Unit # %d, Cell ID: %s, %s \n', iUnit, cellID);

%     curFaceROI=decimate3D(faceROIs(iFR).vol3D, resizeFactor, .25); % turn anat_res ROIs into func_res ROIs
    
    fname = targetMapVol{iMap}.fname;  %  fname = sprintf('%s_ROI_%s+orig.BRIK', nameSubjBOLD, faceROIs(iFR).name(3:end));%
%     curFaceROI=decimate3D(faceROIs(iFR).vol3D, resizeFactor, .25);
    DSP.proc.fncvol_3d = targetMapVol{iMap}.vol;
    vol = single(DSP.proc.fncvol_3d);  %targetMapVol{iMap}.vol;
%     vol = curFaceROI; %single(mapR_Cluster(:,:,:,iK)); % vol = single(DSP.proc.fncvol_3d);  
    
    
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


