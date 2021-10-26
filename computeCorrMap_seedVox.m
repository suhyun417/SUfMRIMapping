% computeCorrMap_seedVox.m
% function [matR_seedvox, matP_seedvox, paramCorr] = computeCorrMap(nameSubjBOLD, flagSaveFile)
%
% get correlation of AF voxels
%


%% Compute correlation for movie 123
%
addpath('/library/matlab_utils/')

nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';

% Load data files
dirDataHome = '/procdata/parksh/';
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);
dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';

filenameBOLD = [nameSubjBOLD, '_movieTS_fMRI_indMov.mat'];

fprintf(1, '\nLoading fMRI data of %s: %s ....\n', nameSubjBOLD, filenameBOLD)
load(fullfile(dirDataBOLD, filenameBOLD))

% Get movie IDs common in two dataset
setMovie = [1 2 3];
paramCorr.setMovie = setMovie;

% 1. fMRI tc
fmritc=[];
% indMovieBOLD = find(ismember(dataBOLD.unimov, setMovie)>0);
for iM = 1:length(setMovie)
    curvoltc = voltcIndMov{setMovie(iM)}; %dataBOLD.mvoltc{iM};
    avgvoltc = repmat(nanmean(curvoltc,4),[1 1 1 size(curvoltc,4)]);
    if ~isempty(find(avgvoltc==0, 1))
        avgvoltc(avgvoltc==0) = realmin; % get rid of zeros because it causes NaNs in percent signals
    end
    pcvoltc = ((curvoltc - avgvoltc)./avgvoltc)*100;
    fmritc = cat(4,fmritc,pcvoltc);
end

% Get the ROI
dirROI = fullfile(dirDataBOLD, 'ROIs');
d_face = dir(fullfile(dirROI, '*faceROIs2.mat'));
load(fullfile(dirROI, d_face.name));
tempName_faceROIs = char(fieldnames(load(fullfile(dirROI, d_face.name))));
eval(['faceROIs=', tempName_faceROIs, ';']) % get face ROI data into "faceROIs"
eval(['clear ' tempName_faceROIs]) % clear up

%     d_vis = dir(fullfile(dirROI, '*VisROIs.mat'));
%     load(fullfile(dirROI, d_vis.name));
%     tempName_visROIs = char(fieldnames(load(fullfile(dirROI, d_vis.name))));
%     eval(['visROIs=', tempName_visROIs, ';']) % get visual ROI data into "visROIs"
%     eval(['clear ' tempName_visROIs])

% Get ROI names & coordinates in EPI coords
resizeFactor = [3 3 3]; % DSP.proc.params3d.res./ROIs(1,1).params.res; % calculate the scaling factor

[tempNameROI{1:length(faceROIs)}] = deal(faceROIs.name);
for iF = 1:length(tempNameROI)
    name_ROI{iF} = tempNameROI{iF}(strfind(tempNameROI{iF}, '_')+1:end);
end
clear tempNameROI

% % Which ROI do you want?
% targetROI = 'ML'; % 'AM'; %'AF';
% indROI = strcmp(targetROI, name_ROI); %

for indROI = [1 2 5 6 8] %[3 7]; %
    targetROI = name_ROI{indROI};
    % switch lower(nameSubjBOLD)
    %     case 'art'
    %         indROI = 4;
    %     case 'ava'
    %         indROI = 7;
    % end
    % nameROI = faceROIs(indROI).name(strfind(faceROIs(indROI).name, '_')+1:end); %faceROIs(iFR).name(strfind(faceROIs(iFR).name, '_')+1:end);
    voxROI=decimate3D(faceROIs(indROI).vol3D, resizeFactor, .25); %decimate3D(faceROIs(iFR).vol3D, resizeFactor, .25); % turn anat_res ROIs into func_res ROIs
    [a, b, c] = ind2sub(size(voxROI), find(voxROI==1)); % Get indices of AF voxels in EPI 3D coords
    indVox_ROI = [a b c];
    indVox_ROI_sub = sub2ind(size(voxROI), a, b, c);
    
    matBOLD_ROI=[];
    for iVox = 1:length(a)
        matBOLD_ROI(:,iVox) = fmritc(a(iVox), b(iVox), c(iVox),:); %matBOLD_shuffle(a(iVox), b(iVox), c(iVox),:); %matBOLD(a(iVox), b(iVox), c(iVox),:);
    end
    meanBOLD_ROI = nanmean(matBOLD_ROI, 2);
    steBOLD_ROI = nanstd(matBOLD_ROI, [], 2)./sqrt(size(matBOLD_ROI,2)-1);
    % matR = corr(matBOLD_ROI, 'rows', 'complete');
    % vecValidR = matR(triu(matR, 1)~=0); % R values above the diagonal
    
    
    % Compute correlation maps
    [nx, ny, nz, nt] = size(fmritc);
    nVox = nx*ny*nz;
    
    [Rvals_meanROI, Pvals] = corr(reshape(fmritc, nVox, nt)', meanBOLD_ROI,...
        'rows','complete', 'type', 'Spearman');
    
    iVox=2; %8;
    curRGR = matBOLD_ROI(:,iVox);% meanBOLD_ROI;
    [Rvals_seedVoxel, Pvals] = corr(reshape(fmritc, nVox, nt)', curRGR,...
        'rows','complete', 'type', 'Spearman');
    
    mapR_meanROI = reshape(Rvals_meanROI, [nx, ny, nz]); % Don't need to multiply -1 because now it's correlation between MION signals %.*(-1);
    mapR_seedVoxel = reshape(Rvals_seedVoxel, [nx, ny, nz]); % Don't need to multiply -1 because now it's correlation between MION signals %.*(-1);
    
    
    % save seed voxel map onto the surface
    % 3) fMRI movie-driven activity mask
    load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp');
    
    
    
    %% fMRI map for each SU
    % average maps for each cluster
    nx = 40; ny = 64; nz = 32;
    
    pname = [dirDataBOLD, '/tempSURF/'];
    dirSPEC = [dirDataBOLD, '/Anatomy/_suma/'];
    
    cd /projects/parksh/_toolbox/BlockAna/ %analysis/BlockAna/
    blockana;
    S_neuralRegressor('Tor', 'Art')
    global STDPATH DSP DATA GH
    
    % convert the map to the surface
    
    for iMap = 1:2
        
        cd /projects/parksh/_toolbox/BlockAna/ %/projects/parksh/NeuralBOLD/analysis/BlockAna/
        
        switch iMap
            case 1 % meanROI
                curMap = mapR_meanROI;
                type = 'meanROI';
            case 2 % one seed voxel
                curMap = mapR_seedVoxel;
                type = 'oneVoxel';
        end
        
        for iMask = 1: 2
            switch iMask
                case 1 %unmasked version
                    fname = sprintf('seedMap_%s_Movie123_%s_%s_unmasked_noFiltering+orig.BRIK', nameSubjBOLD, targetROI, type);%sprintf('new_masked_%s_%s_Movie123_noFiltering+orig.BRIK', cellID, nameSubjNeural);%
                    DSP.proc.fncvol_3d = reshape(curMap, [nx, ny, nz]); %.*movieDrivenAmp.mask_amp1;
                    
                case 2 % masked version
                    fname = sprintf('seedMap_%s_Movie123_%s_%s_noFiltering+orig.BRIK', nameSubjBOLD, targetROI, type);%sprintf('new_masked_LFP_%d%d_%s_Movie123_noFiltering+orig.BRIK', rangeFq(1), rangeFq(2), nameSubjNeural);%
                    DSP.proc.fncvol_3d = reshape(curMap, [nx, ny, nz]).*movieDrivenAmp.mask_amp1;
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
end









