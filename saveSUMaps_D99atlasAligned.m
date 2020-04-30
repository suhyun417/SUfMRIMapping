% saveSUMaps_D99atlasAligned.m
%
% 2018/09/11 SHP
% Generate single-unit fMRI maps in D99-atlas segmentation space 
% using afni function "3dresample" with -master option
%
% 1) For each subject,
%     1-1) Make a new directory under /procdata/parksh/nameSubjNeural to save the BRIK/HEAD files
%     1-2) Copy /procdata/parksh/Art/Anatomy/_D99Atlas/e66_t1_avg_reg2_seg+orig to the new directory
%     1-3) Load the D99 segmentation of Artemis' brain (e66_t1_avg_reg2_seg+orig) and save the 3d volume
% 2) Copy the premade *_refit+orig in /procdata/parksh/nameSubjBOLD/tempSURF to a new directory 
%     -- This is assuming that you are running this code AFTER you made all the surface correlation maps using saveSUMaps_Surface.m
%         because the *_refit+orig.BRIK is made for each neuron while running saveSUMaps_Surface.m 
%         and saved in '/procdata/parksh/nameSubjBOLD/tempSURF
% 3) Perform '3dresample -master ./e66_t1_avg_reg2_seg+orig -prefix new.dset -input ./UNIT_refit+orig
% 4) Read the new resampled BRIK file of correlation maps and save the data matrix (3d volume)
% of all cells in /procdata/parksh/nameSubjNeural

%%
addpath('/library/matlab_utils/')
addpath('/projects/parksh/_toolbox/afni_matlab/')

%% Set directories
% Important input
nameSubjNeural = 'Spi'; % 'Moc'; %'Dex'; %'Tor'; %'Dav'; %'Spi'; %'Mat'; %'Ava'; %'Mat'; %'Spi'; %'Sig'; %'Rho'; % 'Sig'; %'Tor';
nameSubjBOLD = 'Art'; %'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
setMovie = [1 2 3]; %1; %[1 2 3]; %[4 5 6]; %[1 2 3 4 5 6]; %

% Load files
dirDataHome = '/procdata/parksh/_macaque';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

tempS = num2str(setMovie);
MovieStr = tempS(~isspace(tempS));

% 1-1) Make a new directory to save newly made BRIK/HEAD files
dirSaveOutputFiles = fullfile(dirDataNeural, sprintf('_SUmaps%sD99Space', nameSubjBOLD));
if ~exist(dirSaveOutputFiles, 'dir')
    mkdir(dirSaveOutputFiles)
end

% 1-2) Copy the D99 atlas segmentation file to the new directory
filename_seg = 'e66_t1_avg_reg2_seg+orig.BRIK';
if ~exist(fullfile(dirSaveOutputFiles, filename_seg), 'file')    
    copyfile(fullfile(dirDataBOLD, '/Anatomy/_D99Atlas/e66_t1_avg_reg2_seg+orig*'), fullfile(dirSaveOutputFiles, '/'))
end

% 1-3) Load the D99 segmentation file (e66_t1_avg_reg2_seg+orig) and save the 3d volume 
filename_segmat = sprintf('%sD99_ROIsegmentation.mat', nameSubjBOLD);
if ~exist(fullfile(dirDataNeural, filename_segmat), 'file')    
    [err, voltc_seg, info_seg] = BrikLoad(fullfile(dirSaveOutputFiles, filename_seg));
    save(fullfile(dirDataNeural, sprintf('%sD99_ROIsegmentation.mat', nameSubjBOLD)), 'voltc_seg', 'info_seg')
else % if the file already exists
    load(fullfile(dirDataNeural, sprintf('%sD99_ROIsegmentation.mat', nameSubjBOLD)), 'info_seg')
end

% % 1) fMRI correlation maps
% load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie%s_new.mat', nameSubjNeural, nameSubjBOLD, MovieStr)), 'matR_SU', 'paramCorr')
% 
% %  % 2) Clustering results
% %  load(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD))) %
% % % load(fullfile(dirDataNeural, 'Clustering_TorArtAvaMovie123.mat')) % based on two sets of whole brain maps from two different monkeys
% % % Clustering = ClusteringAll;
% 
% % 3) fMRI movie-driven activity mask
% load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp', 'brainMask_BlockAna3D');


%% 2) Copy the premade *_refit+orig in /procdata/parksh/nameSubjBOLD/tempSURF to a new directory
dirOrgData = fullfile(dirDataBOLD, 'tempSURF');
d = dir(fullfile(dirOrgData, sprintf('%s*Movie%s_noFiltering_refit+orig.BRIK', nameSubjNeural, MovieStr))); % only convert non-masked maps

numFile = length(d);
catVoltc_new = NaN(info_seg.DATASET_DIMENSIONS(1), info_seg.DATASET_DIMENSIONS(2), info_seg.DATASET_DIMENSIONS(3), numFile);
% catInfo_new = struct([]);
for iFile = 1:numFile
    fname = d(iFile).name;
    doti = strfind(fname,'+');
    prefix = fname(1:doti-1);
    
    % copy data to new file so original dimensions are saved.
    copyDataCMD = sprintf('3dcopy %s %s', fullfile(dirOrgData, prefix), fullfile(dirSaveOutputFiles, prefix)); % ['3dcopy  ' prefix ' ' prefix '_refit'];
    disp(copyDataCMD);
    system(copyDataCMD);
    
    %% 3) Resample (upsample) the single-unit fMRI map to anatomical resolution to be matched to D99 segmentation file
    resampleDataCMD = sprintf('3dresample -master %s -prefix %s_D99resample -input %s', ...
        fullfile(dirSaveOutputFiles, filename_seg), fullfile(dirSaveOutputFiles, prefix), fullfile(dirSaveOutputFiles, fname));
    disp(resampleDataCMD);
    system(resampleDataCMD);
    
    %% 4-1) Read the new resampled BRIK file 
    filename_new = [prefix '_D99resample+orig.BRIK'];
    [err, voltc_new, info_new] = BrikLoad(fullfile(dirSaveOutputFiles, filename_new));
    
    catVoltc_new(:,:,:,iFile) = voltc_new;
    catInfo_new(iFile) = info_new;
end

 % 4-2) Save the data matrix (3d volume)
 save(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie%s_new_D99resample.mat', nameSubjNeural, nameSubjBOLD, MovieStr)),...
     'catVoltc_new', 'catInfo_new');

% % fileHead{1} = '';               
% % fileHead{2} = 'new_masked_';
%                 
% 
% 
%         
%         
%     
%         fname = sprintf('%s%s_%s_Movie%s_noFiltering+orig.BRIK', fileHead, nameSubjNeural, cellID, MovieStr);%
% 
% % 3) Perform '3dresample -master ./e66_t1_avg_reg2_seg+orig -prefix new.dset -input ./UNIT_refit+orig
% % 4) Read the new resampled BRIK file of correlation maps and save the data matrix (3d volume)
% % of all cells in /procdata/parksh/nameSubjNeural
% 
% filename_corrmapbrik = '/procdata/parksh/Art/tempSURF/Dan_05_Movie123_noFiltering_refit+orig.BRIK';
% [a b c] = find(voltc_seg==21)
% tempLocLv = find(voltc_seg==21);
% tempLocLd = find(voltc_seg==9);
% tempLocLA = cat(1, tempLocLv, tempLocLd);
% tempCorrLA = voltc_cmapnew(tempLocLA);
% ls('/procdata/parksh/Art/tempSURF/new_masked_Dan*')
% 
% %% fMRI map for each SU
% % average maps for each cluster
% nx = 40; ny = 64; nz = 32;
% 
% 
% pname = [dirDataBOLD, '/tempSURF/'];
% dirSPEC = [dirDataBOLD, '/Anatomy/_suma/', nameSubjNeural, '/']; %[dirDataBOLD, '/Anatomy/_suma/'];
% if sum(strcmpi(nameSubjNeural, {'spice', 'spi'}))
%     dirSPEC = [dirDataBOLD, '/Anatomy/_suma/', nameSubjNeural, '/_2018Jan/'];
% end
% if ~exist(dirSPEC, 'dir')
%     mkdir(dirSPEC)
% end
%     
% 
% cd /projects/parksh/_toolbox/BlockAna/ %/projects/parksh/NeuralBOLD/analysis/BlockAna/
% blockana;
% S_neuralRegressor(nameSubjNeural, nameSubjBOLD); %S_neuralRegressor
% global STDPATH DSP DATA GH
% 
% % % convert the map to the surface
% % for iMask = 1:2
% %
% %     for iK = 1:sortTargetK
% %
% %         cd /projects/parksh/NeuralBOLD/analysis/BlockAna/
% %
% %         switch iMask
% %             case 1 % unmasked (everything within brain)
% %                 fileHead = 'new';
% %                 DSP.proc.fncvol_3d = mapR_Cluster(:,:,:,iK).*brainMask_BlockAna3D;
% %             case 2 % masked
% %                 fileHead = 'new_masked';
% %                 DSP.proc.fncvol_3d = mapR_Cluster(:,:,:,iK).*movieDrivenAmp.mask_amp1; %reshape(mapR, [nx, ny, nz]).*movieDrivenAmp.mask_amp1;
% %         end
% %
% % %         fname = sprintf('%s_Cluster%d_%dMeans_Art_AvgCorrMapMovie123_noFiltering+orig.BRIK', fileHead, iK, sortTargetK);%
% %         fname = sprintf('%s_Cluster%d_%dMeansMovieDriven_Art_AvgCorrMapMovie123_noFiltering+orig.BRIK', fileHead, iK, sortTargetK);%
% %
% %         vol = single(DSP.proc.fncvol_3d);
% 
% for iMask = 1:2
%     
%     % convert the map to the surface
%     for iUnit = 1:size(matR_SU,2)
%         
%         cd /projects/parksh/_toolbox/BlockAna/ %/projects/parksh/NeuralBOLD/analysis/BlockAna/
%         cellID = paramCorr.validChanID(iUnit,:);
%         fprintf(1, 'Unit # %d, Cell ID: %s, %s \n', iUnit, cellID);
%         
%         switch iMask
%             case 1 % unmasked (everything within brain)
%                 fileHead = '';
%                 DSP.proc.fncvol_3d =  reshape(matR_SU(:,iUnit), [nx, ny, nz]).*brainMask_BlockAna3D; %mapR_Cluster(:,:,:,iK).*brainMask_BlockAna3D;
%             case 2 % masked
%                 fileHead = 'new_masked_';
%                 DSP.proc.fncvol_3d = reshape(matR_SU(:,iUnit), [nx, ny, nz]).*movieDrivenAmp.mask_amp1; %mapR_Cluster(:,:,:,iK).*movieDrivenAmp.mask_amp1; %reshape(mapR, [nx, ny, nz]).*movieDrivenAmp.mask_amp1;
%         end
%         
% %         fname = sprintf('%s%s_%s_Movie123_noFiltering+orig.BRIK', fileHead, nameSubjNeural, cellID);%
% %         vol = single(DSP.proc.fncvol_3d);
%         
%         fname = sprintf('%s%s_%s_Movie%s_noFiltering+orig.BRIK', fileHead, nameSubjNeural, cellID, MovieStr);%
%         vol = single(DSP.proc.fncvol_3d);
%         
%         % %     % case1: unmasked version
%         % %     fname = sprintf('%s_%s_Movie123_noFiltering_new+orig.BRIK', cellID, nameSubjNeural);%sprintf('new_masked_%s_%s_Movie123_noFiltering+orig.BRIK', cellID, nameSubjNeural);%
%         % %     DSP.proc.fncvol_3d = reshape(matR_SU(:,iUnit), [nx, ny, nz]); %.*movieDrivenAmp.mask_amp1;
%         %
%         % %     % case2: masked version
%         % %     fname = sprintf('new_masked_%s_%s_Movie123_noFiltering+orig.BRIK', cellID, nameSubjNeural);%
%         % %     DSP.proc.fncvol_3d = reshape(matR_SU(:,iUnit), [nx, ny, nz]).*movieDrivenAmp.mask_amp1;
%         %
%         % %     vol = single(DSP.proc.fncvol_3d);  %single(mapR_Cluster(:,:,:,iK)); % vol = single(DSP.proc.fncvol_3d);
%         
%         
%         %% Do dumpFunctionalBrik
%         Info = DATA.Info;
%         Info.DATASET_RANK = DATA.Info.DATASET_RANK;
%         Info.DATASET_DIMENSIONS = DATA.Info.DATASET_DIMENSIONS;
%         Info.DELTA = DATA.Info.DELTA;
%         Info.BYTEORDER_STRING = DATA.Info.BYTEORDER_STRING;
%         Info.TYPESTRING = DATA.Info.TYPESTRING;
%         Info.SCENE_DATA = DATA.Info.SCENE_DATA;
%         Info.ORIGIN = DATA.Info.ORIGIN;
%         Info.TAXIS_FLOATS = DATA.Info.TAXIS_FLOATS;
%         Info.TAXIS_OFFSETS = DATA.Info.TAXIS_OFFSETS;
%         Info.IDCODE_STRING = DATA.Info.IDCODE_STRING;
%         Info.IDCODE_DATE = DATA.Info.IDCODE_DATE;
%         Info.HISTORY_NOTE = DATA.Info.HISTORY_NOTE;
%         
%         Info.DATASET_DIMENSIONS = [size(vol) 0 0];
%         Info.DATASET_RANK(2) = 1;
%         Info.BRICK_TYPES = 3;  %3 = float
%         Info.BRICK_FLOAT_FACS = 0;
%         Info.BRICK_STATS = [min(vol(:)) max(vol(:))];
%         Info.BRICK_LABS = DATA.Info.BRICK_LABS(1:2);
%         Info.BRICK_KEYWORDS = DATA.Info.BRICK_KEYWORDS;
%         Info.TAXIS_NUMS(1) = 1;
%         Info.TypeName = 'float';
%         Info.TypeBytes = 4;
%         Info.ByteOrder = 'ieee-le';
%         Info.Orientation = DATA.Info.Orientation;
%         Info.FileFormat = DATA.Info.FileFormat;
%         Info.RootName = DATA.Info.RootName;
%         Info.Extension_1D = DATA.Info.Extension_1D;
%         
%         code0 = getAFNIOrientationCode('LR');
%         code1 = getAFNIOrientationCode('AP');
%         code2 = getAFNIOrientationCode('DV');
%         
%         Info.ORIENT_SPECIFIC = [code0 code1 code2];
%         
%         doti = strfind(fname,'+');
%         prefix = fname(1:doti-1);
%         Opt.Prefix = prefix;
%         Opt.NoCheck = 0;
%         
%         fprintf('Writing %s as BRIK and HEAD files to %s\n',prefix,pname);
%         curdir = cd(pname);
%         
%         %
%         % Compensating for an AFNI inconsistency
%         %
%         if isfield(Info,'WARP_TYPE'),
%             Info = rmfield(Info,'WARP_TYPE');
%         end
%         
%         WriteBrik(vol,Info,Opt);  % needs to have a modified version to overwrite
%         % an existing file with "w+"
%         
%         %%%% ADDED BY BER 4/5/13 for alignment purposes.
%         % copy data to new file so original dimensions are saved.
%         copyDataCMD = ['3dcopy  ' prefix ' ' prefix '_refit'];
%         disp(copyDataCMD);
%         system(copyDataCMD);
%         
%         % refit the x y z dimensioned based on strange AFNI conventions
%         refitDataCMD = ['3drefit -xdel 1.5 -ydel 1.5 -zdel 1.5 ' prefix '_refit+orig'];
%         disp(refitDataCMD);
%         system(refitDataCMD);
%         
%         monk=DSP.afni_prefix(1:3);
%         % now that data is refit, align the center with the surface anatomy center
%         alignDataCMD = ['@Align_Centers -base ' STDPATH.afni_uw '/Anatomy/' monk '_t1_refit_shft+orig. -dset ' prefix '_refit+orig.'];
%         disp(alignDataCMD);
%         system(alignDataCMD);
%         
%         % finally shift the origin slightly because of strange ANFI stuff
%         shiftDataCMD = ['3drefit -dyorigin -0.5 -dzorigin 0.5 ' prefix '_refit_shft+orig.'];
%         disp(shiftDataCMD);
%         system(shiftDataCMD);
%         %%%%
%         
%         % Copy relevant files to where the surface files are
%         [s,mess,messID]=movefile([pname, '*_refit_shft+orig*'],  dirSPEC);
%         if ~s
%             fprintf(1,' %s \n', mess);
%         end
%     end
% end
% 

