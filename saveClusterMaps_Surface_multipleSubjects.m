% saveClusterMaps_Surface_multipleSubjects.m
%
% 2017/01/22 modified 'saveClusterMaps_Surface.m'
% Clustering results from Toroid and Spice cells together
% 1) load fMRI correlation maps computed based on Toroid's cells and Artemis
% or Avalanche's fMRI data
% 2) load clustering results based on various values (e.g. whole brain corr maps, corr between cells and movie
% regressors, cell spike density function etc.)
% 3) for each cluster, average the corr maps
% 4) save it as AFNI brik/head format so that it can be mapped onto the surface anatomy

clear all;


nameSubjNeural = 'Spi'; % 'Tor'; % 'Sig'; %'Tor';
nameSubjBOLD = 'Art'; %'Ava'; %'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';

%
addpath('/library/matlab_utils/')

% Load files
dirDataHome = '/procdata/parksh/';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% % 1) fMRI correlation maps
% load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'matR_SU', 'paramCorr')

% 1) Clustering results
load(fullfile(dirDataNeural, 'Clustering_SU_allCells_5ROIs_8Means_prob.mat')) % temporary results
% load(fullfile(dirDataNeural, sprintf('Clustering_TorSpi%sMovie123_new_masked.mat', nameSubjBOLD)), 'Clustering_moviemask', 'paramClustering*')

% 2) fMRI masks (movie-driven activity & brain-only mask)
load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp'); %, 'brainMask_BlockAna3D');

% 3) Load the corr map and select valid channel: subject by subject
% Set directories 
setNameSubjNeural = {'Tor', 'Spi'};
numSubject = size(setNameSubjNeural, 2);
matR_SU_all = [];

for iSubj = 1:numSubject
    nameSubjNeural = setNameSubjNeural{iSubj}; %'Spi'; %'Tor';
    dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
    load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'matR_SU');
    
    validChanIndex_clustering = paramClustering(iSubj).validChanIndex;
    
    matR_SU_valid = matR_SU(:, validChanIndex_clustering);
    matR_SU_all = cat(2, matR_SU_all, matR_SU_valid);
    
     clear matR_SU matR_SU_valid
    
%     % 2017/01/19: Excluded some units, in case of Spice's data which hasn't been went through this procedure yet
%     switch lower(nameSubjNeural)
%         case 'spi'
%             excChanIndex = [10 13 22 27 30 49]; % cells were not same acrossd two days
%             validChanIndex_clustering = setdiff(paramCorr.validChanIndex, excChanIndex);
%             validChanID_clustering = cat(2, paramCorr.validChanID(validChanIndex_clustering,:), num2str(ones(size(validChanIndex_clustering)).*iSubj)); %paramCorr.validChanID(validChanIndex_clustering,:);
%         otherwise
%             validChanIndex_clustering = (1:length(paramCorr.validChanIndex))';
%             validChanID_clustering = cat(2, paramCorr.validChanID, num2str(ones(size(paramCorr.validChanIndex)).*iSubj)); %paramCorr.validChanID;
%     end
    
%     paramClustering(iSubj).nameSubj = nameSubjNeural;
%     paramClustering(iSubj).validChanIndex = validChanIndex_clustering;
%     paramClustering(iSubj).validChanID = validChanID_clustering;
        
end

%% 1. Average fMRI maps for each cluster
matIndClust_SU = cat(2, Clustering_moviemask.resultKMeans.SU_indCluster);
% matIndClust_SU = cat(2, Clustering.resultKMeans.SU_indCluster);

for sortTargetK = 2:12 %3:8 %5; %6:8 %4; %:5 %sortTargetK = 4; %8; %6; %7; %6; %7;
[sortedClust, indSortChan]=sort(matIndClust_SU(:,sortTargetK-1));

% average maps for each cluster
nx = 40; ny = 64; nz = 32;
mapR_Cluster=[]; R_avgCluster=[];


for iK = 1:sortTargetK
    R_avgCluster(:,iK) = mean(matR_SU_all(:,indSortChan(sortedClust==iK)), 2);
    tempClustMapR = reshape(R_avgCluster(:,iK), [nx, ny, nz]);
    mapR_Cluster = cat(4, mapR_Cluster, tempClustMapR);
end

% % compute residual
% GRAND average map
% mapR_grandAvg = reshape(mean(matR_SU,2), [nx, ny, nz]);
% R_grandAvg = mean(matR_SU,2);
% R_avgCluster_residual = R_avgCluster - repmat(R_grandAvg, 1, sortTargetK);
% mapR_Cluster_residual = [];
% for iK = 1:sortTargetK
%     tempClustMapR_residual = reshape(R_avgCluster_residual(:,iK), [nx, ny, nz]);
%     mapR_Cluster_residual = cat(4, mapR_Cluster_residual, tempClustMapR_residual);
% end




pname = [dirDataBOLD, '/tempSURF/'];
dirSPEC = [dirDataBOLD, '/Anatomy/_suma/TorSpiClusters'];

cd /projects/parksh/NeuralBOLD/analysis/BlockAna/
blockana;
S_neuralRegressor(nameSubjNeural, nameSubjBOLD);
global STDPATH DSP DATA GH

% convert the map to the surface
for iMask = 1:2
    
    for iK = 1:sortTargetK
        
        cd /projects/parksh/NeuralBOLD/analysis/BlockAna/
        
        switch iMask
            case 1 % unmasked (everything within brain)
                fileHead = sprintf('new_%s', cell2mat(setNameSubjNeural));
                DSP.proc.fncvol_3d = mapR_Cluster(:,:,:,iK).*brainMask_BlockAna3D;
            case 2 % masked
                fileHead = sprintf('new_masked_%s', cell2mat(setNameSubjNeural));
                DSP.proc.fncvol_3d = mapR_Cluster(:,:,:,iK).*movieDrivenAmp.mask_amp1; %reshape(mapR, [nx, ny, nz]).*movieDrivenAmp.mask_amp1;
        end
        
%         fname = sprintf('%s_Cluster%d_%dMeans_Art_AvgCorrMapMovie123_noFiltering+orig.BRIK', fileHead, iK, sortTargetK);%
        fname = sprintf('%s_Cluster%d_%dMeansMovieDriven_Art_AvgCorrMapMovie123_noFiltering+orig.BRIK', fileHead, iK, sortTargetK);%

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
end
end

%% 2. Average time series in each cluster and then compute the map based on the averaged TS
% matIndClust_SU = cat(2, Clustering.resultKMeans.SU_indCluster);
% curK = 7; 
% [sortedClust, indSortChan]=sort(matIndClust_SU(:,curK-1));
% 
% setMovie = [1 2 3];
% validC = paramCorr.validChanIndex; %find(indDataMat*ismember(movieID, setMovie)>0); % valid channel with movie [1 2 3]
% 
% % MION function
% TR=2.4;
% k = gampdf([-40:TR:40],4,2);
% 
% matFR_TR = zeros(375, length(validC));
% for iChan = 1:length(validC) 
% 
%     % Modified 2016/04/05, 2016/04/27 by SHP
%     neuralrgrs=[];
%     for iMov = 1:length(setMovie)
%         curNeuralTC = S(validC(iChan), setMovie(iMov)).mnFR(8:125); %S(validC(iChan), indMovieNeuron(iMov)).mnFR
%         curNeuralTC = curNeuralTC-mean(curNeuralTC); % centering
%         curNeuralTC = doConv(curNeuralTC,k); % convolve MION kernel %conv(neuralrgrs,k,'same');
%         curNeuralTC = cat(2, NaN(1,7), curNeuralTC); %curNeuralTC(1:7) = NaN;
%         
%         neuralrgrs = cat(2, neuralrgrs, curNeuralTC); % concatenation across movies
%     end
%     
%     matFR_TR(:,iChan) = neuralrgrs';
% end
% 
% % Average SDF in each cluster across cells
% matFRCluster=[];steFRCluster=[];
% for iK = 1:curK
%     tempMatFR=[];
%     tempMatFR = matFR(:,indSortChan(sortedClust==iK));
%     matFRCluster(:,iK) = mean(tempMatFR, 2);
%     steFRCluster(:,iK) = std(tempMatFR, 0, 2)./length(indSortChan(sortedClust==indClust));
% end
% 
% % Compute the corr map
% commonSetMovie = intersect(paramBOLD.unimov, paramSDF.setMovIDs);
% dataBOLD.mvoltc = voltcIndMov(commonSetMovie);
% dataBOLD.unimov = commonSetMovie;
% 
% setMovie = [1 2 3];
% paramCorr.setMovie = setMovie;
% 
% % fMRI tc
% fmritc=[];
% indMovieBOLD = find(ismember(dataBOLD.unimov, setMovie)>0);
% for iM = indMovieBOLD %1:length(indMovieBOLD)
%     curvoltc = dataBOLD.mvoltc{iM};
%     avgvoltc = repmat(nanmean(curvoltc,4),[1 1 1 size(curvoltc,4)]);
%     if ~isempty(find(avgvoltc==0, 1))
%         avgvoltc(avgvoltc==0) = realmin; % get rid of zeros because it causes NaNs in percent signals
%     end
%     pcvoltc = ((curvoltc - avgvoltc)./avgvoltc)*100;
%     fmritc = cat(4,fmritc,pcvoltc);
% end
% 
% % Reshape BOLD 4-d data
% [nx, ny, nz, nt] = size(fmritc);
% nVox = nx*ny*nz;
% 
% matR_avgTS = NaN(nVox, curK);
% for iK = 1:curK % compute correlation channel-by-channel
% 
%     neuralrgrs = matFRCluster(:,iK);
%     [Rvals, Pvals] = corr(reshape(fmritc, nVox, nt)', neuralrgrs,...
%         'rows','complete', 'type', 'Spearman');
%     
%     matR_avgTS(:,iK) = Rvals.*(-1); % because of MION        
% end
% 
% 
% pname = [dirDataBOLD, '/tempSURF/'];
% dirSPEC = [dirDataBOLD, '/Anatomy/_suma/'];
% 
% cd /projects/parksh/NeuralBOLD/analysis/BlockAna/
% blockana;
% S_neuralRegressor
% % global STDPATH DSP DATA GH
% 
% % convert the map to the surface
% for iMask = 1:2
%     
%     for iK = 1:curK
%         
%         cd /projects/parksh/NeuralBOLD/analysis/BlockAna/
%         
%         mapR_avgTS = reshape(matR_avgTS(:,iK), [nx, ny, nz]); 
%         switch iMask
%             case 1 % unmasked (everything within brain)
%                 fileHead = 'new';
%                 DSP.proc.fncvol_3d = mapR_avgTS.*brainMask_BlockAna3D;
%             case 2 % masked
%                 fileHead = 'new_masked';
%                 DSP.proc.fncvol_3d = mapR_avgTS.*movieDrivenAmp.mask_amp1; %reshape(mapR, [nx, ny, nz]).*movieDrivenAmp.mask_amp1;
%         end
%         
%         fname = sprintf('%s_Cluster%dAvgTS_%dMeans_Art_CorrMapMovie123_noFiltering+orig.BRIK', fileHead, iK, sortTargetK);%
%         
%         vol = single(DSP.proc.fncvol_3d);
%         
%         %     fname = sprintf('%s_%s_Movie123_noFiltering_%s+orig.BRIK', fileHead, nameSubjNeural, fileTail);%
%         %
%         %     vol = single(DSP.proc.fncvol_3d);  %single(mapR_Cluster(:,:,:,iK)); % vol = single(DSP.proc.fncvol_3d);
%         % %     cellID = validChanID(iUnit,:);
%         % %     fprintf(1, 'Unit # %d, Cell ID: %s, %s \n', iUnit, cellID);       
%         
%         
%         % Do dumpFunctionalBrik
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


%%
%%%%
% %% fMRI map for each SU
% % average maps for each cluster
% nx = 40; ny = 64; nz = 32;
%
%
% pname = [dirDataBOLD, '/tempSURF/'];
% dirSPEC = [dirDataBOLD, '/Anatomy/_suma/'];
%
% cd /projects/parksh/NeuralBOLD/analysis/BlockAna/
% blockana;
% S_neuralRegressor
% global STDPATH DSP DATA GH
%
% % convert the map to the surface
% for iUnit = 1:size(matR_SU,2)
%
%     cd /projects/parksh/NeuralBOLD/analysis/BlockAna/
%     cellID = paramCorr.validChanID(iUnit,:);
%     fprintf(1, 'Unit # %d, Cell ID: %s, %s \n', iUnit, cellID);
%
%     fname = sprintf('new_masked_%s_%s_Movie123_noFiltering+orig.BRIK', cellID, nameSubjNeural);%
%
%     DSP.proc.fncvol_3d = reshape(matR_SU(:,iUnit), [nx, ny, nz]).*movieDrivenAmp.mask_amp1;
%     vol = single(DSP.proc.fncvol_3d);  %single(mapR_Cluster(:,:,:,iK)); % vol = single(DSP.proc.fncvol_3d);
%
%
%     %% Do dumpFunctionalBrik
%     Info = DATA.Info;
%     Info.DATASET_RANK = DATA.Info.DATASET_RANK;
%     Info.DATASET_DIMENSIONS = DATA.Info.DATASET_DIMENSIONS;
%     Info.DELTA = DATA.Info.DELTA;
%     Info.BYTEORDER_STRING = DATA.Info.BYTEORDER_STRING;
%     Info.TYPESTRING = DATA.Info.TYPESTRING;
%     Info.SCENE_DATA = DATA.Info.SCENE_DATA;
%     Info.ORIGIN = DATA.Info.ORIGIN;
%     Info.TAXIS_FLOATS = DATA.Info.TAXIS_FLOATS;
%     Info.TAXIS_OFFSETS = DATA.Info.TAXIS_OFFSETS;
%     Info.IDCODE_STRING = DATA.Info.IDCODE_STRING;
%     Info.IDCODE_DATE = DATA.Info.IDCODE_DATE;
%     Info.HISTORY_NOTE = DATA.Info.HISTORY_NOTE;
%
%     Info.DATASET_DIMENSIONS = [size(vol) 0 0];
%     Info.DATASET_RANK(2) = 1;
%     Info.BRICK_TYPES = 3;  %3 = float
%     Info.BRICK_FLOAT_FACS = 0;
%     Info.BRICK_STATS = [min(vol(:)) max(vol(:))];
%     Info.BRICK_LABS = DATA.Info.BRICK_LABS(1:2);
%     Info.BRICK_KEYWORDS = DATA.Info.BRICK_KEYWORDS;
%     Info.TAXIS_NUMS(1) = 1;
%     Info.TypeName = 'float';
%     Info.TypeBytes = 4;
%     Info.ByteOrder = 'ieee-le';
%     Info.Orientation = DATA.Info.Orientation;
%     Info.FileFormat = DATA.Info.FileFormat;
%     Info.RootName = DATA.Info.RootName;
%     Info.Extension_1D = DATA.Info.Extension_1D;
%
%     code0 = getAFNIOrientationCode('LR');
%     code1 = getAFNIOrientationCode('AP');
%     code2 = getAFNIOrientationCode('DV');
%
%     Info.ORIENT_SPECIFIC = [code0 code1 code2];
%
%     doti = strfind(fname,'+');
%     prefix = fname(1:doti-1);
%     Opt.Prefix = prefix;
%     Opt.NoCheck = 0;
%
%     fprintf('Writing %s as BRIK and HEAD files to %s\n',prefix,pname);
%     curdir = cd(pname);
%
%     %
%     % Compensating for an AFNI inconsistency
%     %
%     if isfield(Info,'WARP_TYPE'),
%         Info = rmfield(Info,'WARP_TYPE');
%     end
%
%     WriteBrik(vol,Info,Opt);  % needs to have a modified version to overwrite
%     % an existing file with "w+"
%
%     %%%% ADDED BY BER 4/5/13 for alignment purposes.
%     % copy data to new file so original dimensions are saved.
%     copyDataCMD = ['3dcopy  ' prefix ' ' prefix '_refit'];
%     disp(copyDataCMD);
%     system(copyDataCMD);
%
%     % refit the x y z dimensioned based on strange AFNI conventions
%     refitDataCMD = ['3drefit -xdel 1.5 -ydel 1.5 -zdel 1.5 ' prefix '_refit+orig'];
%     disp(refitDataCMD);
%     system(refitDataCMD);
%
%     monk=DSP.afni_prefix(1:3);
%     % now that data is refit, align the center with the surface anatomy center
%     alignDataCMD = ['@Align_Centers -base ' STDPATH.afni_uw '/Anatomy/' monk '_t1_refit_shft+orig. -dset ' prefix '_refit+orig.'];
%     disp(alignDataCMD);
%     system(alignDataCMD);
%
%     % finally shift the origin slightly because of strange ANFI stuff
%     shiftDataCMD = ['3drefit -dyorigin -0.5 -dzorigin 0.5 ' prefix '_refit_shft+orig.'];
%     disp(shiftDataCMD);
%     system(shiftDataCMD);
%     %%%%
%
%     % Copy relevant files to where the surface files are
%     [s,mess,messID]=movefile([pname, '*_refit_shft+orig*'],  dirSPEC);
%     if ~s
%         fprintf(1,' %s \n', mess);
%     end
% end


