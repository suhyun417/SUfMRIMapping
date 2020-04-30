% saveAvgMap_Surface.m
%
% 1) Compute the correlation between the average neural response across
% units and fMRI voxels
% 4) Save it as AFNI brik/head format so that it can be mapped onto the surface anatomy

setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi'};
nameSubjBOLD = 'Art'; %'Ava';

%% Compute correlation for movie 123
%
addpath('/library/matlab_utils/')

% Load data files
dirDataHome = '/procdata/parksh/';
% dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
filenameBOLD = [nameSubjBOLD, '_movieTS_fMRI_indMov.mat'];

% fprintf(1, '\nLoading single unit data of %s: %s ....', nameSubjNeural, filenameNeural)
% load(fullfile(dirDataNeural, filenameNeural))
fprintf(1, '\nLoading fMRI data of %s: %s ....\n', nameSubjBOLD, filenameBOLD)
load(fullfile(dirDataBOLD, filenameBOLD))

% Get movie IDs common in two dataset
% commonSetMovie = intersect(paramBOLD.unimov, paramSDF.setMovIDs);
% dataBOLD.mvoltc = voltcIndMov(commonSetMovie);
% dataBOLD.unimov = commonSetMovie;

setMovie = [1 2 3];
paramCorr.setMovie = setMovie;
dataBOLD.mvoltc = voltcIndMov(setMovie);
dataBOLD.unimov = setMovie;

% 1. fMRI tc
fmritc=[];
% indMovieBOLD = find(ismember(dataBOLD.unimov, setMovie)>0);
for iM = setMovie %indMovieBOLD %1:length(indMovieBOLD)
    curvoltc = dataBOLD.mvoltc{iM};
    avgvoltc = repmat(nanmean(curvoltc,4),[1 1 1 size(curvoltc,4)]);
    if ~isempty(find(avgvoltc==0, 1))
        avgvoltc(avgvoltc==0) = realmin; % get rid of zeros because it causes NaNs in percent signals
    end
    pcvoltc = ((curvoltc - avgvoltc)./avgvoltc)*100;
    fmritc = cat(4,fmritc,pcvoltc);
end

% 2. neural regressor: averaged SDF across units
numSubject = size(setNameSubjNeural, 2);
matTS_all = [];

for iSubj = 1:numSubject
    nameSubjNeural = setNameSubjNeural{iSubj}; %'Spi'; %'Tor';
    dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
    filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
    load(fullfile(dirDataNeural, filenameNeural))

    [indDataMat, CellID, movieID] = genDataMatrix_SU(nameSubjNeural, 0); % data matrix
    validC = find(indDataMat*ismember(movieID, setMovie)>0); % valid channel with movie [1 2 3]
    switch lower(nameSubjNeural)
        case 'spi'
            excChanIndex = [10 13 22 27 30 49]; % cells were not same acrossd two days
            validC = setdiff(validC, excChanIndex);
        otherwise
            validC = validC;
    end
    
    indMovieNeuron = find(ismember(paramSDF.setMovIDs, setMovie)>0);
    % MION function
    % % 1. gamma pdf
    TR=2.4;
    k = gampdf([-40:TR:40],4,2);
    % % 2. kernel from AFNI
    % taxis = 0:2.4:50;
    % k = 16.4486 * ( -0.184/ 1.5 * exp(-taxis/ 1.5)...
    % +0.330/ 4.5 * exp(-taxis/ 4.5)...
    % +0.670/13.5 * exp(-taxis/13.5) );
    % k = cat(2, zeros(1, length(k )-1), k );
    
    matTS = NaN(length(validC), 375);
    for iChan = 1:length(validC) % compute correlation channel-by-channel
        
        % Modified 2016/04/05, 2016/04/27 by SHP
        neuralrgrs=[];
        for iMov = 1:length(indMovieNeuron)
            curNeuralTC = S(validC(iChan), indMovieNeuron(iMov)).mnFR(8:125); %S(validC(iChan), indMovieNeuron(iMov)).mnFR
            curNeuralTC = curNeuralTC-mean(curNeuralTC); % centering
            curNeuralTC = doConv(curNeuralTC,k); % convolve MION kernel %conv(neuralrgrs,k,'same');
            curNeuralTC = cat(2, NaN(1,7), curNeuralTC); %curNeuralTC(1:7) = NaN;
            
            neuralrgrs = cat(2, neuralrgrs, curNeuralTC); % concatenation across movies
        end
        
        matTS(iChan,:) = neuralrgrs;
    end
    
    matTS_all = cat(1, matTS_all, matTS);
end

avgTS.meanTS = nanmean(matTS_all);
avgTS.medianTS = nanmedian(matTS_all);
avgTS.steTS = nanstd(matTS_all)./sqrt(size(matTS_all,1));

% Reshape BOLD 4-d data
[nx, ny, nz, nt] = size(fmritc);
nVox = nx*ny*nz;

matR_meanSDF = NaN(nVox, 1);
matR_medianSDF = NaN(nVox, 1);
% Compute correlation
[matR_meanSDF] = corr(reshape(fmritc, nVox, nt)', avgTS.meanTS',...
    'rows','complete', 'type', 'Spearman');
[matR_medianSDF] = corr(reshape(fmritc, nVox, nt)', avgTS.medianTS',...
    'rows','complete', 'type', 'Spearman');


avgCorrMap.matR_meanSDF = matR_meanSDF.*(-1); % because of MION
avgCorrMap.matR_medianSDF = matR_medianSDF.*(-1); % because of MION

paramCorr.validChanIndex = validC;
paramCorr.validChanID = cat(1, S(validC,1).cellID);
paramCorr.MIONkernel = k;

%
if flagSaveFile
    save(fullfile(dirDataNeural, sprintf('CorrMap_SUavg_%s%sMovie123_new.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)), 'avgTS', 'avgCorrMap', 'paramCorr');
    fprintf(1, '\n ...Results are saved in %s as %s \n', dirDataNeural, sprintf('CorrMap_SUavg_%s%sMovie123_new.mat', cell2mat(setNameSubjNeural), nameSubjBOLD))
end

%% fMRI map for each SU
%  fMRI movie-driven activity mask
load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp');

% average maps for each cluster
nx = 40; ny = 64; nz = 32;


pname = [dirDataBOLD, '/tempSURF/'];
dirSPEC = [dirDataBOLD, '/Anatomy/_suma/'];

cd /projects/parksh/NeuralBOLD/analysis/BlockAna/
blockana;
S_neuralRegressor
global STDPATH DSP DATA GH

% convert the map to the surface
for iMap = 1:2    
    
    switch iMap
        case 1 % mean SDF
            mapR = avgCorrMap.matR_meanSDF;
            fileHead = 'avgSDF_mean';
        case 2 % median SDF
            mapR = avgCorrMap.matR_medianSDF;
            fileHead = 'avgSDF_median';
    end
    
    for iMask = 1:2
        
        cd /projects/parksh/NeuralBOLD/analysis/BlockAna/
        
        switch iMask
            case 1 % unmasked
                fileTail = 'new';
                DSP.proc.fncvol_3d = reshape(mapR, [nx, ny, nz]);
            case 2 % masked
                fileTail = 'new_masked';
                DSP.proc.fncvol_3d = reshape(mapR, [nx, ny, nz]).*movieDrivenAmp.mask_amp1;
        end
        
        fname = sprintf('%s_%s_Movie123_noFiltering_%s+orig.BRIK', fileHead, cell2mat(setNameSubjNeural), fileTail);%
        
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
    
    
    
    
    %
    % for iChan=1:size(S,1)
    % mapR_catMov{iChan} = cat(4, S(iChan,:).mapR);
    % meanR_mov{iChan} = mean(mapR_catMov{iChan}, 4);
    % end
    %
    % mapR_catAll = cat(4, meanR_mov{:});
    % meanR_all = mean(mapR_catAll, 4);
    %
    % DSP.proc.fracvarmap_3d = meanR_all; %S(iChan).mapR.*-1; %S(6).mapR.*-1;
    %
    %
    % for iChan=1:size(S,1)
    % mapR_catMov_rect{iChan} = abs(cat(4, S(iChan,:).mapR));
    % meanR_mov_rect{iChan} = mean(mapR_catMov_rect{iChan}, 4);
    % end
    %
    % mapR_catAll_rect = cat(4, meanR_mov_rect{:});
    % meanR_all_rect = mean(mapR_catAll_rect, 4);