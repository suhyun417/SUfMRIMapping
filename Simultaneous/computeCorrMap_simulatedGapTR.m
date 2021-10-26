% computeCorrMap_simulatedGapTR.m
%
% script to compute a whole-brain map based on a simulated 
% single unit time course that contains "gap" period during the half of the TR 
% 2017/10/03 SHP
%

%% 
addpath('/library/matlab_utils/')

nameSubjNeural = 'Tor'; %'Sig'; %'Rho'; %'Tor';
nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';

%
dirDataHome = '/procdata/parksh/';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

%% make new TS for neurons with new TR
TR_sec_org = 2.4;
TR_sec_new = 4; %3.5;
gap_sec = TR_sec_new/2;
setMovie = [1 2 3];

% 2. neural regressor
[indDataMat, CellID, movieID] = genDataMatrix_SU(nameSubjNeural, 0); % data matrix
validC = find(indDataMat*ismember(movieID, setMovie)>0); % valid channel with movie [1 2 3]
indMovieNeuron = find(ismember(paramSDF.setMovIDs, setMovie)>0);

setCellIDs = CellID(validC); % should be cell array of string
FR_dT = createCellRegressor_indMov_discreteTime(dirDataNeural, setCellIDs, setMovie, gap_sec); % ,


% addpath('/library/matlab_utils/')
% 
% nameSubjNeural = 'Tor'; %'Sig'; %'Rho'; %'Tor';
% nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
% 
% % Load data files
% dirDataHome = '/procdata/parksh/';
% dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
% dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
% fileNameNeural_BLP = [nameSubjNeural, '_movieTS_BLPLFP_indMov.mat'];
filenameBOLD = [nameSubjBOLD, '_movieTS_fMRI_indMov.mat'];

fprintf(1, '\nLoading single unit data of %s: %s ....', nameSubjNeural, filenameNeural)
load(fullfile(dirDataNeural, filenameNeural))
% load(fullfile(dirDataNeural, fileNameNeural_BLP))
fprintf(1, '\nLoading fMRI data of %s: %s ....\n', nameSubjBOLD, filenameBOLD)
load(fullfile(dirDataBOLD, filenameBOLD))

% Get movie IDs common in two dataset
commonSetMovie = intersect(paramBOLD.unimov, paramSDF.setMovIDs);
dataBOLD.mvoltc = voltcIndMov(commonSetMovie);
dataBOLD.unimov = commonSetMovie;

setMovie = [1 2 3];
paramCorr.setMovie = setMovie;

% 1. fMRI tc
fmritc_new=[];
indMovieBOLD = find(ismember(dataBOLD.unimov, setMovie)>0);

for iM = indMovieBOLD %1:length(indMovieBOLD)
    curvoltc = dataBOLD.mvoltc{iM};
    avgvoltc = repmat(nanmean(curvoltc,4),[1 1 1 size(curvoltc,4)]);
    if ~isempty(find(avgvoltc==0, 1))
        avgvoltc(avgvoltc==0) = realmin; % get rid of zeros because it causes NaNs in percent signals
    end
    pcvoltc = ((curvoltc - avgvoltc)./avgvoltc)*100;
    
    % Reshape BOLD 4-d data to resample
    [nx, ny, nz, nt] = size(pcvoltc);
    nVox = nx*ny*nz;
    
    curpcvoltc_mat = reshape(pcvoltc, nVox, nt)';
    curpcvoltc_mat_new = resample(curpcvoltc_mat, TR_sec_org*10, TR_sec_new*10); % to match the virtual sampling rate
    
    fmritc_new = cat(1,fmritc_new,curpcvoltc_mat_new);
end

% 2. neural regressor
[indDataMat, CellID, movieID] = genDataMatrix_SU(nameSubjNeural, 0); % data matrix
validC = find(indDataMat*ismember(movieID, setMovie)>0); % valid channel with movie [1 2 3]
indMovieNeuron = find(ismember(paramSDF.setMovIDs, setMovie)>0);
% MION function
typeMION = 1;
switch typeMION
    case 1
        % 1. gamma pdf
        TR=TR_sec_new; %2.4;
        x = -40:TR:40;
        k = gampdf(x, 4, 2);
    case 2
        % % 2. kernel from AFNI
        TR=TR_sec_new; %2.4;
        x = -40:TR:40;
        taxis = x(x>=0); %0:TR:40; %50;
        k = 16.4486 * ( -0.184/ 1.5 * exp(-taxis/ 1.5)...
            +0.330/ 4.5 * exp(-taxis/ 4.5)...
            +0.670/13.5 * exp(-taxis/13.5) );
        k = cat(2, zeros(1, length(k )), k );
    case 3        
        % 3. kernel from Silva et al. (2007) rat somatosensory CBV kernel using two
        % gamma function
        TR=TR_sec_new; %2.4;
        x = -40:TR:40;
        k = gampdf(x, 2, 1) + gampdf(x, 2, 10).*2;
end


matR_SU = NaN(nVox, length(validC)); matP_SU = NaN(nVox, length(validC));
for iChan = 1:length(validC) % compute correlation channel-by-channel
    
    % neuronal tc
    neuralrgrs=[];
    for iMov = 1:length(indMovieNeuron)
        tempTC = FR_dT(iChan, indMovieNeuron(iMov)).mnFR(1:2:end); % (8:125); %
        curNeuralTC = tempTC(15:75); % after NaNs: hard coding for the number of NaNs after resampling based on the new (arbitrary) TR %(8:125); %S(validC(iChan), indMovieNeuron(iMov)).mnFR
        curNeuralTC = curNeuralTC-mean(curNeuralTC); % centering
        curNeuralTC = doConv(curNeuralTC,k); % convolve MION kernel %conv(neuralrgrs,k,'same');
        curNeuralTC = cat(2, NaN(1,14), curNeuralTC); %cat(2, NaN(1,7), curNeuralTC); %curNeuralTC(1:7) = NaN;
        
        neuralrgrs = cat(2, neuralrgrs, curNeuralTC); % concatenation across movies
    end
    
%     % fMRI tc    
%     fmritc_mat = reshape(fmritc, nVox, nt)';
%     fmritc_mat_new = resample(fmritc_mat, TR_sec_org*10, TR_sec_new*10); % to match the virtual sampling rate

%     % Modified 2016/04/05, 2016/04/27 by SHP
%     neuralrgrs=[];
%     for iMov = 1:length(indMovieNeuron)
%         curNeuralTC = S(validC(iChan), indMovieNeuron(iMov)).mnFR(8:125); %S(validC(iChan), indMovieNeuron(iMov)).mnFR
%         curNeuralTC = curNeuralTC-mean(curNeuralTC); % centering
%         curNeuralTC = doConv(curNeuralTC,k); % convolve MION kernel %conv(neuralrgrs,k,'same');
%         curNeuralTC = cat(2, NaN(1,7), curNeuralTC); %curNeuralTC(1:7) = NaN;
%         
%         neuralrgrs = cat(2, neuralrgrs, curNeuralTC); % concatenation across movies
%     end

    
    % Compute correlation    
    fprintf(1, 'Compute correlation using regressors:: \n'); % Put some slider while corr is being calculated..
    fprintf(1, 'Channel: %s \n', S(validC(iChan),1).cellID)
    
        
    [Rvals, Pvals] = corr(fmritc_new, neuralrgrs', ... %reshape(fmritc, nVox, nt)', neuralrgrs',...
        'rows','complete', 'type', 'Spearman');
    
    matR_SU(:,iChan) = Rvals.*(-1); % because of MION
    matP_SU(:,iChan) = Pvals;
    % mapR = reshape(Rvals, [nx, ny, nz]).*(-1); % because of MION
    % mapP = reshape(Pvals, [nx, ny, nz]);
    
end

paramCorr.validChanIndex = validC;
paramCorr.validChanID = cat(1, S(validC,1).cellID);

save(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_simulationGapTR%dsec.mat', nameSubjNeural, nameSubjBOLD, TR_sec_new)), 'matR_SU', 'matP_SU', 'paramCorr');


% if flagSaveFile
%     save(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new_kernel%d.mat', nameSubjNeural, nameSubjBOLD, typeMION)), 'matR_SU', 'matP_SU', 'paramCorr');
%     fprintf(1, '\n ...Results are saved in %s as %s \n', dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new_kernel%d.mat', nameSubjNeural, nameSubjBOLD, typeMION))
% end



