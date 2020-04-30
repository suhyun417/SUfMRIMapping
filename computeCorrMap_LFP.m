 function [matR_LFP, matP_LFP, paramCorr] = computeCorrMap_LFP(nameSubjNeural, nameSubjBOLD, flagSaveFile)
% computeCorrMap_LFP.m
%
% Compute Spearman's rank correlation between each fq BLP and each voxel
% 
% SYNTAX: [matR_SU, matP_SU, paramCorr] = computeCorrMap(nameSubjNeural, nameSubjBOLD, flagSaveFile)
% OUTPUT
%     - matR_SU: matrix of correlation coefficients (voxel x BLP)
%     - matP_SU: matrix of p values provided by "corr" function
%     - paramCorr: some parameters. IDs of movies, IDs and indices of single unit channels 
% INPUT:
%     - nameSubjNeural: Name of the subject for neural data (Tor, Rho, Sig)
%     - nameSubjBOLD: Name of the subject for fMRI data (Art, Ava)
%     - flagSaveFile: if 1, it will save file
%
% 2016/04/05 SHP: modified to do HRF convolution to neural response time course to each
% movie before across-movie-concatenation (so far, convolution was done to
% a long 15-min timecourse obtained by concatnenation)
%

%% Compute correlation for movie 123
% 
addpath('/library/matlab_utils/')

% nameSubjNeural = 'Sig'; %'Rho'; %'Tor';
% nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';

% Load data files
dirDataHome = '/procdata/parksh/';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
filenameNeural_BLP = [nameSubjNeural, '_movieTS_BLPLFP_indMov.mat'];
filenameBOLD = [nameSubjBOLD, '_movieTS_fMRI_indMov.mat'];

% fprintf(1, '\nLoading single unit data of %s: %s ....\n', nameSubjNeural, filenameNeural)
% load(fullfile(dirDataNeural, filenameNeural))
fprintf(1, '\nLoading LFP data of %s: %s ....\n', nameSubjNeural, filenameNeural_BLP)
load(fullfile(dirDataNeural, filenameNeural_BLP))
fprintf(1, '\nLoading fMRI data of %s: %s ....\n', nameSubjBOLD, filenameBOLD)
load(fullfile(dirDataBOLD, filenameBOLD))

% Get movie IDs common in two dataset
% commonSetMovie = intersect(paramBOLD.unimov, paramSDF.setMovIDs);
% dataBOLD.mvoltc = voltcIndMov(commonSetMovie);
% dataBOLD.unimov = commonSetMovie;

setMovie = [1 2 3];
paramCorr.setMovie = setMovie;

% 1. fMRI tc
fmritc=[];
% indMovieBOLD = find(ismember(dataBOLD.unimov, setMovie)>0);
for iM = 1:length(setMovie)% indMovieBOLD %1:length(indMovieBOLD)
    curM = setMovie(iM);
    curvoltc = voltcIndMov{curM};
    avgvoltc = repmat(nanmean(curvoltc,4),[1 1 1 size(curvoltc,4)]);
    if ~isempty(find(avgvoltc==0, 1))
        avgvoltc(avgvoltc==0) = realmin; % get rid of zeros because it causes NaNs in percent signals
    end
    pcvoltc = ((curvoltc - avgvoltc)./avgvoltc)*100;
    fmritc = cat(4,fmritc,pcvoltc);
end

% 2. neural regressor
% [indDataMat, CellID, movieID] = genDataMatrix_SU(nameSubjNeural, 0); % data matrix
% validC = find(indDataMat*ismember(movieID, setMovie)>0); % valid channel with movie [1 2 3]
% indMovieNeuron = find(ismember(paramSDF.setMovIDs, setMovie)>0);
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

% Reshape BOLD 4-d data
[nx, ny, nz, nt] = size(fmritc);
nVox = nx*ny*nz;

matR_LFP = NaN(nVox, length(BLPRGR(1).meanBLP)); matP_LFP = NaN(nVox, length(BLPRGR(1).meanBLP));
for iF = 1:length(BLPRGR(1).meanBLP) % compute correlation frequency-by-frequency

    %     neuralrgrs = cat(1, S(validC(iChan), indMovieNeuron).mnFR);
%     neuralrgrs = neuralrgrs-mean(neuralrgrs); % centering
%     neuralrgrs = doConv(neuralrgrs,k); % convolve MION kernel %conv(neuralrgrs,k,'same');

    % Modified 2016/04/05, 2016/04/27 by SHP
    neuralrgrs=[];
    for iMov = 1:length(setMovie)
        curNeuralTC = BLPRGR(setMovie(iMov)).meanBLP{iF}(8:125); % S(validC(iChan), indMovieNeuron(iMov)).mnFR(8:125); %S(validC(iChan), indMovieNeuron(iMov)).mnFR
        curNeuralTC = curNeuralTC-mean(curNeuralTC); % centering
        curNeuralTC = doConv(curNeuralTC,k); % convolve MION kernel %conv(neuralrgrs,k,'same');
        curNeuralTC = cat(2, NaN(1,7), curNeuralTC); %curNeuralTC(1:7) = NaN;
        
        neuralrgrs = cat(2, neuralrgrs, curNeuralTC); % concatenation across movies
    end
   
        
    % matSDF=[];
    % for iMov = indMovieNeuron %1:length(indMovieNeuron)
    % %     movID = setMovie(iMov);
    %     tempMat=[];
    %     tempMat = cat(2,S(validC,iMov).mnFR); % time x cell
    %     matSDF = cat(1, matSDF, tempMat); % concatenate across movies
    % end
    %     iBP = get(CorrMapHandles.BLPtype, 'Value')-1;
    %     catBLP = cat(1, BLPRGR.meanBLP);
    
    %     neuralrgrs = cat(1,catBLP{selMovID_neural,iBP}); % column vector
    
    % Compute correlation    
    fprintf(1, 'Compute correlation using regressors:: \n'); % Put some slider while corr is being calculated..
    fprintf(1, 'Frequency: %d - %d \n', BLPRGR(1).blpfreq(iF,1), BLPRGR(1).blpfreq(iF,2))
        
    [Rvals, Pvals] = corr(reshape(fmritc, nVox, nt)', neuralrgrs',...
        'rows','complete', 'type', 'Spearman');
    
    matR_LFP(:,iF) = Rvals.*(-1); % because of MION
    matP_LFP(:,iF) = Pvals;
    % mapR = reshape(Rvals, [nx, ny, nz]).*(-1); % because of MION
    % mapP = reshape(Pvals, [nx, ny, nz]);
    
end

paramCorr.blpFreq = BLPRGR(1).blpfreq;
% paramCorr.validChanIndex = validC;
% paramCorr.validChanID = cat(1, S(validC,1).cellID);

if flagSaveFile
    save(fullfile(dirDataNeural, sprintf('CorrMap_LFP_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'matR_LFP', 'matP_LFP', 'paramCorr');
    fprintf(1, '\n ...Results are saved in %s as %s \n', dirDataNeural, sprintf('CorrMap_LFP_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD))
end



