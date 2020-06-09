 function [matR_SU, matP_SU, paramCorr] = computeCorrMap(nameSubjNeural, nameSubjBOLD, setMovie, flagSaveFile, typeMION)
% computeCorrMap.m
%
% Compute Spearman's rank correlation between each neuron and each voxel
% 
% SYNTAX: [matR_SU, matP_SU, paramCorr] = computeCorrMap(nameSubjNeural, nameSubjBOLD, flagSaveFile)
% OUTPUT
%     - matR_SU: matrix of correlation coefficients (voxel x neuron)
%     - matP_SU: matrix of p values provided by "corr" function
%     - paramCorr: some parameters. IDs of movies, IDs and indices of single unit channels 
% INPUT:
%     - nameSubjNeural: Name of the subject for neural data (Tor, Rho, Sig)
%     - nameSubjBOLD: Name of the subject for fMRI data (Art, Ava)
%     - setMovie: Vector of movie ID (e.g. [1 2 3])
%     - flagSaveFile: if 1, it will save file
%
% 2017/11/06 SHP: modified to have set of movie IDs as an input
% 2016/04/05 SHP: modified to do HRF convolution to neural response time course to each
% movie before across-movie-concatenation (so far, convolution was done to
% a long 15-min timecourse obtained by concatnenation)
% 2017/02/27 SHP: add to test different HRF kernels (two more shapes)
%               new input "typeMION" can be one of [1 2 3]: 1 for our
%               original gamma pdf, 2 for kernel from Leite et al. (2002),
%               3 for kernel from Silva et al. (2007)

%% Compute correlation for movie 123
% 
addpath('/library/matlab_utils/')

% nameSubjNeural = 'Sig'; %'Rho'; %'Tor';
% nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';

% Load data files
dirDataHome = '/procdata/parksh/_macaque';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
% fileNameNeural_BLP = [nameSubjNeural, '_movieTS_BLPLFP_indMov.mat'];
filenameBOLD = [nameSubjBOLD,  '_movieTS_fMRI_indMov.mat']; % '_movieTS_fMRI_indMov.mat'];

fprintf(1, '\nLoading single unit data of %s: %s ....', nameSubjNeural, filenameNeural)
load(fullfile(dirDataNeural, filenameNeural))
% load(fullfile(dirDataNeural, fileNameNeural_BLP))
fprintf(1, '\nLoading fMRI data of %s: %s ....\n', nameSubjBOLD, filenameBOLD)
load(fullfile(dirDataBOLD, filenameBOLD))

% Get movie IDs common in two dataset
commonSetMovie = intersect(paramBOLD.unimov, paramSDF.setMovIDs);
dataBOLD.mvoltc = voltcIndMov(commonSetMovie);
dataBOLD.unimov = commonSetMovie;

% setMovie = [4 5 6]; %[1 2 3 4 5 6]; %[4 5 6]; %[1 2 3];
paramCorr.setMovie = setMovie;

% 1. fMRI tc
fmritc=[];
indMovieBOLD = find(ismember(dataBOLD.unimov, setMovie)>0);
for iM = indMovieBOLD %1:length(indMovieBOLD)
    curvoltc = dataBOLD.mvoltc{iM};
    avgvoltc = repmat(nanmean(curvoltc,4),[1 1 1 size(curvoltc,4)]);
    if ~isempty(find(avgvoltc==0, 1))
        avgvoltc(avgvoltc==0) = realmin; % get rid of zeros because it causes NaNs in percent signals
    end
    pcvoltc = ((curvoltc - avgvoltc)./avgvoltc)*100;
    fmritc = cat(4,fmritc,pcvoltc);
end

% 2. neural regressor
[indDataMat, CellID, movieID] = genDataMatrix_SU(nameSubjNeural, 0); % data matrix
validC = find(indDataMat*ismember(movieID, setMovie)>0); % valid channel with the movie set
indMovieNeuron = find(ismember(paramSDF.setMovIDs, setMovie)>0);
% MION function

switch typeMION
    case 1
        % 1. gamma pdf
        TR=2.4;
        x = -40:TR:40;
        k = gampdf(x, 4, 2);
    case 2
        % % 2. kernel from AFNI
        TR=2.4;
        x = -40:TR:40;
        taxis = x(x>=0); %0:TR:40; %50;
        k = 16.4486 * ( -0.184/ 1.5 * exp(-taxis/ 1.5)...
            +0.330/ 4.5 * exp(-taxis/ 4.5)...
            +0.670/13.5 * exp(-taxis/13.5) );
        k = cat(2, zeros(1, length(k )), k );
    case 3        
        % 3. kernel from Silva et al. (2007) rat somatosensory CBV kernel using two
        % gamma function
        TR=2.4;
        x = -40:TR:40;
        k = gampdf(x, 2, 1) + gampdf(x, 2, 10).*2;
end

% Reshape BOLD 4-d data
[nx, ny, nz, nt] = size(fmritc);
nVox = nx*ny*nz;

matR_SU = NaN(nVox, length(validC)); matP_SU = NaN(nVox, length(validC));
for iChan = 1:length(validC) % compute correlation channel-by-channel
    
    %     neuralrgrs = cat(1, S(validC(iChan), indMovieNeuron).mnFR);
%     neuralrgrs = neuralrgrs-mean(neuralrgrs); % centering
%     neuralrgrs = doConv(neuralrgrs,k); % convolve MION kernel %conv(neuralrgrs,k,'same');

    % Modified 2016/04/05, 2016/04/27 by SHP
    neuralrgrs=[];
    for iMov = 1:length(indMovieNeuron)
        curNeuralTC = S(validC(iChan), indMovieNeuron(iMov)).mnFR(8:125); %S(validC(iChan), indMovieNeuron(iMov)).mnFR
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
    fprintf(1, 'Channel: %s \n', S(validC(iChan), indMovieNeuron(iMov)).cellID)
        
    [Rvals, Pvals] = corr(reshape(fmritc, nVox, nt)', neuralrgrs',...
        'rows','complete', 'type', 'Spearman');
    
    matR_SU(:,iChan) = Rvals.*(-1); % because of MION
    matP_SU(:,iChan) = Pvals;
    % mapR = reshape(Rvals, [nx, ny, nz]).*(-1); % because of MION
    % mapP = reshape(Pvals, [nx, ny, nz]);
    
end

paramCorr.validChanIndex = validC;
[catChanID{1:length(validC)}] = deal(S(validC, indMovieNeuron(1)).cellID);
paramCorr.validChanID = char(catChanID); 
% paramCorr.validChanID = cat(1, S(validC,indMovieNeuron(1)).cellID);

tempS = num2str(setMovie);
MovieStr = tempS(~isspace(tempS));

if flagSaveFile
    save(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie%s_new.mat', nameSubjNeural, nameSubjBOLD, MovieStr)), 'matR_SU', 'matP_SU', 'paramCorr');
    fprintf(1, '\n ...Results are saved in %s as %s \n', dirDataNeural, sprintf('CorrMap_SU_%s%sMovie%s_new.mat', nameSubjNeural, nameSubjBOLD, MovieStr))

%     save(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new_kernel%d.mat', nameSubjNeural, nameSubjBOLD, typeMION)), 'matR_SU', 'matP_SU', 'paramCorr');
%     fprintf(1, '\n ...Results are saved in %s as %s \n', dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new_kernel%d.mat', nameSubjNeural, nameSubjBOLD, typeMION))
end



