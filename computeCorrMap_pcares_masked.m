function [matR_SU, matP_SU, paramCorr] = computeCorrMap_pcares_masked(nameSubjNeural, nameSubjBOLD, setMovie, flagSaveFile)
% computeCorrMap_pcares.m
%
% Compute Spearman's rank correlation between each neuron and each voxel
% using residual of fMRI data (1st Principal Component is regressed out)
% 2018/10/27 SHP: for now, setMovie should be within [1 2 3]

%% Compute correlation for movie 123
% 
addpath('/library/matlab_utils/')

typeMION = 1;
% setMovie = 1; %[1 2 3];
% nameSubjBOLD = 'Art';

paramCorr.setMovie = setMovie;

% Load data files
dirDataHome = '/procdata/parksh/_macaque';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
% if sum(strcmpi(nameSubjNeural, {'spice', 'spi'}))
%     dirDataNeural = '/procdata/parksh/Spi/2018Jan_movie/';
% end

filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
load(fullfile(dirDataNeural, filenameNeural))
load('/procdata/parksh/_macaque/Art/Art_movieTS_fMRI_Movie123_PCA.mat', 'resultsPCAres_moviemask')


% 1. fMRI tc
catRes = cat(2, resultsPCAres_moviemask(setMovie).residuals)';
fmritc = [];
for iM = 1:length(setMovie) %3
    fmritc = cat(1, fmritc, NaN(7, size(catRes, 2)), catRes(118*(iM-1)+1:118*iM, :));
end
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

% % Reshape BOLD 4-d data
% % [nx, ny, nz, nt] = size(fmritc);
% % nVox = nx*ny*nz;
nVox = size(fmritc, 2);

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
        
    [Rvals, Pvals] = corr(fmritc, neuralrgrs',... %corr(reshape(fmritc, nVox, nt)', neuralrgrs',...
        'rows','complete', 'type', 'Spearman');
    
    matR_SU(:,iChan) = Rvals.*(-1); % because of MION
    matP_SU(:,iChan) = Pvals;
    % mapR = reshape(Rvals, [nx, ny, nz]).*(-1); % because of MION
    % mapP = reshape(Pvals, [nx, ny, nz]);
    
end

paramCorr.validChanIndex = validC;
[catChanID{1:length(validC)}] = deal(S(validC, indMovieNeuron(1)).cellID);
paramCorr.validChanID = char(catChanID); 
%paramCorr.validChanID = cat(1, S(validC,indMovieNeuron(1)).cellID);

tempS = num2str(setMovie);
MovieStr = tempS(~isspace(tempS));

if flagSaveFile
    save(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie%s_pcares_masked.mat', nameSubjNeural, nameSubjBOLD, MovieStr)), 'matR_SU', 'matP_SU', 'paramCorr');
    fprintf(1, '\n ...Results are saved in %s as %s \n', dirDataNeural, sprintf('CorrMap_SU_%s%sMovie%s_pcares_masked.mat', nameSubjNeural, nameSubjBOLD, MovieStr))

%     save(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new_kernel%d.mat', nameSubjNeural, nameSubjBOLD, typeMION)), 'matR_SU', 'matP_SU', 'paramCorr');
%     fprintf(1, '\n ...Results are saved in %s as %s \n', dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new_kernel%d.mat', nameSubjNeural, nameSubjBOLD, typeMION))
end



