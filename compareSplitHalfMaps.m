% compareSplitHalfMaps.m
%
% 2017/02/27 SHP
% Do a split-half analysis on maps. 
%     
%     1) load fMRI time series. split them in random half. then take the mean.
%     2) load SU time series. split them in random half. then take the mean.
%     3) compute correlations between each of those. 2x2 matrix.
%     4) compare 4 maps by computing correlation between them.

%% Setting
addpath('/library/matlab_utils/')

nameSubjNeural = 'Sig'; %'Rho'; %'Spi'; % 'Tor'; %'Sig'; %'Rho'; %'Tor';
nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';

% Load data files
dirDataHome = '/procdata/parksh/';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
% fileNameNeural_BLP = [nameSubjNeural, '_movieTS_BLPLFP_indMov.mat'];
filenameBOLD = [nameSubjBOLD, '_movieTS_fMRI_indMov.mat'];

load(fullfile(dirDataNeural, filenameNeural))
load(fullfile(dirDataBOLD, filenameBOLD))

%% 1) fMRI time series
% % go and run BlockAna to deal with fMRI time series
% cd /projects/parksh/NeuralBOLD/analysis/BlockAna/
% blockana;
% 
% % 1-1). Load the data
% % Get information from the filelist in .txt
% switch lower(nameSubjBOLD)
%     case 'art'
%         sessionFileList = 'FL_e66_allfiles2.txt';
%     case 'ava'
%         sessionFileList = 'FL_ava_allfiles2.txt';
%     case 'sig'
%         sessionFileList = 'FL_sig_allfiles.txt';
% end
% skip=7; % number of TRs to skip
% SI = SU_createSessInfo(sessionFileList,[],skip); 
% filelist = SU_makeFileList(SI.monkID,SI.sessID,SI.scanID);
% 
% MAX_TR = SI.max_tr;
% skip=SI.skiptr;
% 
% unimov = [1 2 3]; %sort(unique(SI.movID));
% nunimov = length(unimov);
% catimgdat = [];
% rgr = zeros(1,nunimov*MAX_TR);
% procmovs = [];
% 
% for u = 1:nunimov
%     movindxs = find(SI.movID == unimov(u));
%     nmovindxs = length(movindxs);
%     
%     % collect all the data files for this movie
%     filelist = SU_makeFileList(SI.monkID(movindxs),SI.sessID(movindxs),...
%         SI.scanID(movindxs));
%     
%     totfiles = length(filelist)
%     
%     for f=1:totfiles
%         s_sub = filelist{f}.subj;
%         s_ses = filelist{f}.sess;
%         s_sc  = filelist{f}.scan;
%         
%         [fmri_tc,dgz,notes] = SU_loadFile(s_sub,s_ses,s_sc);        
%         if f==1
%             [a,b,c,t] = size(fmri_tc);
%             allvoltc = zeros(a,b,c,MAX_TR,totfiles);
%         end        
%         scanlen = size(fmri_tc,4);
%         if scanlen == 250
%             mov_start_tr = 63;% offset when movie starts in scan
%         elseif scanlen == 125
%             mov_start_tr = 1; % short movie
%         else
%             fprintf('WARNING: bad scan length %d\n',scanlen);
%         end
%         
%         val_tr         = [mov_start_tr:(mov_start_tr+MAX_TR-1)];
%         [clp_fmri_tc,clp_mdgz] = SU_clipMovDat(fmri_tc,dgz,val_tr);
%         allvoltc(:,:,:,:,f) = clp_fmri_tc;
%     end
%     
%     % make the onset response to NaNs
%     allvoltc(:,:,:,1:skip,:) = NaN;
%     
%     % 2-2) & 2-3). take the mean of random half at this point, for each movie
%     if mod(totfiles, 2)>0
%         totfiles = totfiles-1; % make it even
%     end
%     randOrder = reshape(randperm(totfiles), 2, totfiles/2);
%     dataBOLD.mvoltc{1, u} = nanmean(allvoltc(:,:,:,:,randOrder(1,:)), 5); % first half, for each movie
%     dataBOLD.mvoltc{2, u} = nanmean(allvoltc(:,:,:,:,randOrder(2,:)), 5); % second half, for each movie
% end
% 
% % 1. fMRI tc
% fmritc = cell([1 2]);
% indMovieBOLD = [1 2 3]; % find(ismember(dataBOLD.unimov, setMovie)>0);
% for iHalf = 1:2
%     for iM = indMovieBOLD %1:length(indMovieBOLD)
%         curvoltc = dataBOLD.mvoltc{iHalf, iM};
%         avgvoltc = repmat(nanmean(curvoltc,4),[1 1 1 size(curvoltc,4)]);
%         if ~isempty(find(avgvoltc==0, 1))
%             avgvoltc(avgvoltc==0) = realmin; % get rid of zeros because it causes NaNs in percent signals
%         end
%         pcvoltc = ((curvoltc - avgvoltc)./avgvoltc)*100;
%         fmritc{iHalf} = cat(4,fmritc{iHalf},pcvoltc);
%     end
% end
% clear curvoltc avgvoltc pcvoltc dataBOLD fmri_tc
% clear allvoltc clp_fmri_tc clp_mdgz dgz alldgz aa
    
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

%% 2) Single unit time series & 3) Compute correlation
% MION function
typeMION = 1;
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

%
setMovie=[1 2 3];
switch lower(nameSubjNeural)
    case 'spi'
        excChanIndex = [10 13 22 27 30 49]; % cells were not same acrossd two days
        [indDataMat, CellID, movieID] = genDataMatrix_SU(nameSubjNeural, 0);
        validC = find(indDataMat*ismember(movieID, setMovie)>0); %
        validC = setdiff(validC, excChanIndex);
    otherwise
        [indDataMat, CellID, movieID] = genDataMatrix_SU(nameSubjNeural, 0); % data matrix
        validC = find(indDataMat*ismember(movieID, setMovie)>0); % valid channel with movie [1 2 3]
end
length(validC)
indMovieNeuron = find(ismember(paramSDF.setMovIDs, setMovie)>0);
clear splitHalf
for iChan = 1:length(validC) 
    fprintf(1, 'Channel: %d\n', iChan)
    
    % Modified 2016/04/05, 2016/04/27 by SHP
    neuralrgrs = cell([1 2]);
    for iMov = 1:length(indMovieNeuron)
        
        matFR = S(validC(iChan), indMovieNeuron(iMov)).matFR{1};
        nTrial = size(matFR, 2);
        
        if mod(nTrial, 2)>0
            nTrial = nTrial-1; % make it even
        end
        randOrder=[];
        randOrder = reshape(randperm(nTrial), 2, nTrial/2);
        
        for iHalf = 1:2
            curNeuralTC = mean(matFR(8:125, randOrder(iHalf, :)), 2);
            curNeuralTC = curNeuralTC-mean(curNeuralTC); % centering
            curNeuralTC = doConv(curNeuralTC,k); % convolve MION kernel %conv(neuralrgrs,k,'same');
            curNeuralTC = cat(2, NaN(1,7), curNeuralTC); %curNeuralTC(1:7) = NaN;
 
            neuralrgrs{iHalf} = cat(2, neuralrgrs{iHalf}, curNeuralTC); % concatenation across movies
        end
    end
    dataSU(iChan).neuralrgrs = neuralrgrs;

% 3) compute correlation for 2 x 2 combinations of cases

%     for iF = 1:2

        for iS = 1:2
            [Rvals, Pvals] = corr(reshape(fmritc, nVox, nt)', dataSU(iChan).neuralrgrs{iS}',...
                'rows','complete', 'type', 'Spearman');
            
            splitHalf(iS).matR_SU(:,iChan) = Rvals.*(-1); % because of MION
        end
%     end
end

save(fullfile(dirDataNeural, sprintf('CorrMap_%s%s_Movie123_splithalf.mat', nameSubjNeural, nameSubjBOLD)),  'splitHalf')

%     
% save(fullfile(dirDataNeural, 'data4splitHalfCorrComp.mat'), 'fmritc', 'dataSU')

% dirDataHome = '/data/parks20/procdata/NeuroMRI/Tor/'; %biowulf
% load(fullfile(dirDataHome, 'data4splitHalfCorrComp.mat'))
% 
% 
% [nx, ny, nz, nt] = size(fmritc{1});
% nVox = nx*ny*nz;


%% 4) Comparison between maps
% valid voxels
setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi'};
nameSubjBOLD = 'Art';
load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp');
[nx ny nz] = size(movieDrivenAmp.mask_amp1);
nVox = nx*ny*nz;
moviemask_vec = reshape(movieDrivenAmp.mask_amp1, nVox, 1);

matR_valid_1 = [];
matR_valid_2 = [];
for iS = 1:length(setNameSubjNeural)
    nameSubjNeural = setNameSubjNeural{iS};
    dirDataNeural = fullfile('/procdata/parksh/', nameSubjNeural);
    
    load(fullfile(dirDataNeural, sprintf('CorrMap_%s%s_Movie123_splithalf.mat', nameSubjNeural, nameSubjBOLD)))
    
    matR_valid_1 = cat(2, matR_valid_1, splitHalf(1).matR_SU(moviemask_vec, :));
    matR_valid_2 = cat(2, matR_valid_2, splitHalf(2).matR_SU(moviemask_vec, :));
end


dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';

rValSplitHalf=[];
 for iChan = 1:size(matR_valid_1, 2)
     [r, p] = corr(matR_valid_1(:,iChan), matR_valid_2(:,iChan), 'rows', 'complete', 'type', 'Spearman');
     rValSplitHalf(iChan, 1) = r;
 end

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [200 200 375 315])
hist(rValSplitHalf)
medR = median(rValSplitHalf);
hold on
line([medR medR], get(gca, 'YLim'), 'Color', 'r', 'LineWidth', 2)
text(medR, 15, 'median  \rightarrow ', 'HorizontalAlignment', 'Right', 'Color', 'r')
title('Similarity between maps from random half trials')
% title(sprintf('Similarity between maps from random half trials: %s', upper(nameSubjNeural)))
xlabel('Correlation')
ylabel('Frequency')

% difference measure?
diffMatR = matR_valid_1 - matR_valid_2;

% for iK=1:3
%     figure;
%     set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 575 545])
%     imagesc(matR_acrossSUMaps{iK})
%     set(gca, 'CLim', [-1 1])
%     title(sprintf('Correlation between single unit maps: kernel %d', iK))
%     print(gcf, fullfile(dirFig, sprintf('HRFComparison_CorrMatrix_kernel%d', iK)), '-dtiff', '-r120');
% end




