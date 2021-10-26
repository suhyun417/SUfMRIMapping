function [] = computeCorrcoeffCI(nameSubjNeural, nameSubjBOLD, iChan, idVox_start, flagRawSDFShuffle)
% compute confidence interval of correlation values using bootstrap
%
% 2017/02/25 SHP: added the option for shuffling raw SDF before MION
% convolution
%
% 1. generate bootstrapped time series using stationary bootstrap (Politis
% and Romano, 1994)
% 2. compute correlation for bootstrapped sample
% 3. get the confidence interval from the simulated correlation values

%% Settings
% % Add necessary toolbox % Should be done outside of the function in
% advance because of the compile step of this code for Biowulf
% % addpath('/library/matlab_utils/')
% addpath('/projects/parksh/_toolbox/Boot_Time_Series')

% Set directories
dirDataHome =  '/data/parks20/procdata/NeuroMRI'; % Biowulf '/procdata/parksh/'; %nifstorage %'/data/parks20/procdata/NeuroMRI'; % Biowulf
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);
dirSaveFile = fullfile(dirDataNeural, 'corrMap_resultsCI_10hz');
if ~exist(dirSaveFile, 'dir')
    mkdir(dirSaveFile)
end

% % Directory for saving figures as graphic files
% dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';


%% Load the data
filenameNeural = sprintf('neuraltc4computeCorrcoeffCI_%s.mat', nameSubjNeural);
if flagRawSDFShuffle
    filenameNeural = sprintf('neuraltc4computeCorrcoeffCI_%s_10hz.mat', nameSubjNeural);
end
filenameBOLD = sprintf('fmritc4computeCorrcoeffCI_%s.mat', nameSubjBOLD);

load(fullfile(dirDataNeural, filenameNeural))
load(fullfile(dirDataBOLD, filenameBOLD))


%%
% Take care of possible string-double issue from compile
if ischar(iChan); eval(sprintf('iChan = %s;', iChan)); end;
if ischar(idVox_start); eval(sprintf('idVox_start = %s;', idVox_start)); end;

nVox = 1000;
if idVox_start == 81001
    nVox = 920;
end
setVoxID = idVox_start:1:idVox_start+nVox-1;

saveFileName = sprintf('corrcoeffCI_%s%s_fixedBlock_cell%d_vox%d.mat', nameSubjNeural, nameSubjBOLD, iChan, idVox_start);

%% Generate bootstrap samples
% For each movie, use stationary bootstrap with mean size block proportional to 11TR

% Bootstrap parameters
nSample = 10000;
typeBS = 1; % input for "stationaryBB" function: 1 for stationary geometric pdf, 2 for stationary uniform pdf, 3 for circular (non-random blocks)
% for typeBS==3, a fixed block size will be used
avgBlockSize = 11; % MION kernel size in TR (first try just the actual MION kernel size / next, should try 4 times of this value according to Mudelsee (2003))
nOnsetTR = 7; % number of TR excluded beforehands because of the onset time
nTotalTRPerMovie = 125;
avgBlockSize_fineResSDF = 10; % for 100ms-bin (10Hz) resolution SDF, use 1-sec window
nTotalPerMovie_fineResSDF = 3000;
nOnset_fineResSDF = nTotalPerMovie_fineResSDF*(nOnsetTR/nTotalTRPerMovie);


% Matrix allocation (for boostrapped samples)
% matSU_bootstrap = NaN(nTotalTRPerMovie-nOnsetTR, length(validC), nSample);
% matBOLD_bootstrap = NaN(nTotalTRPerMovie-nOnsetTR, nVox, nSample);
matRho_bootstrap = NaN(nSample, nVox);
VoxCI = struct();
% matCI95 = NaN(nVox, 2);
% matCI99 = NaN(nVox, 2);
% matRho_org = NaN(nVox, 1);


if flagRawSDFShuffle % Shuffle the fine res SDF in 10Hz (100ms-bin) resolution, then convolve the MION kernel
    % Original neural time course
    ntc_org=[];
    ntc_org = neuraltc(nOnset_fineResSDF+1:nTotalPerMovie_fineResSDF, :, iChan); % time x movie (118 x 3)
    for iVox = 1:nVox
        
        idVox = setVoxID(iVox);
        
        % original voxel time course for 3 movies
        fmritc_org = [];
        fmritc_org = fmritc(nOnsetTR+1:nTotalTRPerMovie, :, idVox); % time x movie (118 x 3)
        
        matRho_bs = NaN(nSample,1);
        % Bootstrapping
        for iBS = 1:nSample
            
            % Sample the fine res SDF
            ntc_rawbs =[];
            ntc_rawbs = stationaryBB(ntc_org, typeBS, avgBlockSize_fineResSDF);
            
            % convolve MION kernel
            TR=2.4;
            k = gampdf([-40:TR:40],4,2);
            ntc_bs=[];
            ntc_rawbs = ntc_rawbs-repmat(mean(ntc_rawbs), size(ntc_rawbs, 1), 1);
            ntc_bs = doConv(ntc_rawbs, k);
            ntc_bs = ntc_bs(:, nOnsetTR+1:nTotalTRPerMovie)';
            
%             % Sample the voxel response
%             fmritc_bs =[];
%             fmritc_bs = stationaryBB(fmritc_org, typeBS, avgBlockSize);
            
            % Compute correlation for the current voxel original TS x
            % resampled single unit TS
            rho_bs = corr(fmritc_org(:), ntc_bs(:), 'rows','complete', 'type', 'Spearman');
            
            matRho_bs(iBS, 1) = rho_bs;
            
        end
        
        % Get the CI
        sortRho = sort(matRho_bs, 1, 'ascend');
        CI95 = [sortRho(round(nSample*0.025)), sortRho(nSample - round(nSample*0.025))];
        CI99 = [sortRho(round(nSample*0.005)), sortRho(nSample - round(nSample*0.005))];
        VoxCI(iVox).matCI95 = CI95; %matCI95(iVox, :) = CI95;
        VoxCI(iVox).matCI99 = CI99; %matCI99(iVox, :) = CI99;
        
        %     % Original correlation
        %     rho_org = corr(fmritc_org(:), ntc_org(:), 'rows','complete', 'type', 'Spearman');
        %     VoxCI(iVox).rho_org = rho_org; %matRho_org(iVox,1) = rho_org;
        %
        %     matRho_bootstrap(:,iVox) = matRho_bs;
        
        % Print screen
        fprintf(1, 'Vox #%d is done\n', iVox);
        
        % Save
        save(fullfile(dirSaveFile, saveFileName), 'VoxCI', 'matRho_bootstrap');
    end
    
else % bootstrapping of MION-convolved neural regressors
    % Original neural time course
    ntc_org=[];
    ntc_org = neuraltc(nOnsetTR+1:nTotalTRPerMovie, :, iChan); % time x movie (118 x 3)
    
    for iVox = 1:nVox
        
        idVox = setVoxID(iVox);
        
        % original voxel time course for 3 movies
        fmritc_org = [];
        fmritc_org = fmritc(nOnsetTR+1:nTotalTRPerMovie, :, idVox); % time x movie (118 x 3)
        
        matRho_bs = NaN(nSample,1);
        % Bootstrapping
        for iBS = 1:nSample
            
            % Sample the neural response
            ntc_bs =[];
            ntc_bs = stationaryBB(ntc_org, typeBS, avgBlockSize);
            
            % Sample the voxel response
            fmritc_bs =[];
            fmritc_bs = stationaryBB(fmritc_org, typeBS, avgBlockSize);
            
            % Compute correlation for the current voxel x single unit
            rho_bs = corr(fmritc_bs(:), ntc_bs(:), 'rows','complete', 'type', 'Spearman');
            
            matRho_bs(iBS, 1) = rho_bs;
            
        end
        
        % Get the CI
        sortRho = sort(matRho_bs, 1, 'ascend');
        CI95 = [sortRho(round(nSample*0.025)), sortRho(nSample - round(nSample*0.025))];
        CI99 = [sortRho(round(nSample*0.005)), sortRho(nSample - round(nSample*0.005))];
        VoxCI(iVox).matCI95 = CI95; %matCI95(iVox, :) = CI95;
        VoxCI(iVox).matCI99 = CI99; %matCI99(iVox, :) = CI99;
        
        % Original correlation
        rho_org = corr(fmritc_org(:), ntc_org(:), 'rows','complete', 'type', 'Spearman');
        VoxCI(iVox).rho_org = rho_org; %matRho_org(iVox,1) = rho_org;
        
        matRho_bootstrap(:,iVox) = matRho_bs;
        
        % Print screen
        fprintf(1, 'Vox #%d is done\n', iVox);
        
        % Save
        save(fullfile(dirSaveFile, saveFileName), 'VoxCI', 'matRho_bootstrap');
        
        %     flagSig = 0;
        %     if rho_org < CI95(1) || rho_org > CI95(2)
        %         flagSig = 1;
        %     end
        
        %     % Save the results
        %     resultsStationaryBS(iChan).cellID = S(validC(iChan)).cellID;
        %     resultsStationaryBS(iChan).matCI95 = matCI95;
        %     resultsStationaryBS(iChan).matCI99 = matCI99;
        %     resultsStationaryBS(iChan).matRho_org = matRho_org;
    end
    
end


%% Save
save(fullfile(dirSaveFile, saveFileName), 'VoxCI', 'matRho_bootstrap');





%         matR_SU(:,iChan) = Rvals.*(-1); % because of MION
%         matP_SU(:,iChan) = Pvals;
% mapR = reshape(Rvals, [nx, ny, nz]).*(-1); % because of MION
% mapP = reshape(Pvals, [nx, ny, nz]);

%     % fMRI needs to be done through a for-loop for movie-by-movie bootstrap
%     for iM = 1:3
%         fmritc_org = []; fmritc_bs = [];
%         fmritc_org = fmritc(nOnsetTR+1:nTotalTRPerMovie, :, iM);
%         % Get rid of NaNs at the beginning for bootstrap
%         fmritc_bs = stationaryBB(fmritc_org, typeBS, avgBlockSize);


% %% Prepare original time series for individual movie
% % Get movie IDs common in two dataset
% setMovie = [1 2 3];
% indMovieBOLD = find(ismember(paramBOLD.unimov, setMovie)>0);
%
% % 1. fMRI tc in percent signal
% [nx, ny, nz, nt] = size(voltcIndMov{1});
% nVox = nx*ny*nz;
%
% fmritc = [];
% for iM = 1:length(indMovieBOLD)
%     curvoltc =reshape(voltcIndMov{indMovieBOLD(iM)}, nVox, nt)'; %  voltcIndMov{iM};
%     avgvoltc = repmat(nanmean(curvoltc),[size(curvoltc,1), 1]); %repmat(nanmean(curvoltc,4),[1 1 1 size(curvoltc,4)]);
%     if ~isempty(find(avgvoltc==0, 1))
%         avgvoltc(avgvoltc==0) = realmin; % get rid of zeros because it causes NaNs in percent signals
%     end
%     pcvoltc = ((curvoltc - avgvoltc)./avgvoltc)*100;
%     fmritc(:,iM,:) = pcvoltc; % time x movie x voxel
% %     fmritc = cat(3, fmritc, pcvoltc);
% end
% clear voltcIndMov curvoltc avgvoltc pcvoltc
%
% % 2. neural regressor
% [indDataMat, CellID, movieID] = genDataMatrix_SU(nameSubjNeural, 0); % data matrix
% validC = find(indDataMat*ismember(movieID, setMovie)>0); % valid channel with movie [1 2 3]
% indMovieNeuron = find(ismember(paramSDF.setMovIDs, setMovie)>0);
% % MION function (gamma pdf)
% TR=2.4;
% k = gampdf([-40:TR:40],4,2);
%
% neuraltc=[];
% for iChan = 1:length(validC) % for each channel
%     clear tempn neuralrgrs
%     [tempn{1:3}] = deal(S(validC(iChan),indMovieNeuron).mnFR);
%     neuralrgrs = cell2mat(tempn);
%     neuralrgrs = neuralrgrs-repmat(mean(neuralrgrs), size(neuralrgrs,1), 1);
%     neuralrgrs = doConv(neuralrgrs, k);
%
%     neuraltc = cat(3, neuraltc, neuralrgrs');
% end





