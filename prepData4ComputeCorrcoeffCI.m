function []=prepData4ComputeCorrcoeffCI(nameSubjNeural, nameSubjBOLD) %, flagLocal)

%% Settings
ss = pwd;
if ~isempty(strfind(ss, 'Volume')) % if it's local
    dirProjects = '/Volumes/PROJECTS';
    dirProcdata = '/Volumes/PROCDATA';
    dirLibrary = '/Volumes/LIBRARY';
else % on virtual machine
    dirProjects = '/projects';
    dirProcdata = '/procdata';
    dirLibrary = '/library';
end
    
% Add necessary toolbox % Should be 
addpath(fullfile(dirLibrary, 'matlab_utils'))
addpath(fullfile(dirProjects, 'parksh/_toolbox/Boot_Time_Series'))

% Set directories 
% nameSubjNeural = 'Tor';
% nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataHome = fullfile(dirProcdata, 'parksh');
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% Directory for saving figures as graphic files
dirFig = fullfile(dirProjects, 'parksh/NeuralBOLD/_labNote/_figs/');


%% Load the original time series
filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
% fileNameNeural_BLP = [nameSubjNeural, '_movieTS_BLPLFP_indMov.mat'];
filenameBOLD = [nameSubjBOLD, '_movieTS_fMRI_indMov.mat'];

fprintf(1, '\nLoading single unit data of %s: %s ....', nameSubjNeural, filenameNeural)
load(fullfile(dirDataNeural, filenameNeural))
% load(fullfile(dirDataNeural, fileNameNeural_BLP))
fprintf(1, '\nLoading fMRI data of %s: %s ....\n', nameSubjBOLD, filenameBOLD)
load(fullfile(dirDataBOLD, filenameBOLD))


%% Prepare original time series for individual movie
% Get movie IDs common in two dataset
setMovie = [1 2 3];
indMovieBOLD = find(ismember(paramBOLD.unimov, setMovie)>0); 

% 1. fMRI tc in percent signal
[nx, ny, nz, nt] = size(voltcIndMov{1});
nVox = nx*ny*nz;

fmritc = [];
for iM = 1:length(indMovieBOLD)    
    curvoltc =reshape(voltcIndMov{indMovieBOLD(iM)}, nVox, nt)'; %  voltcIndMov{iM}; 
    avgvoltc = repmat(nanmean(curvoltc),[size(curvoltc,1), 1]); %repmat(nanmean(curvoltc,4),[1 1 1 size(curvoltc,4)]);
    if ~isempty(find(avgvoltc==0, 1))
        avgvoltc(avgvoltc==0) = realmin; % get rid of zeros because it causes NaNs in percent signals
    end
    pcvoltc = ((curvoltc - avgvoltc)./avgvoltc)*100;
    fmritc(:,iM,:) = pcvoltc; % time x movie x voxel
%     fmritc = cat(3, fmritc, pcvoltc);
end
clear voltcIndMov curvoltc avgvoltc pcvoltc

% saveFileNameBOLD = sprintf('fmritc4computeCorrcoeffCI_%s.mat', nameSubjBOLD)
% save(fullfile(dirDataBOLD, saveFileNameBOLD), 'fmritc')

% 2. neural regressor
[indDataMat, CellID, movieID] = genDataMatrix_SU(nameSubjNeural, 0); % data matrix
validC = find(indDataMat*ismember(movieID, setMovie)>0); % valid channel with movie [1 2 3]
indMovieNeuron = find(ismember(paramSDF.setMovIDs, setMovie)>0);
% MION function (gamma pdf)
TR=2.4;
k = gampdf([-40:TR:40],4,2);

neuraltc=[];
for iChan = 1:length(validC) % for each channel
    clear tempn neuralrgrs
    [tempn{1:3}] = deal(S(validC(iChan),indMovieNeuron).mnFR);
    neuralrgrs = cell2mat(tempn);
    neuralrgrs = neuralrgrs-repmat(mean(neuralrgrs), size(neuralrgrs,1), 1);
    neuralrgrs = doConv(neuralrgrs, k);
    
    neuraltc = cat(3, neuraltc, neuralrgrs');
end

% saveFileNameNeural = sprintf('neuraltc4computeCorrcoeffCI_%s.mat', nameSubjNeural)
% save(fullfile(dirDataNeural, saveFileNameNeural), 'neuraltc')

% 2-1. Neuronal time course in fine scale
FR_dTfine = createCellRegressor_indMov_discreteTime(dirDataNeural, paramSDF.setCellIDs(validC), ... %cellstr(paramCorr.validChanID),...
    setMovie, 0.1); % in 10Hz (number of spikes for every 100ms)
% then make an array in time x movie x cell
neuraltc=[];
for iChan = 1:length(validC) % for each channel
    clear tempn neuralrgrs
    [tempn{1:3}] = deal(FR_dTfine(iChan, indMovieNeuron).mnFR);
    neuralrgrs = cell2mat(tempn);
    
    neuraltc = cat(3, neuraltc, neuralrgrs);
end

saveFileNameNeural = sprintf('neuraltc4computeCorrcoeffCI_%s_10hz.mat', nameSubjNeural)
save(fullfile(dirDataNeural, saveFileNameNeural), 'neuraltc')

