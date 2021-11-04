% genFig_multipleFP_Revision_correlationCells.m
%
% 2021/11/1 SHP
% Compute correlations across neurons in each area using multiple temporal
% scales
% This code is generated during a revision process of Park et al., 2021 Sci Adv.

%% Settings
flagBiowulf = 0;

if flagBiowulf
    directory.dirDataHome = '/data/parks20/procdata/NeuroMRI/';
    addpath('/data/parks20/analysis/NeuroMRI/'); % to use doConv.m function
else
    ss = pwd;
    if ~isempty(strfind(ss, 'Volume')) % if it's local
        directory.projects = '/Volumes/PROJECTS';
        directory.procdata = '/Volumes/PROCDATA';
        directory.dataHome = fullfile(directory.procdata, 'parksh', '_macaque');
        directory.library = '/Volumes/LIBRARY';
        addpath(fullfile(directory.library, 'matlab_utils'));
    else % on virtual machine
        directory.projects = '/nifvault/NIFVAULT/projects';
        directory.procdata = '/procdata';
        directory.dataHome = fullfile(directory.procdata, 'parksh', '_macaque');
        directory.library = '/library';
        addpath(fullfile(directory.library, 'matlab_utils'));
    end
end

setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi', 'Mat', 'Dan', 'Moc', 'Was', 'Dav'};

load(fullfile(directory.dataHome, sprintf('matRaster_Movie123_allCells.mat')))
load('/procdata/parksh/_macaque/multipleFP_fsi.mat')
locFaceCell =  find(fsi.matFSI(:,1)>0.33); % find(abs(fsi.matFSI(:,1))>0.33);

matRaster = matTS_FP.matRaster(:, locFaceCell); 

setSigma = [5 20 50]; 

% Spike density functions with different sized Guassian
for iSigma = 1:length(setSigma)
    % spike density function
    curSigma  = setSigma(iSigma);   % gaussian kernel SD, in ms
    k = normpdf(-40:40, 0, curSigma)';
    
    cursdf     = conv2(matRaster, k, 'same').*1000;
    
    matSDF{iSigma, 1} = cursdf;
end

    
%% Compute correlation

