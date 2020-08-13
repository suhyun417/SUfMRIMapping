% prepData4ClusteringMultipleFP.m
%
% 2020/08/13 SHP
% Prepare correlation matrices for clustering in Biowulf environment

clear all;

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

nameSubjBOLD = 'Art';

% load correlation matrix for all valid cells from four face patches
load(sprintf('/procdata/parksh/_macaque/CorrMap_SU_AllCells%s_corticalFPMerged.mat', nameSubjBOLD), 'corrMap_merged_FP'); %, 'info*', 'corrMap_Area', 'corrMap_merged');

% load the masks (movie-driven & brain-only mask)
load(sprintf('/procdata/parksh/_macaque/%s/%s_MaskArrays.mat', nameSubjBOLD, nameSubjBOLD), 'movieDrivenAmp') ;%, 'brainMask_BlockAna3D');


brainmask_vec = reshape(movieDrivenAmp.map_sm_brain>0, numel(movieDrivenAmp.map_sm_brain), 1); % change the 3D mask to 1D
matR_brainmask = corrMap_merged_FP.matR(brainmask_vec,:); % 27113 voxels

moviemask_vec = reshape(movieDrivenAmp.mask_amp1, numel(movieDrivenAmp.mask_amp1), 1); % change the 3D mask to 1D
matR_moviemask = corrMap_merged_FP.matR(moviemask_vec,:); % 15495 voxels

infoCells = corrMap_merged_FP;
infoCells = rmfield(infoCells, {'matR', 'meanFractionAcrossArea', 'meanMaxAbs', 'grandMeanR'});

save(sprintf('/procdata/parksh/_macaque/%s/matR4clusteringmultiplepatches.mat', nameSubjBOLD), 'matR*', 'infoCells');