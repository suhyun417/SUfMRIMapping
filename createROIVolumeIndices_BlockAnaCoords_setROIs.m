% createROIVolumeIndices_BlockAnaCoords_setROIs.m
% 2020/06/23 SHP
%
% modified from "createROIVolumeIndices_BlockAnaCoords.m"
% Load manually drawn ROI (converted into AFNI BRIK/HEAD 3d volume from
% surface dataset)
% Get ROI indices
% Make a matrix of ROI indices

clear all;

addpath('/projects/parksh/_toolbox/afni_matlab/')

roifname = '/procdata/parksh/_macaque/Art/Anatomy/_suma/SUmap_ArtRH_ROIs_fromSurf+orig.BRIK';
nx = 40; ny = 64; nz = 32;
nVox = nx*ny*nz;

[err, volroi, info, errM] = BrikLoad(roifname);

roi_vol = permute(volroi, [3 1 2]); % to be matched to the other data that are in BlockAna coords
roi_vec = reshape(roi_vol, nVox, 1); % change the 3D vol to 1D

numROI = max(roi_vec);
matROIIndices = NaN(nVox, numROI);
for iROI = 1:numROI
    matROIIndices(:,iROI) = roi_vec==iROI;
end

paramROI.analCode = '/projects/parksh/NeuralBOLD/analysis/prepROI.m';
paramROI.roifname = roifname;

% Get the name (labels) for each ROI
fileID = fopen('/procdata/parksh/_macaque/Art/Anatomy/_suma/SUmap_ArtRH_ROIs_labels.txt');
C = textscan(fileID, '%f %s', 'delimiter', '\t');
fclose(fileID);

nameROI = C{2};
paramROI.nameROI = nameROI;

save('/procdata/parksh/_macaque/Art/ROIs/Art_ROIs_set00_RH.mat', 'matROIIndices', 'paramROI');
