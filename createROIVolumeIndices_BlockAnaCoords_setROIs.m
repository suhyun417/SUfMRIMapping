% createROIVolumeIndices_BlockAnaCoords_setROIs.m
% 2020/06/23 SHP
% 
% modified from "createROIVolumeIndices_BlockAnaCoords.m"
% Load manually drawn ROI (converted into AFNI BRIK/HEAD 3d volume from
% surface dataset) 
% Get ROI indices from each hemisphere
% Make a matrix of ROI indices for both hemisphere (for central visual
% cortex, motion areas, and face patches)

clear all;

addpath('/projects/parksh/_toolbox/afni_matlab/')

for typeMotionROI = 1:2
    switch typeMotionROI
        case 1
            roifname{1,1} = '/procdata/parksh/Art/Anatomy/_suma/Art_ROIs_CentralVC+Motion1_rh_surfVolFun+orig.BRIK'; % labels are 1 (visual cortex) and 2 (motion areas)
            roifname{1,2} = '/procdata/parksh/Art/Anatomy/_suma/Art_ROIs_CentralVC+Motion1_lh_surfVolFun+orig.BRIK';
            roifname{2,1} = '/procdata/parksh/Art/ROIs/e66_faceROIs2.mat'; % RH face patches
            roifname{2,2} = '/procdata/parksh/Art/ROIs/e66_faceROIs_L.mat'; % LH face patches
            manualROILabel = [1 2]; % for motion1 % [1 3]; % for motion2
            
        case 2
            roifname{1,1} = '/procdata/parksh/Art/Anatomy/_suma/Art_ROIs_CentralVC+Motion2_rh_surfVolFun+orig.BRIK';  % labels are 1 (visual cortex) and 3 (motion areas)
            roifname{1,2} = '/procdata/parksh/Art/Anatomy/_suma/Art_ROIs_CentralVC+Motion2_lh_surfVolFun+orig.BRIK';
            roifname{2,1} = '/procdata/parksh/Art/ROIs/e66_faceROIs2.mat'; % RH face patches
            roifname{2,2} = '/procdata/parksh/Art/ROIs/e66_faceROIs_L.mat'; % LH face patches
            manualROILabel = [1 3]; % for motion2
    end
    
    
    nameROI = {'Central retinotopic cortex', 'Motion areas', 'Face patches'};
    numROI = length(nameROI); %3;
    
    nx = 40; ny = 64; nz = 32;
    nVox = nx*ny*nz;
    
    
    matROIIndices = NaN(nVox, numROI);
    for iROI = 1:numROI
        if iROI < 3
            curROIIndices = [];
            for iH = 1:2
                [err, volroi, info, errM] = BrikLoad(roifname{1, iH});
                roi1_vol = permute(volroi, [3 1 2]); % to be matched to the other data that are in BlockAna coords
                roi1_vec = reshape(roi1_vol, nVox, 1); % change the 3D vol to 1D
                curROIIndices(:,iH) = roi1_vec==manualROILabel(iROI);
            end
            matROIIndices(:,iROI) = sum(curROIIndices, 2);
        else % face patch ROIs
            load(roifname{2, 1}) % RH
            load(roifname{2, 2}) % LH
            faceROIs = struct([]);
            faceROIs = cat(2, e66_faceROIs2, e66_faceROIs_L);
            clear e66*
            
            resizeFactor =[3 3 3]; % calculate the scaling factor
            
            % convert the map to the surface
            % First, get all the map together
            [dvolx, dvoly, dvolz]  = size(faceROIs(1).vol3D); % size(visROIs(1).vol3D); % dimension of volume
            catAllFaceROIs = cat(4, faceROIs.vol3D);
            matAllFaceROIs = reshape(catAllFaceROIs, dvolx*dvoly*dvolz, length(faceROIs));
            
            allFaceROIsInd = sum(matAllFaceROIs,2); %matAllFaceROIs*[1:size(faceROIs, 2)]'; %matAllFaceROIs*[1:8]'; %sum(matAllFaceROIs,2); %matAllFaceROIs*[1:8]';
            allFaceROIsInd_vol =  reshape(allFaceROIsInd, [dvolx, dvoly, dvolz] );
            allFaceROIsInd_fun = decimate3D(allFaceROIsInd_vol, resizeFactor, .25); % to EPI res
            
            matROIIndices(:,iROI) = reshape(allFaceROIsInd_fun, nVox, 1); % change the 3D vol to 1D again
        end
    end
    
    paramROI.analCode = '/projects/parksh/NeuralBOLD/analysis/prepROI.m';
    paramROI.roifname = roifname;
    paramROI.manualROILabel = manualROILabel;
    paramROI.nameROI = nameROI;

    save(sprintf('/procdata/parksh/Art/ROIs/Art_ROIs_RetinoMotion%dFace_bothH.mat', typeMotionROI), 'matROIIndices', 'paramROI')

end
