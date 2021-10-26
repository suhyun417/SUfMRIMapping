function [] = getSigVoxIndex(nameSubjNeural, nameSubjBOLD, indUnit)

% nameSubjNeural = 'Tor';
% nameSubjBOLD = 'Art';

dirDataHome = '/procdata/parksh/';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);
dirDataCI = fullfile(dirDataNeural, 'corrMap_resultsCI/resultsCI_eachCell');


% Original correlation map
load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'paramCorr', 'matR_SU') % get cell IDs from another file


% convert the map to the surface 
% for iUnit = 1:size(matR_SU,2)
    
%     cd /projects/parksh/NeuralBOLD/analysis/BlockAna/
%     cellID = paramCorr.validChanID(indUnit,:);
    fprintf(1, 'Unit # %d, Cell ID: %s \n', indUnit, cellID);
    
%     fname = sprintf('%s_%s_Movie123_CI95_noFiltering+orig.BRIK', cellID, nameSubjNeural);%
    
    
    fileNameCI = sprintf('corrcoeffCI_%s%s_cell%d.mat',...
        nameSubjNeural, nameSubjBOLD, iUnit);
    load(fullfile(dirDataCI, fileNameCI))
    
    catCI95 = cat(1, resultBS(1).VoxCI.matCI95);
    catCI99 = cat(1, resultBS(1).VoxCI.matCI99);
    
    CI95_in = catCI95(mask_amp1_sub, :);
    setCI95_in(iUnit, :) = mean(CI95_in);
    % catRho_org = cat(1, resultBS(1).VoxCI.rho_org).*(-1);
    corrMap_org = matR_SU(:,iUnit);

    aa=cat(2, corrMap_org<catCI95(:,1), corrMap_org>catCI95(:,2));
    sigMask_CI95 = sum(aa,2);
    
    aa2=cat(2, corrMap_org<catCI99(:,1), corrMap_org>catCI99(:,2));
    sigMask_CI99 = sum(aa2,2);


%% 
% 3) fMRI movie-driven activity mask
load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp', 'brainMask_BlockAna3D');

[a, b, c] = ind2sub(size(movieDrivenAmp.mask_amp1), find(movieDrivenAmp.mask_amp1==1)); % Get indices of AF voxels in EPI 3D coords
mask_amp1 = [a b c];
mask_amp1_sub = sub2ind(size(movieDrivenAmp.mask_amp1), a, b, c);
    
    
    
    
    
    
