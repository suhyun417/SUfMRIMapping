% tempVoxCorr
%
% get correlation of AF voxels 

S_neuralRegressor % do this first to get necessary structs such as STDPATH, dataBOLD, etc.

global DSP STDPATH DATA 
% 
% dirROI = fullfile(STDPATH.dataBOLD, 'ROIs');
% d_vis = dir(fullfile(dirROI, '*VisROIs.mat'));
% d_face = dir(fullfile(dirROI, '*faceROIs2.mat'));
% load(fullfile(dirROI, d_vis.name));
% load(fullfile(dirROI, d_face.name));
% 
% % get ROI information into structures with different names
% tempName_visROIs = char(fieldnames(load(fullfile(dirROI, d_vis.name))));
% tempName_faceROIs = char(fieldnames(load(fullfile(dirROI, d_face.name))));
% eval(['visROIs=', tempName_visROIs, ';']) % get visual ROI data into "visROIs"
% eval(['faceROIs=', tempName_faceROIs, ';']) % get face ROI data into "faceROIs"

ROIs = cat(2, visROIs', faceROIs'); % ROIs: 8x2 struct, 1st column visual ROIs, 2nd column face ROIs

% % clear up
% eval(['clear ' tempName_visROIs])
% eval(['clear ' tempName_faceROIs])
% clear visROIs faceROIs


% Get ROI names & coordinates in EPI coords 
resizeFactor = DSP.proc.params3d.res./ROIs(1,1).params.res; % calculate the scaling factor

dirFig = '/projects/parksh/NeuralBold/_labNote/_figs/';

% for each ROI, movie by movie
flagFaceArea = 1; %0;

for iFR=1:8
for selMov=1:9
    
    clear voxROI a b c matBOLD_ROI matR vecValidR

    nameROI = ROIs(iFR,1+flagFaceArea).name(strfind(ROIs(iFR,1+flagFaceArea).name, '_')+1:end); %nameROI =
    voxROI=decimate3D(ROIs(iFR,1+flagFaceArea).vol3D, resizeFactor, .25); % turn anat_res ROIs into func_res ROIs
    
    [a, b, c] = ind2sub(size(voxROI), find(voxROI==1)); % Get indices of AF voxels in EPI 3D coords
    
    for iVox = 1:length(a)
        matBOLD_ROI(:,iVox) = dataBOLD.mvoltc{selMov}(a(iVox), b(iVox), c(iVox),:);
    end
    
    matR = corr(matBOLD_ROI, 'rows', 'complete');
    vecValidR = matR(triu(matR, 1)~=0); % R values above the diagonal
    
    % draw histogram of corr. coef.
    figure(105); clf; set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
    hist(vecValidR, 30);
    hold on;
    line([mean(vecValidR) mean(vecValidR)], get(gca, 'YLim'), 'Color', 'r', 'LineWidth', 2)
    line([median(vecValidR) median(vecValidR)], get(gca, 'YLim'), 'Color', 'm', 'LineWidth', 2)
    xlim([-1 1])
    
    title(sprintf('Histogram of Corr. Coef. between %s voxels', nameROI))
    xlabel('Correlation coefficient')
    ylabel('Frequency')
    input('')
    
%     fileName = sprintf('art%smov%d_corrHist', nameROI, S(1,iMov).movID);
%     print(gcf, '-depsc', fullfile(dirFig, fileName));
end

end

% Correlation maps of AF voxel responses
% first, average AF voxels
clear voxROI a b c matR vecValidR

% iFR=4; % AF
nameROI = ROIs(iFR,1+flagFaceArea).name(strfind(ROIs(iFR,1+flagFaceArea).name, '_')+1:end); %faceROIs(iFR).name(strfind(faceROIs(iFR).name, '_')+1:end);
voxROI=decimate3D(ROIs(iFR,1+flagFaceArea).vol3D, resizeFactor, .25); %decimate3D(faceROIs(iFR).vol3D, resizeFactor, .25); % turn anat_res ROIs into func_res ROIs
[a, b, c] = ind2sub(size(voxROI), find(voxROI==1)); % Get indices of AF voxels in EPI 3D coords
indVox_ROI = [a b c];
indVox_ROI_sub = sub2ind(size(voxROI), a, b, c);

% Get the averaged time series (across runs)
selMov=1:9; %[1 2 3]; %4; %3; %2;

matBOLD = [];
matBOLD_ROI=[];

matBOLD = cat(4, dataBOLD.catmvoltc{selMov});
% shiftSelMov = circshift(selMov', 2)';
% matBOLD_shuffle = cat(4, dataBOLD.catmvoltc{shiftSelMov});
for iVox = 1:length(a)
    matBOLD_ROI(:,iVox) = matBOLD(a(iVox), b(iVox), c(iVox),:); %matBOLD_shuffle(a(iVox), b(iVox), c(iVox),:); %matBOLD(a(iVox), b(iVox), c(iVox),:);
end

matR = corr(matBOLD_ROI, 'rows', 'complete');
vecValidR = matR(triu(matR, 1)~=0); % R values above the diagonal

meanBOLD_ROI = nanmean(matBOLD_ROI, 2);
steBOLD_ROI = nanstd(matBOLD_ROI, [], 2)./sqrt(size(matBOLD_ROI,2)-1);

curRGR = meanBOLD_ROI;
% Compute correlation maps
[nx, ny, nz, nt] = size(matBOLD);
nVox = nx*ny*nz;

[Rvals, Pvals] = corr(reshape(matBOLD, nVox, nt)', curRGR,...
    'rows','complete');

mapR = reshape(Rvals, [nx, ny, nz]); % Don't need to multiply -1 because now it's correlation between MION signals %.*(-1);
DSP.proc.scalarmap_3d = mapR;



% %% Load each run's raw time series 
% addpath ('/projects/parksh/NeuralBold/analysis') %addpath ('/einstein0/USRlab/projects/parksh/analysis')
% 
% % File Details
% nameSubj = 'Art'; %'Ava'; %'Art'; % one of the inputs of this function
% monkID = SU_getMonkIDByName(nameSubj);
% 
% % % set save file names
% % saveFileNameHeader = sprintf('%s_DataFMRI', nameSubj);
% % if optDetrend
% %     saveFileNameHeader = [saveFileNameHeader '_det'];
% % end
% % if optPercentSig
% %     saveFileNameHeader = [saveFileNameHeader '_ps'];
% % end
% 
% global S dataBOLD
% global STDPATH DSP DATA GH
% 
% % First start,open the XL file, and reset the windows
% blockpaths_soohyun;
% STDPATH.xlsfile = 'block_timecourse_files_ber.xls';
% loadXLSFile;
% 
% 
% textfileroot = sprintf('FL_%s_allfiles2', monkID);
% % textfileroot = 'FL_e66_allfiles2'; %'FL_e66_9FILES'; %'FL_e66_allfiles2'; %% MUST be matched to the fMRI data
% 
% %   MCD = loadMeanConcatData(textfileroot,skip,redo);
% 
% skip=0;
% sesstextfile = sprintf('%s.txt',textfileroot);
% SI = SU_createSessInfo(sesstextfile,[],skip); % get infor from text file
% filelist = SU_makeFileList(SI.monkID,SI.sessID,SI.scanID);
% 
% %
% unimov = [1 2 3 10 11 12 13 14 15]; %sort(unique(SI.movID));
% nunimov = length(unimov);
% % catimgdat = [];
% % rgr = zeros(1,nunimov*MAX_TR);
% % procmovs = [];
% 
% for u = 1:nunimov
%     movindxs = find(SI.movID == unimov(u));
%     nmovindxs = length(movindxs);
%     
%     % collect all the data files for this movie
%     %
%     filelist_sub = SU_makeFileList(SI.monkID(movindxs),SI.sessID(movindxs),...
%         SI.scanID(movindxs));
%     
%     totfiles = length(filelist_sub); 
%     matTC_ROI=[];
%     for iFile=1:totfiles
%         
%         s_sub = filelist_sub{iFile}.subj;
%         s_ses = filelist_sub{iFile}.sess;
%         s_sc  = filelist_sub{iFile}.scan;
%         
%         cursession = sprintf('%s.%s',s_sub,s_ses);
%         sessindx   = strmatch(cursession,DSP.filedetails.sess);
%         scanindx   = find(DSP.filedetails.scan == s_sc);
%         xlsindx    = intersect(sessindx,scanindx);
%         
%         dateFiles(iFile,:) = DSP.filedetails.date{xlsindx};
%         
%         [fmri_tc,dgz,notes] = SU_loadFile(s_sub,s_ses,s_sc);  
%         [nx ny nz nt] = size(fmri_tc);
%                 
%         if nt == 250
%             startTR = 63;% offset when movie starts in scan
%         elseif nt == 125
%             startTR = 1; % short movie
%         end      
%         validTR         = [startTR:(startTR+125-1)];
%         fmri_tc = fmri_tc(:,:,:,validTR);
% 
%         
%         % quick & dirty concatenation of a ROI voxel's time series        
%         fmri_tc_reshape = reshape(fmri_tc, nx*ny*nz, 125)'; % time x voxel
%         matTC_ROI{iFile} = fmri_tc_reshape(:,indVox_ROI_sub); %fmri_tc_reshape(:,indVox_ROI_sub);
% 
%     end
%     
%     matTC_concat = cat(1, matTC_ROI{:});
%     
%     for iVox = 1:size(matTC_concat,2)
%     tempTC = reshape(matTC_concat(:,iVox), 125, 36);         
%     tempTC_ps = (tempTC-repmat(mean(tempTC), 125,1))./repmat(mean(tempTC), 125,1).*100;
%     
%     figure; set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1300 300 450 600]);   
%     sp1=subplot(2,1,1);
%     imagesc(tempTC', [0 1600]) 
%     colormap(sp1, hot);     
%     title(sprintf('Raw time series of voxel ID%d, Movie #%d', iVox, unimov(u)))
%     freezeColors 
%     cbfreeze(colorbar)
%     
%     sp2=subplot(2,1,2);
%     imagesc(tempTC_ps(8:end,:)', [-30 30])
%     colormap(sp2, jet); colorbar;
%     title(sprintf('Time series in percent signal of Voxel ID%d, Movie #%d', iVox, unimov(u)))
%     
% %     print(gcf, '-dtiffnocompression', fullfile(dirFig, sprintf('Art_AF_Vox%d_m%d', iVox, unimov(u))))
%     saveas(gcf, fullfile(dirFig, sprintf('Art_V1F_Vox%d_m%d', iVox, unimov(u))), 'tif')
% %     input('')
%     end
%     
% end








