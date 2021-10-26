% genFig_moviefMRIPCA_dataMining.m
%
% 2021/04/15 SHP
% quick & dirty data mining of shared activity during movie fMRI

setSubjName = {'Art', 'Ava'};

for iSubj = 1:length(setSubjName)
    
    nameSubjBOLD = setSubjName{iSubj}; % 'Art';
    load(sprintf('/procdata/parksh/_macaque/%s/%s_movieTS_fMRI_Movie123_PCA.mat', nameSubjBOLD, nameSubjBOLD), '*concat*')
    S(iSubj).resultsPCA_concat_brainmask = resultsPCA_concat_brainmask;
    S(iSubj).resultsCorrPC_concat_brainmask = resultsCorrPC_concat_brainmask;
    
    load(sprintf('/procdata/parksh/_macaque/%s/%s_MaskArrays.mat', nameSubjBOLD, nameSubjBOLD), 'movieDrivenAmp', 'brainMask_BlockAna3D');
    brainmask_vec = reshape(movieDrivenAmp.map_sm_brain>0, nVox, 1); % change the 3D mask to 1D
    S(iSubj).brainmask_vec = brainmask_vec;
    
    load(sprintf('/procdata/parksh/_macaque/%s/%s_movieTS_fMRI_indMov.mat', nameSubjBOLD, nameSubjBOLD))
    indMovieBOLD = [1 2 3];
    catTS = []; catTS_nan = []; catTS_brainmask = []; catTS_brainmask_nan = [];
    for iM = 1:3
        fmritc=[];
        curvoltc = voltcIndMov{iM};
        avgvoltc = repmat(nanmean(curvoltc,4),[1 1 1 size(curvoltc,4)]);
        if ~isempty(find(avgvoltc==0, 1))
            avgvoltc(avgvoltc==0) = realmin; % get rid of zeros because it causes NaNs in percent signals
        end
        pcvoltc = ((curvoltc - avgvoltc)./avgvoltc)*100;
        fmritc = pcvoltc(:,:,:,8:125); %pcvoltc;
        [nx, ny, nz, nt] = size(fmritc);
        nVox = nx*ny*nz;
        matTS = reshape(fmritc, nVox, nt);
        
        catTS = cat(2, catTS, matTS);
        catTS_nan = cat(2, catTS_nan, NaN(size(matTS, 1), 7), matTS);
        
        catTS_brainmask = cat(2, catTS_brainmask, matTS(brainmask_vec, :)); % 27113 voxels
        catTS_brainmask_nan = cat(2, catTS_brainmask_nan, NaN(sum(brainmask_vec), 7), matTS(brainmask_vec, :));
    end
    
    % load(sprintf('/procdata/parksh/_macaque/%s/%s_MaskArrays.mat', nameSubjBOLD, nameSubjBOLD), 'movieDrivenAmp', 'brainMask_BlockAna3D');
    % brainmask_vec = reshape(movieDrivenAmp.map_sm_brain>0, nVox, 1); % change the 3D mask to 1D
    % catTS_brainmask = catTS(brainmask_vec, :); %matR_SU(brainmask_vec,:); % 27113 voxels
    
    S(iSubj).catTS = catTS;
    S(iSubj).catTS_nan = catTS_nan;
    S(iSubj).catTS_brainmask = catTS_brainmask;
    S(iSubj).catTS_brainmask_nan = catTS_brainmask_nan;
    
end

save('/procdata/parksh/_macaque/ArtAva_movieTS_fMRI_concatTS_pca.mat', 'S')


figure
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
subplot(3,1,1)
plot(S(1).resultsPCA_concat_brainmask.coeff(:,1), 'k-', 'LineWidth', 2)
axis tight
axis off
subplot(3,1,2)
plot(S(1).resultsPCA_concat_brainmask.coeff(:,2), 'b-', 'LineWidth', 2)
axis tight
axis off
subplot(3,1,3)
plot(S(1).resultsPCA_concat_brainmask.coeff(:,3), 'm-', 'LineWidth', 2)
axis tight
axis off

clear PC
for iPC = 1:3
for iS = 1:2
    
    tt = reshape(S(iS).resultsPCA_concat_brainmask.coeff(:,iPC), 118, 3);
    tt = cat(1, NaN(7, 3), tt);
    
PC{iPC}(:,iS) = tt(:); %S(iS).resultsPCA_concat_brainmask.coeff(:,iPC);
end
end

figure
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
for iPC = 1:3
SP(iPC) = subplot(3,1,iPC);
P{iPC} = plot(2.4:2.4:900, PC{iPC}, '-', 'LineWidth', 2);
axis tight
end

%%  movie scenes
dirMovie = '/procdata/parksh/Stimulus/Movies/_etc/Rhesus';
videoObj1 = VideoReader(fullfile(dirMovie, 'Movie1.avi'));
videoObj2 = VideoReader(fullfile(dirMovie, 'Movie2.avi'));
videoObj3 = VideoReader(fullfile(dirMovie, 'Movie3.avi'));


%% Compute the highest/lowest time points of the Principal Component time series
% Select the time points where PC time course is high/low
hdelay = 7; % hemodynamic delay in seconds
% crit = 0.1; %0.2; % criterion for high or low
critFrame = 5; % criterion in number of frames

% iPC = 1; % principal component index
% iM = 1; % movie index
axisTime = [8:125].*2.4; % TR = 2.4 sec, first 7 time points are excluded due to Onset effect
indFrame_sec = repmat(1:300, 30, 1);
indFrame_sec = indFrame_sec(:); % vectorize

selectedFrame_PC=[];
for iPC = 1:3 % principal component index
    
    for iM = 1:3 % movie index
        PC_sec = cat(2, axisTime(:), S(1).resultsPCA_concat_brainmask.coeff(118*(iM-1)+1:118*(iM), iPC)); % %resultsPCA(iM).coeff(:,iPC));
        [sortedPC, ind] = sortrows(PC_sec, 2); % default is ascending order
        selectedFrame_PC(:,:,1) = sortedPC(1:critFrame, :); %(1:round(118*crit), :); % highest 20% PC
        selectedFrame_PC(:,:,2) = sortedPC(end-critFrame+1:end, :); %sortedPC(end-round(118*crit)+1:end, :); % lowest 20% PC
        
        indFrame_PC = round(squeeze(selectedFrame_PC(:,1,:))); % selected frames in second (rounded)
        indFrame_PC = indFrame_PC - hdelay; % apply the hemodynamic delay        
        
        switch iM
            case 1
                obj = videoObj1;
            case 2
                obj = videoObj2;
            case 3
                obj = videoObj3;
        end
        % validFrame_range = [(5*(iMovie-1)*60*30)+1, (5*iMovie*60*30)];
        
        for iType = 1:2 % high or low PC
            locValidFrame = find(ismember(indFrame_sec, indFrame_PC(:, iType))>0);
            
            % Preallocate the movie structure for three-movie length
            clear reverseCorrMov
            reverseCorrMov(1:length(locValidFrame)) = struct('cdata',zeros(obj.Height, obj.Width, 3,'uint8'), 'colormap',[]);
            
            % Read one frame at a time.
            for k = 1 : length(locValidFrame)
                reverseCorrMov(k).cdata = read(obj, locValidFrame(k));
            end
            % countFrame = countFrame + length(validLocHigh);
            
            % % Size a figure based on the video's width and height.
            % hf = figure;
            % set(hf, 'position', [150 150 obj.Width  obj.Height])
            % % Play back the movie once at the video's frame rate.
            % movie(hf, reverseCorrMov, 1, obj.FrameRate);
            
            % Write a movie
%             fileName = sprintf('reverseCorrScenes_%s_fMRI_concatTS_PC%d_movie%d_hdelay%d_crit%s_%d.avi', nameSubjBOLD, iPC, iM, ...
%                 hdelay, strrep(sprintf('%0.2f', crit), '.', 'p'), iType);
            fileName = sprintf('reverseCorrScenes_%s_fMRI_concatTS_PC%d_movie%d_hdelay%d_critNumFrame%d_%d.avi', nameSubjBOLD, iPC, iM, ...
                hdelay, critFrame, iType);
            myObj = VideoWriter(fileName);
            myObj.FrameRate = obj.FrameRate; %10;
            %     myObj.Path = dirMovie;
            open(myObj);
            writeVideo(myObj, reverseCorrMov)
            close(myObj);
            
            % Move files to /procdata
            movefile('./reverseCorrScenes_*.avi', dirMovie)
        end
    end
end








