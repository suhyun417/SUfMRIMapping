% analPCA_moviefMRI_movieFeatures.m
% 
% 2019/02/05 SHP
%       - Extract the scenes that drive the shared fMRI responses across
%       the visual cortex 


%% Load the PCA results & movie files
nameSubjBOLD = 'Art';
load(sprintf('/procdata/parksh/%s/%s_movieTS_fMRI_Movie123_PCA.mat', nameSubjBOLD, nameSubjBOLD), 'resultsPCA')

% Movie scenes
% Read videos (takes long time)
dirMovie = '/procdata/parksh/Stimulus/Movies/_etc/Rhesus';
videoObj1 = VideoReader(fullfile(dirMovie, 'Movie1.avi'));
videoObj2 = VideoReader(fullfile(dirMovie, 'Movie2.avi'));
videoObj3 = VideoReader(fullfile(dirMovie, 'Movie3.avi'));

% % Side paragraph to compare the principal components (coefficients) across
% % two monkeys
% setNameSubjBOLD = {'Art', 'Ava'};
% iPC = 1; % ID of the principal component
% for iSubj = 1:length(setNameSubjBOLD)
%     nameSubjBOLD = setNameSubjBOLD{iSubj};
%     load(sprintf('/procdata/parksh/%s/%s_movieTS_fMRI_Movie123_PCA.mat', nameSubjBOLD, nameSubjBOLD), 'resultsPCA')
%     pc_fMRI = [];
%     for iM = 1:3
%         pc_fMRI = cat(1, pc_fMRI, NaN(7,1), resultsPCA(iM).coeff(:,iPC));
%     end
%     catPC(:,iSubj) = pc_fMRI;
% end
% 
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [360 150 530 710])
% for iM = 1:3
%     subplot(3, 1, iM) % for each movie
%     plot([1:125].*2.4, catPC(125*(iM-1)+1:125*iM, :), 'o-') % TR is 2.4 sec
%     set(gca, 'XTick', 30:30:300, 'XTicklabel', 30:30:300) % time in seconds
% end



%% Compute the highest/lowest time points of the Principal Component time series
% Select the time points where PC time course is high/low
hdelay = 7; % hemodynamic delay in seconds
crit = 0.2; % criterion for high or low

% iPC = 1; % principal component index
% iM = 1; % movie index
axisTime = [8:125].*2.4; % TR = 2.4 sec, first 7 time points are excluded due to Onset effect
indFrame_sec = repmat(1:300, 30, 1);
indFrame_sec = indFrame_sec(:); % vectorize


for iPC = 1:3 % principal component index
    
    for iM = 1:3 % movie index
        PC_sec = cat(2, axisTime(:), resultsPCA(iM).coeff(:,iPC));
        [sortedPC, ind] = sortrows(PC_sec, 2); % default is ascending order
        selectedFrame_PC(:,:,1) = sortedPC(1:round(118*crit), :); % highest 20% PC
        selectedFrame_PC(:,:,2) = sortedPC(end-round(118*crit)+1:end, :); % lowest 20% PC
        
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
            fileName = sprintf('reverseCorrScenes_%s_fMRI_movie%d_PC%d_hdelay%d_crit%s_%d.avi', nameSubjBOLD, iM, iPC, ...
                hdelay, strrep(sprintf('%0.2f', crit), '.', 'p'), iType);
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

% %% Make videos
% for iClust = 1:7
%     
%     % Find when activity was high (i.e. exceeds certain criterion in z-score)
%     critHigh = 3; %1.5; %2; % in z-score
%     locHigh = find(meanFRCluster_zscore(:,iClust)>critHigh); 
% %     critHigh = -1; %3; %1.5; %2; % in z-score
% %     locHigh = find(meanFRCluster_zscore(:,iClust)<critHigh);
%     if length(locHigh) < 1
%         continue;
%     end
%     % [i,j] = ind2sub(size(meanFRCluster_zscore), locHigh);
%     % taxis = 1:900; % in second
%     
%     % Preallocate the movie structure for three-movie length
%     clear reverseCorrMov
%     reverseCorrMov(1:length(locHigh)) = struct('cdata',zeros(vidHeight,vidWidth, 3,'uint8'), 'colormap',[]);
%     
%     % Collect the scenes
%     countFrame = 0;
%     for iMovie = 1:3
%         validFrame_range = [(5*(iMovie-1)*60*30)+1, (5*iMovie*60*30)];
%         validLocHigh = [];
%         validLocHigh = locHigh(locHigh>validFrame_range(1) & locHigh<validFrame_range(2));
%         
%         switch iMovie
%             case 1
%                 obj = videoObj1;
%             case 2
%                 obj = videoObj2;
%             case 3
%                 obj = videoObj3;
%         end
%         
%         % Read one frame at a time.
%         for k = 1 : length(validLocHigh)
%             indFrame = validLocHigh(k)-(5*(iMovie-1)*60*30);
%             reverseCorrMov(countFrame+k).cdata = read(obj,indFrame);
%         end
%         countFrame = countFrame + length(validLocHigh);
%     end
%     
% %     % Size a figure based on the video's width and height.
% %     hf = figure;
% %     set(hf, 'position', [150 150 videoObj1.Width  videoObj1.Height])
% %     % Play back the movie once at the video's frame rate.
% %     movie(hf, reverseCorrMov, 1, 10); %videoObj1.FrameRate);
%     
%     % Write a movie
%     newIndCluster = [3 2 7 6 5 1 4];
%     fileName = sprintf('scenesMov123_Cluster%d_meanClustSDF_zCrit%d_OldCluster%d.avi', newIndCluster(iClust), critHigh, iClust);
% %     fileName = sprintf('scenesMov123_Cluster%d_meanClustSDF_zCritLow%d_OldCluster%d.avi', newIndCluster(iClust), critHigh, iClust);
%     myObj = VideoWriter(fileName);
%     myObj.FrameRate = 10;
% %     myObj.Path = dirMovie;
%     open(myObj);
%     writeVideo(myObj, reverseCorrMov)
%     close(myObj);
% end
% 
% % Move files to /procdata
% movefile('./scenesMov123*.avi', dirMovie)

%% Data mining: compare time series & movie scenes, instead of selecting out scenes
% Get the scene info from DM's results
for iMovie = 2:3
    
    switch iMovie
        case 1
            obj = videoObj1;
        case 2
            obj = videoObj2;
        case 3
            obj = videoObj3;
    end
    
    load(sprintf('/procdata/parksh/MovieRegressors/annotationMovie%d.mat', iMovie))
    sceneInfo = cat(2, sta', sto');
    nScene = size(sceneInfo, 1);
    
    oldOrderCluster = [6 2 1 7 5 4 3];
    
    figTS = figure;
    set(figTS, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1200 900])
    SP(1) = subplot('Position', [0.1300    0.65    0.7750    0.3]); %[0.1 0.75 0.8 0.2]); % for time series
    SP(2) = subplot('Position', [0.1 0.05 0.5 0.5]); % for movie image
    
    % % Size a figure based on the video's width and height.
    % hf = figure;
    % set(hf, 'Color', 'w','position', [150 150 obj.Width  obj.Height])
    
    for iScene = 1:nScene
        
        setFrames = sceneInfo(iScene,1):sceneInfo(iScene,2);
        
        
        subplot(SP(1))
        plot(SP(1), setFrames, meanFRCluster_zscore(setFrames+(5*(iMovie-1)*60*30), oldOrderCluster), 'o-')
        LG= legend('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4', 'Cluster 5', 'Cluster 6', 'Cluster 7',...
            'Location', [0.0056    0.8433    0.0804    0.1497]);
        axis tight
        ylim([-1 10])
        xlabel('Frame #')
        ylabel('Normalized response (z)')
        title(sprintf('Movie %d: Scene %d', iMovie, iScene))
        
        frame = read(obj,setFrames(1));
        a(1).cdata = frame;
        a(1).colormap = [];
        imageFrame = frame2im(a);
        subplot(SP(2));
        imagesc(imageFrame);
        
        while 1
            [x, y, but] = ginput(1);
            indFrame = round(x);
            frame = read(obj,indFrame);
            a(1).cdata = frame;
            a(1).colormap = [];
            imageFrame = frame2im(a);
            
            %         figure(hf);
            subplot(SP(2));
            imagesc(imageFrame);
            
            if but ~= 1
                break;
            end
            
        end
    end
end