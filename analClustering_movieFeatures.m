% analClustering_movieFeatures.m
%
% 2016/1/4 SHP
% To compare functions of different clusters, analyze the scenes drive the average activity 
% of each subgroup of cells in terms of visual features (using Brian's
% coding scheme)
% 
% 1) select the time points, 2)

addpath('/library/matlab_utils/')

% Load the data
nameSubjNeural = 'Tor';
nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataHome = '/procdata/parksh/';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% Clustering based on corr maps
load(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new_masked.mat', nameSubjNeural, nameSubjBOLD)));
% % Clustering based on time series
% load(fullfile(dirDataNeural, sprintf('ClusteringSDF_%sMovie123.mat', nameSubjNeural)))

% Directory for saving figures as graphic files
dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';


%% Compare different clustering
% Sort out cells based on particular K-means clustering 
matIndClust_SU = cat(2, Clustering_moviemask.resultKMeans.SU_indCluster); % based on corr map

sortTargetK = 5; %7;
[sortedClust, indSortChan]=sort(matIndClust_SU(:,sortTargetK-1));

%% Single unit time series
setMovie = [1 2 3];
load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'paramCorr') % we only need ID of valid channels

% Prepare neural responses
% First, get the cell response in 30fps time resolution 
FR_dT30 = createCellRegressor_indMov_discreteTime(dirDataNeural, cellstr(paramCorr.validChanID),...
    setMovie, 1/30); 
% concatenate across movies
matFR=[];
for iUnit = 1:size(FR_dT30,1)
    tempFR = cat(1, FR_dT30(iUnit, :).mnFR);
    matFR(:,iUnit) = tempFR;
end
% matFR_zscore = zscore(matFR);

% Prepare Brian's movie features in 30fps (using up-sampling for now)
% Version 1. upsampling from the original 71 features in 4fps (high-level) and 10fps (low-level)
[fullRGR30fps] = createFullMovieRegressors_30fps_indMov(setMovie);
catFullRGR = cat(1, fullRGR30fps.regressors);

% Version 2. upsampling from the 11 features selected features in 4fps
flagSM = 1; % flag for compression and smoothing
fullRGR4fps = createMovieRGR_4fps_indMov(setMovie, flagSM); %createFullMovieRegressors_4fps_indMov(setMovID); %

catRGRfull_30fps=[];matRGRfull_30fps=[];
for iMov=1:length(setMovie)
    m = setMovie(iMov);
    matCurRGR = fullRGR4fps(m).smoRegressors; %fullRGR4fps(iMov).regressors(:,indValidRGR); %fullRGR4fps(iMov).regressors;
    for iRGR = 1:size(matCurRGR, 2)
        tRGR = resample(matCurRGR(:,iRGR), 30, 4); 
        matCurRGR_30fps(:,iRGR) = tRGR;
    end
    catRGRfull_30fps = cat(1, catRGRfull_30fps, matCurRGR_30fps); %matCurRGR); % concatenation across movies
end
% scaleRGR_resampled = resample(scaleRGR, 2.4*100, 0.25*100); %matRGR = resample(catRGR, 0.25*100, 2.4*100);
matRGRfull_30fps = catRGRfull_30fps; %cat(2, catRGRfull, scaleRGR_resampled); %cat(2, matRGR, scaleRGR);
varnamesfull = fullRGR4fps(1).features; % cat(1, fullRGR4fps(1).features, {'Face size'});

% subset of regressors
indValidRGR = [1 2 6 7 3 21 20 32 22 31 25]; %[1, 3, 9, 20, 21, 22, 25]; 
% 1: 'Luminance', 2: 'Contrast', 6: Low spatial Frequency 7: High spatial frequencty 3: 'Motion (Speed)', 
% 21: 'One face', 20: 'Number of faces', 32: 'Face size', 22: 'Body parts', 31: 'Hands', 25: 'Any animal'
matRGRvalid_30fps = matRGRfull_30fps(:,indValidRGR);
varnamesvalid = varnamesfull(indValidRGR);

BRfeatureParams.version1.featureNames = fullRGR30fps(1).features;
BRfeatureParams.version1.notes = 'upsampled 71 regressors (including DM''s 5 size regressors) in 30fps from 4fps (high-level) and 10fps (low-level) using "createFullMovieRegressors_30fps_indMov.m"';
BRfeatureParams.version2.featureNames = varnamesvalid;
BRfeatureParams.version2.notes = 'upsampled 11 regressors in 30fps from 4fps';


% Then collect the scenes that drives average SDF within cluster across cells
meanFRCluster=[];steFRCluster=[];
for iK = 1:sortTargetK
tempMatFR=[];
tempMatFR = matFR(:,indSortChan(sortedClust==iK));
meanFRCluster(:,iK) = mean(matFR(:,indSortChan(sortedClust==iK)), 2);
steFRCluster(:,iK) = std(matFR(:,indSortChan(sortedClust==iK)), 0, 2)./length(indSortChan(sortedClust==iK));
end
meanFRCluster_zscore = zscore(meanFRCluster);

%%
% No Selection of the unique scenes
Clustering_movieFeatures=struct([]);
newIndCluster = [5 4 1 2 3]; % [2 6 4 1 5 3 7]; % what K means cluster 1 is our cluster 2 %[3 2 7 6 5 1 4]; % mapping between original clustering #s to new numbering scheme for the paper (old cluster 1 is new cluster 3)
for iClust = 1:5 %7
    
    % Find when activity was high (i.e. exceeds certain criterion in z-score)
    critHigh = 3; %4; %3; %1.5; %2; % in z-score
    locHigh = find(meanFRCluster_zscore(:,iClust)>critHigh); 
  
    % DM's scene-based features (face, torso, arms, legs)
    % 2016/09/26 Select unique "scenes" and save the feature info +
    % responses at those time points
    avgClusterResponse=[]; indFrameSelected=[]; setIndFrame_valid=[];
    for iMovie = 1:3
        load(sprintf('/procdata/parksh/MovieRegressors/annotationMovie%d.mat', iMovie))
        sceneInfo = cat(2, sta', sto');
        
        validFrame_range = [(5*(iMovie-1)*60*30)+1, (5*iMovie*60*30)];
        validLocHigh = [];
        validLocHigh = locHigh(locHigh>validFrame_range(1) & locHigh<validFrame_range(2));
        
        % Extract scene information frame-by-frame from DM's scene analysis
        indSceneList=[]; tempS=[];
        for iFrame = 1:length(validLocHigh)
            indFrame = validLocHigh(iFrame)-(5*(iMovie-1)*60*30);
            indScene = find(sum([sta<=indFrame; sto>=indFrame], 1)>1);
            
            if isempty(indScene)
                indSceneList(iFrame,1) = NaN;
                continue;
            end
            
            indSceneList(iFrame,1) = indScene;
            
            tempS(iFrame, 1) = epoch(indScene).notes.face.A;
            tempS(iFrame, 2) = epoch(indScene).notes.torso.A;
            tempS(iFrame, 3)= epoch(indScene).notes.arms.A;
            tempS(iFrame, 4) = epoch(indScene).notes.legs.A;
            tempS(iFrame, 5) = epoch(indScene).notes.viewAngle;
            
        end
        
        DMfeatures(iMovie).sceneIndex = indSceneList;
        DMfeatures(iMovie).features = tempS; 
        
        setIndFrame_valid = cat(1, setIndFrame_valid, validLocHigh(~isnan(indSceneList)));
        
        
        %         setFrame=[]; matResp=[]; tempS = [];
        %         for iScene = 1:length(validUniqueIndScene)
        %             tempFrames = validLocHigh(ic == iScene); % in case that there are multiple frames from one scene
        %             [tempResp, indMax] = max(meanFRCluster_zscore(tempFrames, iClust));
        %
        %             setFrame(iScene, 1) = tempFrames(indMax);
        %             matResp(iScene, 1) = max(tempResp);
        %
        %             tempS(iScene, 1) = epoch(validUniqueIndScene(iScene)).notes.face.A;
        %             tempS(iScene, 2) = epoch(validUniqueIndScene(iScene)).notes.torso.A;
        %             tempS(iScene, 3)= epoch(validUniqueIndScene(iScene)).notes.arms.A;
        %             tempS(iScene, 4) = epoch(validUniqueIndScene(iScene)).notes.legs.A;
        %             tempS(iScene, 5) = epoch(validUniqueIndScene(iScene)).notes.viewAngle;
        
    end
                
        
%         paramFeatureAnalysis(iMovie).sceneIndex = uniqueIndScene;
        
%         indFrameSelected = cat(1, indFrameSelected, setFrame);
%         avgClusterResponse = cat(1, avgClusterResponse, matResp);
        
        
    
    Clustering_movieFeatures(newIndCluster(iClust)).setIndFrame = setIndFrame_valid; %indFrameSelected = indFrameSelected;
    Clustering_movieFeatures(newIndCluster(iClust)).avgClusterResponse = meanFRCluster_zscore(setIndFrame_valid,iClust); % avgClusterResponse;
    
    Clustering_movieFeatures(newIndCluster(iClust)).DMfeatures = cat(1, DMfeatures.features);
%     Clustering_movieFeatures(newIndCluster(iClust)).sceneIndex = DMfeatures.sceneIndex;
    
    % Brian's features
    setBRFeatureVal_version1 = catFullRGR(setIndFrame_valid,:); %catFullRGR(indFrameSelected, :);     %catFullRGR(locHigh, :);  
    setBRFeatureVal_version2 = matRGRvalid_30fps(setIndFrame_valid, :); %matRGRvalid_30fps(indFrameSelected, :); 
    
    Clustering_movieFeatures(newIndCluster(iClust)).clusterID = newIndCluster(iClust);
    Clustering_movieFeatures(newIndCluster(iClust)).clusterID_old = iClust;
    Clustering_movieFeatures(newIndCluster(iClust)).BRfeatureParams = BRfeatureParams;    
    
    Clustering_movieFeatures(newIndCluster(iClust)).BRfeatures_version1.matRawFeatureValue_zCrit3 = setBRFeatureVal_version1;
    Clustering_movieFeatures(newIndCluster(iClust)).BRfeatures_version1.meanFeatureValue_zCrit3 = nanmean(setBRFeatureVal_version1);
    Clustering_movieFeatures(newIndCluster(iClust)).BRfeatures_version1.steFeatureValue_zCrit3 = nanstd(setBRFeatureVal_version1)./sqrt(size(setBRFeatureVal_version1, 1));

    Clustering_movieFeatures(newIndCluster(iClust)).BRfeatures_version2.matRawFeatureValue_zCrit3 = setBRFeatureVal_version2;
    Clustering_movieFeatures(newIndCluster(iClust)).BRfeatures_version2.meanFeatureValue_zCrit3 = nanmean(setBRFeatureVal_version2);
    Clustering_movieFeatures(newIndCluster(iClust)).BRfeatures_version2.steFeatureValue_zCrit3 = nanstd(setBRFeatureVal_version2)./sqrt(size(setBRFeatureVal_version2, 1));

end    
% save(fullfile(dirDataNeural, 'ClusterMovieFeatures_new5Means.mat'), 'Clustering_movieFeatures')
save(fullfile(dirDataNeural, 'ClusterMovieFeatures_new5Means_revised.mat'), 'Clustering_movieFeatures')







%%
% 2016/09/26: select unique scenes
Clustering_movieFeatures=struct([]);
newIndCluster = [5 4 1 2 3]; % [2 6 4 1 5 3 7]; % what K means cluster 1 is our cluster 2 %[3 2 7 6 5 1 4]; % mapping between original clustering #s to new numbering scheme for the paper (old cluster 1 is new cluster 3)
for iClust = 1:5 %7
    
    % Find when activity was high (i.e. exceeds certain criterion in z-score)
    critHigh = 3; %4; %3; %1.5; %2; % in z-score
    locHigh = find(meanFRCluster_zscore(:,iClust)>critHigh); 
  
    % DM's scene-based features (face, torso, arms, legs)
    % 2016/09/26 Select unique "scenes" and save the feature info +
    % responses at those time points
    avgClusterResponse=[]; indFrameSelected=[]; setIndFrame_valid=[];
    for iMovie = 1:3
        load(sprintf('/procdata/parksh/MovieRegressors/annotationMovie%d.mat', iMovie))
        sceneInfo = cat(2, sta', sto');
        
        validFrame_range = [(5*(iMovie-1)*60*30)+1, (5*iMovie*60*30)];
        validLocHigh = [];
        validLocHigh = locHigh(locHigh>validFrame_range(1) & locHigh<validFrame_range(2));
        
        % Extract scene information frame-by-frame from DM's scene analysis
        indSceneList=[];
        for iFrame = 1:length(validLocHigh)
            indFrame = validLocHigh(iFrame)-(5*(iMovie-1)*60*30);
            indScene = find(sum([sta<=indFrame; sto>=indFrame], 1)>1);                    

            if isempty(indScene)
                indSceneList(iFrame,1) = NaN;
                continue;
            end
            
            indSceneList(iFrame,1) = indScene;

        end
        
%         setIndFrame_valid = cat(1, setIndFrame_valid, validLocHigh(~isnan(indSceneList)));
        
        
        % Select unique scenes and save the feature info
        [uniqueIndScene, ia, ic] = unique(indSceneList, 'stable'); 
        % uniqueIndScene will have "unique scene indices" that I can use in
        % 'epoch' struct
        % ic will have the row numbers of unique scene indices as in
        % uniqueIndScene (so, length(uniqueIndScene) = max(ic) )
        
        validUniqueIndScene = uniqueIndScene(~isnan(uniqueIndScene));
        
        setFrame=[]; matResp=[]; tempS = [];
        for iScene = 1:length(validUniqueIndScene)
            tempFrames = validLocHigh(ic == iScene); % in case that there are multiple frames from one scene
            [tempResp, indMax] = max(meanFRCluster_zscore(tempFrames, iClust));
            
            setFrame(iScene, 1) = tempFrames(indMax);
            matResp(iScene, 1) = max(tempResp);
            
            tempS(iScene, 1) = epoch(validUniqueIndScene(iScene)).notes.face.A;
            tempS(iScene, 2) = epoch(validUniqueIndScene(iScene)).notes.torso.A;
            tempS(iScene, 3)= epoch(validUniqueIndScene(iScene)).notes.arms.A;
            tempS(iScene, 4) = epoch(validUniqueIndScene(iScene)).notes.legs.A;
            tempS(iScene, 5) = epoch(validUniqueIndScene(iScene)).notes.viewAngle;
            
        end
                
        DMfeatures(iMovie).sceneIndex = validUniqueIndScene;
        DMfeatures(iMovie).features = tempS; 
%         paramFeatureAnalysis(iMovie).sceneIndex = uniqueIndScene;
        
        indFrameSelected = cat(1, indFrameSelected, setFrame);
        avgClusterResponse = cat(1, avgClusterResponse, matResp);
        
    end    
    
    Clustering_movieFeatures(newIndCluster(iClust)).indFrameSelected = indFrameSelected;
    Clustering_movieFeatures(newIndCluster(iClust)).avgClusterResponse = avgClusterResponse;
    
    Clustering_movieFeatures(newIndCluster(iClust)).DMfeatures = cat(1, DMfeatures.features);
%     Clustering_movieFeatures(newIndCluster(iClust)).sceneIndex = DMfeatures.sceneIndex;
    
    % Brian's features
    setBRFeatureVal_version1 = catFullRGR(indFrameSelected, :);     %catFullRGR(locHigh, :);  
    setBRFeatureVal_version2 = matRGRvalid_30fps(indFrameSelected, :); 
    
    Clustering_movieFeatures(newIndCluster(iClust)).clusterID = newIndCluster(iClust);
    Clustering_movieFeatures(newIndCluster(iClust)).clusterID_old = iClust;
    Clustering_movieFeatures(newIndCluster(iClust)).BRfeatureParams = BRfeatureParams;    
    
    Clustering_movieFeatures(newIndCluster(iClust)).BRfeatures_version1.matRawFeatureValue_zCrit3 = setBRFeatureVal_version1;
    Clustering_movieFeatures(newIndCluster(iClust)).BRfeatures_version1.meanFeatureValue_zCrit3 = mean(setBRFeatureVal_version1);
    Clustering_movieFeatures(newIndCluster(iClust)).BRfeatures_version1.steFeatureValue_zCrit3 = std(setBRFeatureVal_version1)./sqrt(size(setBRFeatureVal_version1, 1));

    Clustering_movieFeatures(newIndCluster(iClust)).BRfeatures_version2.matRawFeatureValue_zCrit3 = setBRFeatureVal_version2;
    Clustering_movieFeatures(newIndCluster(iClust)).BRfeatures_version2.meanFeatureValue_zCrit3 = nanmean(setBRFeatureVal_version2);
    Clustering_movieFeatures(newIndCluster(iClust)).BRfeatures_version2.steFeatureValue_zCrit3 = nanstd(setBRFeatureVal_version2)./sqrt(size(setBRFeatureVal_version2, 1));

end    

save(fullfile(dirDataNeural, 'ClusterMovieFeatures_new5Means_sceneSelected.mat'), 'Clustering_movieFeatures')


% Plot
%     for iPart = 1:4 %face, torso, arms, legs
for iClust = 1:5 %7
    propPart(iClust,:) = sum(~isnan(Clustering_movieFeatures(iClust).DMfeatures))...
        ./size(Clustering_movieFeatures(iClust).DMfeatures, 1);
    sizePart_mean(iClust,:) = nanmean(Clustering_movieFeatures(iClust).DMfeatures);
    sizePart_ste(iClust,:) = nanstd(Clustering_movieFeatures(iClust).DMfeatures)...
        ./sqrt(sum(~isnan(Clustering_movieFeatures(iClust).DMfeatures)));
end
catDMfeatures = cat(1, Clustering_movieFeatures.DMfeatures);

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% Scene type
subplot(2,1,1);
bar(propPart(:,1:4))
legend('face', 'torso', 'arms', 'legs')
title('Scene type')
xlabel('Cluster #')
ylabel('Proportion')
% Size
subplot(2,1,2);
bar(sizePart_mean(:,1:4))
legend('face', 'torso', 'arms', 'legs')
title('Size')
xlabel('Cluster #')
ylabel('Size (deg^2)')

% Plot of Brian's features
matMeanFeatureVal = cat(1, Clustering_movieFeatures.meanFeatureValue_zCrit3);
% [matMeanFeatureVal_sort, indFeature]  = sortrows(matMeanFeatureVal');
matSteFeatureVal = cat(1, Clustering_movieFeatures.steFeatureValue_zCrit3);

featureNames = Clustering_movieFeatures(1).featureParams.featureNames;

% setRGR = [1 2 3 6 7]; %[1, 3, 9, 20, 21, 22, 25]; 
% 1: 'Luminance', 2: 'Contrast', 3: 'Motion (Speed)', 6: Low spatial Frequency 7: High spatial frequencty 3: 'Motion (Speed)', 

% Plot selected features for selected clusters
figClusterFeatures_compare = figure;
set(figClusterFeatures_compare, 'Color', 'w', 'PaperPositionMode', 'auto');%, 'Position', [100 100 780 195])
setCluster = [1 2 3 4 5]; %[1 2 3 4 5 6 7]; %[3 1 2 6]; %[1 2 3 4 5 6 7]; % [3 5 1 4];
barColor = [0.5 0.5 0.5];
% motion
SP(1) = subplot(1,3,1); 
hb(1) = bar(matMeanFeatureVal(setCluster, 3), 'FaceColor', barColor);
line([get(hb(1), 'XData'); get(hb(1), 'XData')],...
    [get(hb(1), 'YData')-matSteFeatureVal(setCluster,3)';get(hb(1), 'YData')+matSteFeatureVal(setCluster,3)'],...
    'Color', 'k', 'LineWidth', 2)
ylim([0 2.2])
ylabel('Motion (a.u.)')
% face size
SP(2) = subplot(1,3,2); 
hb(2) = bar(sizePart_mean(setCluster,1), 'FaceColor', barColor);
line([get(hb(2), 'XData'); get(hb(2), 'XData')],...
    [get(hb(2), 'YData')-sizePart_ste(setCluster,1)';get(hb(2), 'YData')+sizePart_ste(setCluster,1)'],...
    'Color', 'k', 'LineWidth', 2)
ylim([0 41])
ylabel('Face size (deg^2)')
% contrast
SP(3) =subplot(1,3,3); 
hb(3) = bar(matMeanFeatureVal(setCluster, 2), 'FaceColor', barColor);
line([get(hb(3), 'XData'); get(hb(3), 'XData')],...
    [get(hb(3), 'YData')-matSteFeatureVal(setCluster,2)';get(hb(3), 'YData')+matSteFeatureVal(setCluster,2)'],...
    'Color', 'k', 'LineWidth', 2)
ylabel('Contrast (a.u.)')
set(SP, 'Box', 'off', 'XLim', [0 length(setCluster)+1], 'LineWidth', 2, 'XTickLabel', {'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7'})
set(SP, 'TickDir', 'out', 'FontSize', 15, 'TickLength', [0.03 0.0250])

ylabel(SP(1), 'Motion (a.u.)')
ylabel(SP(2), 'Face size (deg^2)')
ylabel(SP(3), 'Contrast (a.u.)')
    

%% Correlation between Cluster time series and Movie regressors
% First, get the cell response in 4fps time resolution to compare it with movie regressors
FR_dT4fps = createCellRegressor_indMov_discreteTime(dirDataNeural, cellstr(paramCorr.validChanID),...
    setMovie, 0.25); % in 
% concatenate across movies
matFR4fps=[];
for iUnit = 1:size(FR_dT4fps,1)
    tempFR = cat(1, FR_dT4fps(iUnit, :).mnFR);
    matFR4fps(:,iUnit) = tempFR;
end

flagSM = 1; % flag for compression and smoothing

fullRGR4fps = createMovieRGR_4fps_indMov(setMovie, flagSM); %createFullMovieRegressors_4fps_indMov(setMovID); %
% ttt=load('/procdata/parksh/MovieRegressors/dbtmMriReg.mat'); % Face scale regressor (in TR unit)
% scaleRGR = ttt.reg.xx(7,:)';

% full regressors
catRGRfull_30fps=[];matRGRfull_30fps=[];
for iMov=1:length(setMovie)
    m = setMovie(iMov);
    matCurRGR = fullRGR4fps(m).smoRegressors; %fullRGR4fps(iMov).regressors(:,indValidRGR); %fullRGR4fps(iMov).regressors;
    catRGRfull_30fps = cat(1, catRGRfull_30fps, matCurRGR); % concatenation across movies
end
% scaleRGR_resampled = resample(scaleRGR, 2.4*100, 0.25*100); %matRGR = resample(catRGR, 0.25*100, 2.4*100);
matRGRfull_30fps = catRGRfull_30fps; %cat(2, catRGRfull, scaleRGR_resampled); %cat(2, matRGR, scaleRGR);
varnamesfull = fullRGR4fps(1).features; % cat(1, fullRGR4fps(1).features, {'Face size'});


% % Compute correlation
% [R_SUmovieRGRvalid] = corr(matRGRvalid, matFR4fps, 'rows', 'complete', 'type', 'Spearman');

% % Compute the average correlation for each cluster
% meanCorrCluster=[];steCorrCluster=[];
% oldIndCluster = [3 4 5 2 1]; %[4 1 6 3 5 2 7]; % [6 2 1 7 5 4 3]; % mapping between original clustering #s to new numbering scheme for the paper (new cluster 1 is old cluster 6)
% for iK = 1:sortTargetK
%     indClust = oldIndCluster(iK);
%     tempMatR=[];
%     tempMatR = R_SUmovieRGRvalid(:,indSortChan(sortedClust==indClust));
%     meanCorrCluster(:,iK) = mean(tempMatR, 2);
%     steCorrCluster(:,iK) = std(tempMatR, 0, 2)./length(indSortChan(sortedClust==indClust));
% end

% OR compute correlation between averaged TS in each cluster
% Compute the average correlation for each cluster
meanFRCluster4fps=[]; %steFRCluster4fps=[];
oldIndCluster =  [3 4 5 2 1]; %[4 1 6 3 5 2 7]; %[6 2 1 7 5 4 3];
for iK = 1:sortTargetK
    indClust = oldIndCluster(iK);
    tempMatFR=[];
    tempMatFR = matFR4fps(:,indSortChan(sortedClust==indClust));
    meanFRCluster4fps(:,iK) = mean(tempMatFR, 2);
%     steFRCluster4fps(:,iK) = std(matFR(:,tempMatFR, 0, 2)./length(indSortChan(sortedClust==iK));
end
% [R_ClusterMovieRGRvalid] = corr(matRGRvalid, meanFRCluster4fps, 'rows', 'complete', 'type', 'Spearman');


% Compute correlation between neural TS and feature TS
R_ClusterMovieRGRfull=NaN(size(matRGRfull_30fps,2), size(meanFRCluster4fps,2));
R_SUmovieRGRfull=NaN(size(matRGRfull_30fps,2), size(matFR4fps,2));
for iRGR = 1:size(matRGRfull_30fps,2)
    r_c=[]; r_su=[];
    
    % averaged TS in each cluster
    r_c = corr(matRGRfull_30fps(:,iRGR), meanFRCluster4fps, 'rows', 'complete', 'type', 'Spearman');
    % single unit TS
    r_su = corr(matRGRfull_30fps(:,iRGR), matFR4fps, 'rows', 'complete', 'type', 'Spearman');
    
    R_ClusterMovieRGRfull(iRGR, :) = r_c;
    R_SUmovieRGRfull(iRGR, :) = r_su;    
    
end

oldIndCluster = [3 4 5 2 1]; % [4 1 6 3 5 2 7]; 
indSortChan_new = [];
for iC = 1:sortTargetK %7
    curC = oldIndCluster(iC);
    tempind = indSortChan(sortedClust==curC);
    indSortChan_new = cat(1, indSortChan_new, tempind);
end

% subset of regressors
indValidRGR = [1 2 6 7 3 21 20 32 22 31 25]; %[1, 3, 9, 20, 21, 22, 25]; 
% 1: 'Luminance', 2: 'Contrast', 6: Low spatial Frequency 7: High spatial frequencty 3: 'Motion (Speed)', 
% 21: 'One face', 20: 'Number of faces', 32: 'Face size', 22: 'Body parts', 31: 'Hands', 25: 'Any animal'
% matRGRvalid = matRGRfull(:,indValidRGR);
varnamesvalid = varnamesfull(indValidRGR);

% R_SUmovieRGRvalid = R_SUmovieRGRfull(indValidRGR,:);
% R_ClusterMovieRGRvalid = R_ClusterMovieRGRfull(indValidRGR,:);

% correlation matrix between SUs and regressors
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 400 800])
imagesc(R_SUmovieRGRfull(indValidRGR, indSortChan_new)');
set(gca, 'XTick', 1:length(indValidRGR), 'XTickLabel', []) %varnamesfull)
set(gca, 'YTick', 1:length(paramCorr.validChanIndex), 'YTickLabel', cellstr(paramCorr.validChanID(indSortChan_new,:)))

% make blue-white-red colorbar
cval = 0.7;
cmin = -cval; cmax = cval;
colornum = 256;
colorInput = [1 0 0; 1 1 1; 0 0 1];
oldSteps = linspace(-1, 1, length(colorInput));
newSteps = linspace(-1, 1, colornum);
for j=1:3 % RGB
    newmap_all(:,j) = min(max(transpose(interp1(oldSteps, colorInput(:,j), newSteps)), 0), 1); 
end
endPoint = round((cmax-cmin)/2/abs(cmin)*colornum);
newmap = squeeze(newmap_all(1:endPoint, :));
figure(gcf)
set(gca, 'CLim', [cmin cmax])
colormap(flipud(newmap))
set(gca, 'TickDir', 'out')
box off
c=colorbar;
set(gca, 'FontSize', 12)
print(gcf, fullfile(dirFig, 'Corr_SUMovieRGR_movie123Tor_orderedClusterNew_5clusters'), '-dtiff', '-r150')
% print(gcf, fullfile(dirFig, 'Corr_SUMovieRGR_movie123Tor_rgrSubset'), '-dtiff', '-r150')




% Plot
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1120 320])
hBars = bar(R_ClusterMovieRGRvalid');
% legend(varnamesvalid); %, 'location', 'Best')
ylim([-.6 .6])
set(gca, 'Box', 'off', 'LineWidth', 2, 'FontSize', 15, 'TickDir', 'out')
set(gca, 'YTick', -0.6:0.2:0.6, 'YTickLabel', -0.6:0.2:0.6)
ylabel('Correlation (r)')
xlabel('Cluster')

% Sort and get the idea of what is strongly correlated for each cluster
featureCluster=struct([]);
for iC = 1:size(R_ClusterMovieRGRvalid,2)
    [aa, ii]=sort(abs(R_ClusterMovieRGRvalid(:,iC)), 'descend');
    
    featureCluster(iC).nameFeature_descendAbsR = varnamesvalid(ii);
    featureCluster(iC).valCorr_descendAbsR = R_ClusterMovieRGRvalid(ii, iC);
end
    
    

% plot
cMap = [0 0 0; 230 159 0; 86 180 233; 0 158 115; 240 228 66; 0 114 178; 213 94 0; 204 121 167]./255;
marker = {'o', '*', 'x', 's', 'd', '+', '^'};
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto') %, 'Position', [100 100 200 320])
for iClust = 1:7 
% tempX=get(hPol(iClust), 'XData');
% tempY=get(hPol(iClust), 'YData');
    plot(R_ClusterMovieRGRvalid(:, iClust), 'o-', 'LineWidth', 2, 'Color', cMap(iClust,:),...
        'MarkerSize', 8, 'MarkerEdgeColor', cMap(iClust,:), 'MarkerFaceColor', cMap(iClust,:)); %'Marker', marker{iClust}, 
    hold on;
end
line([0 12], [0 0], 'LineStyle', '--', 'Color', 'k')
ylim([-.6 .6])
set(gca, 'TickDir', 'out', 'LineWidth', 2)
box off
set(gca, 'FontSize', 15)
set(gca, 'XTick', 1:11)

% polar plot
testPolar_rho = abs(R_ClusterMovieRGRvalid);
testPolar_rho = cat(1, testPolar_rho, testPolar_rho(1,:));   % make it circular
tempTheta = 0:2*pi./size(R_ClusterMovieRGRvalid, 1):2*pi;
testPolar_theta = repmat(tempTheta', 1, 7);
% just for plot
limit = (pi/180)*2;
randroom_theta = linspace(-limit, limit, 7);
testPolar_theta = testPolar_theta + repmat(randroom_theta, 12, 1);

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100, 300, 600, 600])
hPol = polar(testPolar_theta, testPolar_rho); hold on;
hLine = findall(gcf, 'type', 'line');
delete(hLine(8:13))
t = findall(gcf, 'type', 'text');
delete(t)
set(hPol(:), 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 20)
set(hLine(14), 'LineWidth', 1.5)
legend('', '', '', '', '', '', '', '', 'Location', 'northeastoutside', 'boxoff')
legend boxoff
% mark negatvie correlation with open symbol
locNegCorr = R_ClusterMovieRGRvalid<0;
for iClust = 1:7 
tempX=get(hPol(iClust), 'XData');
tempY=get(hPol(iClust), 'YData');
plot(tempX(locNegCorr(:,iClust)), tempY(locNegCorr(:,iClust)), 'o', 'LineWidth', 2,...
    'MarkerEdgeColor', get(hPol(iClust), 'Color'), 'MarkerFaceColor', 'w')
end
% print(gcf, fullfile(dirFig, 'Cluster_avgClusterCorrMovieRGR_abs_polar'), '-depsc')

%% Another polar plot
testPolar_rho = R_ClusterMovieRGRvalid+1;  %abs(R_ClusterMovieRGRvalid);
testPolar_rho = cat(1, testPolar_rho, testPolar_rho(1,:));   % make it circular
tempTheta = 0:2*pi/size(R_ClusterMovieRGRvalid, 1):2*pi;
testPolar_theta = repmat(tempTheta', 1, 7);
% just for plot (to make the points at same axis of angle not overlap to
% each other, add some randomness in theta)
limit = (pi/180)*2;
randroom_theta = linspace(-limit, limit, 7);
testPolar_theta = testPolar_theta + repmat(randroom_theta, 12, 1);

fig5e=figure;
set(fig5e, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100, 300, 600, 600])
set(gca, 'NextPlot', 'replacechildren')
set(gca, 'ColorOrder', cMap)
hPol = polar(testPolar_theta, testPolar_rho); hold on;
hLine = findall(gcf, 'type', 'line');
delete(hLine(8:end))
t = findall(gcf, 'type', 'text');
delete(t)
set(hPol(:), 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 20)

% draw axes
rholim = 0.6;
rhocenter = 1;
% lower axis
axis0_rho = ones(1,200).*rhocenter;
axis0_theta = linspace(0, 2*pi, 200);
hAxis(1) = polar(axis0_theta, axis0_rho, ':');

axis1_rho = ones(1,200).*(rhocenter-rholim);
axis1_theta = linspace(0, 2*pi, 200);
hAxis(2) = polar(axis1_theta, axis1_rho, '-');

axis2_rho = ones(1,200).*(rhocenter+rholim);
axis2_theta = linspace(0, 2*pi, 200);
hAxis(3) = polar(axis2_theta, axis2_rho, '-');

set(hAxis, 'Color', 'k', 'LineWidth', 1.5)

legend('', '', '', '', '', '', '', 'Location', 'northeastoutside')
legend boxoff

%% Another polar plot: hedgehog style
rhocenter = .45; %.6; %1;

testPolar_rho = R_ClusterMovieRGRvalid+rhocenter;  %abs(R_ClusterMovieRGRvalid);
testPolar_rho = cat(1, testPolar_rho, testPolar_rho(1,:));   % make it circular
tempTheta = 0:2*pi/size(R_ClusterMovieRGRvalid, 1):2*pi;
testPolar_theta = repmat(tempTheta', 1, 7);
% % just for plot (to make the points at same axis of angle not overlap to
% % each other, add some randomness in theta)
% limit = (pi/180)*2;
% randroom_theta = linspace(-limit, limit, 7);
% testPolar_theta = testPolar_theta + repmat(randroom_theta, 12, 1);

cMap_rgrs = jet(11).*0.8;

posStart_rho = ones(1,12).*rhocenter;
for iC=1:7;
    curPoint_rho = testPolar_rho(:,iC)';
    curPoint_theta = testPolar_theta(:,iC)';
    posStart_theta = curPoint_theta;
    
    figtemp=figure;
    set(figtemp, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100, 300, 450, 450])
    set(gca, 'NextPlot', 'replacechildren')
    set(gca, 'ColorOrder', cMap_rgrs)
    
    [xc, yc] = pol2cart([posStart_theta; curPoint_theta], [posStart_rho; curPoint_rho]);
    L = line(xc, yc); hold on;
    set(L, 'LineWidth', 3);
    
%     hAxis=[];
    axis0_rho = ones(1,200).*rhocenter;
    axis0_theta = linspace(0, 2*pi, 200);
    [xAxis0, yAxis0] = pol2cart(axis0_theta, axis0_rho);
    hAxis(1) = line(xAxis0, yAxis0);
%     hAxis(1) = polar(axis0_theta, axis0_rho, ':');
    
%     axis1_rho = ones(1,200).*(rhocenter-rholim);
%     axis1_theta = linspace(0, 2*pi, 200);
%     [xAxis1, yAxis1] = pol2cart(axis1_theta, axis1_rho);
%     hAxis(2) = line(xAxis1, yAxis1);
% %     hAxis(2) = polar(axis1_theta, axis1_rho, '-');
    
    axis2_rho = ones(1,200).*(rhocenter+rholim);
    axis2_theta = linspace(0, 2*pi, 200);
    [xAxis2, yAxis2] = pol2cart(axis2_theta, axis2_rho);
    hAxis(2) = line(xAxis2, yAxis2);
%     hAxis(3) = polar(axis2_theta, axis2_rho, '-');

    set(hAxis, 'Color', ones(1,3).*0.6, 'LineWidth', 1.5)    

%     hout = polar(curPoint_theta, curPoint_rho);
    [xout, yout] = pol2cart(curPoint_theta, curPoint_rho);
    hout = line(xout, yout);
    set(hout, 'Color', 'k', 'LineWidth', 1.5)
    xlim([-1 1].*(rhocenter+rholim))
    ylim([-1 1].*(rhocenter+rholim))
    
    axis square
    axis off
    
    print(gcf, fullfile(dirFig, sprintf('Cluster%d_avgClusterCorrMovieRGR_polarHedgeHog_outerline_tight', iC)), '-depsc')
    
end
%     hPol = polar([posStart_theta; curPoint_theta], [posStart_rho; curPoint_rho]);
%     hold on;
%     hLine = findall(gcf, 'type', 'line');
%     delete(hLine(13:end))
%     t = findall(gcf, 'type', 'text');
%     delete(t)
%     set(hPol(:), 'LineWidth', 3) %, 'Marker', '.', 'MarkerSize', 20)
    
%     % zero axis
%     axis0_rho = ones(1,200).*rhocenter;
%     axis0_theta = linspace(0, 2*pi, 200);
%     hAxis(1) = polar(axis0_theta, axis0_rho, ':');
%     
%     axis1_rho = ones(1,200).*(rhocenter-rholim);
%     axis1_theta = linspace(0, 2*pi, 200);
%     hAxis(2) = polar(axis1_theta, axis1_rho, '-');
%     
%     axis2_rho = ones(1,200).*(rhocenter+rholim);
%     axis2_theta = linspace(0, 2*pi, 200);
%     hAxis(3) = polar(axis2_theta, axis2_rho, '-');
    
%     set(hAxis, 'Color', ones(1,3).*0.6, 'LineWidth', 1.5)
%     xlim([-2 2])
%     
%     hold on
%     hout = polar(curPoint_theta, curPoint_rho);
%     set(hout, 'Color', 'k', 'LineWidth', 1.5)
% %     ylim([-2 2])
        
% end


% hBaseline = get(hBars,'BaseLine');
% set(hBaseline,'LineStyle',':',...
%               'Color','k',...
%               'LineWidth',2);

% figCorrFeatures = figure;
% set(figCorrFeatures, 'Color', 'w', 'PaperPositionMode', 'auto')
% setClust = [3 5 1 4]; % for comparison
%  for iK = 1:length(setClust)
% subplot(length(setClust), 1, iK);
% bar(R_ClusterMovieRGRvalid(:,setClust(iK)))
% ylim([-0.4 0.4])
% hold on;
% line([0 12], [0 0], 'Color', 'k', 'LineStyle', '--')
% box off
% set(gca, 'XTickLabel', [], 'TickDir', 'out')
% end




            