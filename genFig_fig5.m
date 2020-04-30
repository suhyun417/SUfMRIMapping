% genFig_fig5.m
%
% 5a: methods of scene reverse correlation
% 5bcd: bar graph of different movie features
% 5e: comprehensive correlations of different features

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

% Set directories 
nameSubjNeural = 'Tor';
nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataHome = fullfile(dirProcdata, 'parksh');
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
matIndClust_SU = cat(2, Clustering_moviemask.resultKMeans.SU_indCluster); % based on movie-driven-masked corr map 
% matIndClust_SU = cat(2, Clustering.resultKMeans.SU_indCluster); % based on corr map

sortTargetK = 5; %4; %7;
[sortedClust, indSortChan]=sort(matIndClust_SU(:,sortTargetK-1));


%% Prepare dataset: Single unit time series, movie feature time series, 
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
% Average SDF in each cluster across cells
meanFRCluster=[];steFRCluster=[];
oldIndCluster =  [3 4 5 2 1]; %[1 2 3 4 5 6 7]; %[1 3 2 5 6 7 4]; %[4 1 6 3 5 2 7]; %
for iK = 1:sortTargetK
    indClust = oldIndCluster(iK);
    tempMatFR=[];
    tempMatFR = matFR(:,indSortChan(sortedClust==indClust));
    meanFRCluster(:,iK) = mean(tempMatFR, 2);
    steFRCluster(:,iK) = std(tempMatFR, 0, 2)./length(indSortChan(sortedClust==indClust));
end
meanFRCluster_zscore = zscore(meanFRCluster);
% matFR_zscore = zscore(matFR);

% In 4fps time resolution to compare it with movie regressors
FR_dT4fps = createCellRegressor_indMov_discreteTime(dirDataNeural, cellstr(paramCorr.validChanID),...
    setMovie, 0.25); % in 
% concatenate across movies
matFR4fps=[];
for iUnit = 1:size(FR_dT4fps,1)
    tempFR = cat(1, FR_dT4fps(iUnit, :).mnFR);
    matFR4fps(:,iUnit) = tempFR;
end
% Compute the average correlation for each cluster
meanFRCluster4fps=[]; %steFRCluster4fps=[];
oldIndCluster = [3 4 5 2 1]; %[1 2 3 4 5 6 7]; %[1 3 2 5 6 7 4]; %  [4 1 6 3 5 2 7]; %[6 2 1 7 5 4 3];
for iK = 1:sortTargetK
    indClust = oldIndCluster(iK);
    tempMatFR=[];
    tempMatFR = matFR4fps(:,indSortChan(sortedClust==indClust));
    meanFRCluster4fps(:,iK) = mean(tempMatFR, 2);
%     steFRCluster4fps(:,iK) = std(matFR(:,tempMatFR, 0, 2)./length(indSortChan(sortedClust==iK));
end
meanFRCluster4fps_zscore = zscore(meanFRCluster4fps);

% Movie regressors
flagSM = 1; % flag for compression and smoothing

fullRGR4fps = createMovieRGR_4fps_indMov(setMovie, flagSM); %createFullMovieRegressors_4fps_indMov(setMovID); %
% ttt=load('/procdata/parksh/MovieRegressors/dbtmMriReg.mat'); % Face scale regressor (in TR unit)
% scaleRGR = ttt.reg.xx(7,:)';

% full regressors
catRGRfull=[];matRGRfull=[];
for iMov=1:length(setMovie)
    m = setMovie(iMov);
    matCurRGR = fullRGR4fps(m).smoRegressors; %fullRGR4fps(iMov).regressors(:,indValidRGR); %fullRGR4fps(iMov).regressors;
    catRGRfull = cat(1, catRGRfull, matCurRGR); % concatenation across movies
end
% scaleRGR_resampled = resample(scaleRGR, 2.4*100, 0.25*100); %matRGR = resample(catRGR, 0.25*100, 2.4*100);
matRGRfull = catRGRfull; %cat(2, catRGRfull, scaleRGR_resampled); %cat(2, matRGR, scaleRGR);
varnamesfull = fullRGR4fps(1).features; % cat(1, fullRGR4fps(1).features, {'Face size'});

% subset of regressors
indValidRGR = [1 2 6 7 3 21 20 32 22 31 25]; %[1, 3, 9, 20, 21, 22, 25]; 
% 1: 'Luminance', 2: 'Contrast', 6: Low spatial Frequency 7: High spatial frequencty 3: 'Motion (Speed)', 
% 21: 'One face', 20: 'Number of faces', 32: 'Face size', 22: 'Body parts', 31: 'Hands', 25: 'Any animal'
% matRGRvalid = matRGRfull(:,indValidRGR);
varnamesvalid = varnamesfull(indValidRGR);

% Compute correlation between neural TS and feature TS
R_ClusterMovieRGRfull=NaN(size(matRGRfull,2), size(meanFRCluster4fps,2));
R_SUmovieRGRfull=NaN(size(matRGRfull,2), size(matFR4fps,2));
for iRGR = 1:size(matRGRfull,2)
    r_c=[]; r_su=[];
    
    % averaged TS in each cluster
    r_c = corr(matRGRfull(:,iRGR), meanFRCluster4fps, 'rows', 'complete', 'type', 'Spearman');
    % single unit TS
    r_su = corr(matRGRfull(:,iRGR), matFR4fps, 'rows', 'complete', 'type', 'Spearman');
    
    R_ClusterMovieRGRfull(iRGR, :) = r_c;
    R_SUmovieRGRfull(iRGR, :) = r_su;    
    
end

R_SUmovieRGRvalid = R_SUmovieRGRfull(indValidRGR,:);
R_ClusterMovieRGRvalid = R_ClusterMovieRGRfull(indValidRGR,:);

% [R_ClusterMovieRGRfull] = corr(matRGRfull, meanFRCluster4fps, 'rows', 'complete', 'type', 'Spearman');
% [R_SUmovieRGRfull] = corr(matRGRfull, matFR4fps, 'rows', 'complete', 'type', 'Spearman');

%% Fig 5a: Feature correlation method
iClust = 1; %1; %oldIndCluster(1);
cMap_rgrs = jet(11).*0.8;

% Averaged time series of one example cell group
fig5a1 = figure;
set(fig5a1, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [683 606 1000 200])
plot(meanFRCluster4fps_zscore(:, iClust), 'k-', 'LineWidth', 2); hold on;
% plot(meanFRCluster_zscore(1:9000, iClust)-steFRCluster(1:9000, iClust), '-', 'LineWidth', 1, 'Color', ones(1,3).*0.5)
% line([0 3600], [critHigh critHigh], 'LineStyle', ':', 'Color', 'k')

set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 15)
set(gca, 'XLim', [0 3600], 'XTick', 0:240:3600, 'XTickLabel', [])
set(gca, 'YLim', [-2 7], 'YTick', [-2 0 7], 'YTickLabel', [])

print(fig5a1, fullfile(dirFig, 'fig5a_1'), '-depsc')


% Movie feature time courses
matRGRvalid = matRGRfull(:,indValidRGR);
indFeature = [1 2 5 7 8];

for iR = 1:length(indFeature)
    id = indFeature(iR);
    
    fig5a2 = figure;
    set(fig5a2, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [683 606 1000 150])
    plot(matRGRvalid(:,id), 'Color', cMap_rgrs(id,:), 'LineWidth', 2)
    axis tight
    set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 15)
    set(gca, 'XLim', [0 3600], 'XTick', 0:240:3600, 'XTickLabel', [])
    set(gca, 'YTick', [])
    print(fig5a2, sprintf(fullfile(dirFig, 'fig5a_2_%s'), char(varnamesvalid(id))), '-depsc')

%     set(gca, 'YTick', [],'XTick', 1:400:3600, 'XTickLabel', [], 'TickDir', 'out')
%     print(gcf, sprintf(fullfile(dirFig, 'movieRGR_%s'), char(varnames(id))), '-depsc')
%     
% 
end


%% FIg 5b: Correlation between Cluster time series and Movie regressors
cMap_rgrs = jet(11).*0.8;
rhocenter = .45; %.6; %1;
rholim = 0.7;

testPolar_rho = R_ClusterMovieRGRvalid+rhocenter;  %abs(R_ClusterMovieRGRvalid);
testPolar_rho = cat(1, testPolar_rho, testPolar_rho(1,:));   % make it circular
tempTheta = 0:2*pi/size(R_ClusterMovieRGRvalid, 1):2*pi;
testPolar_theta = repmat(tempTheta', 1, 7);
% % just for plot (to make the points at same axis of angle not overlap to
% % each other, add some randomness in theta)
% limit = (pi/180)*2;
% randroom_theta = linspace(-limit, limit, 7);
% testPolar_theta = testPolar_theta + repmat(randroom_theta, 12, 1);

posStart_rho = ones(1,12).*rhocenter;
for iC=1:sortTargetK
    curPoint_rho = testPolar_rho(:,iC)';
    curPoint_theta = testPolar_theta(:,iC)';
    posStart_theta = curPoint_theta;
    
    fig5e=figure;
    set(fig5e, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100, 300, 450, 450])
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

    set(hAxis(1), 'Color', ones(1,3).*0.6, 'LineWidth', 1.5)  
    set(hAxis(2), 'Color', ones(1,3), 'LineWidth', 1.5)   

%     hout = polar(curPoint_theta, curPoint_rho);
    [xout, yout] = pol2cart(curPoint_theta, curPoint_rho);
    hout = line(xout, yout);
    set(hout, 'Color', 'k', 'LineWidth', 2)
    xlim([-1 1].*(rhocenter+rholim))
    ylim([-1 1].*(rhocenter+rholim))
    
    axis square
    axis off
    
    print(gcf, fullfile(dirFig, sprintf('fig5b_Cluster%d_5Means_usingValidVox', iC)), '-depsc')
    
end


%% Fig 5c: Reverse correlation method

% Averaged time series of one example cell group
critHigh = 3; %4; %3; %1.5; %2; % in z-score
for iClust = [1 4]; %[1 6]; %1; %oldIndCluster(1);   
    
    locHigh = meanFRCluster_zscore(:,iClust)>critHigh; % find(meanFRCluster_zscore(:,iClust)>critHigh);
    tempLocs = ones(size(locHigh)).*critHigh;
    tempLocs(locHigh) = meanFRCluster_zscore(locHigh, iClust);
%     tempLocs = [meanFRCluster_zscore(locHigh, iClust)'; ones(size(locHigh')).*critHigh];
%     tempLocs = tempLocs(:);
    
    % Averaged time series of one example cell group
    fig5c1 = figure;
    set(fig5c1, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [683 606 1000 200])
    plot(meanFRCluster_zscore(:, iClust), 'Color', ones(1,3).*0.7, 'LineWidth', 2); hold on;
    plot(tempLocs, 'k-', 'LineWidth', 2)
%     plot(meanFRCluster4fps_zscore(:, iClust), 'Color', ones(1,3).*0.7, 'LineWidth', 2); hold on;
%     plot(meanFRCluster4fps_zscore(meanFRCluster4fps_zscore(:,iClust)>critHigh, iClust), 'k-', 'LineWidth', 2)
    line([0 27000], [critHigh critHigh], 'LineStyle', ':', 'Color', 'w', 'LineWidth', 2)

    set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 15)
    set(gca, 'XLim', [0 27000], 'XTick', 0:3000:27000, 'XTickLabel', [])
    set(gca, 'YLim', [-2 7], 'YTick', [-2 0 7], 'YTickLabel', [])
    if iClust==4 %6
        set(gca, 'YLim', [-1.5 11], 'YTick', [-1.5 0 11], 'YTickLabel', [])
    end
    
    print(fig5c1, fullfile(dirFig, sprintf('fig5c_TScluster%d_5Means_usingValidVox', iClust)), '-depsc')
    
end



% Example image frames
dirMovie = '/procdata/parksh/Movies';
videoObj1 = VideoReader(fullfile(dirMovie, 'Movie1.avi'));
videoObj2 = VideoReader(fullfile(dirMovie, 'Movie2.avi'));
videoObj3 = VideoReader(fullfile(dirMovie, 'Movie3.avi'));

iClust = 1; %6; %1;
locHigh = find(meanFRCluster_zscore(:,iClust)>critHigh);
for iMovie = 1:3
    load(sprintf('/procdata/parksh/MovieRegressors/annotationMovie%d.mat', iMovie))
    sceneInfo = cat(2, sta', sto');
    validFrame_range = [(5*(iMovie-1)*60*30)+1, (5*iMovie*60*30)];
    validLocHigh = [];
    validLocHigh = locHigh(locHigh>validFrame_range(1) & locHigh<validFrame_range(2));
    indFrame = validLocHigh-(5*(iMovie-1)*60*30);
%     setIndFrame = indFrame(find(indFrame>0 & indFrame< 1800));
    
    % dirMovie = '/procdata/parksh/Movies';
    % obj = VideoReader(fullfile(dirMovie, sprintf('Movie%d.avi', iMovie)));
    switch iMovie
        case 1
            obj = videoObj1;
        case 2
            obj = videoObj2;
        case 3
            obj = videoObj3;
    end
    
    for iFrame = 1:length(indFrame)
        curFrame = indFrame(iFrame);
        frame = read(obj, curFrame);
        a(1).cdata = frame;
        a(1).colormap = [];
        imageFrame = frame2im(a);
        figure(100);
        set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
        imagesc(imageFrame);
        axis off
        title(sprintf('iFrame: %d, curFrame: %d', iFrame, curFrame))
%         input('')
            print(gcf, fullfile(dirFig, sprintf('exampleFrames_Cluster%d_movie%d_frame%d', iClust, iMovie, curFrame)), '-dtiff', '-r150')
    end
    
end


%% Fig 5d
load(fullfile(dirDataNeural, 'ClusterMovieFeatures_new5Means_revised.mat'), 'Clustering_movieFeatures')
% load(fullfile(dirDataNeural, 'ClusterMovieFeatures_new5Means_sceneSelected.mat'), 'Clustering_movieFeatures')
% load(fullfile(dirDataNeural, 'ClusterMovieFeatures_new5Means.mat'), 'Clustering_movieFeatures')
%load(fullfile(dirDataNeural, 'ClusterMovieFeatures.mat'), 'Clustering_movieFeatures')

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

% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% % Scene type
% subplot(2,1,1);
% bar(propPart)
% legend('face', 'torso', 'arms', 'legs')
% title('Scene type')
% xlabel('Cluster #')
% ylabel('Proportion')
% % Size
% subplot(2,1,2);
% bar(sizePart_mean)
% legend('face', 'torso', 'arms', 'legs')
% title('Size')
% xlabel('Cluster #')
% ylabel('Size (deg^2)')

% Plot of Brian's features
BRfeatures_version1 = cat(1, Clustering_movieFeatures.BRfeatures_version1);
matMeanFeatureVal_version1 = cat(1, BRfeatures_version1.meanFeatureValue_zCrit3);
matSteFeatureVal_version1 = cat(1, BRfeatures_version1.steFeatureValue_zCrit3);
% 
% BRfeatures_version2 = cat(1, Clustering_movieFeatures.BRfeatures_version2);
% matMeanFeatureVal_version2 = cat(1, BRfeatures_version2.meanFeatureValue_zCrit3);
% matSteFeatureVal_version2 = cat(1, BRfeatures_version2.steFeatureValue_zCrit3);

% matMeanFeatureVal = cat(1, Clustering_movieFeatures.BRfeatures_version1.meanFeatureValue_zCrit3);
% % [matMeanFeatureVal_sort, indFeature]  = sortrows(matMeanFeatureVal');
% matSteFeatureVal = cat(1, Clustering_movieFeatures.steFeatureValue_zCrit3);

featureNames_version1 = Clustering_movieFeatures(1).BRfeatureParams.version1.featureNames;
% featureNames_version2 = Clustering_movieFeatures(1).BRfeatureParams.version2.featureNames;

meanFeatureComp=[];
for iC = 1:5
    [aa, sortedInd] = sort(Clustering_movieFeatures(iC).avgClusterResponse, 'descend');
    meanFeatureComp(iC).setMotionVal_top300 = Clustering_movieFeatures(iC).BRfeatures_version1.matRawFeatureValue_zCrit3(sortedInd(1:300), indMotion_version1);
    meanFeatureComp(iC).setFaceSizeVal_top300 = Clustering_movieFeatures(iC).BRfeatures_version1.matRawFeatureValue_zCrit3(sortedInd(1:300), indFaceSize_version1);
    meanFeatureComp(iC).setContrastVal_top300 = Clustering_movieFeatures(iC).BRfeatures_version1.matRawFeatureValue_zCrit3(sortedInd(1:300), indContrast_version1);

    meanFeatureComp(iC).set11Features_top300 = Clustering_movieFeatures(iC).BRfeatures_version2.matRawFeatureValue_zCrit3(sortedInd(1:300), :);

end

matMotionVal = cat(2, meanFeatureComp.setMotionVal_top300);
matFaceSizeVal = cat(2, meanFeatureComp.setFaceSizeVal_top300);
matContrastVal = cat(2, meanFeatureComp.setContrastVal_top300);


% Plot selected features for selected clusters
setCluster = [1 2 3 4 5]; %[1 2 3 4 5 6 7]; % [3 1 2 6]; %
% barColor = [0.5 0.5 0.5];
cMap_rgrs = jet(11).*0.8;
indRgrs = [5 8 2];
barw = 0.7;

% Fig 5d1: Motion
fig5d1 = figure;
set(fig5d1, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [600 600 350 125])
plot(nanmedian(matMotionVal), 'o', 'MarkerEdgeColor', cMap_rgrs(indRgrs(1),:),...
    'MarkerFaceColor', 'w',  'MarkerSize', 10, 'LineWidth', 2);
hold on
% line([1:5; 1:5], ...
%     [matMeanFeatureVal_version1(:, indMotion_version1)-matSteFeatureVal_version1(:, indMotion_version1) matMeanFeatureVal_version1(:, indMotion_version1)+matSteFeatureVal_version1(:, indMotion_version1)]', ...
%     'Color', cMap_rgrs(indRgrs(1),:), 'LineWidth', 2);
ylim([0 2])
ylabel('Motion (a.u.)')
% set(h_5b, 'EdgeColor', 'none', 'BarWidth', barw)
set(gca, 'Box', 'off', 'XLim', [0.5 length(setCluster)+.5], 'LineWidth', 2, 'XTick', 1:5, 'XTickLabel', []) % {'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7'})
set(gca, 'YTickLabel', [])
set(gca, 'TickDir', 'out', 'FontSize', 15, 'TickLength', [0.02 0.02])
% print(fig5d1, fullfile(dirFig, 'fig5d_1'), '-depsc')

% Fig 5d2: Face size
fig5d2 = figure;
set(fig5d2, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [600 600 350 125])
plot(nanmedian(matFaceSizeVal), 'o', 'MarkerEdgeColor', cMap_rgrs(indRgrs(2),:),...
    'MarkerFaceColor', 'w',  'MarkerSize', 10, 'LineWidth', 2);
hold on
% line([1:5; 1:5], ...
%     [matMeanFeatureVal_version1(:, indFaceSize_version1)-matSteFeatureVal_version1(:, indFaceSize_version1) matMeanFeatureVal_version1(:, indFaceSize_version1)+matSteFeatureVal_version1(:, indFaceSize_version1)]', ...
%     'Color', cMap_rgrs(indRgrs(2),:), 'LineWidth', 2);
ylim([0 50])
ylabel('Face size (deg^2)')
% set(h_5c, 'EdgeColor', 'none', 'BarWidth', barw)
set(gca, 'Box', 'off', 'XLim', [0.5 length(setCluster)+.5], 'LineWidth', 2, 'XTick',  1:5, 'XTickLabel', []) % {'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7'})
set(gca, 'YTickLabel', [])
set(gca, 'TickDir', 'out', 'FontSize', 15, 'TickLength', [0.02 0.02])
% print(fig5d2, fullfile(dirFig, 'fig5d_2'), '-depsc')

% Fig 5d3: Contrast
fig5d3 = figure;
set(fig5d3, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [600 600 350 125])
plot(nanmedian(matContrastVal), 'o', 'MarkerEdgeColor', cMap_rgrs(indRgrs(3),:),...
    'MarkerFaceColor', 'w',  'MarkerSize', 10, 'LineWidth', 2);
hold on
% line([1:5; 1:5], ...
%     [matMeanFeatureVal_version1(:, indContrast_version1)-matSteFeatureVal_version1(:, indContrast_version1) matMeanFeatureVal_version1(:, indContrast_version1)+matSteFeatureVal_version1(:, indContrast_version1)]', ...
%     'Color', cMap_rgrs(indRgrs(3),:), 'LineWidth', 2);
ylim([0 43])
ylabel('Contrast (a.u.)')
% set(h_5d, 'EdgeColor', 'none', 'BarWidth', barw)
set(gca, 'Box', 'off', 'XLim', [0.5 length(setCluster)+.5], 'LineWidth', 2, 'XTick',  1:5, 'XTickLabel', []) % {'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7'})
set(gca, 'YTickLabel', [])
set(gca, 'TickDir', 'out', 'FontSize', 15, 'TickLength', [0.02 0.02])
% print(fig5d3, fullfile(dirFig, 'fig5d_3'), '-depsc')

% % Fig 5d1: Motion
% indMotion_version1 = 3;
% fig5d1 = figure;
% set(fig5d1, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [600 600 350 125])
% plot(matMeanFeatureVal_version1(setCluster, indMotion_version1), 'o', 'MarkerEdgeColor', cMap_rgrs(indRgrs(1),:),...
%     'MarkerFaceColor', 'w',  'MarkerSize', 8, 'LineWidth', 2);
% hold on
% line([1:5; 1:5], ...
%     [matMeanFeatureVal_version1(:, indMotion_version1)-matSteFeatureVal_version1(:, indMotion_version1) matMeanFeatureVal_version1(:, indMotion_version1)+matSteFeatureVal_version1(:, indMotion_version1)]', ...
%     'Color', cMap_rgrs(indRgrs(1),:), 'LineWidth', 2);
% % h_5b = bar(matMeanFeatureVal_version1(setCluster, 3), 'FaceColor', cMap_rgrs(indRgrs(1), :));
% % line([get(h_5b, 'XData'); get(h_5b, 'XData')],...
% %     [get(h_5b, 'YData')-matSteFeatureVal_version1(setCluster,3)';get(h_5b, 'YData')+matSteFeatureVal_version1(setCluster,3)'],...
% %     'Color', 'k', 'LineWidth', 2)
% ylim([0 2.5])
% ylabel('Motion (a.u.)')
% % set(h_5b, 'EdgeColor', 'none', 'BarWidth', barw)
% set(gca, 'Box', 'off', 'XLim', [0 length(setCluster)+1], 'LineWidth', 2, 'XTick', []) % {'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7'})
% set(gca, 'TickDir', 'out', 'FontSize', 15, 'TickLength', [0.02 0.02])
% % print(fig5d1, fullfile(dirFig, 'fig5d_1'), '-depsc')
% 
% % Fig 5d2: Face size
% indFaceSize_version1 = 67;
% fig5d2 = figure;
% set(fig5d2, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [600 600 350 125])
% plot(matMeanFeatureVal_version1(setCluster, indFaceSize_version1), 'o', 'MarkerEdgeColor', cMap_rgrs(indRgrs(2),:),...
%     'MarkerFaceColor', 'w',  'MarkerSize', 8, 'LineWidth', 2);
% hold on
% line([1:5; 1:5], ...
%     [matMeanFeatureVal_version1(:, indFaceSize_version1)-matSteFeatureVal_version1(:, indFaceSize_version1) matMeanFeatureVal_version1(:, indFaceSize_version1)+matSteFeatureVal_version1(:, indFaceSize_version1)]', ...
%     'Color', cMap_rgrs(indRgrs(2),:), 'LineWidth', 2);
% % h_5c = bar(sizePart_mean(setCluster, 1), 'FaceColor', cMap_rgrs(indRgrs(2), :));
% % line([get(h_5c, 'XData'); get(h_5c, 'XData')],...
% %     [get(h_5c, 'YData')-sizePart_ste(setCluster,1)';get(h_5c, 'YData')+sizePart_ste(setCluster,1)'],...
% %     'Color', 'k', 'LineWidth', 2)
% ylim([0 40])
% ylabel('Face size (deg^2)')
% % set(h_5c, 'EdgeColor', 'none', 'BarWidth', barw)
% set(gca, 'Box', 'off', 'XLim', [0.5 length(setCluster)+.5], 'LineWidth', 2, 'XTick', []) % {'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7'})
% set(gca, 'TickDir', 'out', 'FontSize', 15, 'TickLength', [0.02 0.02])
% % print(fig5d2, fullfile(dirFig, 'fig5d_2'), '-depsc')
% 
% % Fig 5d3: Contrast
% indContrast_version1 = 2;
% fig5d3 = figure;
% set(fig5d3, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [600 600 350 125])
% plot(matMeanFeatureVal_version1(setCluster, indContrast_version1), 'o', 'MarkerEdgeColor', cMap_rgrs(indRgrs(3),:),...
%     'MarkerFaceColor', 'w',  'MarkerSize', 8, 'LineWidth', 2);
% hold on
% line([1:5; 1:5], ...
%     [matMeanFeatureVal_version1(:, indContrast_version1)-matSteFeatureVal_version1(:, indContrast_version1) matMeanFeatureVal_version1(:, indContrast_version1)+matSteFeatureVal_version1(:, indContrast_version1)]', ...
%     'Color', cMap_rgrs(indRgrs(3),:), 'LineWidth', 2);
% % h_5d = bar(matMeanFeatureVal(setCluster, 2), 'FaceColor', cMap_rgrs(indRgrs(3), :));
% % line([get(h_5d, 'XData'); get(h_5d, 'XData')],...
% %     [get(h_5d, 'YData')-matSteFeatureVal(setCluster,2)';get(h_5d, 'YData')+matSteFeatureVal(setCluster,2)'],...
% %     'Color', 'k', 'LineWidth', 2)
% ylim([0 43])
% ylabel('Contrast (a.u.)')
% % set(h_5d, 'EdgeColor', 'none', 'BarWidth', barw)
% set(gca, 'Box', 'off', 'XLim', [0 length(setCluster)+1], 'LineWidth', 2, 'XTick', []) % {'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7'})
% set(gca, 'TickDir', 'out', 'FontSize', 15, 'TickLength', [0.02 0.02])
% % print(fig5d3, fullfile(dirFig, 'fig5d_3'), '-depsc')


%% Comparison of mean feature values
% Compare everything!
indMotion_version1 = 3;
indFaceSize_version1 = 67;
indContrast_version1 = 2;

meanFeatureComp=[];
for iC = 1:5
    meanFeatureComp(iC).setMotionVal = Clustering_movieFeatures(iC).BRfeatures_version1.matRawFeatureValue_zCrit3(:, indMotion_version1);
    meanFeatureComp(iC).setFaceSizeVal = Clustering_movieFeatures(iC).BRfeatures_version1.matRawFeatureValue_zCrit3(:, indFaceSize_version1);
    meanFeatureComp(iC).setContrastVal = Clustering_movieFeatures(iC).BRfeatures_version1.matRawFeatureValue_zCrit3(:, indContrast_version1);

    meanFeatureComp(iC).set11Features = Clustering_movieFeatures(iC).BRfeatures_version2.matRawFeatureValue_zCrit3(:, :);
end
matMotionVal = cat(1, meanFeatureComp.setMotionVal);
matFaceSizeVal = cat(1, meanFeatureComp.setFaceSizeVal);
matContrastVal = cat(1, meanFeatureComp.setContrastVal);

groupVal = cat(1, ones(length(Clustering_movieFeatures(1).setIndFrame), 1), ones(length(Clustering_movieFeatures(2).setIndFrame), 1).*2,... 
    ones(length(Clustering_movieFeatures(3).setIndFrame), 1).*3, ...
    ones(length(Clustering_movieFeatures(4).setIndFrame), 1).*4, ones(length(Clustering_movieFeatures(5).setIndFrame), 1).*5);

[p_motion, tbl_motion, stats_motion] = kruskalwallis(matMotionVal, groupVal);
expVar_motion = cell2mat(tbl_motion(2, 2))/cell2mat(tbl_motion(4, 2));

[p_facesize, tbl_facesize, stats_facesize] = kruskalwallis(matFaceSizeVal, groupVal);
expVar_facesize = cell2mat(tbl_facesize(2, 2))/cell2mat(tbl_facesize(4, 2));

[p_contrast, tbl_contrast, stats_contrast] = kruskalwallis(matContrastVal, groupVal);
expVar_contrast = cell2mat(tbl_contrast(2, 2))/cell2mat(tbl_contrast(4, 2));


% % First, select 300 frames according to the response amplitude
% meanFeatureComp=[];
% for iC = 1:5
%     [aa, sortedInd] = sort(Clustering_movieFeatures(iC).avgClusterResponse, 'descend');
%     meanFeatureComp(iC).setMotionVal_top300 = Clustering_movieFeatures(iC).BRfeatures_version1.matRawFeatureValue_zCrit3(sortedInd(1:300), indMotion_version1);
%     meanFeatureComp(iC).setFaceSizeVal_top300 = Clustering_movieFeatures(iC).BRfeatures_version1.matRawFeatureValue_zCrit3(sortedInd(1:300), indFaceSize_version1);
%     meanFeatureComp(iC).setContrastVal_top300 = Clustering_movieFeatures(iC).BRfeatures_version1.matRawFeatureValue_zCrit3(sortedInd(1:300), indContrast_version1);
% 
%     meanFeatureComp(iC).set11Features_top300 = Clustering_movieFeatures(iC).BRfeatures_version2.matRawFeatureValue_zCrit3(sortedInd(1:300), :);
% 
% end
% 
% matMotionVal = cat(2, meanFeatureComp.setMotionVal_top300);
% matFaceSizeVal = cat(2, meanFeatureComp.setFaceSizeVal_top300);
% matContrastVal = cat(2, meanFeatureComp.setContrastVal_top300);
% 
% [p_motion, tbl_motion, stats_motion] = kruskalwallis(matMotionVal);
% expVar_motion = cell2mat(tbl_motion(2, 2))/cell2mat(tbl_motion(4, 2));
% 
% [p_facesize, tbl_facesize, stats_facesize] = kruskalwallis(matFaceSizeVal);
% expVar_facesize = cell2mat(tbl_facesize(2, 2))/cell2mat(tbl_facesize(4, 2));
% 
% [p_contrast, tbl_contrast, stats_contrast] = kruskalwallis(matContrastVal);
% expVar_contrast = cell2mat(tbl_contrast(2, 2))/cell2mat(tbl_contrast(4, 2));
% 
% % [p_motion, anovatab_motion, stats_motion] = anova1(matMotionVal); % 10.48 percent explain
% % [p_facesize, anovatab_facesize, stats_facesize] = anova1(matFaceSizeVal); % 37.13 percent explain
% % [p_contrast, anovatab_contrast, stats_contrast] = anova1(matContrastVal); % 3.38 percent explain
% 
% for iR = 1:11
%     [p, tbl, stats] = kruskalwallis(squeeze(tempcat(:,iR,:)));
%     ss(iR,1) = tbl(2, 2); % ss_columns
%     ss(iR,2) = tbl(4, 2); % ss_total
%     pVal(iR,1) = p;
%     setTBL(iR).tbl = tbl;
%     % input('')
end






% SP(3) =subplot(1,3,3); 
% hb(3) = bar(matMeanFeatureVal(setCluster, 2), 'FaceColor', barColor);
% line([get(hb(3), 'XData'); get(hb(3), 'XData')],...
%     [get(hb(3), 'YData')-matSteFeatureVal(setCluster,2)';get(hb(3), 'YData')+matSteFeatureVal(setCluster,2)'],...
%     'Color', 'k', 'LineWidth', 2)
% ylabel('Contrast (a.u.)')
% set(SP, 'Box', 'off', 'XLim', [0 length(setCluster)+1], 'LineWidth', 2, 'XTickLabel',  {'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7'})
% set(SP, 'TickDir', 'out', 'FontSize', 15, 'TickLength', [0.03 0.0250])
% 
% ylabel(SP(1), 'Motion (a.u.)')
% ylabel(SP(2), 'Face size (deg^2)')
% ylabel(SP(3), 'Contrast (a.u.)')
    




%% 
% % correlation matrix between SUs and regressors
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 400 800])
% imagesc(R_SUmovieRGRfull(indValidRGR, indSortChan_new)');
% set(gca, 'XTick', 1:size(R_SUmovieRGRvalid,1), 'XTickLabel', []) %varnamesfull)
% set(gca, 'YTick', 1:length(paramCorr.validChanIndex), 'YTickLabel', cellstr(paramCorr.validChanID(indSortChan_new,:)))
% 
% % make blue-white-red colorbar
% cval = 0.7;
% cmin = -cval; cmax = cval;
% colornum = 256;
% colorInput = [1 0 0; 1 1 1; 0 0 1];
% oldSteps = linspace(-1, 1, length(colorInput));
% newSteps = linspace(-1, 1, colornum);
% for j=1:3 % RGB
%     newmap_all(:,j) = min(max(transpose(interp1(oldSteps, colorInput(:,j), newSteps)), 0), 1); 
% end
% endPoint = round((cmax-cmin)/2/abs(cmin)*colornum);
% newmap = squeeze(newmap_all(1:endPoint, :));
% figure(gcf)
% set(gca, 'CLim', [cmin cmax])
% colormap(flipud(newmap))
% set(gca, 'TickDir', 'out')
% box off
% c=colorbar;
% set(gca, 'FontSize', 12)



% % Compute correlation
% [R_SUmovieRGRvalid] = corr(matRGRvalid, matFR4fps, 'rows', 'complete', 'type', 'Spearman');
% 
% % Compute the average correlation for each cluster
% meanCorrCluster=[];steCorrCluster=[];
% oldIndCluster = [4 1 6 3 5 2 7]; % [6 2 1 7 5 4 3]; % mapping between original clustering #s to new numbering scheme for the paper (new cluster 1 is old cluster 6)
% for iK = 1:sortTargetK
%     indClust = oldIndCluster(iK);
%     tempMatR=[];
%     tempMatR = R_SUmovieRGRvalid(:,indSortChan(sortedClust==indClust));
%     meanCorrCluster(:,iK) = mean(tempMatR, 2);
%     steCorrCluster(:,iK) = std(tempMatR, 0, 2)./length(indSortChan(sortedClust==indClust));
% end


% % Plot
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1120 320])
% hBars = bar(R_ClusterMovieRGRvalid');
% % legend(varnamesvalid); %, 'location', 'Best')
% ylim([-.6 .6])
% set(gca, 'Box', 'off', 'LineWidth', 2, 'FontSize', 15, 'TickDir', 'out')
% set(gca, 'YTick', -0.6:0.2:0.6, 'YTickLabel', -0.6:0.2:0.6)
% ylabel('Correlation (r)')
% xlabel('Cluster')
% 
% % plot
% cMap = [0 0 0; 230 159 0; 86 180 233; 0 158 115; 240 228 66; 0 114 178; 213 94 0; 204 121 167]./255;
% marker = {'o', '*', 'x', 's', 'd', '+', '^'};
% 
% 
% % 1. Polar plot
% testPolar_rho = R_ClusterMovieRGRvalid+1;  %abs(R_ClusterMovieRGRvalid);
% testPolar_rho = cat(1, testPolar_rho, testPolar_rho(1,:));   % make it circular
% tempTheta = 0:2*pi/size(R_ClusterMovieRGRvalid, 1):2*pi;
% testPolar_theta = repmat(tempTheta', 1, 7);
% % just for plot (to make the points at same axis of angle not overlap to
% % each other, add some randomness in theta)
% limit = (pi/180)*2;
% randroom_theta = linspace(-limit, limit, 7);
% testPolar_theta = testPolar_theta + repmat(randroom_theta, 12, 1);
% 
% fig5e=figure;
% set(fig5e, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100, 300, 600, 600])
% set(gca, 'NextPlot', 'replacechildren')
% set(gca, 'ColorOrder', cMap)
% hPol = polar(testPolar_theta, testPolar_rho); hold on;
% hLine = findall(gcf, 'type', 'line');
% delete(hLine(8:end))
% t = findall(gcf, 'type', 'text');
% delete(t)
% set(hPol(:), 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 20)
% 
% % draw axes
% rholim = 0.6;
% rhocenter = 1;
% % lower axis
% axis0_rho = ones(1,200).*rhocenter;
% axis0_theta = linspace(0, 2*pi, 200);
% hAxis(1) = polar(axis0_theta, axis0_rho, ':');
% 
% axis1_rho = ones(1,200).*(rhocenter-rholim);
% axis1_theta = linspace(0, 2*pi, 200);
% hAxis(2) = polar(axis1_theta, axis1_rho, '-');
% 
% axis2_rho = ones(1,200).*(rhocenter+rholim);
% axis2_theta = linspace(0, 2*pi, 200);
% hAxis(3) = polar(axis2_theta, axis2_rho, '-');
% 
% set(hAxis, 'Color', 'k', 'LineWidth', 1.5)
% 
% legend('', '', '', '', '', '', '', 'Location', 'northeastoutside')
% legend boxoff
% 
% % print(fig5e, fullfile(dirFig, 'Cluster_avgClusterCorrMovieRGR_polar'), '-depsc')
% 
% % 2. Overlapped plot
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto') %, 'Position', [100 100 200 320])
% for iClust = 1:7 
% % tempX=get(hPol(iClust), 'XData');
% % tempY=get(hPol(iClust), 'YData');
%     plot(R_ClusterMovieRGRvalid(:, iClust), 'o-', 'LineWidth', 2, 'Color', cMap(iClust,:),...
%         'MarkerSize', 8, 'MarkerEdgeColor', cMap(iClust,:), 'MarkerFaceColor', cMap(iClust,:)); %'Marker', marker{iClust}, 
%     hold on;
% end
% line([0 12], [0 0], 'LineStyle', '--', 'Color', 'k')
% ylim([-.6 .6])
% set(gca, 'TickDir', 'out', 'LineWidth', 2)
% box off
% set(gca, 'FontSize', 15)
% set(gca, 'XTick', 1:11)
% 
% % % 1-2. Polar plot with absolute corr: Mark negatvie correlation with open symbol
% % figure;
% % set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100, 300, 600, 600])
% % hPol = polar(testPolar_theta, testPolar_rho); hold on;
% % hLine = findall(gcf, 'type', 'line');
% % delete(hLine(8:13))
% % t = findall(gcf, 'type', 'text');
% % delete(t)
% % set(hPol(:), 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 20)
% % set(hLine(14), 'LineWidth', 1.5)
% % legend('', '', '', '', '', '', '', '', 'Location', 'northeastoutside', 'boxoff')
% % legend boxoff
% 
% % locNegCorr = R_ClusterMovieRGRvalid<0;
% % for iClust = 1:7 
% % tempX=get(hPol(iClust), 'XData');
% % tempY=get(hPol(iClust), 'YData');
% % plot(tempX(locNegCorr(:,iClust)), tempY(locNegCorr(:,iClust)), 'o', 'LineWidth', 2,...
% %     'MarkerEdgeColor', get(hPol(iClust), 'Color'), 'MarkerFaceColor', 'w')
% % end
% % print(gcf, fullfile(dirFig, 'Cluster_avgClusterCorrMovieRGR_abs_polar'), '-depsc')
% 


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

