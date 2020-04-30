function [fullRGR4fps] = createMovieRGR_4fps_indMov(setMovID, flagSM)
%
% Creates 31 movie regressors (i.e. feature time series), including 19 low-level and a subset of high-level
% regressors
%
% input
%   -- setMovID: set of movie IDs in vector
%   -- flagSM: flag for smoothing & compression
%
% First created on 10/24/14 by Soo Hyun Park
% 11/3/14 Compression and smoothing feature added by SHP
% 5/19/16 Added DM's scene based size regressors

% cd /projects/parksh/NeuralBOLD/analysis
% setMovID = [1 2 3]; %

flagfig = 0; %1; %0; % 1 to check the compression and smoothing 

dataDir ='/procdata/parksh/MovieRegressors/'; %'/Volumes/PROCDATA/parksh/MovieRegressors/'; %'/procdata/parksh/MovieRegressors/';
load(fullfile(dataDir, 'HighLevelRGR.mat'))


%% Regressor storage in fullRGR4fps(movID).regressors variable for each movie
nFramePerMovie = 4*300;
fullRGR4fps = struct([]);

for iMov = 1:length(setMovID);
    movID = setMovID(iMov);
    % movID=1; %:3; % loop
    
    fullRGR4fps(iMov).movieID = movID;
    
    % Low level RGRs
    fName_lowLevRGR = sprintf('Movie%d_10fps_rgr.mat', movID);
    tempStruct = load(fullfile(dataDir, fName_lowLevRGR));
    tempLowRGR = resample(tempStruct.RGR10fps.regressors', 4, 10);  % resample from 10hz to 4hz sampling rate
    
    fullRGR4fps(iMov).regressors = tempLowRGR; % frame x rgrs
    fullRGR4fps(iMov).features = tempStruct.RGR10fps.features';
    nLowLevRGR = size(fullRGR4fps(iMov).regressors,2);
    
    % % smoothing
    % indRGR_log = [1, 2, 3, 6:17, 19];
    % % tempInd = zeros(nFramePerMovie, length(tempStruct.RGR10fps.features));
    % % tempInd(:, indRGR_log) = 1;
    % fullRGR4fps(movID).regressors4fps(iMov).regressors(:,indRGR_log) = real(log(tempLowRGR(:,indRGR_log)));
    
    % High level RGRs
    %20 Faces (sqrt)
    %21 one face (sqrt)
    %22 body parts (sqrt)
    %23 conspecifics (sqrt)
    %24 humans (sqrt)
    %25 any animal (sqrt)
    %26 dyadic interaction
    %27 One Face (full) (sqrt)
    %28 One Face (profile) (sqrt)
    %29 Faces (full) (sqrt)
    %30 Faces (profile/side view) (sqrt)
    %31 Hands(sqrt_MION) 11
    highLevFeatures = {'Faces', 'One face', 'Body parts', 'Conspecifics', 'Humans', 'Any animal', 'Dyadic interaction',...
        'One Face (full)', 'One Face (side view)', 'Faces (full)', 'Faces (side view)', 'Hands'}';
    nHighLevRGR = length(highLevFeatures);
    highLevRGR4fps = zeros(nFramePerMovie, nHighLevRGR);
    
    for iRgr = 1:+nHighLevRGR
        
        switch iRgr
            case 1 %face
                nowdata = RAWRGR.rgrs(:,8) + RAWRGR.rgrs(:,9);
                highLevRGR4fps(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
                
            case 2 %one-face
                nowdata = RAWRGR.rgrs(:,46) + RAWRGR.rgrs(:,47); %1x3600 entries
                highLevRGR4fps(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
                
            case 3 %body parts
                nowdata = RAWRGR.rgrs(:,10) + RAWRGR.rgrs(:,11) + RAWRGR.rgrs(:,12);
                highLevRGR4fps(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
                
            case 4 %conspecifics
                nowdata = RAWRGR.rgrs(:,2);
                highLevRGR4fps(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
                
            case 5 %humans
                nowdata = RAWRGR.rgrs(:,6);
                highLevRGR4fps(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
                
            case 6 %any animal
                nowdata = RAWRGR.rgrs(:,2) + RAWRGR.rgrs(:,3) + RAWRGR.rgrs(:,4) + RAWRGR.rgrs(:,5) + RAWRGR.rgrs(:,6);
                highLevRGR4fps(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
                
            case 7 %dyadic interaction
                nowdata = RAWRGR.rgrs(:,33);
                highLevRGR4fps(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
                
            case 8 %One Face (full)
                nowdata = RAWRGR.rgrs(:,46);
                highLevRGR4fps(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
                
            case 9 %One Face (profile)
                nowdata = RAWRGR.rgrs(:,47);
                highLevRGR4fps(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
                
            case 10 %Faces (full)
                nowdata = RAWRGR.rgrs(:,8);
                highLevRGR4fps(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
                
            case 11 %Faces (profile)
                nowdata = RAWRGR.rgrs(:,9);
                highLevRGR4fps(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
                
            case 12 %Hands
                nowdata = RAWRGR.rgrs(:,11);
                highLevRGR4fps(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
                
        end
        
        fullRGR4fps(movID).regressors(:,nLowLevRGR+1:nLowLevRGR+nHighLevRGR) = highLevRGR4fps;
        fullRGR4fps(movID).features(nLowLevRGR+1:nLowLevRGR+nHighLevRGR) = highLevFeatures;
        
    end
    
    % DM's scene-based features (face, torso, arms, legs)
    load(sprintf('/procdata/parksh/MovieRegressors/annotationMovie%d.mat', movID))
    sceneInfo = cat(2, sta', sto');
    nScene = length(epoch);
    validFrame_range = [1 9000];  % in 30fps
    
    matSizeRGR = [];
    for iS = 1:nScene
        matSizeRGR(sta(iS):sto(iS),1) = epoch(iS).notes.face.A;
        matSizeRGR(sta(iS):sto(iS),2) = epoch(iS).notes.torso.A;
        matSizeRGR(sta(iS):sto(iS),3) = epoch(iS).notes.arms.A;
        matSizeRGR(sta(iS):sto(iS),4) = epoch(iS).notes.legs.A;
        matSizeRGR(sta(iS):sto(iS),5) = epoch(iS).notes.viewAngle;
    end
    
    if sto(nScene) < validFrame_range(2) % add NaN to the end of the regressor
        matPad = NaN(validFrame_range(2) - sto(nScene), 5);
        matSizeRGR = cat(1, matSizeRGR, matPad);
    elseif sto(nScene) > validFrame_range(2) % exclude extra frames
        matSizeRGR = matSizeRGR(validFrame_range(1):validFrame_range(2), :);        
    end
    
    tempDMRGR = resample(matSizeRGR, 4, 30);  % down sample in 4 hz
    
    varnames_DM = {'Face size', 'Torso size', 'Arm size', 'Leg size', 'View angle'};
    
    fullRGR4fps(movID).regressors(:,nLowLevRGR+nHighLevRGR+1:nLowLevRGR+nHighLevRGR+5) = tempDMRGR;
    fullRGR4fps(movID).features(nLowLevRGR+nHighLevRGR+1:nLowLevRGR+nHighLevRGR+5) = varnames_DM;
    
    %% Regressor storage in fullRGR4fps(movID).regressors variable
    if flagSM
        % [4, 5, 19]: binary
        % [1, 2, 3, 6:17, 19]: take log
        % [20:30]: sqrt
        %1 'Luminance' (log)
        %2 'Contrast' (log)
        %3 'Speed' (log)
        %4 'Rightward Motion' (binary)
        %5 'Leftward Motion' (binary)
        %6 'SF (0.0 to 0.2)' (log)
        %7 'SF (3.0 to 100.0)' (log)
        %8 'SF Ratio' (log)
        %9 'Beta Contrast' (log)
        %10 'Gamma Contrast' (log)
        %11 'Motion Contrast'
        %12 'Motion Dir Beta'
        %13 'Motion Dir Gamma'
        %14 'Motion Div Mean'
        %15 'Motion Div Rectify' (log)
        %16 'Motion Div STD' (log)
        %17 'Motion Div Beta' (log)
        %18 'Motion Div Gamma' (sqrt_sm2)
        %19 'Scene Cuts' (binary) (log_sm2)
        %20 Faces (sqrt)
        %21 one face (sqrt)
        %22 body parts (sqrt)
        %23 conspecifics (sqrt)
        %24 humans (sqrt)
        %25 any animal (sqrt)
        %26 dyadic interaction
        %27 One Face (full) (sqrt)
        %28 One Face (profile) (sqrt)
        %29 Faces (full) (sqrt)
        %30 Faces (profile/side view) (sqrt)
        %31 Hands(sqrt_MION)
        %32 Face size
        %33 Torso size
        %34 Arm size
        %35 Leg size
        %36 View angle
        
        % index for different ways to compress
        %         indBin = [4, 5, 19]; % binary
        indLog = [1, 2, 3, 6:17]; % use log to compress (usually in values)
        indSqrt = [18, 20:31]; % use sqrt to compress (usually in numbers)
        indDM = [32:36]; %DM's regressors
        
        % Set a Gaussian smoothing kernel
        smooth_k = normpdf([-50:50],0,2); % normpdf([-10:10], 0, 2);
        smooth_k = smooth_k/sum(smooth_k);        
        
        % do it one by one
        rgrs=[];
        for i = 1:length(fullRGR4fps(movID).features)
            rgr_org = fullRGR4fps(movID).regressors(:,i);
            
            % Compression
            if ismember(i, indLog)
                rgr_cp = real(log(rgr_org));
                rgr_cp(isinf(rgr_cp))=0;
            elseif ismember(i, indSqrt)
                rgr_cp = sqrt(rgr_org);                 % compress with a sqrt function (and make sure we have no imaginary numbers
                tmp=imag(conj(rgr_cp));              % and make sure we have no imaginary numbers
                rgr_cp(tmp~=0)=tmp(tmp~=0);          % fill the now negative numbers into the vector
            else
                rgr_cp = rgr_org;
            end
            
            % Gaussian smoothing            
            kernelLen = length(smooth_k);
            rgr_sm =  conv(rgr_cp, smooth_k, 'same');
            rgr_sm(1:floor(kernelLen./2)) = rgr_cp(1:floor(kernelLen./2));
            rgr_sm(end- floor(kernelLen./2) : end) = rgr_cp( end-floor(kernelLen./2) : end);          
            
            rgrs(:,i) = rgr_sm;
            if ismember(i, indDM) % DM's size regressor: don't smooth
                rgrs(:,i) = rgr_org;
            end
            
            if flagfig
                % Check
                figure(200);
                set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [10 411 1029 306])
                subplot(2,1,1); cla
                plot(rgr_org, 'bo-')
                title([fullRGR4fps(1).features{i}, ': original'])
                
                subplot(2,1,2); cla;
                plot(rgr_cp, 'mo-')
                hold on;
                plot(rgr_sm, 'k-', 'LineWidth', 2)
                title([fullRGR4fps(1).features{i}, ': compressed and smoothed'])
                
                input('')
            end
        end
        
        
        fullRGR4fps(movID).smoRegressors = rgrs;
        
    end
    
end


% function cl = convList(l,k)
% hlen    = floor(length(k)/2);
% llen    = length(l);
% begmean = mean(l(1:hlen))*ones(1,hlen);
% endmean = mean(l(llen-hlen:llen))*ones(1,hlen);
% ltmp    = [begmean l endmean];
% tmp    = conv(ltmp,k);
% start  = round(length(k)-1);
% finish = start+length(l)-1;
% cl     = tmp([start:finish]);





% RAWRGR
% 'TIME'
%     'conspecifics'
%     'heterospecific macaques'
%     'heterospecific mammals'
%     'heterospecific other'
%     'humans'6
%     'Heads'
%     'Faces (full)'8
%     'Faces (side view)'9
%     'Butts'10
%     'Hands'11
%     'Feet'12
%     'Grasslands'
%     'Jungle'
%     'Rocky'
%     'Tundra'
%     'City'
%     'Water'
%     'Manmade structures (non city)'
%     'Food (being eaten)'
%     'Food (being manipulated)'
%     'Food (not manipulated)'
%     'Grooming (active)'
%     'Grooming (sleeping)'
%     'Copulations (displaying)'
%     'Copulations (mounting)'
%     'Play (juvenile juvenile)'
%     'Play (juvenile adult)'
%     'Agression (submissive)'
%     'Agression (threats)'
%     'Agression (chasing)'
%     'Agression (fighting)'
%     'Dyadic interaction'
%     'Sleeping'
%     'Inactive'
%     'Feeding'
%     'Animation'
%     'Sunny'
%     'Raining'
%     'Snowing'
%     'Nighttime '
%     'Fear'
%     'Swimming'
%     'Play (adult adult)'
%     'Yawning'
%     'One Face (full)'
%     'One Face (side view)'






% %smooth the data
%  boxCar = ones(1,50) ./ 50;
%  s_gInfoEvenAll =  conv(gInfoEvenAll, boxCar, 'valid');
%  s_gInfoOddAll = conv(gInfoOddAll, boxCar, 'valid');
%  s_gInfoControl = conv(gInfoControl, boxCar, 'valid');
%  for reg = 1 : 22
%      s_fullRGR4fps(movID).regressors(reg,:) = conv(fullRGR4fps(movID).regressors(reg,:), boxCar, 'valid');
%  end
% clear corrMap;
% for reg = 1 : 22
%    corrMap(reg,1) = corr(s_gInfoEvenAll, s_fullRGR4fps(movID).regressors(reg,:)');
%    corrMap(reg,2) = corr(s_gInfoOddAll, s_fullRGR4fps(movID).regressors(reg,:)');
%    corrMap(reg,3) = corr(s_gInfoControl, s_fullRGR4fps(movID).regressors(reg,:)');
% end
% figure; imagesc(corrMap);
% caxis([-0.4 0.4]);
%
% figure;
% hold on;
% plot(s_gInfoEvenAll);
% plot(s_gInfoOddAll, 'r-');
% plot(s_fullRGR4fps(movID).regressors(9,:), 'k-');
% hold off;



