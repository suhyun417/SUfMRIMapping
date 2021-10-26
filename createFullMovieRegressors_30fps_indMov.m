function [fullRGR30fps] = createFullMovieRegressors_30fps_indMov(unimov)
%
% 2016/05/19 Add DM's scene based size regressors in 30fps

% unimov = [1 2 3];


% Get the raw highlevel movie features
highlevpath = '/procdata/russbe/CodingXLS';
lowlevpath = '/procdata/parksh/MovieRegressors'; %/procdata/russbe/LowlevelMovRegressors';
%   eyempath   = '/procdata/russbe/EyeMovements';

nFramePerMovie = 30*300; % 4 frames per second for 5-min (300 sec) movie
fullRGR30fps = struct([]);

% First, get all regressors available
for m=1:length(unimov)
    
    rgrs = [];
    
    highlevfile = fullfile(highlevpath, sprintf('CodingSheetmov%d_inv.xls',unimov(m)));
    lowlevfile = fullfile(lowlevpath, sprintf('Movie%d_10fps_rgr.mat',unimov(m)));
    
    % Low level features    
    lowRGR = load(lowlevfile);  % loads RGR variable
    rgrs = resample(lowRGR.RGR10fps.regressors', 30, 10);  % resample from 10hz to 30hz sampling rate    
    nLowVars = length(lowRGR.RGR10fps.features);  
        
    % High level features
    [num,str] = xlsread(highlevfile);                
    nHighVars = length(str(1,:));
    for i=1:nHighVars
        tmp   = num(:,i);
        rtmp  = resample(tmp, 30, 4); % resample(tmp,100*Fs,100*TR);
        rtmp(nFramePerMovie+1:end) = [];
        rgrs(:,nLowVars+i) = rtmp;
    end
    
    % DM's scene-based features (face, torso, arms, legs)
    load(sprintf('/procdata/parksh/MovieRegressors/annotationMovie%d.mat', unimov(m)))
    sceneInfo = cat(2, sta', sto');
    nScene = length(epoch);
    validFrame_range = [1 nFramePerMovie]; 
    
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
    
    rgrs = cat(2, rgrs, matSizeRGR);
    
%     cregressors = [cregressors rgrs'];
    
    if m==1
        varnames_low = lowRGR.RGR10fps.features;
        varnames_high = str(1,:);
        varnames_DM = {'Face size', 'Torso size', 'Arm size', 'Leg size', 'View angle'};
    end
    
    fullRGR30fps(m).movieID = unimov(m);
    fullRGR30fps(m).regressors = rgrs;
    
end

[fullRGR30fps.features] = deal([varnames_low varnames_high varnames_DM]');



% highLevFeatures = {'Faces', 'One face', 'Body parts', 'Conspecifics', 'Humans', 'Any animal', 'Dyadic interaction',...
%         'One Face (full)', 'One Face (side view)', 'Faces (full)', 'Faces (side view)', 'Hands'}';



