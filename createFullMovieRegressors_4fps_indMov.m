function [fullRGR4fps] = createFullMovieRegressors_4fps_indMov(unimov)


% unimov = [1 2 3];


% Get the raw highlevel movie features
highlevpath = '/procdata/russbe/CodingXLS';
lowlevpath = '/procdata/parksh/MovieRegressors'; %/procdata/russbe/LowlevelMovRegressors';
%   eyempath   = '/procdata/russbe/EyeMovements';

nFramePerMovie = 4*300; % 4 frames per second for 5-min (300 sec) movie
fullRGR4fps = struct([]);

% First, get all regressors available
for m=1:length(unimov)
    
    rgrs = [];
    
    highlevfile = fullfile(highlevpath, sprintf('CodingSheetmov%d_inv.xls',unimov(m)));
    lowlevfile = fullfile(lowlevpath, sprintf('Movie%d_10fps_rgr.mat',unimov(m)));
    
    % Low level RGRs    
    lowRGR = load(lowlevfile);  % loads RGR variable
    rgrs = resample(lowRGR.RGR10fps.regressors', 4, 10);  % resample from 10hz to 4hz sampling rate    
    nLowVars = length(lowRGR.RGR10fps.features);  
        
    % High level features in 4fps
    [num,str] = xlsread(highlevfile);                
    nHighVars = length(str(1,:));
    for i=1:nHighVars
%         tmp   = num(1:nFramePerMovie,i);
        % 	  rtmp  = resample(tmp,100*Fs,100*TR);
        % 	  rtmp(max_tr+1:end) = [];
        rgrs(:,nLowVars+i) = num(1:nFramePerMovie,i); %rtmp;
    end
    
%     cregressors = [cregressors rgrs'];
    
    if m==1
        varnames_low = lowRGR.RGR10fps.features;
        varnames_high = str(1,:);
    end
    
    fullRGR4fps(m).movieID = unimov(m);
    fullRGR4fps(m).regressors = rgrs;
    
end

[fullRGR4fps.features] = deal([varnames_low varnames_high]');



% highLevFeatures = {'Faces', 'One face', 'Body parts', 'Conspecifics', 'Humans', 'Any animal', 'Dyadic interaction',...
%         'One Face (full)', 'One Face (side view)', 'Faces (full)', 'Faces (side view)', 'Hands'}';



