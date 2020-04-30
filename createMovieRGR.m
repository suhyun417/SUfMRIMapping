% function [smRGR] = smoothCompressMovieRegressors(fullRGR4fps, indexRgr).m



% NAMING SCHEME
%   sum     Sum the categories
%   log     Take the log of the values
%   z       Zscore the values
%   sqrt    Take the square root of the values
%   sm#     smooth the Regressor with a # kernel
%   c#      compress the values by # (divide by #)
%   lp#     Low Pass with max being #
%   hp#     High Pass with min being #
%   bp#_#   Band Pass with min_max
%   ave     Average the categories
%   ms#     Split the data at # becomes 0's and 1's
%   msMN    Split the data at the mean of the variable
%   msMD    Split the data at the median of the variable

%  ORDER OF FEATURES
%  1   =  Speed(log_MION)
%  2   =  Mot_Diverg_Rect(log_MION)
%  3   =  Mot_Diverg_STD(log_MION)
%  4   =  Mot_Diverg_Beta(log_MION)
%  5   =  Mot_Diverg_Gamma(MION_sqrt_sm2)
%  6   =  Scene_Cuts(log_MION_sm2)
%  7   =  Faces_Full(sqrt_MION)
%  8   =  Faces_SV(sqrt_MION)
%  9   =  1_Face_Full(sqrt_MION)
%  10  =  1_Face_SV(sqrt_MION)
%  11  =  Heads(sqrt_MION)
%  12  =  Faces_all(sqrt_MION)
%  13  =  1_Face_comb(sqrt_MION)
%  14  =  Luminance(log_MION)
%  15  =  Contrast_STD(log_MION)
%  16  =  Contrast_Beta(log_MION)
%  17  =  Contrast_Gamma(log_MION)
%  18  =  SF_low(log_MION)
%  19  =  SF_high(log_MION)
%  20  =  SF_ratio(log_MION)
%  21  =  Butts(MION_sqrt)
%  22  =  Bodies(sqrt_MION)
%  23  =  Extremities(sqrt_MION)
%  24  =  Hands(sqrt_MION)
%  25  =  Animals(MION_sqrt)
%  26  =  Macaque(MION_sqrt)
%  27  =  Rhesus(MION)
%  28  =  Human(MION)
%  29  =  Mot_Local_STD(MION)
%  30  =  Mot_Bartel_Local(MION)
%  31  =  Mot_Bartel_Global(MION)
%  32  =  Mot_Bartel_Resid(MION)
%  33  =  Mot_Bartel_Total(MION)
%  34  =  Aggression(sqrt_MION)
%  35  =  Affiliative(sqrt_MION)
%  36  =  Food(sqrt_MION)
%  37  =  Aggr_Play(sqrt_MION)
%  38  =  Dyadic(sqrt_MION)
%  39  =  Saccades(sqrt_MION)

%   function RAWRGR = getAllRawRegressors_(unimov,TR,Fs,max_tr)


unimov = [1 2 3];

% Get the raw highlevel movie features
highlevpath = '/procdata/russbe/CodingXLS';
lowlevpath = '/procdata/parksh/MovieRegressors'; %/procdata/russbe/LowlevelMovRegressors';
%   eyempath   = '/procdata/russbe/EyeMovements';

nFramePerMovie = 4*300; % 4 frames per second for 5-min (300 sec) movie
fullRGR4fps = struct([]);

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


rgrlist=[1 29 13 12 23 25 16 17 14 39]; % [1:Global_Speed 2:Local_Motion 3:1_Face_comb 4:Faces_all ....
%  5:Extremities 6:Animals 7:Contrast_Beta 8:Contrast_Gamma 9:Luminance]







%   lowlevfile = sprintf('Movie1_rgr');
%   lowlevfullpath = sprintf('%s/%s.mat',lowlevpath,lowlevfile);
%   load(lowlevfullpath);  % loads RGR variable
%
%   eyefile = sprintf('e66SacRGR');
%   eyefullpath = sprintf('%s/%s.mat',eyempath,eyefile);
%   load(eyefullpath);  % loads e66 RGR variable
%
%   varnames = [varnames RGR.features];
%   varnames = [varnames {'NumSacs'}];

%   if nargin < 2  %% ADDED skipTR variable 10/04/12 BER
%       skipTRs=0; % don't skip any TRs if no skip is sent
%   end

%   cregressors = [];
%   rgrs = [];
%   for m=1:length(unimov)
% 	codingfile   = sprintf('CodingSheetmov%d_inv',unimov(m));
% 	lowlevfile = sprintf('Movie%d_rgr',unimov(m));
%
% 	xlsfullpath = sprintf('%s/%s.xls',codingpath,codingfile);
% 	[num,str] = xlsread(xlsfullpath);
% % 	lowlevfullpath = sprintf('%s/%s.mat',lowlevpath,lowlevfile);
% % 	load(lowlevfullpath);  % loads RGR variable
%
%
% 	nvars = length(str(1,:));
% 	for i=1:nvars
% 	  tmp   = num(1:1200,i);
%
% % 	  rtmp  = resample(tmp,100*Fs,100*TR);
% % 	  rtmp(max_tr+1:end) = [];
% 	  rgrs(:,i) = tmp; %rtmp;
%     end
%
% %     % low levels
% % 	nrgr = size(RGR.regressors,1);
% % 	for i=1:nrgr
% % 	  rgrs(:,nvars+i)=RGR.regressors(i,:);
% %     end
%
%     rgrs(:,nvars+i+1) = SacRGR(unimov(m)).rgr';
% % 	tvals = [(-20*TR):TR:(20*TR)];
% % 	% This is just made-up
% % 	MION_k = gampdf(tvals,TR,2*TR);
% % 	MION_k = MION_k./sum(MION_k);
% % 	crgrs = doConv(rgrs,MION_k);
% % 	cregressors = [cregressors crgrs];
%     cregressors = [cregressors rgrs'];
%   end
%
%   RAWRGR.names = varnames;
%   RAWRGR.rgrs  = cregressors';


%
% Here we replace the initial list of regressors
% with some that we tailor (e.g. combining, compressing, etc)
%
%   ALLRGR = SU_addSpecificRegressors_Motion(RAWRGR,MCD.TR);
%   ALLRGR = SU_addSpecificRegressors_Motion2(RAWRGR,MCD.TR);
%   ALLRGR = SU_addSpecificRegressors_Faces(RAWRGR,MCD.TR);
%   ALLRGR = SU_addSpecificRegressors_LowVision(RAWRGR,MCD.TR);
%   ALLRGR = SU_addSpecificRegressors_Bodies(RAWRGR,MCD.TR);
%   ALLRGR = SU_addSpecificRegressors_Behav(RAWRGR,MCD.TR);
ALLRGR = SU_addSpecificRegressors_Compare(RAWRGR,MCD.TR);




% % setMovID = [1 2 3];
%
% cd /projects/parksh/NeuralBold/analysis
%
% dataDir = '/procdata/parksh/MovieRegressors/';
% load(fullfile(dateDir, 'HighLevelRGR.mat'))
%
% % Create Movie RGRs: rgrs for each movie in separate struct (so that it can
% % be concatenated if you want)
%
% nFramePerMovie = 4*300;
%
% iMov = 1; %:length(setMovID);
% movID = setMovID(iMov);
% % movID=1; %:3; % loop
%
% fullRGR4fps(movID).movieID = movID;
%
% % Low level RGRs
% fName_lowLevRGR = sprintf('Movie%d_10fps_rgr.mat', movID);
% tempStruct = load(fullfile(dataDir, fName_lowLevRGR));
% nLowLevRGR = length(tempStruct.RGR10fps.features);
%
% tempLowRGR = zeros(nFramePerMovie, nLowLevRGR); % frame x rgrs
% tempLowRGR = resample(tempStruct.RGR10fps.regressors', 4, 10);  %resample from 10hz to 4hz sampling rate
%
% fullRGR4fps(movID).regressors4fps(movID).regressors = tempLowRGR;
% fullRGR4fps(movID).regressors4fps(movID).features = tempStruct.RGR10fps.features';
%
% % % smoothing
% % indRGR_log = [1, 2, 3, 6:17, 19];
% % % tempInd = zeros(nFramePerMovie, length(tempStruct.RGR10fps.features));
% % % tempInd(:, indRGR_log) = 1;
% % fullRGR4fps(movID).regressors4fps(iMov).regressors(:,indRGR_log) = real(log(tempLowRGR(:,indRGR_log)));
%
% % High level RGRs
% %20 Faces (sqrt)
% %21 one face (sqrt)
% %22 body parts (sqrt)
% %23 conspecifics (sqrt)
% %24 humans (sqrt)
% %25 any animal (sqrt)
% %26 dyadic interaction
% %27 One Face (full) (sqrt)
% %28 One Face (profile) (sqrt)
% %29 Faces (full) (sqrt)
% %30 Faces (profile/side view) (sqrt)
% %31 Hands(sqrt_MION) 11
% highLevFeatures = {'Faces', 'One face', 'Body parts', 'Conspecifics', 'Humans', 'Any animal', 'Dyadic interaction',...
%     'One Face (full)', 'One Face (side view)', 'Faces (full)', 'Faces (side view)', 'Hands'}';
% nHighLevRGR = length(highLevFeatures);
% % highLevRGR
%
% for iRgr = nLowLevRGR+1:nLowLevRGR+nHighLevRGR
%
%
%     fullRGR4fps(movID).regressors
%
%     if iRgr == 20 %face
%         nowdata = RAWRGR.rgrs(:,8) + RAWRGR.rgrs(:,9); %1x 3600entries
%         fullRGR4fps(movID).regressors(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
%         fullRGR4fps(movID).features{iRgr} = 'Faces';
%     elseif iRgr == 21 %one-face
%         nowdata = RAWRGR.rgrs(:,46) + RAWRGR.rgrs(:,47); %1x3600 entries
%         fullRGR4fps(movID).regressors(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
%         fullRGR4fps(movID).features{iRgr} = 'One face';
%     elseif iRgr == 22 %body parts
%         nowdata = RAWRGR.rgrs(:,10) + RAWRGR.rgrs(:,11) + RAWRGR.rgrs(:,12);
%         fullRGR4fps(movID).regressors(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
%         fullRGR4fps(movID).features{iRgr} = 'Body parts';
%     elseif iRgr == 23 %conspecifics
%         nowdata = RAWRGR.rgrs(:,2);
%         fullRGR4fps(movID).regressors(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
%         fullRGR4fps(movID).features{iRgr} = 'Conspecifics';
%     elseif iRgr == 24 %humans
%         nowdata = RAWRGR.rgrs(:,6);
%         fullRGR4fps(movID).regressors(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
%         fullRGR4fps(movID).features{iRgr} = 'Humans';
%     elseif iRgr == 25 %any animal
%         nowdata = RAWRGR.rgrs(:,2) + RAWRGR.rgrs(:,3) + RAWRGR.rgrs(:,4) + RAWRGR.rgrs(:,5) + RAWRGR.rgrs(:,6);
%         fullRGR4fps(movID).regressors(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
%         fullRGR4fps(movID).features{iRgr} = 'Any animal';
%     elseif iRgr == 26 %dyadic interaction
%         nowdata = RAWRGR.rgrs(:,33);
%         fullRGR4fps(movID).regressors(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
%         fullRGR4fps(movID).features{iRgr} = 'Dyadic interaction';
%     elseif iRgr == 27 %One Face (full)
%         nowdata = RAWRGR.rgrs(:,46);
%         fullRGR4fps(movID).regressors(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
%         fullRGR4fps(movID).features{iRgr} = 'One face (full)';
%     elseif iRgr == 28 %One Face (profile)
%         nowdata = RAWRGR.rgrs(:,47);
%         fullRGR4fps(movID).regressors(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
%         fullRGR4fps(movID).features{iRgr} = 'One face (side view)';
%     elseif iRgr == 29 %Faces (full)
%         nowdata = RAWRGR.rgrs(:,8);
%         fullRGR4fps(movID).regressors(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
%         fullRGR4fps(movID).features{iRgr} = 'Faces (full)';
%     elseif iRgr == 30 %Faces (profile)
%         nowdata = RAWRGR.rgrs(:,9);
%         fullRGR4fps(movID).regressors(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
%         fullRGR4fps(movID).features{iRgr} = 'Faces (side view)';
%     elseif iRgr == 31 %Hands
%         nowdata = RAWRGR.rgrs(:,11);
%         fullRGR4fps(movID).regressors(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
%         fullRGR4fps(movID).features{iRgr} = 'Hands';
%     end
%
% end
%
%
%
%
%
%
%
%
%
% %construct the 56 sec x 15 clips long regressor
% fullRGR4fps(movID).regressors = zeros(30, 224*15);
% for iRgr = 1 : 30
%     counter = 1;
%     for clip = 1 : 15
%
%             nowNdx = 9+(clip-1)* 240;
%            % nowNdx
%             if iRgr <= 19
%                 fullRGR4fps(movID).regressors(iRgr,counter: counter+224-1) =  lowLevelRGR_4Hz(iRgr, nowNdx: nowNdx+224-1);
%
%             elseif iRgr == 20 %face
%                 nowdata = RAWRGR.rgrs(:,8) + RAWRGR.rgrs(:,9); %1x 3600entries
%                 fullRGR4fps(movID).regressors(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
%             elseif iRgr == 21 %one-face
%                 nowdata = RAWRGR.rgrs(:,46) + RAWRGR.rgrs(:,47); %1x3600 entries
%                 fullRGR4fps(movID).regressors(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
%
%             elseif iRgr == 22 %body parts
%                 nowdata = RAWRGR.rgrs(:,10) + RAWRGR.rgrs(:,11) + RAWRGR.rgrs(:,12);
%                 fullRGR4fps(movID).regressors(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
%             elseif iRgr == 23 %conspecifics
%                 nowdata = RAWRGR.rgrs(:,2);
%                 fullRGR4fps(movID).regressors(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
%             elseif iRgr == 24 %humans
%                 nowdata = RAWRGR.rgrs(:,6);
%                 fullRGR4fps(movID).regressors(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
%             elseif iRgr == 25 %any animal
%                 nowdata = RAWRGR.rgrs(:,2) + RAWRGR.rgrs(:,3) + RAWRGR.rgrs(:,4) + RAWRGR.rgrs(:,5) + RAWRGR.rgrs(:,6);
%                 fullRGR4fps(movID).regressors(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
%             elseif iRgr == 26 %dyadic interaction
%                 nowdata = RAWRGR.rgrs(:,33);
%                 fullRGR4fps(movID).regressors(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
%             elseif iRgr == 27 %One Face (full)
%                 nowdata = RAWRGR.rgrs(:,46);
%                 fullRGR4fps(movID).regressors(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
%             elseif iRgr == 28 %One Face (profile)
%                 nowdata = RAWRGR.rgrs(:,47);
%                 fullRGR4fps(movID).regressors(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
%             elseif iRgr == 29 %Faces (full)
%                 nowdata = RAWRGR.rgrs(:,8);
%                 fullRGR4fps(movID).regressors(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
%             elseif iRgr == 30 %Faces (profile)
%                 nowdata = RAWRGR.rgrs(:,9);
%                 fullRGR4fps(movID).regressors(:, iRgr) = nowdata((movID-1)*nFramePerMovie+1:movID*nFramePerMovie,1);
%             end
%             counter = counter+224;
%
%
%     end % for clip
% end %for iRgr
%
%
% % for mov = 1  : 3
% %     fName = ['Movie' num2str(mov) '_10fps_rgr.mat'];
% %
% %     rgrInfo(mov)= load(fullfile(dataDir, fName));
% % end
% % clear mov
% %
% % lowLevelRGR_4Hz = zeros(19, 3600); % 4 frame / sec * 300 sec * 3 movies
% % for mov = 1 : 3
% %     for iRgr = 1 : 19
% %         nowRAW = rgrInfo(mov).RGR10fps.regressors(iRgr,:); %1x3000
% %         nowInterp = resample( nowRAW, 2, 5); %resmaple from 10hz to 4hz sampling rate
% %
% %         lowLevelRGR_4Hz( iRgrfrahco, (1 + (mov-1).*1200) : ( mov .*1200) ) = nowInterp;
% %
% %
% %     end
% % end
%
% %% Regressor storage in fullRGR4fps(movID).regressors variable
%
% % [4, 5, 19]: binary
% % [1, 2, 3, 6:17, 19]: take log
% % [20:30]: sqrt
% % Gaussian smoothing
%
% %1 'Luminance' (log)
% %2 'Contrast' (log)
% %3 'Speed' (log)
% %4 'Rightward Motion' (binary)
% %5 'Leftward Motion' (binary)
% %6 'SF (0.0 to 0.2)' (log)
% %7 'SF (3.0 to 100.0)' (log)
% %8 'SF Ratio' (log)
% %9 'Beta Contrast' (log)
% %10 'Gamma Contrast' (log)
% %11 'Motion Contrast'
% %12 'Motion Dir Beta'
% %13 'Motion Dir Gamma'
% %14 'Motion Div Mean'
% %15 'Motion Div Rectify' (log)
% %16 'Motion Div STD' (log)
% %17 'Motion Div Beta' (log)
% %18 'Motion Div Gamma' (sqrt_sm2)
% %19 'Scene Cuts' (binary) (log_sm2)
% %20 Faces (sqrt)
% %21 one face (sqrt)
% %22 body parts (sqrt)
% %23 conspecifics (sqrt)
% %24 humans (sqrt)
% %25 any animal (sqrt)
% %26 dyadic interaction
% %27 One Face (full) (sqrt)
% %28 One Face (profile) (sqrt)
% %29 Faces (full) (sqrt)
% %30 Faces (profile/side view) (sqrt)
% %31 Hands(sqrt_MION) 11
%
%
%
% %31 %  21  =  Butts(MION_sqrt)
%
% %33 %  34  =  Aggression(sqrt_MION)
% %34 %  35  =  Affiliative(sqrt_MION)
% %35 %  36  =  Food(sqrt_MION)
%
% % 'TIME'
% %     'conspecifics'
% %     'heterospecific macaques'
% %     'heterospecific mammals'
% %     'heterospecific other'
% %     'humans'6
% %     'Heads'
% %     'Faces (full)'8
% %     'Faces (side view)'9
% %     'Butts'10
% %     'Hands'11
% %     'Feet'12
% %     'Grasslands'
% %     'Jungle'
% %     'Rocky'
% %     'Tundra'
% %     'City'
% %     'Water'
% %     'Manmade structures (non city)'
% %     'Food (being eaten)'
% %     'Food (being manipulated)'
% %     'Food (not manipulated)'
% %     'Grooming (active)'
% %     'Grooming (sleeping)'
% %     'Copulations (displaying)'
% %     'Copulations (mounting)'
% %     'Play (juvenile juvenile)'
% %     'Play (juvenile adult)'
% %     'Agression (submissive)'
% %     'Agression (threats)'
% %     'Agression (chasing)'
% %     'Agression (fighting)'
% %     'Dyadic interaction'
% %     'Sleeping'
% %     'Inactive'
% %     'Feeding'
% %     'Animation'
% %     'Sunny'
% %     'Raining'
% %     'Snowing'
% %     'Nighttime '
% %     'Fear'
% %     'Swimming'
% %     'Play (adult adult)'
% %     'Yawning'
% %     'One Face (full)'
% %     'One Face (side view)'
%
%
%
%
%
%
% % %smooth the data
% %  boxCar = ones(1,50) ./ 50;
% %  s_gInfoEvenAll =  conv(gInfoEvenAll, boxCar, 'valid');
% %  s_gInfoOddAll = conv(gInfoOddAll, boxCar, 'valid');
% %  s_gInfoControl = conv(gInfoControl, boxCar, 'valid');
% %  for reg = 1 : 22
% %      s_fullRGR4fps(movID).regressors(reg,:) = conv(fullRGR4fps(movID).regressors(reg,:), boxCar, 'valid');
% %  end
% % clear corrMap;
% % for reg = 1 : 22
% %    corrMap(reg,1) = corr(s_gInfoEvenAll, s_fullRGR4fps(movID).regressors(reg,:)');
% %    corrMap(reg,2) = corr(s_gInfoOddAll, s_fullRGR4fps(movID).regressors(reg,:)');
% %    corrMap(reg,3) = corr(s_gInfoControl, s_fullRGR4fps(movID).regressors(reg,:)');
% % end
% % figure; imagesc(corrMap);
% % caxis([-0.4 0.4]);
% %
% % figure;
% % hold on;
% % plot(s_gInfoEvenAll);
% % plot(s_gInfoOddAll, 'r-');
% % plot(s_fullRGR4fps(movID).regressors(9,:), 'k-');
% % hold off;
%
%
%
