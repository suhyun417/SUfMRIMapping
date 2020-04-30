% quick check on McMahon's scenes and frame

% Read videos (takes long time)
dirMovie = '/procdata/parksh/Movies';
videoObj1 = VideoReader(fullfile(dirMovie, 'Movie1.avi'));
videoObj2 = VideoReader(fullfile(dirMovie, 'Movie2.avi'));
videoObj3 = VideoReader(fullfile(dirMovie, 'Movie3.avi'));

iMovie = 1;

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

% Size a figure based on the video's width and height.
hf = figure;
set(hf, 'position', [150 150 obj.Width  obj.Height])
for iScene = 1:nScene
    
    setFrames = sceneInfo(iScene,1):sceneInfo(iScene,2);
    nFrame = length(setFrames);
    clear checkmov
    checkmov(1:nFrame) = struct('cdata',zeros(obj.Height,obj.Width, 3,'uint8'), 'colormap',[]);
    
    
    % Read one frame at a time.
    for k = 1 : nFrame
        indFrame = k; %setFrames(k);
        checkmov(k).cdata = read(obj, indFrame);
    end
    
    
    % Play back the movie once at the video's frame rate.
    figure(hf); clf;
    movie(hf, checkmov, 1, 10); %obj.FrameRate);
    
    input('')
end




% % Collect the scenes
% countFrame = 0;
% for iMovie = 1:3
%     validFrame_range = [(5*(iMovie-1)*60*30)+1, (5*iMovie*60*30)];
%     validLocHigh = [];
%     validLocHigh = locHigh(locHigh>validFrame_range(1) & locHigh<validFrame_range(2));
%     
%     switch iMovie
%         case 1
%             obj = videoObj1;
%         case 2
%             obj = videoObj2;
%         case 3
%             obj = videoObj3;
%     end
%     
%     % Read one frame at a time.
%     for k = 1 : length(validLocHigh)
%         indFrame = validLocHigh(k)-(5*(iMovie-1)*60*30);
%         checkmov(countFrame+k).cdata = read(obj,indFrame);
%     end
%     countFrame = countFrame + length(validLocHigh);
% end
    
    

