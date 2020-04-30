
% Read videos (takes long time)
dirMovie = '/procdata/parksh/Stimulus/Movies/Rhesus';
videoObj1 = VideoReader(fullfile(dirMovie, 'Movie1.avi'));
% videoObj2 = VideoReader(fullfile(dirMovie, 'Movie2.avi'));
% videoObj3 = VideoReader(fullfile(dirMovie, 'Movie3.avi'));

iMovie = 1;

nameSubjNeural = 'Tor';
dirDataNeural = sprintf('/procdata/parksh/%s', nameSubjNeural);
load(fullfile(dirDataNeural, sprintf('%s_mov%d_eyeSignal.mat', nameSubjNeural, iMovie)))
%     switch iMovie
%         case 1
obj = videoObj1;
%         case 2
%             obj = videoObj2;
%         case 3
%             obj = videoObj3;
%     end

load(sprintf('/procdata/parksh/MovieRegressors/annotationMovie%d.mat', iMovie))
sceneInfo = cat(2, sta', sto');
nScene = size(sceneInfo, 1);

 for iScene = 1:nScene
        
        setFrames = sceneInfo(iScene,1):sceneInfo(iScene,2);
        
        a(1).cdata = read(obj, epoch(iScene).frame);
        a(1).colormap = [];
        imageFrame = frame2im(a);
        
        figure; % draw the scene
        imagesc(imageFrame)
        hold on
        plot(epoch(iScene).notes.face.x, epoch(iScene).notes.face.y, 'r')
        
        figure;
        [sta(iScene)*30 sto(iScene)*30]
        plot(gazePos.x(1:3, sta(iScene)*30:sto(iScene)*30)', gazePos.y(1:3, sta(iScene)*30:sto(iScene)*30)')
        min(time)
        max(time)
        sto(iS)*30
        
 end
 
 
 
 % eyetraceOverlayEtude.m
%
% Development script for overlaying eye traces on movie frames.
%
% last modified 2013-jul-16
% dbtm
%

movdir = '/einstein0/USRlab/projects/heba/analysis/movies/';
cd(movdir)
pp = setpathsMovies('a');
movImidDir = [ pp.stim 'Movie1imageSequence/'];


trial = 3; % first presentation of movie 1.
mov = 1;
chan = 122;
totMovImid = dir([movImidDir '*.bmp']);


fname = ['tmov' num2str(mov) 'sig' num2str(chan) 'a.mat'];
load([pp.mas fname]);
load(['/archive0/USRlab/projects/heba/analysis/movies/'  'fixEMOffsets/fixed735037_' num2str(trial)]); % ([pp.analysis 'fixEMOffsets/fixed735037_' num2str(trial)])
dat = movieTimeScale(dat,'ms');
unpack;

subdat = selectTrials(dat,TRID,trial); %subdat is different from dat in 
% it only includes one trial, this is for the eyetrace purposes as opposed
% to dat which includes all days all trials. 



% movie specs
FrameDur = 1000/30; % ms per frame.
movieFrameX = [-5.2 +5.2];
movieFrameY = [-3.8 +3.8]; % dva
Xdva = movieFrameX(2)-movieFrameX(1);
Ydva = movieFrameY(2)-movieFrameY(1);


fname = ['tmov' num2str(mov) 'sig' num2str(chan) 'a.mat'];
load([pp.mas fname]);
load([movdir 'fixEMOffsets/fixed735037_' num2str(trial)])
dat = movieTimeScale(dat,'ms');
unpack;

dat = selectTrials(dat,TRID,trial);

% Nframes = length(totMovImid.name);
for f = 1:10
    imgFile = totMovImid(f).name;
    frame = imread([movImidDir imgFile]);
    
    w1 = (f-1)*FrameDur;
    w2 = w1+FrameDur;
    snip = find(dat.t>=w1 & dat.t<=w2);
%     xx = dat.eye(1,snip,1);
%     yy = dat.eye(1,snip,2);
    xx = UNWARP.unwarped_ems(snip,1);
    yy = UNWARP.unwarped_ems(snip,2);
    if f==1
        % compute pixels per degree (ppd)
        [Ypxl Xpxl rgb] = size(frame);
    end
    [cc rr] = eyeDva2Pixels(xx,yy,[Xpxl Ypxl],[Xdva Ydva]);
    figure(f);
    imshow(frame);
    hold on; 
    plot(cc,rr,'m','LineWidth',2);
    %plot(rr,cc,'m','LineWidth',2);
%     fname = ['tmov' num2str(mov) 'UNWARPeye' num2str(f)];
%     saveFig = ['/einstein0/USRlab/projects/heba/results/movies/eyetraceOverlay/' fname];
%     disp(['image' num2str(f)])
%     saveas(gcf, saveFig, 'png');
end



        
        
%         subplot(SP(1))
%         plot(SP(1), setFrames, meanFRCluster_zscore(setFrames+(5*(iMovie-1)*60*30), oldOrderCluster), 'o-')
%         LG= legend('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4', 'Cluster 5', 'Cluster 6', 'Cluster 7',...
%             'Location', [0.0056    0.8433    0.0804    0.1497]);
%         axis tight
%         ylim([-1 10])
%         xlabel('Frame #')
%         ylabel('Normalized response (z)')
%         title(sprintf('Movie %d: Scene %d', iMovie, iScene))
        
%         frame = read(obj,setFrames(1));
%         a(1).cdata = frame;
%         a(1).colormap = [];
%         imageFrame = frame2im(a);
%         subplot(SP(2));
%         imagesc(imageFrame);