% eyetraceOverlayEtude
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