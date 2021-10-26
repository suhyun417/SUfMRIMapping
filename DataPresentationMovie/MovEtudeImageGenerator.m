% MovEtudeImageGenerator
%
% script that combines eye traces as well as rasters to create movie
% images.
%
% last modified 2013-jul-19
% hde
%

cd /einstein0/USRlab/projects/heba/analysis/movies;
pp = setpathsMovies('a');
movImidDir = [ pp.stim 'Movie1imageSequence/'];

trial = 3; % first presentation of movie 1.
mov = 1;
chan = 122;
totMovImid = dir([movImidDir '*.bmp']);
fullwin = [0 300];
kernel = 0.1;


fname = ['tmov' num2str(mov) 'sig' num2str(chan) 'a.mat'];
load([pp.mas fname]);
load([movdir 'fixEMOffsets/fixed735037_' num2str(trial)])
dat = movieTimeScale(dat,'ms');
unpack;

m = (fullwin(2)/length(totMovImid));

% for imid = 1:length(totMovImid)
for imid = 1
    n = ((imid*m)-m);
    win = [(n-1.5) (n+1.5)];
    figure(imid)
    subplot(3,2,1:2);
    eyetraceoverlay (mov,imid,totMovImid,movImidDir)

    subplot(3,2,3:4);
    hold on;
    part2ult(dat,imid,win);
    yline(n);

    subplot(3,2,5:6);
    hold on;
imagesc(sdf);colormap(hot);
set(gca,'XLim',win)
k = win(1):0.5:win(2);
set(gca, 'XTickLabel', k)
 yline(n);
%     time = (time*1000);
%     plot(time,sdf,'Color','r');
%     ymax = ceil(max(max(sdf)));
%     axis([win 0 ymax])
%     set(gca,'YTick', [0 ymax/4 ymax/2 3*ymax/4 ymax]);
%     yline(n);

%     fname = ['tmov' num2str(mov) 'sig' num2str(chan) '_' num2str(imid)];
%     saveFig = [pp.results 'ultMovRas2/' totMovImid.name{imid}];
    disp(['image' num2str(imid)])
%     saveas(gcf, saveFig, 'png');
end


home = '/einstein0/USRlab/projects/heba/analysis/movies/DataPresentationMovie/';
cd(home)






dat = selectTrials(dat,TRID,trial);

% Nframes = length(totMovImid.name);
for f = 101:600
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
    fname = ['tmov' num2str(mov) 'UNWARPeye' num2str(f)];
    saveFig = ['/einstein0/USRlab/projects/heba/results/movies/eyetraceOverlay/' fname];
    disp(['image' num2str(f)])
    saveas(gcf, saveFig, 'png');
end


