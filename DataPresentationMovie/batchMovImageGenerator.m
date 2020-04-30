% batchMovImageGenerator
% copyt this from 

mov = 1;
chan = 122;
movImidDir = '/archive0/USRlab/projects/heba/stimuli/movies/Movie1imageSequence/'; %'/einstein0/USRlab/projects/heba/stimuli/Movie1imageSequence/';
totMovImid = dir([movImidDir '*.bmp']);
fullwin = [0 300];
kernel = 0.1;
pp = setpathsMovies('a');
pp.mas = '/procdata/parksh/Tor/';

fname = ['tmov' num2str(mov) 'sig' num2str(chan) 'a.mat'];
load([pp.mas fname]);
unpack;
% [sdf time] = getSdf(dat, fullwin, kernel);
m = (fullwin(2)/length(totMovImid));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MovCorrEtudeImageGenerator;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for imid = 1:length(totMovImid)
for imid = 1
    n = ((imid*m)-m);
    win = [(n-1.5) (n+1.5)];
    figure(imid)
    subplot(3,2,1:2);
    plotMovieFrame(mov,imid,totMovImid,movImidDir)
    
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