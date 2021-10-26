%MovCorrEtudeImageGenerator
%
% Exploratory script for looking at correlations.
%
% last modified 2013-apr-26
% dbtm

detail = [0 300];
tickInterval = 30;

win = [0 300];
kernel = 0.5;
pp = setpathsMovies;
pp.results = '/einstein0/USRlab/projects/heba/results/movies/posterFigures';

if ~exist('sdf')
    day1 = 735037;
    flist = dir([pp.mas 'tmov1*']);
    for f=1:length(flist)
        fname = flist(f).name;
        load([pp.mas fname]);
        dat = movieTimeScale(dat,'sec');
        unpack;

        dat = selectTrials(dat,DATE,day1);
        [resp time] = getSdf(dat,win,kernel);
        resp = resp-min(resp);
        resp = resp./max(resp);
        sdf(f,:) = resp;
        disp(f);
    end
    missing = find(isnan(mean(sdf')));
    sdf(missing,:) = [];
    resam = 1:10:length(time);
    sdf = sdf(:,resam);
    time = time(resam);
    
    
end

% figure(1);clf;
% imagesc(sdf);colormap(hot);
% set(gca,'XLim',win)
% n = win(1):50:win(2);
% set(gca, 'XTickLabel', n/100)
