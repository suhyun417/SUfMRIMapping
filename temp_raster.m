dirData = '/einstein0/USRlab/data/parksh/';
W = what(dirData);

% indFileNeural = strfind(W.mat, 'smov');

iFileNeural = 10;
fNameNeural = W.mat{iFileNeural+1};
load([dirData, fNameNeural])

nTrial = length(dat.s);


figure(101);
clf

w1 = 0;
n = w1/1000;
prop = 'Color';
val = 'k';

win = [0 300]; % seconds

hold on;
for t=1:length(dat.s)
    set(gca,'YLim',[t-1 t]);
    ts = dat.s{t};
    spikes = ts(find((ts>=win(1)) & ts<=win(2)));
    yline(spikes,prop,val);
end
set(gca,'YLim',[0 length(dat.s)]);
set(gca,'XLim', win);
