function showRasters(filename)

global dat

% dirData = '/procdata/parksh/Tor'; %/rmov/';
addpath('/library/matlab_utils/');
load(filename) %load([dirData, filename])

nTrial = length(dat.s);



prop = 'Color';
val = 'k';

win = [0 300]; % seconds

spikes = {};
for t=1:length( dat.s)
    ts = dat.s{t};
    spikes{t} = ts(find((ts>=win(1)) & ts<=win(2)));
end

figure(101);
clf
subplot(2,1,1);
hold on;

for i=1:length(spikes)
    set(gca,'YLim',[i-1 i]);
    yline(spikes{i},prop,val);
end

set(gca,'YLim',[0 length(dat.s)]);
set(gca,'XLim', win);
title(filename)

subplot(2,1,2);

% sdf sample period
sdf_samp = 100;

sdf = [];
for i=1:length(spikes)
  sdf(:,i) = getSDF(spikes{i},win);
end

plot(sdf)
hold on;
plot(mean(sdf,2),'k','LineWidth',4)

function sdf_out = getSDF(ts,win)

% note: window is in seconds
  upsamp = 1000;
  timeline = zeros(upsamp*win(2)-upsamp*win(1),1);
  its = (round(ts*upsamp));
  timeline(its) = 1;
  
  sd = 2000;
  k  = normpdf(-3*sd:+3*sd,0,sd);
  sum(k)
  k = k./sum(k);
  
  % convolve, scale to spikes/sec, and decimate
  %
  sdf_out = doConv(timeline,k).*upsamp;
