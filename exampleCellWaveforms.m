% exampleCellWaveforms
%
% Spike stability figure for movie paper
%
% last modified 2014-may-08
% dbtm

%sigs 6a 12a 13a

% pp = setpathsMovies('a');
pp.rare = '/archive0/USRlab/data/mcmahond/moviesRhombus/';
flist = dir([pp.rare '*-sorted.plx']);

wftime = [1/25:1/25:[1/25]*32];
 
chan = [6 12 13];
%chan = [9 12 13];

unit = 1;
bins = [0:101];
for d=1:length(flist)
%for d=1
    plxfile = [pp.rare flist(d).name];
    for c=1:length(chan)
        %for c=1
        [n, npw, ts, wave] = plx_waves_v(plxfile, chan(c), unit);
        spk{c}(d,:) = 1000*mean(wave);
        ts = ts*1000;
        deltats = ts(2:end)-ts(1:end-1);
        isihist = hist(deltats,bins);
        isihist(end) = [];
        isihist = isihist./max(isihist);
        isi{c}(d,:) = isihist;
        fname{n} = flist(d).name;
    end
end
% 
% cc = 'kgrb';
% left = [1 3 5];
% right = [2 4 6];
% figure(5);clf;
% for c=1:length(chan)
%     
%     for d=1:size(spk{c},1)
%         subplot(4,2,left(c));
%         hold on;
%         offset = 5;
%         plot(wftime,spk{c}(d,:)-offset,'Color',cc(d));
%         set(gca,'XLim',[wftime(1) wftime(end)]);
%         set(gca,'YLim',[-80 50]);
%         axis off
%         subplot(4,2,right(c));
%         hold on;
%         offset = 5;
%         plot(isi{c}(d,:)-offset,'Color',cc(d));
%         plot(0,isi{c}(d,1)-offset,'.','Color',cc(d));
%         set(gca,'XLim',[-20 100]);
%         axis off        
%     end
%     title(chan(c));
% end

%cc = 'kgrb';
left = [1 3 5];
right = [2 4 6];
figure(5);clf;
for c=1:length(chan)
    
    for d=1:size(spk{c},1)
        subplot(4,2,left(c));
        hold on;
        offset = 7*(d-1);
        %plot(wftime,spk{c}(d,:)-offset,'Color',cc(d));
        plot(wftime,spk{c}(d,:)-offset,'Color','b');
        set(gca,'XLim',[wftime(1) wftime(end)]);
        set(gca,'YLim',[-80 50]);
        axis off
        subplot(4,2,right(c));
        hold on;
        offset = 0.1*(d-1);
% $$$         plot(isi{c}(d,:)-offset,'Color',cc(d));
% $$$         plot(0,isi{c}(d,1)-offset,'.','Color',cc(d));
        plot(isi{c}(d,:)-offset,'Color','b');
        plot(0,isi{c}(d,1)-offset,'.','Color','b');

        set(gca,'XLim',[-20 100]);
        axis off        
    end
    %    title(chan(c));
end


%scalebars
subplot(4,2,7);
xx = [0.5 1 1];
yy = [-50 -50 0];      
line(xx,yy,'Color','k');
set(gca,'XLim',[wftime(1) wftime(end)]);
set(gca,'YLim',[-80 50]);
axis off        

subplot(4,2,8);
line([0 100],[0.5 0.5],'Color','k');
set(gca,'YLim',[0 100]);
set(gcf,'Color','w');   
axis off

%scalebars
subplot(4,2,7);
xx = [0.5 1 1];
yy = [-50 -50 0];      
line(xx,yy,'Color','k');
set(gca,'XLim',[wftime(1) wftime(end)]);
set(gca,'YLim',[-80 50]);
axis off        

subplot(4,2,8);
line([0 100],[0.5 0.5],'Color','k');
set(gca,'YLim',[0 100]);
set(gcf,'Color','w');   
axis off

saveName = ['exampleCellsWaveforms.eps']
saveas(gcf,saveName,'epsc');
