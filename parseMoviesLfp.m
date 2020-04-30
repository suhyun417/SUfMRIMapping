function dat = parseMoviesLfp(header,pp)
%
% dat = parseMoviesLfp(header,pp)
%
% Builds data structure for movies project containing LFP and BLP fields.
%
% Processing stream revisions in 2012-dec:
% 1. Spikes extracted from .plx files rather than .nex
% 2. Event codes and timestamps based on merged dgz and tdt files.
% 3. Standardized naming conventions for data files, scripts, etc.
%
% last modified 2013-dec-20
% dbtm

% build header
namestr = header.fname;
%if isempty(header.date)
%   header.date = str2num(namestr(find(ismember(namestr,'0123456789'))));
%end
header.ver = 0.6;
dat.h = header;
dat.h.parseDate = datestr(now);
dat.c = [];
%dat.s = [];
%dat.wf = [];
dat.lfp = [];
dat.blp = [];
dat.eye = [];
dat.t = [];

lfpFile = [pp.rare 'lfp/' dat.h.fname '_' num2str(dat.h.lfpchan) '.lfp'];
fid = fopen(lfpFile);
lowpass = fread(fid,'int16');
fclose(fid);
% convert to microVolts
lowpass = 1000*double(lowpass)./[2^15];
msTime = 1:length(lowpass);


%dgzfile = [pp.rare dat.h.
%dgz = dg_read([pp.rare 'toroid20120619a1.dgz']);
%[tt cc AbsoluteTimes tCurr] = dgz2events(fname);


% % get spikes
% SR = 24414.4; % sample rate in Hz
% chan = dat.h.chan;
% abc = ['a':'z'];
% sig = 0;
% waveform = [];
% wftime = [];
% for c=1:length(chan)
%     plxfile = [pp.rare dat.h.src.tank dat.h.src.plx '_' num2str(chan(c)) '.plx'];
%     try
%         tsc = plx_info(plxfile,1); %fullread
%     catch
%         tsc = [];
%         disp([plxfile ' not read.']);
%     end
%     if ~isempty(tsc)
%         Nspikes = max(find(tsc(:,2)~=0))-1;
%         for s = 1:Nspikes
%             [n npw ts wave] = plx_waves(plxfile,1,s); % nominally ch1
%             if n>0
%                 sig=sig+1;
%                 TS{sig} = ts*1000; % in milliseconds
%               %  trace = nan(1,max(length(wave),size(waveform,2)));
%               %  trace(1:size(wave,2)) = 1000*mean(wave)./[2^15]; % in microvolts
%               %  waveform(sig,:) = trace;
%               %  wftime(sig,:) = 1000*[[1/SR]:[1/SR]:[1/SR]*npw]; % in milliseconds
%                 snames(sig,:) = ['sig' sprintf('%03d',chan(c)) abc(s)];
%             end
%         end
%     end
% end
% dat.h.snames = snames;
% dat.wf.v = waveform;
% dat.wf.t = wftime;
% dat.wf.ts = TS;

% get tdt tank info
load ([pp.rare dat.h.fname '-tdt.mat']);
dat.h.date = tdt.dateS;
vec = datevec(dat.h.date);
yearMonthDay = datenum(vec(1:3));
hourMinuteSecond = datenum([0 0 0 vec(4:6)]);
tt = tdt.times;
ee = tdt.codes;

% create contents of dat.c and dat.s
movieid = [4001:4030];
movieidShift = 4000;
GoodTrials = find(ismember(ee,movieid));
eval(unpack_header(dat));
t=0;
for m = 1:length(GoodTrials)
   % for sig = 1:length(TS)
        t=t+1;
        dat.c(t,TRID) = NaN;
        dat.c(t,SESSION) = dat.h.session;
        dat.c(t,DATE) = yearMonthDay;
        trel = tt(GoodTrials(m));
        dat.c(t,TREL) = trel;
        dat.c(t,TR) = m;
        dat.c(t,MOV) = ee(GoodTrials(m))-movieidShift;
        %         dat.c(t,SIG) = sig;
        %         objetTrouve = find((TS{sig}>=(trel+dat.h.movieWin(1))) & (TS{sig}<=(trel+dat.h.movieWin(2))));
        %         spikes = TS{sig}(objetTrouve);
        %         spikes = spikes-trel;
        %         dat.s{t} = spikes;
         tsnip = find((msTime>=(trel+dat.h.movieWin(1))) ...
                       & (msTime<=(trel+dat.h.movieWin(2))));
            vtrace = lowpass(tsnip);
            lfp(t,:) = vtrace;
        
    %end
end

dat.lfp = lfp;
dat = build_blpDB(dat,dat.h.bands);

% smooth with adjustable Gaussian kernel
for b=1:size(dat.blp,3)
    for t=1:size(dat.blp,1)
        blp = dat.blp(t,:,b);
        band = dat.h.bands(b,:);
        sigma = 2*1000/band(2);
        mirror = fliplr(blp);
        padded = [mirror blp mirror];
        kwidth = -sigma*3:3*sigma;
        kernel = normpdf(kwidth,0,sigma);
        s = conv(padded,kernel);
        s(1:length(mirror)) = [];
        s(length(s)-length(mirror)+1:length(s)) = [];
        s(1:floor(length(kernel)/2)) = [];
        s(length(blp)+1:length(s)) = [];               
        dat.blp(t,:,b) = s;
    end
end

for t=1:size(dat.lfp,1)
    lfpTrace = dat.lfp(t,:);
    
    sigma = 10;
    mirror = fliplr(lfpTrace);
    padded = [mirror lfpTrace mirror];
    kwidth = -sigma*3:3*sigma;
    kernel = normpdf(kwidth,0,sigma);
    s = conv(padded,kernel);
    s(1:length(mirror)) = [];
    s(length(s)-length(mirror)+1:length(s)) = [];
    s(1:floor(length(kernel)/2)) = [];
    s(length(lfpTrace)+1:length(s)) = [];
    dat.lfp(t,:) = s;
end


% dat = align_timebase(dat,T_STIMON);

dat.h.units = 'ms';
dat = parseEyepos(dat,pp);

% downsample dat.lfp and dat.lfp to match dat.eye and dat.t
snipTime = dat.h.movieWin(1):dat.h.movieWin(2);
dsTime = find(ismember(snipTime,dat.t));
dat.lfp = dat.lfp(:,dsTime);
dat.blp = dat.blp(:,dsTime,:);











