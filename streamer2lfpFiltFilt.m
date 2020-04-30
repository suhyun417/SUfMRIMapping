function lfp = streamer2lfpFiltFilt(raw);
%
% lfp = streamer2lfpFiltFilt(raw);
%
% Input argument:
% RAW is a vector of floating point values in Volts (i.e. the raw data from streamer).
% 
% Output argument:
% LFP is a low-pass (<300 Hz) filtered vector of int16 voltage values 
% downsampled to 1000 Hz.
%
% Example of usage:
%    sev = readStreamerFile(sevFile);
%    lfp = streamer2lfpFiltFilt(sev);
%    writebinFile(lfp,fname);
%
% To convert lfp to microVolts:
% uV = 1000*double(bin)./[2^15];
%
% last modified 2013-apr-06
% dbtm

waveSampF = 24414;
% hp = 300 ;  % Hz
% lp = 5000 ; % Hz
% WnHP = hp / ((waveSampF)/2) ; % cuttoff frequency, between 0.0 < Wn < 1.0, with 1.0 = 1/2 sample rate
% WnLP = lp / ((waveSampF)/2) ; % cuttoff frequency, between 0.0 < Wn < 1.0, with 1.0 = 1/2 sample rate
% N = 2 ; % Nth order filter
% [bHP, aHP] = butter(N, WnHP,'high');
% [bLP, aLP] = butter(N, WnLP,'low');
% output = filtfilt(bHP,aHP,raw);    % match TDT 1st biquad
% output = filtfilt(bLP,aLP,output);     % match TDT 3rd biquad
% 
% output = int16(output*1000*2^15);


% filter lfp band
lp = 300 ; % Hz
%WnHP = hp / ((waveSampF)/2) ; % cuttoff frequency, between 0.0 < Wn < 1.0, with 1.0 = 1/2 sample rate
WnLP = lp / ((waveSampF)/2) ; % cuttoff frequency, between 0.0 < Wn < 1.0, with 1.0 = 1/2 sample rate
N = 2; % Nth order filter
[bLP, aLP] = butter(N, WnLP,'low');
lfp = filtfilt(bLP,aLP,raw);
lfp = int16(lfp*1000*2^15);

% downsample lfp to 1000 Hz
index = 1:[waveSampF/1000]:length(lfp);
lfp = lfp(round(index));
%mstime = [0:1:length(lfp)];

