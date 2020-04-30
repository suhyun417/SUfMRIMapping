%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% COMPLICATED SOUNDS IN MATLAB %%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% amplitude modulated (AM) sound
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs = 44100;
dur = 1;
t = [0 : 1/Fs : dur];

mrate = 1.5;	% modulation frequency
freq = 500;	    % carrier frequency
mindex = 1;     % modulation index
amp = 2;        % amplitude of carrier

% carrier frequency sound vector
f_c = amp * sin ( 2*pi * freq * t );

% modulation frequency sound vector
f_m = mindex * sin ( 2*pi * mrate * t );

% amplitude modulated sound vector
f_am = [f_c .* f_m];

sound(f_am,Fs)

figure(1)
plot(t,f_c,'b') % plot carrier frequency
hold on
plot(t,f_am, 'g')  % plot amplitude modulated sound
hold on
plot(t,f_m,'r', 'LineWidth', 2)  % plot modulation frequency



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frequency modulated (FM) sound
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs = 44100;
dur = 1;
fc = 400;		% carrier frequency

mrate = 4;		% modulation rate
mindex = 100;	% modulation index (for fm = max_freq_change/modulation rate)
amp = 0.5;

t = [ 0 : 1/Fs : dur ];

f_fm  = amp * sin((2 * pi * fc * t) + (-mindex * sin(2 * mrate * pi * t))) ;

sound(f_fm,Fs)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% random noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs = 44100;
amp = .1;

y = amp * randn(1,Fs);

sound(y,Fs)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% chirp -- sound frequency varies with time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fs = 44100;
t = 0:1/Fs:1;

% linear chirp (f(t) = f0 + kt)
f0 = 200; % freq at time 0
k = 300; % rate of freq change (chirp rate)

chirp = sin(2*pi*(f0 + k/2*t).*t);

sound(chirp,Fs)

subplot(4,1,1);
plot(t,chirp);


% exponential chirp (f(t) = f0*k^t)
f0 = 200;
k = 10;

exp_chirp = sin(2*pi*f0/log(k)*(k.^t-1));

sound(exp_chirp,Fs)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a gating window -- fade in and out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fs = 44100;
t = [0:1/Fs:1-1/Fs];
amp = 1;   % amplitude
wdm = 200;  % gating window in milliseconds

% noise
y = amp * randn(1,Fs);

% sine wave
freq = 440;
y = amp * sin ( 2*pi * freq * t );

% gate sound
y1 = wind(wdm,Fs,y);

sound(y1,Fs)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% play a scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs = 44100;
dur = .2;
t = [0 : 1/Fs : dur-1/Fs];

% general formula for n-split equitempered octave
% f(m) = sin ( 2*pi * freq * 2^(m/n) * t)

f0 = 400;   % fundamental frequency
n_split = 12;   % number of equidistant tones that octave is split into
n = [0:n_split];
for m = 1:length(n)
	f(m,:) = sin ( 2*pi * f0 * 2^(n(m)/n_split) * t );
end

y = [];  % prepare sound vector
for k = 1:13
	y = [y f(k,:)];
    y = wind(Fs,20,y);
end
y = .4 * y;
sound(y,Fs)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% play a tune
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Script to create tunes in the form of pure tones
% It really is only designed for simple tunes spanning less than 1 octave in 4/4 time

clear all;
% make the pitches for the 12 notes from a start pitch using equal temperament tuning
f0 = 440;
n_split = 12;
num_pitches = 25; % number of different pitches (2 octaves in this case)

for p = 1:num_pitches
    pitch(p) = f0 * 2^((p-1)/n_split);
end

% pitch sequence for tune
pitch_seq = [1 3 5 1 1 3 5 1 5 6 8 5 6 8 ];

% note/rest duration sequence
dur_seq = [.5 .5 .5 .5 .5 .5 .5 .5 .5 .5 1 .5 .5 1];

note = 500;
srate = 44100;  % sample rate
wdms = 20;      % gating window in ms

sig1 = [];  % prepare sound vector
for q = 1:length(pitch_seq)     % for number of pitches

    t = [0: 1/srate : dur_seq(q)*note/1000];    % define time vector

    if pitch_seq(q) == 0    % if rest
        wav = zeros(1,length(t));
    else
        wav = sin ( 2*pi * pitch(pitch_seq(q)) * t );   % frequency pitch(pitch_seq(q)) vector
        wav = 0.2 * wav/std(wav);       % fix rms level noise
        wav = wind(srate,wdms,wav);     % window (need wind.m function)
    end
    sig1 = [sig1 wav];  % add
end

sig1 = [sig1', sig1'];  % make stereo signal
wavwrite(sig1, srate, 'sequence');
sound(sig1,srate)
% end play tune
 

Published with MATLAB® 7.1