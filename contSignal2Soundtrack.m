function [] = contSignal2Soundtrack(inputSig, Fs_org, flagAM, saveFileName)
%
% Change a continuous signal to a soundtrack
% -- input: 1. inputSig: continuous vector in certain length
%           2. Fs_org: Sampling rate of original vector (in samples per sec)
%           3. flagAM: 1 for amplitude modulation, 0 for frequency
%           modulation
%           4. saveFileName: file name of the soundtrack (it can include
%           directories: if the directory doesn't exist, it will make the
%           directory)
%
% -- output: soundtrack made out of the vector
%            Amplitude changes in the input vector can be reflected as
%           either frequeny changes or intensity changes of the soundtrack
%           (should be chosen by user)
%


Fs_out = 44100; % 44kHz
dur_sec = length(inputSig)/Fs_org; % duration of input
t = 1/Fs_out:1/Fs_out:dur_sec; % sampling time

% Resampling with padding in advance, to prevent an arbitrary up and drop
% at the beginning and the end of the signal
numPadSec = 10;
inputSig_pad = cat(1, repmat(inputSig(1), numPadSec*Fs_org, 1), inputSig, repmat(inputSig(end), numPadSec*Fs_org, 1));
resmpSig_pad = resample(inputSig_pad, Fs_out, Fs_org);
resmpSig = resmpSig_pad(numPadSec*Fs_out+1:end-(numPadSec*Fs_out));



if flagAM
    % Case 1: AM        
    
%     % Make a tail for saveFileName to indicate it's AM
%     saveFileName_tail = '_am';
    
    % Carrier signal: use random noise in case of AM
    amp = 0.01;
    Sc = amp*randn(1, length(t)); % random white noise
    
%     Fc=100; % carrier frequency
%     amp = 0.5;
% %     Sc = amp*randn(1, length(t)); % random white noise
%     Sc = amp*sin(2*pi*Fc*t); % carrier signal
    
    
    % Generate pink (1/f) noise
%     % number of samples
%     N = fs*t;
%     
%     if rem(N,2)
%         M = N+1;
%     else
%         M = N;
%     end
%     
%     % generation of the signal
%     % generate white noise
%     x = randn(M, 1);
%     
%     % FFT
%     X = fft(x);
%     
%     % prepare a vector for 1/f multiplication
%     NumUniquePts = M/2 + 1;
%     n = 1:NumUniquePts;
%     n = sqrt(n');
%     
%     % multiplicate the left half of the spectrum so the power spectral density
%     % is proportional to the frequency by factor 1/f, i.e. the
%     % amplitudes are proportional to 1/sqrt(f)
%     X(1:NumUniquePts) = X(1:NumUniquePts)./n;
%     
%     % prepare a right half of the spectrum - a copy of the left one,
%     % except the DC component and Nyquist frequency - they are unique
%     X(NumUniquePts+1:M) = real(X(M/2:-1:2)) -1i*imag(X(M/2:-1:2));
%     
%     % IFFT
%     x = ifft(X);
%     
%     % prepare output vector y
%     x = real(x(1:N, 1));
%     x = x./max(abs(x));
    % Sc = amp*sin(2*pi*Fc*t); % carrier signal
    
    
    % Modulator (message) signal
    % smarter/more ideal way would be...
    %       1) if the signal includes both positive and negative values,
    %           handle each sign separately (either by using different types of carrier noise or 
    %           by making it stereo)
    % after make the modulatory signal, normalize it as maximum < 1
    
    if min(resmpSig)<0 
        % Should be positive numbers
        % : try to take exponential -> it exaggerates differences of
        % originally positive part of the signal
        % : try to simply shift the distribution to positivie        
        convertedSig = abs(resmpSig - min(resmpSig)); %exp(resmpSig);
    elseif max(resmpSig) > 1
        convertedSig = resmpSig;
    end
    
    indModulation = 20/max(convertedSig); %3; % modulation index
    Sm = indModulation*convertedSig'; % modulation (message) signal
    
    % Amplitude modulation
    S_out = Sm.*Sc;%AM Signal, Amplitude of Carrier changes to (A+Message) %ampSignal.*Sc;
    
else
    % Case 2: FM
    
%     % Make a tail for saveFileName to indicate it's FM
%     saveFileName_tail = '_fm';
    
%     % Carrier signal
    Fc=100; % carrier frequency
    amp = 0.5;
%     Sc = amp*randn(1, length(t)); % random white noise
    Sc = amp*sin(2*pi*Fc*t); % carrier signal
%     
%     % Modulator (message) signal
%     % since it's z-scored value, just multiply 10
%     % maybe in the future I can 1) first conver inputSig to z-score 2) then
%     % multiply 10 to maintain the difference
%     convertedSig = resmpSig*10; %abs(resmpSig - min(resmpSig)); %exp(resmpSig);
% %     f0 = 100;
% %     Sm = sin(2*pi*(f0 + convertedSig').*t);
%     indModulation = 4; % modulation index
%     Sm = indModulation*cos(2*pi*(Fc+convertedSig').*t); % modulation (message) signal
%     
%     % Frequency modulation
%     S_out = 0.5*sin(Sc + Sm); %FM Signal
    
    convertedSig = resmpSig*10; %abs(resmpSig - min(resmpSig)); %exp(resmpSig);
    Fc = 100; % 
    S_out = sin(2*pi*(Fc + convertedSig').*t);
            
end


% final save
tloc = strfind(saveFileName, '/');
saveDir =saveFileName(1:tloc(end)); 
if ~exist(saveDir,'dir') % if the desired directory doesn't exist, make it
    mkdir(saveDir)
    audiowrite(saveFileName, S_out, Fs_out)
else
    audiowrite(saveFileName, S_out, Fs_out)
end



