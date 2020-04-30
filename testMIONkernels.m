% testMIONkernels.m
%
% short script for testing different kernels of MION which is convolved to
% mean SDF to create neural regressor

%% Set parameters
% movie
setMovie = [1 2 3];
% fmri time series
fmritc = cat(4, dataBOLD.catmvoltc{setMovie});
[nx, ny, nz, nt] = size(fmritc); 
nVox = nx*ny*nz;
reshapedfmritc = reshape(fmritc, nVox, nt)';

% test voxel
a = 36; b = 29; c = 16;
ind = sub2ind([40 64 32], a, b, c);
% test cell and spike density function
iC = 1; % 065a
mnsdf = cat(1,S(iC,setMovie).mnsdf); % gaussian kernel convolved

% non-smoothed firing rate
FR = createCellRegressor_indMov_discreteTime(STDPATH.dataNeural, {'065a'}, [1 2 3], dataBOLD.TR);
mnfr = cat(1,FR(iC,setMovie).mnFR);

% some parameters
fMRI_TR_sec = dataBOLD.TR; %2.4;
timeResNeural_sec = 0.001; % 1000 Hz

% Figure for check
figMIONtest = figure;
set(gcf, 'Color','w', 'PaperPositionMode', 'auto')


%% 1. original MION regressor
% MION function from AFNI
lengthMovie_sec = 300; %length(catmnsdf)/1000/9; %3;
t_org = 0:fMRI_TR_sec:lengthMovie_sec;
k_org = 16.4486 * ( -0.184/ 1.5 * exp(-t_org/ 1.5)...
+0.330/ 4.5 * exp(-t_org/ 4.5)...
+0.670/13.5 * exp(-t_org/13.5) );
rgr_org = doConv(mnsdf, k_org); % no difference from original mnsdf

% let's change the time scale
% -- no difference with t=0:1:50
% -- with t=0:1:50*1000, looks weird: everything's actually advanced, not
% delayed
t1 = 0:1:50; %0:fMRI_TR_sec:lengthMovie_sec;
k1 = 16.4486 * ( -0.184/ 1.5 * exp(-t1/ 1.5)...
+0.330/ 4.5 * exp(-t1/ 4.5)...
+0.670/13.5 * exp(-t1/13.5) );


% let's use non-smoothed firing rate (discrete time bin) in fMRI time unit
% -- not right with doConv(mnfr, k_org);
% -- not right with doConv(mnfr, k1); % when t1=0:1:50;
%       again, it was advanced, not delayed


% let's use centered version of AFNI MION kernel
%  -- should be matched to the timescale (even if I make the time axis of 
%       this kernel longer, the shape of kernel kept same) 
t2 = 0:1:5000;
k2 = 16.4486 * ( -0.184/ 1.5 * exp(-t2/ 1.5)...
+0.330/ 4.5 * exp(-t2/ 4.5)...
+0.670/13.5 * exp(-t2/13.5) );
k2_zeropadded = cat(2, zeros(1, length(k2)-1), k2);
rgr_zeropadded = doConv(mnsdf, k2_zeropadded);

t3_tr = 0:2.4:50;
k3_tr = 16.4486 * ( -0.184/ 1.5 * exp(-t3_tr/ 1.5)...
+0.330/ 4.5 * exp(-t3_tr/ 4.5)...
+0.670/13.5 * exp(-t3_tr/13.5) );
MION_k_afni_tr = cat(2, zeros(1, length(k3_tr)-1), k3_tr);

rgr_afni_tr = doConv(mnfr, MION_k_afni_tr);


% use what BR used
% set up MION function for convolution
tvals = [(-20*dataBOLD.TR):dataBOLD.TR:(20*dataBOLD.TR)];
% This is just made-up
MION_k = gampdf(tvals,dataBOLD.TR,2*dataBOLD.TR);
MION_k = MION_k./sum(MION_k);
% -- nothing happend with convolution to mnsdf
% -- mnfr: looks good but too much smoothed

k = gampdf([-40:2.4:40],4,2);
rgr_gamma = doConv(mnfr, k);




%%% Summary
%%% 1. Use a kernel which x-axis is centered to zero (because the
%%% conv function applies kernel centered to the time series
%%% 2. Make the structure of function be matched to the current time scale
%%% : for example, what BR used is in TR unit, so it doesn't do much
%%% smoothing on spike density function data








%%
% figure;
% plot(t, k1, 'o-')

% if a boxcar function is already convolved: convolve only mion function
% BUT it should be in second unit



% 'MION(d)'     
% = 1 parameter block stimulus of duration 'd',
% intended to model the response of MION.
% The zero-duration impulse response 'MION(0)' is
% h(t) = 16.4486 * ( -0.184/ 1.5 * exp(-t/ 1.5)
% +0.330/ 4.5 * exp(-t/ 4.5)
% +0.670/13.5 * exp(-t/13.5) )
% which is adapted from the paper
% FP Leite, et al.  NeuroImage 16:283-294 (2002)
% http://dx.doi.org/10.1006/nimg.2002.1110
% ** Note that this is a positive function, but MION
% produces a negative response to activation, so the
% beta and t-statistic for MION are usually negative.
% ***** If you want a negative MION function (so you get
% a positive beta), use the name 'MIONN' instead.
% ** After convolution with a square wave 'd' seconds
% long, the resulting single-trial waveform is
% scaled to have magnitude 1.  For example, try
%     this fun command to compare BLOCK and MION:
%     3dDeconvolve -nodata 300 1 -polort -1 -num_stimts 2   \
%     -stim_times 1 '1D: 10 150' 'MION(70)'    \
%     -stim_times 2 '1D: 10 150' 'BLOCK(70,1)' \
%     -x1D stdout: | 1dplot -stdin -one -thick
%     You will see that the MION curve rises and falls
%     much more slowly than the BLOCK curve.
%     ==> ** Note that 'MION(d)' is already convolved with a
%     square wave of duration 'd' seconds.  Do not
%     convolve it again by putting in multiple closely
%     spaced stimulus times (this mistake has been made)!
%     ** Scaling the single-trial waveform to have magnitude
%     1 means that trials with different durations 'd'
%     will have the same magnitude for their regression
%     models.


