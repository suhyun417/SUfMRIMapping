function updateXCorr()

  global DATA GH DSP
  
  
  if isfield(GH,'shiftpos') | ~isempty(GH.shiftpos),
	set(GH.shiftpos,'String',sprintf('%d',DATA.xcorr.shift));
  end

  if ~isfield(DSP,'timecourse') | isempty(DSP.timecourse)
	DSP.timecourse = zeros(1,length(DATA.ftimes));
  end
  
  vox = DSP.timecourse;
  rgr = DATA.stim.curstimmodel;
  lags = DATA.xcorr.lags;
  
  
  if ~isempty(DSP.roi.tcarray)
	voxmean   = mean(DSP.roi.tcarray,2);
	m         = repmat(voxmean,[1 size(DSP.roi.tcarray,2)]);
	pct       = 100*(DSP.roi.tcarray-m)./m;
	tcroimean = mean(pct,1);
  end

  if ~isempty(DSP.sel.tcarray)
	voxmean   = mean(DSP.sel.tcarray,2);
	m         = repmat(voxmean,[1 size(DSP.sel.tcarray,2)]);
	pct       = 100*(DSP.sel.tcarray-m)./m;
	tcselmean = mean(pct,1);
  end


  DATA.valtp = DATA.valtp_dork & DATA.valtp_outlier & DATA.valtp_rgr...
	  & DATA.valtp_epoch;
  
  if sum(DATA.valtp)
	%	fprintf('Warning: selecting %d of %d time points\n',...
	%			sum(DATA.valtp_epoch),length(DATA.valtp_epoch));
	vox(~DATA.valtp) = NaN;
	rgr(~DATA.valtp) = NaN;
	tcroimean(~DATA.valtp) = NaN;
	tcselmean(~DATA.valtp) = NaN;
  end
  

  
  if DSP.highpassFilter_X
	hpc     = DSP.highpassValue_X;
	Fs      = 1/DATA.params.tr;
	vox = filterXCorrTraces(vox,hpc,Fs);
  end
  xc = xcorr_nan(rgr,vox,lags);
  set(GH.voxelXCorr,'YData',xc);

  if nansum(tcselmean)
	mvox = tcselmean;
	if DSP.highpassFilter_X
	  mvox = filterXCorrTraces(mvox,hpc,Fs);
	end
	mxc = xcorr_nan(rgr,mvox,lags);
	set(GH.selectedXCorr,'YData',mxc);
  else
%	disp('No NAN vals');
  end

  if nansum(tcroimean)
	mvox = tcroimean;
	if DSP.highpassFilter_X
	  mvox = filterXCorrTraces(mvox,hpc,Fs);
	end
	mxc = xcorr_nan(rgr,mvox,lags);
	set(GH.roiXCorr,'YData',mxc);
  else
%	disp('No NAN vals');
  end
  
  
%%%%%%%%%
%
% cross correlate, but only for the valid data points
%
function nresult = xcorr_nan(x,y,lags)
nresult = [];

  tot = 0;
  for l=-lags:lags
	tot=tot+1;
	tmp = corrcoef(x([(lags+l+1):(end-lags+l)]),y([(lags+1):(end-lags)]),...
				   'rows', 'pairwise');
	nresult(tot) = tmp(1,2);
  end

%---------------------------
function stat = xcorr_nan_marieke(x,y) % x=LFP regressor, y=MRI timeseries

%% TESTING FOR STATISTICAL SIGNIFICANCE

% Cublic spline interpolation of NaNs in regressor data
nans = find(isnan(x)==1);
xTmp = spline(1:length(x),x,nans);
x(nans) = xTmp;

% Cross-correlation of x and y at 0 lags
P = xcorr(x,y,0); 

% Randomise phase of fft(x) and compute x-corr 100 times
p = zeros(1,100); % peaks
X = fft(x);
R = abs(X); % magnitude of X 
theta = angle(X); % phase of X 
for run = 1:100
    % scramble phase of X from 2:half
    theta_shuff = shuffle(theta(2:(ceil(length(x)/2-1)+1)));
    % mirror-reverse around midpoint
    if rem(length(x),2)
        theta_scramble = [theta(1) theta_shuff -theta_shuff(end:-1:1)];
    else
        theta_scramble = [theta(1) theta_shuff  theta(length(x)/2+1)  -theta_shuff(end:-1:1)];
    end
    % convert back to x and compute cross-correlation again
    X = R.*exp(i*theta_scramble); 
    x = ifft(X);
    x = real(x); % imaginary parts are almost 0; make them completely 0
    p(run) = xcorr(x,y,0);
end
% See how often P occurs by chance
stat = sum(logical(p>=P))/length(p);
%------------------------------------------------------------------



function nresult = xcorr_nan_long(x,y,lags)
  
 
  n1 = length(x);
  n2 = length(y);
  mx = nanmean(x);  % already centered data
  my = nanmean(y);
  
  for l=-lags:lags
	xysum = 0;
	xytot = 0;
	for i=lags+1:n1-lags
	  if ~isnan(x(i)) & ~isnan(y(i-l))
		corrpt = (x(i)-mx)*(y(i-l)-my);
		xysum = xysum+corrpt;
		xytot = xytot+1;
	  end
	end
	result(lags+1+l) = xysum/xytot;
  end
  
  % autocorrs
  for l=-lags:lags
	xxsum = 0;
	xxtot = 0;
	for i=lags+1:n1-lags
	  if ~isnan(x(i))
		xxsum = xxsum+(x(i)-mx)^2;
		xxtot = xxtot+1;
	  end
	end
	acorr1(lags+1+l) = sqrt(xxsum/xxtot);
  end
  for l=-lags:lags
	yysum = 0;
	yytot = 0;
	for i=lags+1:n2-lags
	  if ~isnan(y(i-l))
		yysum = yysum+(y(i-l)-my)^2;
		yytot = yytot+1;
	  end
	end
	acorr2(lags+1+l) = sqrt(yysum/yytot);
  end

  normxc = acorr1.*acorr2;
  nresult = result./normxc;




