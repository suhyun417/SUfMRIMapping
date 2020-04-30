function stat = stats(x,y)

%% TESTING FOR STATISTICAL SIGNIFICANCE

% Cublic spline interpolation of NaNs in regressor data
M = length(x);
nans = find(isnan(x)==1);
xTmp = spline(1:M,x,nans);
x(nans) = xTmp;

% Cross-correlation of x and y at 0 lags
P = xcorr(x,y,0); 

% Randomise phase of X and compute x-corr 100 times
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