function sdf_out = getSDF_SHP(ts, win, sd)

% sd = 1000; %2000; %1000; %2000;
k  = normpdf(-3*sd:+3*sd,0,sd);
k = k./sum(k);
upsamp = 1000;
%   TR = 2.4;
%   MION_k = gampdf(-3*TR*upsamp:3*TR*upsamp,TR,2*TR);
%   MION_k = MION_k./sum(MION_k);

if max(win)<1000 % if window is in seconds, we need to upsample it
    % it's not the exact way to do this, but just for now..
    timeline = zeros(upsamp*win(2)-upsamp*win(1),1);
    its = (ceil(ts*upsamp));
    timeline(its) = 1;    
    
else % units are milliseconds
    timeline = zeros(win(2)-win(1),1);
    its = ceil(ts);
    timeline(its) = 1;
    
end
    
% convolve, scale to spikes/sec, and decimate
sdf_out = doConv(timeline,k).*upsamp; %doConv(timeline, MION_k).*upsamp; %  