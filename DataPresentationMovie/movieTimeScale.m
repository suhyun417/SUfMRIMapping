function dat = movieTimeScale(dat,tscale);
%
%  dat = movieTimeScale(dat,tscale);
% 
% Sets movie time scale to seconds or milliseconds.
%
% Input argments:
% DAT movie data structure
% TSCALE character string: 'sec' or 'ms'
%
% last modified 2013-apr-13
% dbtm

tscale = lower(tscale(1));
if isequal(tscale,'m')
    units = 'ms';
elseif isequal(tscale,'s')
    units = 'sec';
else
    error('Invalid unit input argument');
    help movieTimeScale;
end

currentUnits = dat.h.units;
if isequal(units,currentUnits)
    return;
elseif isequal(units,'sec') & isequal(currentUnits,'ms')
    for s=1:length(dat.s)
       dat.s{s} = dat.s{s}/1000;
    end
    dat.t = dat.t/1000;
    %dat.lfp/1000;
    %dat.blp/1000;
    dat.h.units = 'sec';
elseif isequal(units,'ms') & isequal(currentUnits,'sec')
    for s=1:length(dat.s)
       dat.s{s} = dat.s{s}*1000;
    end
    dat.t = dat.t*1000;
    %dat.lfp*1000;
    %dat.blp*1000;
    dat.h.units = 'ms';
else
    error('Invalid header field dat.h.units! Must be "ms" or "sec".');
end