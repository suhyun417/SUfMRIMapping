function dat = selectTrials(dat,varargin)

% subdat = selectTrials(dat,key1,val1,key2,val2,...)
%
% Select trials from a movie data structure.
% KEY refers to a column in dat.c
% VALUE refers to the contents within KEY to be searched for.
%
% last modified 2013-apr-13
% dbtm



if mod(length(varargin),2)~=0 % if not iseven
    error('KEY and VALUE input arguments must be paired.');
end

Npairs = length(varargin)/2;
for n=1:Npairs
    key = varargin{2*(n-1)+1};
    value = varargin{2*(n-1)+2};
    dat = binary_select_trials(dat,key,value);
end


function dat = binary_select_trials(dat,key,value)

unpack;
trials = find(ismember(dat.c(:,key),value));
if isfield(dat.h,'analog')
   analogTrid = unique(dat.c(trials,TRID));
   anaInd = find(ismember(dat.h.analog.index,analogTrid));
end

dat.c = dat.c(trials,:);

if key==SIG
    dat.h.snames = dat.h.snames(value,:);
end

if isfield(dat,'s')
    dat.s = dat.s(trials);    
end

if isfield(dat,'blp');
    %dat.blp = dat.blp(trials,:,:);
    dat.blp = dat.blp(anaInd,:,:);
    
end
%
if isfield(dat,'lfp')
    %dat.lfp = dat.lfp(trials,:,:); 
    dat.lfp = dat.lfp(anaInd,:,:);
end

if isfield(dat,'eye');
       dat.eye = dat.eye(anaInd,:,:);
end

dat.h.analog.index = analogTrid;

