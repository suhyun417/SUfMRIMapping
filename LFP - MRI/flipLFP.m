function flipLFP
% This function flips an LFP matrix upside-down, i.e. reverses the channel 
% index. Its output is an LFP matrix of the same size as the original one, 
% i.e. the first dimension represents timepoints, the second epoch (segment) 
% and the third channel. MLS July 2008

date = '14-07-08';  date2 = '140708';  sess = '_0001'; monkey = 'Aspen'; 

cd(['/einstein0/USRlab/projects/scholvinckm/data/' monkey '/inside scanner/' date '/Matlab']);
load([date2 sess '_LFP_epochs']);

LFP = flipdim(LFP,3);

eval(['save ' date2 sess '_LFP_epochs LFP remove']);
