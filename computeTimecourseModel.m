function computeTimecourseModel()

  global DATA GH DSP

  if ~isfield(DATA,'stim'), return; end
  upat = unique(DATA.stim.stimtype);

  [TR,NR] = getTRandNR();
  
  if TR >=8, 
	DATA.stim.sparse = 1; 
  else
	DATA.stim.sparse = 0;
  end
  nupat = length(upat);  % number of unique patterns
    
  %  hrf = gampdf([0:TR:30],3,3);  % hemodynamic response function
  %  hrf = gampdf([0:TR:30],3,2);  % hemodynamic response function
  hrf = gampdf([0:TR:15],2,3);  % hemodynamic response function

  DATA.stim.hrf      = hrf/sum(hrf);
  modeldefault = 1;

  DATA.stim.stimtc = zeros(100, NR); %set it too 100, default # is too waste of space
  DATA.stim.stimmodel = [];
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Here we make the two different kinds of models.
  %


  if strcmp(DATA.stim.curstattype,'t-stat') |...
		strcmp(DATA.stim.curstattype,'pct') | ...
		strcmp(DATA.stim.curstattype,'diff')
	
	% The FIRST model is the "levels" model, which considers time points
	% for which the model, following HRF convolution, is > 60% of the way
	% to its maximum value
	%
	
	maxmodels = 20;
	stimmodel = nan(maxmodels,NR);
	trace = 0;
	for i=1:nupat
	  for j=1:nupat
		%
		% This is for computing T-tests when comparing two conditions out
		% of many
		%
		if j<i
		  trace = trace+1;
		  if isfield(DATA.stim, 'sparse') & DATA.stim.sparse
			% make model without reference to hrf
			%
			DATA.stim.stimtc = [];  % no need to compute analog model
			posi = find(DATA.stim.alltrstim == upat(i));
			negi = find(DATA.stim.alltrstim == upat(j));
			stimmodel(trace,posi) = 0.5;
			stimmodel(trace,negi) = -0.5;
		  else
			% number of instances that stim type i turns on
			nrep = length(find(DATA.stim.stimtype == upat(i)));
			% what times (in sec) does stim i turn on
			valont = DATA.stim.stimont(DATA.stim.stimtype == upat(i));  
			% what times (in sec) does stim i turn off
			valofft = DATA.stim.stimofft(DATA.stim.stimtype == upat(i));
			
			% add 1's for those TR's that lie between on and off
			for k=1:nrep
			  trl = (DATA.ftimes >= valont(k)) & (DATA.ftimes <= valofft(k));
			  DATA.stim.stimtc(trace,:) = DATA.stim.stimtc(trace,:) + trl;
			end
			
			% number of instances that stim type j turns on
			nrep = length(find(DATA.stim.stimtype == upat(j)));
			% what times (in sec) does stim j turn on
			valont = DATA.stim.stimont(DATA.stim.stimtype == upat(j));
			% what times (in sec) does stim j turn off
			valofft = DATA.stim.stimofft(DATA.stim.stimtype == upat(j));
			
			% add -1's for those TR's that lie between on and off
			for k=1:nrep
			  trl = (DATA.ftimes >= valont(k)) & (DATA.ftimes <= valofft(k));
			  DATA.stim.stimtc(trace,:) = DATA.stim.stimtc(trace,:) - trl;
			end
			
			model= conv(DATA.stim.stimtc(trace,:),DATA.stim.hrf);
			
			% change range from -.5 to +.5
			tmp = model([1:NR])/2; 
			
			% excluse values within middle 60%
			tmp(abs(tmp) < 0.3) = NaN;
			
			% set to max
			tmp(tmp >= 0.3) = 0.5;
			tmp(tmp <= -0.3) = -0.5;
			
			% populate models with 0.5, -0.5, and NaN
			stimmodel(trace,:) = tmp;
		  end
		end
	  end
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Then set the labels for the popup menu
	%
	
	
	txt = {};
	tot = 0;
	for i=1:nupat
	  for j=1:nupat
		if j < i
		  tot = tot+1;
		  %		  txt{tot} = sprintf('S%dvS%d',i-1,j-1);
		  txt{tot} = sprintf('S%dvS%d',upat(i),upat(j));
		end
	  end
	end
	
	%%%%
	%  Finally add the "All except zero" model
	%
	trace = trace+1;
	
	if isfield(DATA.stim, 'sparse') & DATA.stim.sparse
	  
	  posi = find(DATA.stim.alltrstim > 0);
	  negi = find(DATA.stim.alltrstim == 0);
	  stimmodel(trace,posi) = 0.5;
	  stimmodel(trace,negi) = -0.5;
	else
	  % for all stim types that are nonzero
	  nrep = length(find(DATA.stim.stimtype > 0));
	  
	  % designate on and off times as before
	  valont = DATA.stim.stimont(DATA.stim.stimtype > 0 );
	  valofft = DATA.stim.stimofft(DATA.stim.stimtype > 0);
	  
	  % as before, now for all nonzero stim combined
	  for k=1:nrep
		trl = (DATA.ftimes >= valont(k)) & (DATA.ftimes <= valofft(k));
		DATA.stim.stimtc(trace,:) = DATA.stim.stimtc(trace,:) + trl;
	  end
	  
	  % now for all zero stim
	  nrep = length(find(DATA.stim.stimtype == 0));
	  valont = DATA.stim.stimont(DATA.stim.stimtype == 0);
	  valofft = DATA.stim.stimofft(DATA.stim.stimtype == 0);
	  for k=1:nrep
		trl = (DATA.ftimes >= valont(k)) & (DATA.ftimes <= valofft(k));
		DATA.stim.stimtc(trace,:) = DATA.stim.stimtc(trace,:) - trl;
	  end
	  
	  model= conv(DATA.stim.stimtc(trace,:),DATA.stim.hrf);
	  tmp = model([1:NR])/2;
	  tmp(abs(tmp) < 0.3) = NaN;
	  tmp(tmp >= 0.3) = 0.5;
	  tmp(tmp <= -0.3) = -0.5;
	  stimmodel(trace,:) = tmp;
	end
	tot = tot+1;
	txt{tot} = 'ALLv0';
	stimmodel([trace+1:end],:) = [];
	DATA.stim.stimtc([trace+1:end],:) = [];
  else
	% for future work
  end	
  DATA.stim.stimmodel = stimmodel;
  DATA.stim.modelnames = txt;
  DATA.stim.curstimmodel = DATA.stim.stimmodel(modeldefault,:);
  DATA.stim.curmodelname =DATA.stim.modelnames{modeldefault};
  DATA.stim.curmodelindx = modeldefault;
  set(GH.main.menu.selectmodel,'String',DATA.stim.modelnames);  
  
  DSP.recomputeMap = 1;
  
