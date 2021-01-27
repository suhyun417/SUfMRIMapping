function calculateFunctionalMap()

  global DATA DSP

  
  if ~isfield(DATA.stim,'curmodelname') | isempty(DATA.stim.curmodelname), return; end

  
  voltc = DSP.proc.fmri_tc_3d;
  stattype   = DATA.stim.curstattype;

  
  if strcmp(DATA.stim.curmodelname,'REGRESSOR')
	
	DATA.stim.curstimmodel = DATA.stim.stimmodel;
	
	%
	% If the movie is shown starting in the middle of the scan, adjust
    % the offset 
	%
	
	% Note that comparing the 
	if length(DATA.ftimes) > length(DATA.stim.curstimmodel)
	  off = DATA.stim.stimon_offsett_TR;	  
	  ntp = length(DATA.stim.curstimmodel);
	  voltc = DSP.proc.fmri_tc_3d(:,:,:,[off:(off+ntp-1)]);
	  mnvoltc = repmat(squeeze(mean(voltc,4)),[1 1 1 size(voltc,4)]);
	  voltc = voltc-mnvoltc;
	  dsc = 1;
	else
	  voltc = DSP.proc.fmri_tc_3d;
	  dsc = DSP.TR_Discard;
	end
	
	wb = waitbar(0, 'Computing Correlation with Regressor');
	
	for i=1:size(voltc,1)
	  for j=1:size(voltc,2)
		for k=1:size(voltc,3)
		  tmp = corrcoef(squeeze(voltc(i,j,k,[dsc:end])),...
						 DATA.stim.curstimmodel([dsc:end]));
		  vals(i,j,k) = tmp(1,2);
		end
	  end
	  waitbar(i/size(voltc,1),wb);
	end
	close(wb);

  else

	% For computing t-scores, and various otherthings etc.
	%
	
	modelname  = DATA.stim.curmodelname;
	
	nmodels = size(DATA.stim.stimmodel,1);
	
	if nmodels > 1
	  for i=1:length(DATA.stim.modelnames)
		if strcmp(DATA.stim.modelnames{i},modelname), val = i; break; end
	  end	  
	  DATA.stim.curstimmodel = DATA.stim.stimmodel(val,:);
	else
	  DATA.stim.curstimmodel = DATA.stim.stimmodel;
	end
	
	model = DATA.stim.curstimmodel; 
	
	ndatatp = size(DSP.proc.fmri_tc_3d,4);
	if ~isfield(DATA,'valtp') | length(DATA.valtp)~=ndatatp
	  DATA.valtp = ones(1,size(DSP.proc.fmri_tc_3d,4)); 
	end

	if ndatatp < length(model) 
	  model(ndatatp+1:end) = [];
	end

	
	d = size(voltc);
	x0 = d(2);
	y0 = d(1);
	z0 = d(3);
	nt = d(4);
	
	vals = zeros(y0,x0,z0);
	
	if strcmp(stattype,'t-stat')
	  
	  % When there are more conditions, it computes a t-score between any
	  % two conditions
	  
	  % First group together the lightbox time points into two groups,
	  % corresponding to the desired stimulus conditions	
	  apts = find(model == 0.5 & DATA.valtp);
	  bpts = find(model == -0.5 & DATA.valtp);
	  
	  
	  disp('Computing T-Statistic map');
	  
	  wb = waitbar(0, 'Computing T-map');
	  for z=1:z0
		cmatrx = squeeze(voltc(:,:,z,([apts bpts])));
		s   = std(cmatrx,0,3);
		cmatrx = [];
		
		amatrx = squeeze(voltc(:,:,z,apts));
		amn = mean(amatrx,3);
		al  = size(amatrx,3);
		amatrx = [];
		
	  bmatrx = squeeze(voltc(:,:,z,bpts));
	  bmn = mean(bmatrx,3);
	  bl  = size(bmatrx,3);
	  bmatrx = [];
	  
	  % Here computing the t-stat directly (embarassingly,..
	  % from the matlab ttest2 help).
	  %
	  %  T = (amn-bmn)/(s(sqrt((1/al)+(1/bl))))
	  %
	  % where amn = mean value for a voxel from amatrx
	  %       bmn = mean value for a voxel from bmatrx
	  %       s = stddev for a voxel across both conditions
	  %       bl = number of pts in b condition
	  %       al = number off pts in a condition
	  %
	  vals(:,:,z) = (amn-bmn)./(s*sqrt((1/al)+(1/bl)));
	  waitbar(z/d(3),wb);
	  end
	  close(wb);
	elseif strcmp(stattype,'pct')
	  
	  apts = find(model == 0.5 & DATA.valtp);
	  bpts = find(model == -0.5 & DATA.valtp);
	  
	  
	  disp('PCT map');
	  
	  wb = waitbar(0, 'Computing PCT map');
	  for z=1:z0

		amatrx = squeeze(voltc(:,:,z,apts));
		amn = mean(amatrx,3);
		al  = size(amatrx,3);
		amatrx = [];
		
		bmatrx = squeeze(voltc(:,:,z,bpts));
		bmn = mean(bmatrx,3);
		bl  = size(bmatrx,3);
		bmatrx = [];
	  
		vals(:,:,z) = 100*2*(amn-bmn)./(amn+bmn);
		waitbar(z/d(3),wb);
	  end
	  close(wb);
	elseif strcmp(stattype,'diff')
	  
	  apts = find(model == 0.5 & DATA.valtp);
	  bpts = find(model == -0.5 & DATA.valtp);
	  
	  
	  disp('DIFF map');
	  
	  wb = waitbar(0, 'Computing DIFF map');
	  for z=1:z0

		amatrx = squeeze(voltc(:,:,z,apts));
		amn = mean(amatrx,3);
		al  = size(amatrx,3);
		amatrx = [];
		
		bmatrx = squeeze(voltc(:,:,z,bpts));
		bmn = mean(bmatrx,3);
		bl  = size(bmatrx,3);
		bmatrx = [];
	  
		vals(:,:,z) = amn-bmn;
		waitbar(z/d(3),wb);
	  end
	  close(wb);
	else
	  disp(sprintf('WARNING: Unknown levels stat type %s',stattype));
	end
  end
  
  DSP.proc.modmap_3d = vals;
