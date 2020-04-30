function ALLRGR = SU_addSpecificRegressors_Compare(r,TR)
  
%% All REGRESSORS ALL FORMATS!!!
  
% NAMING SCHEME
%   sum     Sum the categories
%   log     Take the log of the values
%   z       Zscore the values
%   sqrt    Take the square root of the values
%   sm#     smooth the Regressor with a # kernel
%   c#      compress the values by # (divide by #)
%   lp#     Low Pass with max being #
%   hp#     High Pass with min being #
%   bp#_#   Band Pass with min_max 
%   ave     Average the categories
%   ms#     Split the data at # becomes 0's and 1's
%   msMN    Split the data at the mean of the variable
%   msMD    Split the data at the median of the variable
  
%  ORDER OF FEATURES
%  1   =  Speed(log_MION)
%  2   =  Mot_Diverg_Rect(log_MION)
%  3   =  Mot_Diverg_STD(log_MION)
%  4   =  Mot_Diverg_Beta(log_MION)
%  5   =  Mot_Diverg_Gamma(MION_sqrt_sm2)
%  6   =  Scene_Cuts(log_MION_sm2)
%  7   =  Faces_Full(sqrt_MION)
%  8   =  Faces_SV(sqrt_MION)
%  9   =  1_Face_Full(sqrt_MION)
%  10  =  1_Face_SV(sqrt_MION)
%  11  =  Heads(sqrt_MION)
%  12  =  Faces_all(sqrt_MION)
%  13  =  1_Face_comb(sqrt_MION)
%  14  =  Luminance(log_MION)
%  15  =  Contrast_STD(log_MION)
%  16  =  Contrast_Beta(log_MION)
%  17  =  Contrast_Gamma(log_MION)
%  18  =  SF_low(log_MION)
%  19  =  SF_high(log_MION)
%  20  =  SF_ratio(log_MION)
%  21  =  Butts(MION_sqrt)
%  22  =  Bodies(sqrt_MION)
%  23  =  Extremities(sqrt_MION)
%  24  =  Hands(sqrt_MION)
%  25  =  Animals(MION_sqrt)
%  26  =  Macaque(MION_sqrt)
%  27  =  Rhesus(MION)
%  28  =  Human(MION)
%  29  =  Mot_Local_STD(MION)
%  30  =  Mot_Bartel_Local(MION)
%  31  =  Mot_Bartel_Global(MION)
%  32  =  Mot_Bartel_Resid(MION)
%  33  =  Mot_Bartel_Total(MION)
%  34  =  Aggression(sqrt_MION)
%  35  =  Affiliative(sqrt_MION)
%  36  =  Food(sqrt_MION)
%  37  =  Aggr_Play(sqrt_MION)
%  38  =  Dyadic(sqrt_MION)
%  39  =  Saccades(sqrt_MION)
    
  %%
  NSTD = 5;
  R = [];
  
  % set up MION function for convolution
  tvals = [(-20*TR):TR:(20*TR)];
	% This is just made-up
  MION_k = gampdf(tvals,TR,2*TR);
  MION_k = MION_k./sum(MION_k);
%   crgrs = doConv(rgrs,MION_k);
  
  % set up standard smoothing kernel
  smooth_k = normpdf([-10:10],0,2);  
  smooth_k = smooth_k/sum(smooth_k);
  
  %%%%%%%%%%%%%%%%
  %% GLOBAL MOTION
  % Add the speed regressor
  speed_indx = find(strcmp('Speed',r.names));
  speed_rgr=r.rgrs(:,speed_indx);
  
  speed4=real(log(speed_rgr));          % compress with a log function 
  speed4(isinf(speed4))=0;              % make sure we have no INF numbers
  speed5=doConv(speed4,MION_k)';        % apply MION function to raw data

  R = addRegressor(R,speed5,'Speed(log_MION)');
  
  clear speed*

  %%%%%%%%%%%
  %% Motion Divergence
  dmot_indx = find(strcmp('Motion Div Rectify',r.names));
  dmot_rgr=r.rgrs(:,dmot_indx);

  dmot4=real(log(dmot_rgr));          % compress with a log function 
  dmot4(isinf(dmot4))=0;              % make sure we have no INF numbers
  dmot5=doConv(dmot4,MION_k)';        % apply MION function to raw data

  R = addRegressor(R,dmot5,'Mot_Diverg_Rect(log_MION)');

  clear dmot*
                    
  %%%%%%%%%%%
  %% Motion Divergence
  dmot_indx = find(strcmp('Motion Div STD',r.names));
  dmot_rgr=r.rgrs(:,dmot_indx);

  dmot4=real(log(dmot_rgr));          % compress with a log function 
  dmot4(isinf(dmot4))=0;              % make sure we have no INF numbers
  dmot5=doConv(dmot4,MION_k)';        % apply MION function to raw data

  R = addRegressor(R,dmot5,'Mot_Diverg_STD(log_MION)');

  clear dmot*
  
  %%%%%%%%%%%
  %% Motion Divergence
  dmot_indx = find(strcmp('Motion Div Beta',r.names));
  dmot_rgr=r.rgrs(:,dmot_indx);

  dmot4=real(log(dmot_rgr));          % compress with a log function 
  dmot4(isinf(dmot4))=0;              % make sure we have no INF numbers
  dmot5=doConv(dmot4,MION_k)';        % apply MION function to raw data
  dmot6=doConv(dmot5,smooth_k)';      % apply a 2 standard deviation smoothing kernel

  R = addRegressor(R,dmot5,'Mot_Diverg_Beta(log_MION)');

  clear dmot*
  
  %%%%%%%%%%%
  %% Motion Divergence
  dmot_indx = find(strcmp('Motion Div Gamma',r.names));
  dmot_rgr=r.rgrs(:,dmot_indx);
  
  dmot7=doConv(dmot_rgr,MION_k)';     % apply MION function to raw data
  dmot7=sqrt(dmot7);                  % compress with a sqrt function (and make sure we have no imaginary numbers
  tmp=imag(conj(dmot7));              % and make sure we have no imaginary numbers
  dmot7(tmp~=0)=tmp(tmp~=0);          % fill the now negative numbers into the vector
  dmot7=doConv(dmot7,smooth_k)';      % apply a 2 standard deviation smoothing kernel

  R = addRegressor(R,dmot7,'Mot_Diverg_Gamma(MION_sqrt_sm2)');
  
  clear dmot*
  
  %%%%%%%%%%%
  %% Scene Cuts
  SCmot_indx = find(strcmp('Scene Cuts',r.names));
  SCmot_rgr=r.rgrs(:,SCmot_indx);
  
  SCmot4=real(log(SCmot_rgr));          % compress with a log function (and make sure we have no imaginary numbers
  SCmot4(isinf(SCmot4))=0;
  SCmot5=doConv(SCmot4,MION_k)';        % apply MION function to raw data
%   SCmot6=doConv(SCmot5,smooth_k)';      % apply a 2 standard deviation smoothing kernel
 
  R = addRegressor(R,SCmot5,'Scene_Cuts(log_MION)');
  
  clear SCmot*
  ALLRGR = R;

  %%%%%%%%%%%%%%%%
  %% Fullview Face 
  faces_indx = find(strcmp('Faces (full)', r.names));
  faces_rgr  = r.rgrs(:,faces_indx);
  
  faces4=sqrt(faces_rgr);               % compress with a log function
  tmp=imag(conj(faces4));               % and make sure we have no imaginary numbers
  faces4(tmp~=0)=tmp(tmp~=0);           % fill the now negative numbers into the vector
  faces5=doConv(faces4,MION_k)';        % apply MION function to raw data
 
  R = addRegressor(R,faces5,'Faces_Full(sqrt_MION)');
 
  ALLRGR = R;
  clear faces*
  
  %%%%%%%%%%%%%%%%
  %% Sideview Face 
  faces_indx = find(strcmp('Faces (side view)', r.names));
  faces_rgr  = r.rgrs(:,faces_indx);

  faces4=sqrt(faces_rgr);               % compress with a square root function
  tmp=imag(conj(faces4));               % and make sure we have no imaginary numbers
  faces4(tmp~=0)=tmp(tmp~=0);           % fill the now negative numbers into the vector
  faces5=doConv(faces4,MION_k)';        % apply MION function to raw data

  R = addRegressor(R,faces5,'Faces_SV(sqrt_MION)');

  ALLRGR = R;
  clear faces*
  
  %%%%%%%%%%%%%%%%
  %% 1 Fullview Face 
  faces_indx = find(strcmp('One Face (full)', r.names));
  faces_rgr  = r.rgrs(:,faces_indx);

  faces4=sqrt(faces_rgr);               % compress with a log function
  tmp=imag(conj(faces4));               % and make sure we have no imaginary numbers
  faces4(tmp~=0)=tmp(tmp~=0);           % fill the now negative numbers into the vector
  faces5=doConv(faces4,MION_k)';        % apply MION function to raw data

  R = addRegressor(R,faces5,'1_Face_Full(sqrt_MION)');

  ALLRGR = R;
  clear faces*
  
  %%%%%%%%%%%%%%%%
  %% 1 Sideview Face 
  faces_indx = find(strcmp('One Face (side view)', r.names));
  faces_rgr  = r.rgrs(:,faces_indx);
 
  faces4=sqrt(faces_rgr);               % compress with a log function
  tmp=imag(conj(faces4));               % and make sure we have no imaginary numbers
  faces4(tmp~=0)=tmp(tmp~=0);           % fill the now negative numbers into the vector
  faces5=doConv(faces4,MION_k)';        % apply MION function to raw data
 
  R = addRegressor(R,faces5,'1_Face_SV(sqrt_MION)');
 
  ALLRGR = R;
  clear faces*
  
  %%%%%%%%%%%%%%%%
  %% HEADS  
  heads_indx = find(strcmp('Heads',r.names));
  heads_rgr  = r.rgrs(:,heads_indx);

  heads4=sqrt(heads_rgr);               % compress with a square root function
  tmp=imag(conj(heads4));               % and make sure we have no imaginary numbers
  heads4(tmp~=0)=tmp(tmp~=0);           % fill the now negative numbers into the vector
  heads5=doConv(heads4,MION_k)';        % apply MION function to raw data

  R = addRegressor(R,heads5,'Heads(sqrt_MION)');

  ALLRGR = R;
  clear heads*
  
  %%%%%%%%%%%%%%%%
  %% all Faces
  faces_indx = find(strncmp('Faces', r.names,5));
  faces_rgr  = sum(r.rgrs(:,faces_indx),2);

  faces4=sqrt(faces_rgr);               % compress with a square root function
  tmp=imag(conj(faces4));               % and make sure we have no imaginary numbers
  faces4(tmp~=0)=tmp(tmp~=0);           % fill the now negative numbers into the vector
  faces5=doConv(faces4,MION_k)';        % apply MION function to raw data

  R = addRegressor(R,faces5,'Faces_all(sqrt_MION)');

  ALLRGR = R;
  clear faces*
  
  %%%%%%%%%%%%%%%%
  %% One Face Fullview or Profile
  faces_indx = find(strncmp('One Face', r.names,8));
  faces_rgr  = sum(r.rgrs(:,faces_indx),2);

  faces4=sqrt(faces_rgr);               % compress with a square root function
  tmp=imag(conj(faces4));               % and make sure we have no imaginary numbers
  faces4(tmp~=0)=tmp(tmp~=0);           % fill the now negative numbers into the vector
  faces5=doConv(faces4,MION_k)';        % apply MION function to raw data

  R = addRegressor(R,faces5,'1_Face_comb(sqrt_MION)');
  
  ALLRGR = R;
  clear faces*
  %%%%%%%%%%%
  %% Luminance
  lum_indx = find(strcmp('Luminance',r.names));
  lum_rgr=r.rgrs(:,lum_indx);
 
  lum4=real(log(lum_rgr));          % compress with a log function 
  lum4(isinf(lum4))=0;              % make sure we have no INF numbers
  lum5=doConv(lum4,MION_k)';        % apply MION function to raw data
 
  R = addRegressor(R,lum5,'Luminance(log_MION)');
 
  clear lum*
  
  %%%%%%%%%%%
  %% Contrast 
  cont_indx = find(strcmp('Contrast',r.names));
  cont_rgr=r.rgrs(:,cont_indx);
  
  cont4=real(log(cont_rgr));          % compress with a log function 
  cont4(isinf(cont4))=0;              % make sure we have no INF numbers
  cont5=doConv(cont4,MION_k)';        % apply MION function to raw data
  
  R = addRegressor(R,cont5,'Contrast_STD(log_MION)');
 
  clear cont*
  
  %%%%%%%%%%%
  %% Contrast Beta
  cont_indx = find(strcmp('Beta Contrast',r.names));
  cont_rgr=r.rgrs(:,cont_indx);
  
  cont4=real(log(cont_rgr));          % compress with a log function 
  cont4(isinf(cont4))=0;              % make sure we have no INF numbers
  cont5=doConv(cont4,MION_k)';        % apply MION function to raw data
  
  R = addRegressor(R,cont5,'Contrast_Beta(log_MION)');
 
  clear cont*
  
  %%%%%%%%%%%
  %% Contrast Gamma
  cont_indx = find(strcmp('Gamma Contrast',r.names));
  cont_rgr=r.rgrs(:,cont_indx);
 
  cont4=real(log(cont_rgr));          % compress with a log function 
  cont4(isinf(cont4))=0;              % make sure we have no INF numbers
  cont5=doConv(cont4,MION_k)';        % apply MION function to raw data
  
  R = addRegressor(R,cont5,'Contrast_Gamma(log_MION)');
 
  clear cont*
  
  %%%%%%%%%%%
  %% Spatial Frequency (low)
  low_range  = [0 0.2];
  low_sf_feat  = sprintf('SF (%.1f to %.1f)',low_range(1), low_range(2));
  
  sff_indx = find(strcmp(low_sf_feat,r.names));
  sff_rgr=r.rgrs(:,sff_indx);
 
  sff4=real(log(sff_rgr));          % compress with a log function 
  sff4(isinf(sff4))=0;              % make sure we have no INF numbers
  sff5=doConv(sff4,MION_k)';        % apply MION function to raw data

  R = addRegressor(R,sff5,'SF_low(log_MION)');
  
  clear sff*
  
  %%%%%%%%%%%
  %% Spatial Frequency (high)
  high_range = [3 100];
  high_sf_feat = sprintf('SF (%.1f to %.1f)',high_range(1), high_range(2));
  
  sff_indx = find(strcmp(high_sf_feat,r.names));
  sff_rgr=r.rgrs(:,sff_indx);
  
  sff4=real(log(sff_rgr));          % compress with a log function 
  sff4(isinf(sff4))=0;              % make sure we have no INF numbers
  sff5=doConv(sff4,MION_k)';        % apply MION function to raw data
  
  R = addRegressor(R,sff5,'SF_high(log_MION)');
  
  clear sff*
  
  %% Spatial Frequency (ratio)
  sff_indx = find(strcmp('SF Ratio',r.names));
  sff_rgr=r.rgrs(:,sff_indx);

  sff4=real(log(sff_rgr));          % compress with a log function 
  sff4(isinf(sff4))=0;              % make sure we have no INF numbers
  sff5=doConv(sff4,MION_k)';        % apply MION function to raw data

  R = addRegressor(R,sff5,'SF_ratio(log_MION)');
 
  clear sff*
  
  %%%%%%%%%%%%%%%%
  %% Butts
  butts_indx = find(strcmp('Butts',r.names));
  butts_rgr  = r.rgrs(:,butts_indx);
  
  butts1=doConv(butts_rgr,MION_k)';     % apply MION function to raw data
  butts2=sqrt(butts1);                  % compress with a log function
  tmp=imag(conj(butts2));               % and make sure we have no imaginary numbers
  butts2(tmp~=0)=tmp(tmp~=0);           % fill the now negative numbers into the vector
  
  R = addRegressor(R,butts2,'Butts(MION_sqrt)');

  clear butts*
  
  %%%%%%%%%%%%%%%%
  %% Bodies 
  body_indx = find(strcmp('Hands',r.names) | strcmp('Feet',r.names) | strcmp('Butts',r.names));
  body_rgr  = sum(r.rgrs(:,body_indx),2);
 
  body4=sqrt(body_rgr);                % compress with a log function
  tmp=imag(conj(body4));               % and make sure we have no imaginary numbers
  body4(tmp~=0)=tmp(tmp~=0);           % fill the now negative numbers into the vector
  body5=doConv(body4,MION_k)';         % apply MION function to raw data
 
  R = addRegressor(R,body5,'Bodies(sqrt_MION)');
 
  clear body*
  
  %%%%%%%%%%%%%%%%
  %% Extremities (Hands & Feet)
  body_indx = find(strcmp('Hands',r.names) | strcmp('Feet',r.names));
  body_rgr  = sum(r.rgrs(:,body_indx),2);
  
  body4=sqrt(body_rgr);                % compress with a log function
  tmp=imag(conj(body4));               % and make sure we have no imaginary numbers
  body4(tmp~=0)=tmp(tmp~=0);           % fill the now negative numbers into the vector
  body5=doConv(body4,MION_k)';         % apply MION function to raw data
  
  R = addRegressor(R,body5,'Extremities(sqrt_MION)');
 
  clear body*
  
  %%%%%%%%%%%%%%%%
  %% Hands 
  body_indx = find(strcmp('Hands',r.names));
  body_rgr  = r.rgrs(:,body_indx);
  
  body4=sqrt(body_rgr);                % compress with a log function
  tmp=imag(conj(body4));               % and make sure we have no imaginary numbers
  body4(tmp~=0)=tmp(tmp~=0);           % fill the now negative numbers into the vector
  body5=doConv(body4,MION_k)';         % apply MION function to raw data
  
  R = addRegressor(R,body5,'Hands(sqrt_MION)');
  
  clear body*
  
  %%%%%%%%%%%%
  %% Animals
  ani_indx = find(strcmp('conspecifics',r.names) |... 
                 strncmp('heterospecific', r.names,14) |...
                 strcmp('humans',r.names));
  ani_rgr  = sum(r.rgrs(:,ani_indx),2);
  
  ani1=doConv(ani_rgr,MION_k)';       % apply MION function to raw data
  ani2=sqrt(ani1);                    % compress with a log function
  tmp=imag(conj(ani2));               % and make sure we have no imaginary numbers
  ani2(tmp~=0)=tmp(tmp~=0);           % fill the now negative numbers into the vector
  
  R = addRegressor(R,ani2,'Animals(MION_sqrt)');
 
  clear ani*
  
  %%%%%%%%%%%%
  %% Macaques
  ani_indx = find(strcmp('conspecifics',r.names) | strcmp('heterospecific macaques', r.names));
  ani_rgr  = sum(r.rgrs(:,ani_indx),2);
  
  ani1=doConv(ani_rgr,MION_k)';       % apply MION function to raw data
  ani2=sqrt(ani1);                    % compress with a log function
  tmp=imag(conj(ani2));               % and make sure we have no imaginary numbers
  ani2(tmp~=0)=tmp(tmp~=0);           % fill the now negative numbers into the vector
  
  R = addRegressor(R,ani2,'Macaque(MION_sqrt)');
  
  clear ani*
  
  %%%%%%%%%%%%
  %% Rhesus
  ani_indx = find(strcmp('conspecifics',r.names));
  ani_rgr  = r.rgrs(:,ani_indx);
  
  ani1=doConv(ani_rgr,MION_k)';       % apply MION function to raw data
  
  R = addRegressor(R,ani1,'Rhesus(MION)');
  
  clear ani*
  
  %%%%%%%%%%%%
  %% Human
  ani_indx = find(strcmp('humans',r.names));
  ani_rgr  = r.rgrs(:,ani_indx);
  
  ani1=doConv(ani_rgr,MION_k)';       % apply MION function to raw data

  R = addRegressor(R,ani1,'Human(MION)');

  clear ani*
  
  %%%%%%%%%%%%
  %% Local Motion (ITTI)
  dmot_indx = find(strcmp('Local Motion STD',r.names));
  dmot_rgr=r.rgrs(:,dmot_indx);
  
  dmot1=doConv(dmot_rgr,MION_k)';     % apply MION function to raw data
  
  R = addRegressor(R,dmot1,'Mot_Local_STD(MION)');
  
  clear dmot*
  
    %%%%%%%%%%%%
  %% BARTEL Local Motion 
  dmot_indx = find(strcmp('Bartel Local',r.names));
  dmot_rgr=r.rgrs(:,dmot_indx);
  
  dmot1=doConv(dmot_rgr,MION_k)';     % apply MION function to raw data
  
  R = addRegressor(R,dmot1,'Mot_Bartel_Local(MION)');
  
  clear dmot*
  
    %%%%%%%%%%%%
  %% BARTEL Global Motion 
  dmot_indx = find(strcmp('Bartel Global',r.names));
  dmot_rgr=r.rgrs(:,dmot_indx);
  
  dmot1=doConv(dmot_rgr,MION_k)';     % apply MION function to raw data
  
  R = addRegressor(R,dmot1,'Mot_Bartel_Global(MION)');
  
  clear dmot*
  
    %%%%%%%%%%%%
  %% BARTEL Residual Motion 
  dmot_indx = find(strcmp('Bartel Resid',r.names));
  dmot_rgr=r.rgrs(:,dmot_indx);
  
  dmot1=doConv(dmot_rgr,MION_k)';     % apply MION function to raw data
  
  R = addRegressor(R,dmot1,'Mot_Bartel_Resid(MION)');
  
  clear dmot*
  
    %%%%%%%%%%%%
  %% BARTEL Total Motion
  dmot_indx = find(strcmp('Bartel Total',r.names));
  dmot_rgr=r.rgrs(:,dmot_indx);
  
  dmot1=doConv(dmot_rgr,MION_k)';     % apply MION function to raw data
  
  R = addRegressor(R,dmot1,'Mot_Bartel_Total(MION)');
  
  clear dmot*

  %%%%%%%%%%%
  %% Aggression
  % Add the agression regressors together
  aggr_indx = find(strncmp('Agression',r.names,9));

  aggr_rgr  = sum(sqrt(r.rgrs(:,aggr_indx)),2); % compress with a square root function
  tmp=imag(conj(aggr_rgr));                     % and make sure we have no imaginary numbers
  aggr_rgr2 = aggr_rgr;
  aggr_rgr2(tmp~=0)=tmp(tmp~=0);                % fill the now negative numbers into the vector
  aggr_rgr2 = doConv(aggr_rgr2,MION_k)';             % apply MION function to raw data
  R = addRegressor(R,aggr_rgr2,'Aggression(sqrt_MION)');
  
  clear aggr* tmp
  
  %%%%%%%%%%%
  %% Affliative 
  % Add the agression regressors together
  affil_indx = find(strncmp('Play',r.names,4) | strncmp('Groom',r.names,5) | strncmp('Copulations',r.names,11));

  affil_rgr  = sum(sqrt(r.rgrs(:,affil_indx)),2); % compress with a square root function
  tmp=imag(conj(affil_rgr));                     % and make sure we have no imaginary numbers
  affil_rgr2 = affil_rgr;
  affil_rgr2(tmp~=0)=tmp(tmp~=0);                % fill the now negative numbers into the vector
  affil_rgr2 = doConv(affil_rgr2,MION_k)';             % apply MION function to raw data
  R = addRegressor(R,affil_rgr2,'Affiliative(sqrt_MION)');

  clear affil* tmp
  %%%%%%%%%%%
  %% Food
  % Add the agression regressors together
  food_indx = find(strncmp('Food',r.names,4));

  food_rgr  = sum(sqrt(r.rgrs(:,food_indx)),2); % compress with a square root function
  tmp=imag(conj(food_rgr));                     % and make sure we have no imaginary numbers
  food_rgr2 = food_rgr;
  food_rgr2(tmp~=0)=tmp(tmp~=0);                % fill the now negative numbers into the vector
  food_rgr2 = doConv(food_rgr2,MION_k)';             % apply MION function to raw data
  R = addRegressor(R,food_rgr2,'Food(sqrt_MION)');

  clear food* tmp
  %%%%%%%%%%%
  %% Aggression & Play
  % Add the agression regressors together
  aggr_indx = find(strncmp('Agression',r.names,9) | strncmp('Play',r.names,4));

  aggr_rgr  = sum(sqrt(r.rgrs(:,aggr_indx)),2); % compress with a square root function
  tmp=imag(conj(aggr_rgr));                     % and make sure we have no imaginary numbers
  aggr_rgr2 = aggr_rgr;
  aggr_rgr2(tmp~=0)=tmp(tmp~=0);                % fill the now negative numbers into the vector
  aggr_rgr2 = doConv(aggr_rgr2,MION_k)';             % apply MION function to raw data
  R = addRegressor(R,aggr_rgr2,'Aggr_Play(sqrt_MION)');
  
  clear aggr* tmp
  
  %% Dyadic Interactions
  % Add the agression regressors together
  dyad_indx = find(strncmp('Dyadic',r.names,6));
  
  dyad_rgr=sum(r.rgrs(:,dyad_indx),2);
  tmp=imag(conj(dyad_rgr));                     % and make sure we have no imaginary numbers
  dyad_rgr2 = dyad_rgr;
  dyad_rgr2(tmp~=0)=tmp(tmp~=0);                % fill the now negative numbers into the vector
  dyad_rgr2 = doConv(dyad_rgr2,MION_k)';             % apply MION function to raw data
  R = addRegressor(R,dyad_rgr2,'Dyadic(sqrt_MION)');
  
  %% Number of Saccades 
  % Add the agression regressors together
  sac_indx = find(strncmp('NumSacs',r.names,7));
  sac_rgr = r.rgrs(:,sac_indx);
  
  sac_rgr2=sqrt(sac_rgr);
  tmp=imag(conj(sac_rgr2));                     % and make sure we have no imaginary numbers
  sac_rgr2(tmp~=0)=tmp(tmp~=0);                % fill the now negative numbers into the vector
  sac_rgr2 = doConv(sac_rgr2,MION_k)';             % apply MION function to raw data
  R = addRegressor(R,sac_rgr2,'Saccades(sqrt_MION)');
  
  ALLRGR = R;
  
function R = addRegressor(R,rgr,name)
  
  if isempty(R)
	R.names = {};
	R.rgrs = [];
  end
  nrgr = length(R.names);
  R.names{nrgr+1} = name;
  R.rgrs(:,nrgr+1) = rgr;
  