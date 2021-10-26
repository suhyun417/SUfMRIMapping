
% parameters
setMovID = [1 2 3];


%% Create movie regressors
flagSM = 1; % flag for compression and smoothing

fullRGR4fps = createMovieRGR_4fps_indMov(setMovID, flagSM); %createFullMovieRegressors_4fps_indMov(setMovID); %
indValidRGR = [1, 3, 9, 20, 21, 22, 25]; 
% 1: 'Luminance', 3: 'Speed', 9: 'Beta Contrast', 20: 'Faces', 21: 'One face', 22: 'Body parts', 25: 'Any animal'

% Face scale regressor (in TR unit)
ttt=load('/Volumes/PROCDATA/parksh/MovieRegressors/dbtmMriReg.mat'); %load('/procdata/parksh/MovieRegressors/dbtmMriReg.mat');
scaleRGR = ttt.reg.xx(7,:)';
% [Rscale]=corr(resample(pcaMovCat.coeff(:,1), 0.1*10, 2.4*10), scaleRGR', 'rows', 'complete');



%% Apply regressors to PCs
iCase = 1; % 2; % 1 for discrete FR, 2 for Gaussian SDF
catRGR=[];
for iMov=1:length(setMovID)
    m = setMovID(iMov);
    matCurRGR = fullRGR4fps(iMov).regressors; %fullRGR4fps(m).smoRegressors(:,indValidRGR); %fullRGR4fps(iMov).regressors(:,indValidRGR); %fullRGR4fps(iMov).regressors;
    catRGR = cat(1, catRGR, matCurRGR); % concatenation across movies
end
matRGR = resample(catRGR, 0.25*100, 2.4*100);
matRGR = cat(2, matRGR, scaleRGR);

switch iCase
    case 1
        
        % load PC
        load('/procdata/parksh/Tor/eigen/pcaCatMovie123_FR_tor.mat')
        % setMovID = [1 2 3];
        
        
        setPC = 1:5;
        matPC = [];
        matPC = pcaMovCat.coeff(:,setPC); %resample(pcaMovCat.coeff(:,setPC), 4, 10);
        
    case 2
        
        % load PC
        load('/procdata/parksh/Tor/eigen/pcaCatMovie123_SDF_tor.mat')
%         
%         matRGR = resample(catRGR, 0.25*100, 2.4*100);
%         matRGR = cat(2, matRGR, scaleRGR);
        
        setPC = 1:5;
        matPC = [];
        matPC = resample(pcaMovCat.coeff(:,setPC), 0.1*10, 2.4*10);
        
end

[Rpca_validRGR] = corr(matRGR, matPC, 'rows', 'complete');


% figure
varnames_validRGR = cat(1, fullRGR4fps(1).features(indValidRGR), {'Face size'});
varnames_fullRGR = cat(1, fullRGR4fps(1).features, {'Face size'}); %cat(1, fullRGR4fps(1).features(indValidRGR), {'Face size'});

figPCARGR_catMov = figure;
set(figPCARGR_catMov, 'Color', 'w', 'PaperPositionMode', 'auto')
imagesc(Rpca_validRGR);
set(gca, 'YTick', 1:size(matRGR,2), 'YTickLabel', varnames_validRGR);
set(gca, 'FontSize', 15)
% set(gca, 'YTick', 1:length(indValidRGR), 'YTickLabel', fullRGR4fps(1).features(indValidRGR));
xlabel('PC #')
ylabel('Features')

% make blue-white-red colorbar
cval = 0.5;
cmin = -cval; cmax = cval;
colornum = 256;
colorInput = [1 0 0; 1 1 1; 0 0 1];
oldSteps = linspace(-1, 1, length(colorInput));
newSteps = linspace(-1, 1, colornum);
for j=1:3 % RGB
    newmap_all(:,j) = min(max(transpose(interp1(oldSteps, colorInput(:,j), newSteps)), 0), 1); 
end

endPoint = round((cmax-cmin)/2/abs(cmin)*colornum);
newmap = squeeze(newmap_all(1:endPoint, :));

figure(figPCARGR_catMov)
set(gca, 'CLim', [cmin cmax])
colormap(flipud(newmap))
set(gca, 'TickDir', 'out')
box off
c=colorbar;

% for poster
set(gca, 'XTick', [], 'YTick', [])
set(figPCARGR_catMov, 'Position', [100 100 360 525])
xlabel(''); ylabel('');
box on

print(gcf, fullfile(dirFig, 'pca_catMovie_FR_movieRGR'), '-dtiff', '-r600')

% figBarMovieRGR = figure;
% set(figBarMovieRGR, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1300 450 560 420])
% for i=1:5
%     
%     figure(figBarMovieRGR); clf;
%     bar(Rpca_validRGR(:,i))
%     ylim([-.5 .5])
%     box off    
%     set(gca, 'LineWidth', 2, 'FontSize', 15)
%     set(gca, 'XTickLabel', [])
%     
%     fName = sprintf(fullfile(dirFig, 'pca_catMovie_movieRGR_PC%d'), i);
%     print(figBarMovieRGR, fName, '-depsc')
%     
% end


% %% Apply movie regressors to single unit data
% % load the spikes 
% if ~ exist('S') || isempty(S)
%     load('/procdata/parksh/Tor/Tor_movieTS_SU_indMov.mat')
% %     load(fullfile(dataDir, saveFileName), 'S')
% end
% 
% % Make a response matrix
% unimov = cat(1, S(1,:).movID);
% mnSDF_org = struct([]);
% for iMov = 1:length(setMovID)
%     tempInd = find(unimov==setMovID(iMov));
%     matSDF = cat(2, S(:,tempInd).mnsdf); % get spikes of all the cells for that movie
%     mnSDF_org(iMov).matSDF = resample(matSDF, 4, 1000); % matSDF;
%     mnSDF_org(iMov).cellIDs = cat(1,S(:,tempInd).cellID);
% end
% 
% % Compute correlations
% catRGR = []; catResp =[];
% % figCellRGR_indMov = figure;
% % set(figCellRGR_indMov, 'Color', 'w', 'PaperPositionMode', 'auto')
% for iMov=1:3
% 
% matCurRGR = fullRGR4fps(iMov).smoRegressors(:,indValidRGR); %fullRGR4fps(iMov).regressors(:,indValidRGR); %sqrt_RGR; %fullRGR4fps(iMov).regressors;
% matResp = mnSDF_org(iMov).matSDF;
% 
% catRGR = cat(1, catRGR, matCurRGR); % concatenation across movies
% catResp = cat(1, catResp, matResp);
%     
% [R, P] = corr(matCurRGR, matResp); % all cells for one movie
% 
% % figure(figCellRGR_indMov);
% % subplot(1, length(setMovID), iMov);
% % imagesc(R);
% % set(gca, 'CLim', [-1 1])
% % set(gca, 'YTick', 1:length(indValidRGR), 'YTickLabel', fullRGR4fps(1).features(indValidRGR));
% % xlabel('Cell ID')
% % ylabel('Features')
% % title(sprintf('Movie #%d', fullRGR4fps(iMov).movieID))
% 
% end
% 
% [Rcat, Pcat] = corr(catRGR, catResp); % all cells for one movie
% 
% figCellRGR_catMov = figure;
% set(figCellRGR_catMov, 'Color', 'w', 'PaperPositionMode', 'auto')
% imagesc(Rcat')
% set(gca, 'YTick', 1:length(mnSDF_org(1).cellIDs), 'YTickLabel', mnSDF_org(1).cellIDs)
% set(gca, 'XTick', 1:length(indValidRGR), 'XTickLabel', fullRGR4fps(1).features(indValidRGR));
% set(gca, 'CLim', [-1 1].*.5)
% ylabel('Cell ID')
% xlabel('Features')
% title('Concatenated Movies')
% 
% % imagesc(Rcat);
% % set(gca, 'YTick', 1:length(indValidRGR), 'YTickLabel', fullRGR4fps(1).features(indValidRGR));
% % xlabel('Cell ID')
% % ylabel('Features')
% % title('Concatenated Movies')

    






% for iMov = 1:3
%     figure(100);
%        set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
%     subplot(1, length(setMovID), iMov); cla;
%     
%  iPC = 1:3;
%    matCurRGR = cat(1,fullRGR4fps(iMov).regressors);
%    catPC = cat(1, pcaMov(iMov).coeff);
%    matResp = resample(catPC(:,iPC), 4, 10);
%    
%    [Rpca] = corr(matCurRGR, matResp);
%    
% %    figure
%    imagesc(Rpca)
%    set(gca, 'CLim', [-1 1].*.5)
% %    set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
%    if iMov==1
%    set(gca, 'YTick', 1:length(fullRGR4fps(1).features), 'YTickLabel', fullRGR4fps(1).features);
%    end
%    
% end
% end


        
%     figure;
%     hold on;
%     %plot(rawData, 'k-');
%     %plot(smoData, 'r-');
%     plot(sqrt_raw, 'b-');
%     plot(sqrt_smoData, 'g-');
%     plot(log_raw, 'c-');
%     plot(log_smoData, 'm-');
%    
%     hold off;
%     
%     zScore_raw = zscore(smoData);
%     zScore_sqrt = zscore(sqrt_smoData);
%     zScore_log = zscore(log_smoData);
%     figure;
%     hold on;
%     plot(zScore_raw, 'k-');
%     plot(zScore_log, 'm-');
%     plot(zScore_sqrt, 'b-');
%     hold off;
%     
%     
%      
%     
      

    
    
% end





% curRGR = meanBOLD_ROI;
% % Compute correlation maps
% [nx, ny, nz, nt] = size(matBOLD);
% nVox = nx*ny*nz;
% 
% [Rvals, Pvals] = corr(reshape(matBOLD, nVox, nt)', curRGR,...
%     'rows','complete');
% 
% mapR = reshape(Rvals, [nx, ny, nz]); % Don't need to multiply -1 because now it's correlation between MION signals %.*(-1);
% DSP.proc.scalarmap_3d = mapR;



% function S_correlate_ALLR_e66()
%   
%   global S
%   %
%   % examination of regressors combined across days
%   %
%   
%   global STDPATH DSP DATA GH
%   
%   
%   % first start,open the XL file, and reset the windows
%   %
%   blockpaths_ber;
%   STDPATH.xlsfile = 'block_timecourse_files_ber.xls';
%   loadXLSFile;
% 
%   resetBlockanaWindows;
% 
% 
%   % save as data file in session directory 
%   %
%   saveflag = 0;  loadS = 1; %10/16/14
%   S_filename = 'e66_S_rgr_corr_lots';
%   
%   redo = 0;
%   skip=7; % number of TRs to skip  %% ADDED 10/04/12 BER
% 
%   textfileroot = 'FL_e66_allfiles2';
% %   textfileroot = 'e66_2file_allmovies';
% 
%   MCD  = loadMeanConcatData(textfileroot,skip,redo); % load or
% 
%   % 
%   % here read in some pre-computed regressors and read some from the xls
%   % file (all in raw form)
%   %
%   RAWRGR = getAllRawRegressors(MCD.unimov,MCD.TR,MCD.Fs,MCD.max_tr);
% 
%   %
%   % Here we replace the initial list of regressors 
%   % with some that we tailor (e.g. combining, compressing, etc)
%   %
% %   ALLRGR = SU_addSpecificRegressors_Motion(RAWRGR,MCD.TR); 
% %   ALLRGR = SU_addSpecificRegressors_Motion2(RAWRGR,MCD.TR); 
% %   ALLRGR = SU_addSpecificRegressors_Faces(RAWRGR,MCD.TR); 
% %   ALLRGR = SU_addSpecificRegressors_LowVision(RAWRGR,MCD.TR);
% %   ALLRGR = SU_addSpecificRegressors_Bodies(RAWRGR,MCD.TR);
% %   ALLRGR = SU_addSpecificRegressors_Behav(RAWRGR,MCD.TR);
%   ALLRGR = SU_addSpecificRegressors_Compare(RAWRGR,MCD.TR);
% %   ALLRGR = SU_addSpecificRegressors_ShortCompare(RAWRGR,MCD.TR);
%   
%   %
%   % compute the correlation maps (single variable regression)
%   %
%   %
%   %  ORDER OF FEATURES
%     %  1   =  Speed(log_MION)
%     %  2   =  Mot_Diverg_Rect(log_MION)
%     %  3   =  Mot_Diverg_STD(log_MION)
%     %  4   =  Mot_Diverg_Beta(log_MION)
%     %  5   =  Mot_Diverg_Gamma(MION_sqrt_sm2)
%     %  6   =  Scene_Cuts(log_MION_sm2)
%     %  7   =  Faces_Full(sqrt_MION)
%     %  8   =  Faces_SV(sqrt_MION)
%     %  9   =  1_Face_Full(sqrt_MION)
%     %  10  =  1_Face_SV(sqrt_MION)
%     %  11  =  Heads(sqrt_MION)
%     %  12  =  Faces_all(sqrt_MION)
%     %  13  =  1_Face_comb(sqrt_MION)
%     %  14  =  Luminance(log_MION)
%     %  15  =  Contrast_STD(log_MION)
%     %  16  =  Contrast_Beta(log_MION)
%     %  17  =  Contrast_Gamma(log_MION)
%     %  18  =  SF_low(log_MION)
%     %  19  =  SF_high(log_MION)
%     %  20  =  SF_ratio(log_MION)
%     %  21  =  Butts(MION_sqrt)
%     %  22  =  Bodies(sqrt_MION)
%     %  23  =  Extremities(sqrt_MION)
%     %  24  =  Hands(sqrt_MION)
%     %  25  =  Animals(MION_sqrt)
%     %  26  =  Macaque(MION_sqrt)
%     %  27  =  Rhesus(MION)
%     %  28  =  Human(MION)
%     %  29  =  Mot_Local_STD(MION)
%     %  30  =  Mot_Bartel_Local(MION)
%     %  31  =  Mot_Bartel_Global(MION)
%     %  32  =  Mot_Bartel_Resid(MION)
%     %  33  =  Mot_Bartel_Total(MION)
%     %  34  =  Aggression(sqrt_MION)
%     %  35  =  Affiliative(sqrt_MION)
%     %  36  =  Food(sqrt_MION)
%     %  37  =  Aggr_Play(sqrt_MION)
%     
%     
%     function RAWRGR = getAllRawRegressors(unimov,TR,Fs,max_tr)
%   codingpath = '/procdata/russbe/CodingXLS';
%   lowlevpath = '/procdata/russbe/LowlevelMovRegressors/FPS10';
%   eyempath   = '/procdata/russbe/EyeMovements';
%   
%   codingfile = sprintf('CodingSheetmov1_inv');
%   xlsfullpath = sprintf('%s/%s.xls',codingpath,codingfile);
%   [num,str] = xlsread(xlsfullpath);
%   varnames = str(1,:);
%  
%   lowlevfile = sprintf('Movie1_10fps_rgr'); %sprintf('Movie1_rgr'); 
%   lowlevfullpath = sprintf('%s/%s.mat',lowlevpath,lowlevfile);
%   load(lowlevfullpath);  % loads RGR10fps variable
%   
%   eyefile = sprintf('e66SacRGR');
%   eyefullpath = sprintf('%s/%s.mat',eyempath,eyefile);
%   load(eyefullpath);  % loads e66 RGR variable
%   
%   varnames = [varnames RGR10fps.features]; %[varnames RGR.features];
%   varnames = [varnames {'NumSacs'}];
%   
%   if nargin < 2  %% ADDED skipTR variable 10/04/12 BER
%       skipTRs=0; % don't skip any TRs if no skip is sent
%   end
%   
%   cregressors = [];
%   rgrs = [];
%   for m=1:length(unimov)
% 	codingfile   = sprintf('CodingSheetmov%d_inv',unimov(m));
% 	lowlevfile = sprintf('Movie%d_10fps_rgr',unimov(m)); 
% 	
% 	xlsfullpath = sprintf('%s/%s.xls',codingpath,codingfile);
% 	[num,str] = xlsread(xlsfullpath);
% 	lowlevfullpath = sprintf('%s/%s.mat',lowlevpath,lowlevfile);
% 	load(lowlevfullpath);  % loads RGR variable
% 	
% 	
% 	nvars = length(str(1,:));
% 	for i=1:nvars
% 	  tmp   = num(:,i);
% 	  rtmp  = resample(tmp,100*Fs,100*TR);
% 	  rtmp(max_tr+1:end) = [];
% 	  rgrs(:,i) = rtmp;
% 	end
% 	
% 	nrgr = size(RGR.regressors,1);
% 	for i=1:nrgr
% 	  rgrs(:,nvars+i)=RGR.regressors(i,:);
%     end
% 	
%     rgrs(:,nvars+i+1) = SacRGR(unimov(m)).rgr';
% % 	tvals = [(-20*TR):TR:(20*TR)];
% % 	% This is just made-up
% % 	MION_k = gampdf(tvals,TR,2*TR);
% % 	MION_k = MION_k./sum(MION_k);
% % 	crgrs = doConv(rgrs,MION_k);
% % 	cregressors = [cregressors crgrs];  
%     cregressors = [cregressors rgrs'];
%   end
%   
%   RAWRGR.names = varnames;
%   RAWRGR.rgrs  = cregressors';