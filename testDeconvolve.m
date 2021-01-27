
% Load data
dirData = '/Volumes/USRlab/data/parksh/'; % '/einstein0/USRlab/data/parksh/'; %  '/Volumes/USRlab/data/parksh/'; %/rmov/';
nameSubjNeural = 'Tor';
nameSubjBOLD = 'Art';
dirDataNeural = [dirData, nameSubjNeural, '/'];
dirDataBOLD = [dirData, nameSubjBOLD, '/'];

loadFileName = sprintf('dataNeuralBOLD_%s%s_indMov.mat', nameSubjNeural, nameSubjBOLD);
load(fullfile(dirData, 'dataNeuralBOLD_TorArt_indMov.mat'))

%%
% Get the movie list from single unit data
d_n = dir(fullfile(dirDataNeural, '*sig*.mat'));
[listMatSUFile{1:length(d_n), 1}] = deal(d_n.name); % Get the list of single unit files 
listMovSU = unique(regexp(cat(2,d_n.name), '\d*(?=sig)', 'match')); % Movie IDs of the single units in cellstr


% Define movie set
setMovIDs = [1 2 3 10 11 12 13 14 15];

% Get the list of single units & their movie
listSU_all = regexp(cat(2,d_n.name), '(?<=sig)\w*', 'match')';
listSUchannelID = unique(listSU_all);

setCellIDs = listSUchannelID;

% get discrete-time firing rate
sizeTimeBin_sec = dataBOLD.TR; % Size of time bin in firing rate
FR = createCellRegressor_indMov_discreteTime(dirDataNeural, setCellIDs, setMovIDs, sizeTimeBin_sec);
% --FR(iChan, iMov): Struct for each cell (FR.cellID) & each movie (FR.movID)
%       1) matFR{1}: firing rate for all trials (time x trial)
%       2) mnFR: mean firing rate across trials (time x 1), i.e. average in
%       second dimension of matFR


% 2. linear algebra
selChan = [1 5 6 11 18 21 28 34 42 47 51];
for iC=1:length(selChan)
    
    concatMatStim = [];
    concatMatMION = [];
    
    iChan = selChan(iC);
    
    for iMov=1:9
        if ~isempty(FR(iChan, iMov).mnFR)
            fprintf('Channel: %s, Movie: %d\n', FR(iChan,iMov).cellID, FR(iChan,iMov).movID)
            
            %     % using spike density function
            %     sdf_test  = S(iChan,iMov).mnsdf;
            %     neuralRgrs = resample(sdf_test, 0.001*1000, dataBOLD.TR*1000); % sdf as an event sequence
            
            % using dicrete-time firing rate
            neuralRgrs = FR(iChan,iMov).mnFR; %>mean(FR(iChan,iMov).mnFR);
            
            MION_test = dataBOLD.catmvoltc{iMov};
            locValid = squeeze(~isnan(MION_test(1,1,1,:))); % exclude NaNs
            
            validMION = MION_test(:,:,:,locValid);
            [nx ny nz nt] = size(validMION);
            matMION = reshape(validMION, nx*ny*nz, nt)';
            
            nhrf = 10; % length of the hrf in TR
            xAxis = [1:nhrf].*dataBOLD.TR;
            
            % make a design (event) matrix
            matStim = [];
            matStim = zeros(nt, nhrf);
            tempSeq = neuralRgrs(locValid);
            for iC = 1:nhrf
                matStim(:,iC) = tempSeq;
                tempSeq = [0;tempSeq(1:end-1)];
            end            
            
            PmatStim = pinv(matStim); %inv(matStim'*matStim)*matStim';
            hest = PmatStim*matMION;
            
            figure(100);
            clf
            subplot(1,2,1)
            imagesc(matStim)
            axis equal
            axis off
            colormap(gray)

            subplot(1,2,2)
            plot(xAxis, hest(:,51935:51955), 'o-')
            xlim([xAxis(1) xAxis(end)])
            xlabel('Time (s)')
            title(sprintf('Channel %s, Movie %d', S(iChan, iMov).cellID, S(iChan, iMov).movID))
%             input('')
            
            indMovDeconv(iChan,iMov).matStim = matStim;
            indMovDeconv(iChan,iMov).matMION = matMION;
            indMovDeconv(iChan,iMov).hest = hest;
            
            
        end
    end
end

iChan=1;
iMov=1;
amp = max(indMovDeconv(iChan,iMov).hest);
mapAmp = reshape(amp, [nx, ny, nz]);
DSP.proc.scalarmap_3d = mapAmp;



% selMovieIndex = 4:9;
for iC=1:length(selChan)
    iChan = selChan(iC);
    
    concatMatStim=[]; concatMatMION=[]; PConcatMatStim=[]; concatHEst=[];
    
    concatMatStim = cat(1, indMovDeconv(iChan,:).matStim);
    concatMatMION = cat(1, indMovDeconv(iChan,:).matMION);
    
    PConcatMatStim = pinv(concatMatStim);
    concatHEst = PConcatMatStim*concatMatMION;
    
    concatDeconv(iC).cellID = FR(iChan,iMov).cellID;
    concatDeconv(iC).matStim = concatMatStim;
    concatDeconv(iC).matMION = concatMatMION;
    concatDeconv(iC).hest = concatHEst;
    
    figure; clf;
    plot(xAxis, concatHEst(:,51935:51955), 'o-')
    xlim([xAxis(1) xAxis(end)])
    xlabel('Time (s)')
    title(sprintf('Channel %s, 9-movie Concatenated', S(iChan, iMov).cellID))
       
    
%     input('')
    
end

% iChan=1;
% iMov=1;


amp = max(concatDeconv(iC).hest); %max(indMovDeconv(iChan,iMov).hest);
mapAmp = reshape(amp, [nx, ny, nz]);
DSP.proc.scalarmap_3d = mapAmp;


% Use JG's method to get estimated hdr and r2 for each voxel


selMovInd = 4:9;
for iC=1:11 % among 11 selected cells

    matStim=[]; matMION=[];
    
    iChan = selChan(iC);
    matStim = cat(1, indMovDeconv(iChan, selMovInd).matStim); % concatDeconv(iC).matStim
    matMION = cat(1, indMovDeconv(iChan, selMovInd).matMION); % concatDeconv(iC).matMION
    
    % % replicate and tile a single movie data
    % matStim = repmat(indMovDeconv(iChan, iMov).matStim, 9, 1); % concatDeconv(iC).matStim
    % matMION = repmat(indMovDeconv(iChan, iMov).matMION, 9, 1); % concatDeconv(iC).matMION
    
    d.scm = matStim;
    d.dim = [nx ny nz];
    d.nhdr = 1;
    d.hdrlen = 10;
    d.volumes = 1:size(matMION,1);
    d.data = reshape(matMION', [nx, ny, nz, size(matMION,1)]);
    
    
    verbose = 0;
    d = getr2(d, verbose);
    fprintf('cell #%d/%d\n', iC, length(selChan))
    Decon(iC) = d;
end

% init some variables
ehdr=[];r2 = [];

% precalculate the normal equation (this dramatically speeds up things)
covarianceMatrix = d.scm'*d.scm;
covarianceMatrixRank = rank(covarianceMatrix);
% use normal equations if we have a full ranked covariance matrix
if covarianceMatrixRank == size(covarianceMatrix,1)
  inverseCovarianceMatrix = covarianceMatrix^-1;
  precalcmatrix = inverseCovarianceMatrix*d.scm';
  diagOfInverseCovariance = diag(inverseCovarianceMatrix);
% otherwise use pinv
else
  % note that if we need to use the pseudo inverse it means that there is ambiguity in the design
  % such that there are an infinite number of possible solutions. The psuedo-inverse solution
  % chosses the solution with the minimum length (i.e. Euclidian norm)
  if verbose,disp(sprintf('(getr2) Design covariance matrix (%ix%i) is rank %i. Using pseudo-inverse to invert.',size(covarianceMatrix,1),size(covarianceMatrix,2),covarianceMatrixRank));end
  precalcmatrix = pinv(d.scm);
  % get the diagonal of the inverse of the design covariance matrix (used for estimating standard errors)
  diagOfInverseCovariance = diag(pinv(d.scm'*d.scm));
end

% check roi
slices = 1:d.dim(3);slicen = length(slices);
xvals = 1:d.dim(1);xvaln = length(xvals);
yvals = 1:d.dim(2);yvaln = length(yvals);
  
% preallocate memory
d.ehdr = zeros(d.dim(1),d.dim(2),d.dim(3),d.nhdr,d.hdrlen);
d.ehdrste = zeros(d.dim(1),d.dim(2),d.dim(3),d.nhdr,d.hdrlen);
d.r2 = zeros(d.dim(1),d.dim(2),d.dim(3));

% turn off warnings to avoid divide by zero warning
warning('off','MATLAB:divideByZero');

% display string
if verbose,disppercent(-inf,'(getr2) Calculating r2');end
% cycle through images calculating the estimated hdr and r^2s of the 
% estimate.
%
% this following section has been optimized to run faster by
% eliminating one of the loops. Various different methods were
% tested eliminating all the loops and doing one big calculation
% which thrashed memory too much, or eliminating different
% dimensions and it was found that eliminating the first dimension
% was by far the faster by a factor of about 2-3. 
onesmatrix = ones(length(d.volumes),1);
for j = yvals
  for k = slices
    % get the time series we are working on
    % this includes all the rows of one column from one slice
    % and all data points for each of these
    % thus the time series is a nxm matrix where each of the m columns
    % contains the n time points recording for that voxel
    timeseries = squeeze(d.data(:,j,k,d.volumes))';
    % subtract off column means
    colmeans = mean(timeseries,1);
    timeseries = timeseries - onesmatrix*colmeans;
    % convert to percent signal change
    timeseries = 100*timeseries./(onesmatrix*colmeans);
    % get hdr for the each voxel
    ehdr{j,k} = precalcmatrix*timeseries;
    % calculate error bars, first get sum-of-squares of residual
    % (in percent signal change)
    sumOfSquaresResidual = sum((timeseries-d.scm*ehdr{j,k}).^2);
    % now calculate the sum-of-squares of that error
    % and divide by the degrees of freedom (n-k where n
    % is the number of timepoints in the scan and k is 
    % the number of timepoints in all the estimated hdr)
    S2 = sumOfSquaresResidual/(length(d.volumes)-size(d.scm,2));
    % now distribute that error to each one of the points
    % in the hemodynamic response according to the inverse
    % of the covariance of the stimulus convolution matrix.
    ehdrste{j,k} = sqrt(diagOfInverseCovariance*S2);
    % calculate variance accounted for by the estimated hdr
    r2{j,k} = (1-sumOfSquaresResidual./sum(timeseries.^2));
  end
  if verbose,disppercent(max((j-min(yvals))/yvaln,0.1));end
end
if verbose,disppercent(inf);end

% reshape matrix. this also seems the fastest way to do things. we
% could have made a matrix in the above code and then reshaped here
% but the reallocs needed to continually add space to the matrix
% seems to be slower than the loops needed here to reconstruct
% the matrix from the {} arrays.
if verbose,disppercent(-inf,'(getr2) Reshaping matrices');end
for i = xvals
  for j = yvals
    for k = slices
      % get the ehdr
      d.ehdr(i,j,k,1:d.nhdr,:) = reshape(squeeze(ehdr{j,k}(:,i)),[d.hdrlen d.nhdr])';
      % and the stderror of that
      d.ehdrste(i,j,k,1:d.nhdr,:) = reshape(squeeze(ehdrste{j,k}(:,i)),[d.hdrlen d.nhdr])';
      % now reshape r2 into a matrix
      d.r2(i,j,k) = r2{j,k}(i);
    end
  end
  if verbose,disppercent((i-min(xvals))/xvaln);end
end

% display time took
if verbose,disppercent(inf);end

warning('on','MATLAB:divideByZero');



% % 1. use matlab deconv function
% sdf_test  = S(1,1).mnsdf;
% MION_test = dataBOLD.catmvoltc{1};
% TR = dataBOLD.TR;
% Fs = 1000;
% sample_vox = squeeze(MION_test(31,19,21,:));
% 
% sdf_ts = [1:length(sdf_test)]/Fs;
% MION_ts = [1:length(MION_test)]*TR;
% figure(1000); clf
% subplot(4,1,1);
% plot(sdf_ts,sdf_test);
% subplot(4,1,2);
% plot(MION_ts,sample_vox,'r')
% 
% % get rid of nans
% % resample_vox = interp(sample_vox,Fs*TR);
% % val = ~isnan(resample_vox);
% 
% % sdf_testval = sdf_test(val);
% % sdf_test_short = sdf_testval(15000:end-15000);
% % 
% % [K,R] = deconv(resample_vox(val),sdf_test_short);
% 
% 
% resample_sdf = resample(sdf_test, 1, 2400);
% val = ~isnan(sample_vox);
% 
% [K,R] = deconv(sample_vox(val),resample_sdf(val));
% 
% subplot(4,1,3);
% plot(resample_vox(val),'k');
% hold on
% plot(sdf_test(val),'r')
% 
% subplot(4,1,4);
% plot(K);





% X = [];
% for j=1:4
%     Xj = zeros(m,n);
%     temp = s==j;
% 
%     for i=1:n
%         Xj(:,i) = temp;
%         temp = [0;temp(1:end-1)];
%     end
%     X = [X,Xj];
% end
% 
% figure(1)
% clf
% imagesc(X)
% axis equal
% axis off
% colormap(gray)
% 
% r = X*h(:);
% 
% %add noise
% noiseSD = .25;
% rNoise = r+noiseSD*randn(size(r));
% 
% hest = pinv(X)*rNoise;