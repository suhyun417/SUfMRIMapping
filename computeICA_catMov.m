% computeICA_catMov

%%%%%%%%%%%%%%%%%%%%%%%%% NEED TO MODIFY TO USE 1kHz resolution spikes,
%%%%%%%%%%%%%%%%%%%%%%%%% instead of using TR-resolution SDF

% ICA toolbox
addpath('/projects/parksh/_toolbox/RobustICA')

dirData = '/procdata/parksh/';
dirDataNeural = '/procdata/parksh/Tor/';
dirDataEig = '/procdata/parksh/Tor/eigen/';

setMovID = [1 2 3]; %[1 2 3 10 11 12 13 14 15];
validC = [1:8, 11:34, 36:51]; % valid channel with movie [1 2 3]
% validC = [1:8, 11, 13, 15:18, 20:21, 23, 25:28, 30:34, 36:43, 46:49, 51]; % valid channel with all 9 movie data

% load the mean FR
load(fullfile(dirDataNeural,'Tor_movieTS_SU_indMov.mat'));
matSDF=[];
for iMov = 1:length(setMovID)
    movID = setMovID(iMov);
    tempMat=[];
    tempMat = cat(2,S(validC,movID).mnFR); %cat(2,S(validC,movID).mnsdf);
    matSDF = [matSDF ; tempMat]; % concatenate across movies
end
cellData = S;
% clear S;

% Robust ICA
arguments = {};
% 2) Regression-based deflation:
%
% arguments = {'deftype', 'regression'};
% 3) Aim at two sub-Gaussian sources first, running a fixed number of
% iterations per independent component:
%
% arguments = {'deftype', 'regression', 'kurtsign', [-1, -1, zeros(1, K-2)], 'tol', -1, 'maxiter', 10}
% 4) No prewhitening, regression-based deflation, dimensionality reduction
%
% arguments = {'deftype', 'regression', 'prewhi', false, 'dimred', true}
arguments = {arguments{:}, 'verbose', true}; % execute in verbose mode

X = matSDF';
[SourceMat, H, iter, W] = robustica(X, arguments);

figure;
for isub=1:8
    subplot(8,1,isub)
    plot(SourceMat(isub,:))
    axis tight
end


% z-score normalization
matSDF_norm = matSDF-repmat(mean(matSDF), size(matSDF, 1), 1);
matSDF_norm = matSDF_norm./repmat(std(matSDF), size(matSDF,1), 1);

X_norm = matSDF_norm';
[SourceMat_norm, H_norm, iter_norm, W_norm] = robustica(X_norm, arguments);

figure;
for isub=1:8
subplot(8,1,isub)
plot(SourceMat_norm(isub,:))
axis tight
end

%% Compute correlation with movie regressors
setMovID = [1 2 3];
% Create movie regressors
flagSM = 1; % flag for compression and smoothing
fullRGR4fps = createMovieRGR_4fps_indMov(setMovID, flagSM); %createFullMovieRegressors_4fps_indMov(setMovID); %
indValidRGR = [1, 3, 9, 20, 21, 22, 25];
% 1: 'Luminance', 3: 'Speed', 9: 'Beta Contrast', 20: 'Faces', 21: 'One face', 22: 'Body parts', 25: 'Any animal'
% Face scale regressor (in TR unit)
ttt=load('/procdata/parksh/MovieRegressors/dbtmMriReg.mat');
scaleRGR = ttt.reg.xx(7,:)';
% [Rscale]=corr(resample(pcaMovCat.coeff(:,1), 0.1*10, 2.4*10), scaleRGR', 'rows', 'complete');
%% Apply regressors to PCs
catRGR=[];
for iMov=1:length(setMovID)
m = setMovID(iMov);
matCurRGR = fullRGR4fps(m).smoRegressors(:,indValidRGR); %fullRGR4fps(iMov).regressors(:,indValidRGR); %fullRGR4fps(iMov).regressors;
catRGR = cat(1, catRGR, matCurRGR); % concatenation across movies
end
matRGR = resample(catRGR, 0.25*100, 2.4*100);
matRGR = cat(2, matRGR, scaleRGR);


matIC = SourceMat';
[Rica_validRGR] = corr(matRGR, matIC, 'rows', 'complete');

max(max(Rica_validRGR))
min(min(Rica_validRGR))

% Plot the correlation
varnames = cat(1, fullRGR4fps(1).features(indValidRGR), {'Face size'});
figPCARGR_catMov = figure;
set(figPCARGR_catMov, 'Color', 'w', 'PaperPositionMode', 'auto')
imagesc(Rica_validRGR)
set(gca, 'YTick', 1:size(matRGR,2), 'YTickLabel', varnames);
set(gca, 'FontSize', 15)
% set(gca, 'YTick', 1:length(indValidRGR), 'YTickLabel', fullRGR4fps(1).features(indValidRGR));
xlabel('IC #')
ylabel('Features')
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


figure
imagesc(H)
figure
[coeff,pcascore,latent,tsquared,explained] = pca(matSDF_norm');

[R_icaPC1] = corr(coeff(:,1), SourceMat', 'rows', 'complete');
size(R_icaPC1)

% Reconstruct 
estSDF = H*SourceMat;


figure;
plot(sum(abs(H)), 'ro-')
grid on
xlabel('Source')
ylabel('Sum of H')

[a,i]=sort(sum(abs(H)), 'descend');
% S_sort = S(i,:);
% matIC_sort = S_sort';
Rica_validRGR_sort = Rica_validRGR(:,i);

figure
varnames = cat(1, fullRGR4fps(1).features(indValidRGR), {'Face size'});
figPCARGR_catMov = figure;
set(figPCARGR_catMov, 'Color', 'w', 'PaperPositionMode', 'auto')
imagesc(Rica_validRGR_sort)
set(gca, 'YTick', 1:size(matRGR,2), 'YTickLabel', varnames);
set(gca, 'XTick', 1:48, 'XTickLabel', i)
% set(gca, 'FontSize', 15)
% set(gca, 'YTick', 1:length(indValidRGR), 'YTickLabel', fullRGR4fps(1).features(indValidRGR));
xlabel('IC #')
ylabel('Features')
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


% relationship with firing rate
fr = mean(matSDF)./2.4; % in spikes/s
[b, indFrSort] = sort(fr, 'descend');

figure
imagesc(H(indFrSort, :))

selCell_lowFR = indFrSort(6:end);
X_sel = matSDF(:,selCell_lowFR)';
[SourceMat_sel, H_sel, iter_sel, W_sel] = robustica(X_sel, arguments);

% matSDF_norm2 = matSDF-repmat(mean(matSDF), size(matSDF, 1), 1);
% X_norm2 = matSDF_norm2';
% [S_norm2, H_norm2, iter_norm2, W_norm2] = robustica(X_norm2, arguments);
% 
% figure
% imagesc(H_norm2)
