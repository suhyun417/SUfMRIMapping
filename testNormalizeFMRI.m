% testNormalizeFMRI.m
%
% short script for testing different ways of normalization of fMRI
% timeseries


% movie
setMovie = [1 2 3];
% regressor (1st PC)
rgr = resample(pcaMovCat.coeff(:,1), 0.1*10, 2.4*10); % time x 1st PC
% test voxel
a = 36; b = 29; c = 16;
ind = sub2ind([40 64 32], a, b, c);

% 1. original way
fmritc = cat(4, dataBOLD.catmvoltc{setMovie});
[nx, ny, nz, nt] = size(fmritc); 
nVox = nx*ny*nz;
reshapedfmritc = reshape(fmritc, nVox, nt)';

R_org = corr(reshapedfmritc(:, ind), rgr, 'rows', 'complete').*(-1); % -0.1787


% 2. normalize before concatenation (i.e. change to percent signal for each
% movie)
reshapedfmritc2 = [];
for m = 1:length(selMovID_fmri)
    curtc = dataBOLD.catmvoltc{selMovID_fmri(m)};
    [nx, ny, nz, nt] = size(curtc);
    nVox = nx*ny*nz;
    tempTC = reshape(curtc, nVox, nt)';z`
    tempTC_sub=tempTC-repmat(nanmean(tempTC), size(tempTC,1),1);
    tempPS = 100.*tempTC_sub./repmat(nanmean(tempTC), size(tempTC,1),1);

    reshapedfmritc2 = cat(1, reshapedfmritc2, tempPS);
    tempPS=[]; tempTC=[]; tempTC_sub=[];
end

R_2 = corr(reshapedfmritc2(:, ind), rgr, 'rows', 'complete').*(-1); %-0.1787

% 3. normalize after concatenation 
tempTC = reshapedfmritc;
tempTC_sub=tempTC-repmat(nanmean(tempTC), size(tempTC,1),1);
tempPS = 100.*tempTC_sub./repmat(nanmean(tempTC), size(tempTC,1),1);
reshapedfmritc3 = tempPS;

R_3 = corr(reshapedfmritc3(:, ind), rgr, 'rows', 'complete').*(-1); % -0.1787

%%% Summary
%%% 1. Correlation wasn't affected by the way of normalization 
%%% 2. Normalization before concatenation (i.e. for each movie) and after
%%% concatenation (i.e. for concatenated movie) were actually the same.






