% computeEigenneuron_indMov.m

% compute eigenneuron for each individual movie 
% and then concatenate across different movies
% and then use it as a regressor

homeDir_eig = '/procdata/parksh/Tor/eigen/';

setMovID = [1 2 3 10 11 12 13 14 15];
validC = [1:8, 11, 13, 15:18, 20:21, 23, 25:28, 30:34, 36:43, 46:49, 51]; % valid channel with all 9 movie data

paramPcaMov.validCell = validC;
paramPcaMov.monkey = 'Toroid';

for iMov=1:length(setMovID)
    
    matSDF=[];matSDF_norm=[];
    matSDF = cat(2,S(validC,iMov).mnsdf)';
    matSDF_norm = matSDF;
    
%     matSDF_norm = matSDF-repmat(mean(matSDF,2), 1, size(matSDF, 2));
%     matSDF_norm = matSDF_norm./repmat(std(matSDF')', 1, size(matSDF,2));
    
    [coeff,pcascore,latent,tsquared,explained] = pca(matSDF_norm);
    
    coeff = coeff(1:100:end, :);
    
    pcaMov(iMov).movID = setMovID(iMov);
    pcaMov(iMov).coeff = coeff;
    pcaMov(iMov).pcascore = pcascore;
    pcaMov(iMov).expVar = explained;
    
    
%     % test the other way
%     matSDF=[];matSDF_norm=[];
%     matSDF = cat(2,S(validC,iMov).mnsdf);
%     
%     matSDF_norm = matSDF-repmat(mean(matSDF), size(matSDF, 1),1);
%     matSDF_norm = matSDF_norm./repmat(std(matSDF), size(matSDF,1),1);
%     
%     [coeff,pcascore,latent,tsquared,explained] = pca(matSDF_norm);

end

% for 
iMov=1:3;
    
    matSDF=[];matSDF_norm=[];
    matSDF = cat(2,S(validC,iMov).mnsdf)';
    matSDF_norm = matSDF;
    
%     matSDF_norm = matSDF-repmat(mean(matSDF,2), 1, size(matSDF, 2));
%     matSDF_norm = matSDF_norm./repmat(std(matSDF')', 1, size(matSDF,2));
    
    [coeff,pcascore,latent,tsquared,explained] = pca(matSDF_norm);
    
    coeff = coeff(1:100:end, :);
    
    pcaMov123.movID = setMovID(iMov);
    pcaMov123.coeff = coeff;
    pcaMov123.pcascore = pcascore;
    pcaMov123.expVar = explained;
    
    
    % test the other way
    matSDF=[];matSDF_norm=[];
    matSDF = cat(2,S(validC,iMov).mnsdf);
    
    matSDF_norm = matSDF-repmat(mean(matSDF), size(matSDF, 1),1);
    matSDF_norm = matSDF_norm./repmat(std(matSDF), size(matSDF,1),1);
    
    [coeff,pcascore,latent,tsquared,explained] = pca(matSDF_norm);

% end

figure
plot(coeff(:,1), 'o-')

pcaDimen = 1;
pcaMov
size(pcaMov(1).pcascore)
[dum, si] = sort(pcaMov(1).pcascore(:, pcaDimen));
si = flipud(si);
size(si)
si
size(mnSDF_org(1).cellIDs)
paramPcaMov.validCell
sname = paramPcaMov.validCell(si)
S
ssdf = cat(2,S(sname,1).mnsdf);
size(ssdF)
size(ssdf)
ssdf = ssdf(1:100:300000, :);
figure
plot(ssdf)
ssdf_norm = zscore(ssdf);
imagesc(ssdf)
imagesc(ssdf')
imagesc(ssdf_norm')
ls
edit icicleHeatmapFigure.m

saveFileName = 'pcaIndMovie_tor.mat';
save(saveFileName, 'pcaMov', 'paramPcaMov')


% % quick & dirty way to compute and see correlation maps
% catPC_mov = cat(1,pcaMov.coeff);
% matBOLD = cat(dataBOLD.catmvoltc
% 
% 
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



% cmap = jet;
% win = [0 300];
% kernel = 0.5;
% pp = setpathsMovies;
% subj = 'tsr';
% mov = '123';
% flist1 = dir([pp.mas 'tmov1sig*']);
% flist2 = dir([pp.mas 'tmov2sig*']);
% flist3 = dir([pp.mas 'tmov3sig*']);
% snames = [];
% subjResp = [];
% for s=1:3
%     disp(['monk ' subj(s)]);
%     movResp = [];
%     for m=1:3
%         disp(['mov ' mov(m)]);
%         fstem = [subj(s) 'mov' mov(m) 'sig*.mat'];
%         flist = dir([pp.mas fstem]);
%         spkResp = [];
%         for f=1:length(flist)
%             tmp = load([pp.mas flist(f).name]);
%             dat = tmp.dat;
%             dat = movieTimeScale(dat,'sec');
%             unpack;
%             [spk time] = getSdf(dat,win,kernel);
%             spkResp = [spkResp ; spk];
%             if m==1
%                 snames = [snames ; [subj(s) dat.h.snames]];
%             end
%         end
%         movResp = [movResp spkResp];
%     end
%     breaks(s) = size(snames,1);
%     subjResp = [subjResp ; movResp];
% end
% 
% % %  normalize, z-transform
% sdf = subjResp;
% sdf = sdf-repmat(mean(sdf,2),[1 size(sdf,2)]);
% sdf = sdf./repmat(std(sdf')',[1 size(sdf,2)]);
% 
% 
% % sort the signal in different ways
% [coeff,pcascore,latent,tsquared,explained] = pca(sdf);
