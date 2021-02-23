% genFig_multipleFP_fig2_FSI.m
%
% 2021/02/22 SHP
% Compute face selectivity for example neurons shown in Figure 2

clear all;

%% Figure 2 example cells
% setExampleCellIDs = {'Dav_33', 'Moc_130AF', 'Mat_097a', 'Dan_10', ...
%     'Dav_27', 'Tor_065a', 'Was_39AM', 'Moc_117AM', ...
%     'Dav_25', 'Spi_022b', 'Was_51AM', 'Moc_109AM', ...
%     'Dav_16', 'Spi_045a', 'Was_33AM', 'Dan_05', ...
%     'Dav_06', 'Moc_122AF', 'Was_06AM', 'Moc_115AM'};

setExampleCellIDs = {'33Dav', '130AFMoc', '097aMat', '10Dan'; ...
    '27Dav', '065aTor', '39AMWas', '117AMMoc'; ...
    '25Dav', '022bSpi', '51AMWas', '109AMMoc'; ...
    '16Dav', '045aSpi', '33AMWas', '05Dan'; ...
    '06Dav', '122AFMoc', '06AMWas', '115AMMoc'};

%% Read the spreadsheet to load corresponding fingerprinting data file name
filename_xls = '/procdata/parksh/_macaque/multipleFP_4FPneurons_CellIDFingerPrinting.xls';
C = readcell(filename_xls); % 1st col: cell ID in movie data, 2st col: fingerprinting results directory, 3rd col: cell file name

% load Matcha's data
load('/procdata/parksh/_macaque/Mat/_orgData/MatMov1to6_CatData.mat');

%% Cell-by-cell computation of FSI etc.
matFaceSelective = NaN(size(setExampleCellIDs));
for iCell = 1:numel(setExampleCellIDs)
    curCellID = setExampleCellIDs{iCell};
    nameSubjNeural = char(curCellID(end-2:end));
    
    % get the index from the spreadsheet 
    indCell = find(strcmp(C(:,1), curCellID)>0);
    
    if ~cellfun(@ischar, C(indCell,2)) % has zero (not character) if there's no fingerprinting results        
        matFaceSelective(iCell) = 0.5; % NaN for gray in the furture plot
        continue;
    end
    
    % load fingerprinting file
    load(sprintf('/procdata/parksh/_macaque/%s%s/%s.mat', nameSubjNeural, C{indCell, 2}, C{indCell, 3}))
    
    % compute fsi & other things
    fprintf(1, ':::::%s:::::\n', curCellID)
    curFSI = calc_fr_from_multidays_for_fsi(multiday);
    
    if abs(curFSI) > 0.33
        matFaceSelective(iCell) = 1; %face-selective
    else
        matFaceSelective(iCell) = 0; % not face-selective
    end
    
%     input('')

end

% Matcha's neurons recorded by Brian
% load('/procdata/parksh/_macaque/Mat/_orgData/MatMov1to6_CatData.mat');
response_win_Ma = [0 200]; % calculation window, ms
baseline_win_Ma = [-50 0]; % calculation window, ms
iCell_Ma = 14;
fr_resp = cellfun(@mean, MatMov1to6_CatData(iCell_Ma).cat_resp); 
fr_base= mean(cellfun(@mean, MatMov1to6_CatData(iCell_Ma).cat_base)); 


%% face-selective index (fsi)
fsi_fun = @(x,y) (x-y)/(abs(x)+abs(y))*sign(double(sign(x)>0|sign(y)>0)-0.5);
response_face    = nanmean(fr_resp(1:2)) -fr_base;  % human face, monkey face and whole monkey
response_nonface = nanmean(fr_resp(4))-fr_base;  % object 
fsi_ma = fsi_fun(response_face,response_nonface);

fprintf(1, 'Cell #%d: %03d%s: faces:%2.2f  / objects:%2.2f :: FSI: %s\n', iCell_Ma, MatMov1to6_CatData(iCell_Ma).chan(1), ...
    abc(MatMov1to6_CatData(iCell_Ma).chan(2)),response_face, response_nonface, num2str(fsi_ma))

if abs(fsi_ma) > 0.33
    matFaceSelective(11) = 1;
else
    matFaceSelective(11) = 0;
end

%% Quick overview of Fig 2b
figure;
imagesc(matFaceSelective);
colormap(gray)
hold on;

xC = repmat(1:4, 5, 1);
yC = repmat([1:5]', 1, 4);
T = text(xC(:), yC(:), setExampleCellIDs(:), 'HorizontalAlignment', 'center');
set(T(matFaceSelective<0.5), 'Color', 'w')
set(gca, 'XTick', [], 'YTick', [])


%% Fig 2a example: heatmap?
curCellID = '33Dav';
nameSubjNeural = char(curCellID(end-2:end));
    
% get the index from the spreadsheet
filename_xls = '/procdata/parksh/_macaque/multipleFP_4FPneurons_CellIDFingerPrinting.xls';
C = readcell(filename_xls); % 1st col: cell ID in movie data, 2st col: fingerprinting results directory, 3rd col: cell file name
indCell = find(strcmp(C(:,1), curCellID)>0);
    
% load fingerprinting file
load(sprintf('/procdata/parksh/_macaque/%s%s/%s.mat', nameSubjNeural, C{indCell, 2}, C{indCell, 3}))
    
response_time_win = [0 200]; %[50 250];    % calculation window, ms
baseline_time_win = [-50 0]; %[-100 50];    % calculation window, ms

% calc firing rate
time_raster = multiday.rasterwin(1):multiday.rasterwin(2);
indx_time   = time_raster>=response_time_win(1) & time_raster<=response_time_win(2);

for dayno =1:length(multiday.rasterdata)
    raster = multiday.rasterdata{dayno};
    stim   = multiday.stim{dayno};
    for stimno =1:60
        indx_stim = stim==stimno;
        raster_w_the_stim = raster(indx_time,indx_stim);
        spikes_w_the_stim = sum(raster_w_the_stim(:))/size(raster_w_the_stim,2);    % spike/trial
        fr_w_the_stim     = spikes_w_the_stim/diff(response_time_win)*1000; % Hz
        fr(dayno,stimno)  = fr_w_the_stim;
    end
end

fr = nanmean(fr,1);
baseline_fr = nanmean(multiday.fr.base(:));

max_fr = nanmax(abs(fr-baseline_fr));
norm_fr = (fr - baseline_fr)/(max_fr);

figure
plot(ones(1,20), norm_fr(1:20), 'ko', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'w');
hold on;
plot(ones(1,10).*2, norm_fr(31:40), 'ko', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'w')
xlim([0.5 2.5])
hold on
line([0.7 1.3], [mean(norm_fr(1:20)) mean(norm_fr(1:20))], 'LineWidth', 3, 'Color', 'k')
line([1.7 2.3], [mean(norm_fr(31:40)) mean(norm_fr(31:40))], 'LineWidth', 3, 'Color', 'k')

%% entire population from Matcha
for iCC = 1:30
fprintf(1, 'Cell #%d: %03d%s: baseline: %2.2f / faces:%2.2f  / objects:%2.2f\n', iCC, MatMov1to6_CatData(iCC).chan(1), abc(MatMov1to6_CatData(iCC).chan(2)), ...
    mean(cellfun(@mean, MatMov1to6_CatData(iCC).cat_base)), mean(cellfun(@mean, MatMov1to6_CatData(iCC).cat_resp([1 2]))), mean(MatMov1to6_CatData(iCC).cat_resp{4}));
end



%% Population heatmap
% 
% C = readcell(filename_xls); % 1st col: cell ID in movie data, 2st col: fingerprinting results directory, 3rd col: cell file name
% indCell = find(strcmp(C(:,1), curCellID)>0);
%     
% % load fingerprinting file
% load(sprintf('/procdata/parksh/_macaque/%s%s/%s.mat', nameSubjNeural, C{indCell, 2}, C{indCell, 3}))
%     
% response_time_win = [0 200]; %[50 250];    % calculation window, ms
% baseline_time_win = [-50 0]; %[-100 50];    % calculation window, ms
% 
% % calc firing rate
% time_raster = multiday.rasterwin(1):multiday.rasterwin(2);
% indx_time   = time_raster>=response_time_win(1) & time_raster<=response_time_win(2);
% 
% for dayno =1:length(multiday.rasterdata)
%     raster = multiday.rasterdata{dayno};
%     stim   = multiday.stim{dayno};
%     for stimno =1:60
%         indx_stim = stim==stimno;
%         raster_w_the_stim = raster(indx_time,indx_stim);
%         spikes_w_the_stim = sum(raster_w_the_stim(:))/size(raster_w_the_stim,2);    % spike/trial
%         fr_w_the_stim     = spikes_w_the_stim/diff(response_time_win)*1000; % Hz
%         fr(dayno,stimno)  = fr_w_the_stim;
%     end
% end
% 
% fr = nanmean(fr,1);
% baseline_fr = nanmean(multiday.fr.base(:));
% 
% max_fr = nanmax(abs(fr-baseline_fr));
%     norm_fr(:,cellno) = (fr - baseline_fr)/(max_fr);
%     
%     %% calculate face-selective index
%     fsi_fun = @(x,y) (x-y)/(abs(x)+abs(y))*sign(double(sign(x)>0|sign(y)>0)-0.5);
%     response_face    = nanmean(fr(1:30))-baseline_fr;    % Human face, monkey face and monkey body
% %     response_nonface = nanmean(fr(40:60))-baseline_fr;   % Object, scene and bird
% %     response_face    = nanmean(fr(1:20))-baseline_fr;    % Human face and monkey face
%     response_nonface = nanmean(fr(40:50))-baseline_fr;   % Object and scene
%     fsi(cellno,1) = fsi_fun(response_face,response_nonface);
% end
% 
% %% save data
% fsidata.fsi = fsi;
% fsidata.norm_fr = norm_fr;
% fsidata.info.cellname = cellnamelist;
% fsidata.info.duration.response = response_win;
% fsidata.info.duration.baseline = baseline_win;



