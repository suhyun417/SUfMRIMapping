function [imageResp] = saveFingerPrintingData_multipleFP_ForHarish(flagSave)
% modified from "computeFSI.m"
% 
% 2021/08/17 SHP
%       face-selective index using responses to fingerprinting
%       stimulus set for all multiple face patch neurons 
%       (except Tor, Rho, Sig, Spi 2016 data that don't have fingerprinting
%       results)

%% Parameters
response_win_Ma = [0 200]; % calculation window, ms
baseline_win_Ma = [-50 0]; % calculation window, ms
cond_face = 1:2; % human face & monkey face
cond_obj = 4; % object


%% FSI
%% Read the spreadsheet to load corresponding fingerprinting data file name
filename_xls = '/nifvault/procdata/parksh/_macaque/multipleFP_4FPneurons_CellIDFingerPrinting.xls';
C = readcell(filename_xls); % 1st col: cell ID in movie data, 2st col: fingerprinting results directory, 3rd col: cell file name

% load Matcha's data
load('/nifvault/procdata/parksh/_macaque/Mat/_orgData/MatMov1to6_CatData.mat');
% colHeadings = {'Category', 'StimID', 'SpikeTime', 'SDF', 'etc'};


%% Cell-by-cell computation of FSI etc.
% matFaceSelective = NaN(size(C, 1), 2);
imageResp = struct();

for iCell = 1:size(C, 1) %numel(setExampleCellIDs)
    curCellID = C{iCell, 1};
    nameSubjNeural = char(curCellID(end-2:end));
    imageResp(iCell).curCellID = curCellID;
    imageResp(iCell).nameSubj = nameSubjNeural;
    
    % get the index from the spreadsheet 
    indCell = find(strcmp(C(:,1), curCellID)>0);
    
    if ~cellfun(@ischar, C(indCell, 2)) % has zero (not character) if there's no fingerprinting results   
        imageResp(iCell).fr_resp = [];
        imageResp(iCell).fr_base = [];
%         matFaceSelective(iCell, 1) = NaN; %
%         matFaceSelective(iCell, 2) = 0; %0.5; % NaN for gray in the furture plot
        continue;
    end
    
        
    fr_base = [];
    fr_resp = [];
    
    if ~strcmpi(nameSubjNeural, 'mat') % everyone except Matcha   
        % load fingerprinting file
        load(sprintf('/nifvault/procdata/parksh/_macaque/%s%s/%s.mat', nameSubjNeural, C{indCell, 2}, C{indCell, 3}))
        
        % compute fsi & other things
        fprintf(1, ':::::%s:::::\n', curCellID)
        [fr_resp, fr_base] = calc_fr_from_multidays(multiday, cond_face, cond_obj);
        fr_resp_norm = fr_resp - fr_base;
    else % in case of Matcha's data
        iCell_Ma = C{indCell, 3};
        
        tempC = MatMov1to6_CatData(iCell_Ma).cat_spikes;
%         tempS = cell2struct(tempC, colHeadings, 2);
        catStimID = cell2mat(tempC(:,2)); %cat(1, tempS.StimID);
        
        for ii = 1:6 % 6 categories
            aa = ismember(catStimID, 40*(ii-1)+1:40*ii);
            
            for jj = 1:40
                indStim = catStimID(aa) == 40*(ii-1)+jj;
                fr_base(40*(ii-1)+jj) = mean(MatMov1to6_CatData(iCell_Ma).cat_base{ii}(indStim));
                fr_resp(40*(ii-1)+jj) = mean(MatMov1to6_CatData(iCell_Ma).cat_resp{ii}(indStim));
            end
        end
        fr_base = mean(fr_base);
        
        fr_resp_norm = fr_resp - fr_base;
        
%         fr_resp = cellfun(@mean, MatMov1to6_CatData(iCell_Ma).cat_resp);
%         fr_base= mean(cellfun(@mean, MatMov1to6_CatData(iCell_Ma).cat_base));
%         
%         % face-selective index (fsi)
%         fsi_fun = @(x,y) (x-y)/(abs(x)+abs(y))*sign(double(sign(x)>0|sign(y)>0)-0.5);
%         response_face    = nanmean(fr_resp(cond_face)) -fr_base;  % human face, monkey face and whole monkey
%         response_nonface = nanmean(fr_resp(cond_obj))-fr_base;  % object
%         curFSI = fsi_fun(response_face,response_nonface);
    end    
    
    
    imageResp(iCell).fr_base = fr_base;
    imageResp(iCell).fr_resp = fr_resp;
    imageResp(iCell).fr_resp_norm = fr_resp_norm;
        
%     matFaceSelective(iCell, 1) = curFSI;
%     
%     if curFSI > 0.33
%         matFaceSelective(iCell, 2) = 1; %face-selective
%     else
%         matFaceSelective(iCell, 2) = -1; % not face-selective
%     end
    
end

% fsi.matFSI = matFaceSelective;
% fsi.catChanID = C(:,1);

if flagSave
    save('/nifvault/procdata/parksh/_macaque/_Harish/multipleFP_ImageResp.mat', 'imageResp')
end

end
 
function [fr_resp, fr_base] = calc_fr_from_multidays(multiday, cond_face, cond_obj)

% 2021/05/03 SHP
%       - modify to incorporate indices for face & object conditions 
%       - fixed the final check property: face selectivity criterion 
%           (before: cur_fsi < -0.33 : face selective neuron (negative))


% %% load data
% clear all
% % load /procdata/parksh/_macaque/Dav/_orgData/Davida180723/FPrint/6_037_2.mat
% load /procdata/parksh/_macaque/Dav/_orgData/Davida180723/FPrint/22_109_1.mat

%% parameters
response_time_win =  [0 200]; %[50 250];    % [0 200]; %[50 250];    % calculation window, ms
baseline_time_win = [-50 0]; %[-100 50];    % calculation window, ms

min_mean_fr = 1;  % minimum mean firing rate required
min_max_fr = 3;  % minimum max  firing rate required

%% calc firing rate
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
        fr_resp(dayno,stimno)  = fr_w_the_stim;
    end
end

fr_resp = nanmean(fr_resp,1);
fr_base = nanmean(multiday.fr.base(:));

% 
% %% face-selective index (fsi)
% % stimulus index for face & object category
% ind_face = [];
% for icond_face = 1:length(cond_face)
%      ind_face = cat(2, ind_face, (cond_face(icond_face)-1)*10+1:cond_face(icond_face)*10);
% end
% ind_nonface = [];
% for icond_obj = 1:length(cond_obj)
%      ind_nonface = cat(2, ind_nonface, (cond_obj(icond_obj)-1)*10+1:cond_obj(icond_obj)*10);
% end
% 
% fsi_fun = @(x,y) (x-y)/(abs(x)+abs(y))*sign(double(sign(x)>0|sign(y)>0)-0.5);
% response_face    = nanmean(fr(ind_face)) -baseline_fr;  % face category
% response_nonface = nanmean(fr(ind_nonface))-baseline_fr;  % non-face category 
% fsi = fsi_fun(response_face,response_nonface);

end

% %% loop for cell
% norm_fr =[];    % normalized firing rate
% for cellno = 1:length(celllist)
%     [~,cellname,~] = fileparts(celllist{cellno});
%     disp(cellname);
%     
%     cellnamelist{1,cellno} = cellname;
%     
%     %% load data and format raster data into an array
%     load([directory.session cellname '.mat']);
%     
%     %% calculate firing rate of each trial
%     response_calculate_win = response_win - raster.win(1);
%     fr_trial = sum(raster.data(response_calculate_win(1):response_calculate_win(2),:),1)/diff(response_win)*1000;  % Hz
%         
%     %% loop for stimulus
%     fr = nan(length(code.event),1);
%     for stimno = 1:length(code.event)
%         stimcode = code.event(stimno);
%         stim_indx = raster.stimcode==stimcode;
%         fr(stimno) = nanmean(fr_trial(stim_indx));
%     end
%     fr = fr(plot_order);    % sort stimuli
%         
%     %% normalize firing rate
%     baseline_calculate_win = baseline_win - raster.win(1);
%     baseline_fr_trial = sum(raster.data(baseline_calculate_win(1):baseline_calculate_win(2),:),1)/diff(baseline_win)*1000;  % Hz
%     baseline_fr = nanmean(baseline_fr_trial);
% %     max_fr = nanmax(fr);
%     max_fr = nanmax(abs(fr-baseline_fr));
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
% 
% if ~exist(directory.save)
%     mkdir(directory.save);
% end
% save([directory.save 'fsidata_AF' TankName '.mat'],'fsidata');
% 
% 
% %% sort normalized firing rate matrix
% sort_indx = [fsi,(1:length(fsi))'];
% sort_indx = sortrows(sort_indx);
% sort_indx = sort_indx(:,2);
% 
% %% plot population heat map, sorted by fsi
% figure('Position',[200 200 600 400],'Color','w');
% hold on;
% imagesc(norm_fr(:,sort_indx)');
% for ii = 1:length(category_boundary)
%     plot([category_boundary(ii),category_boundary(ii)]+0.5,[0, size(norm_fr,2)]+0.5,'w');
% end
% set(gca,'XTick',(category_boundary(1:end-1)+diff(category_boundary)/2)+0.5);
% set(gca,'XTickLabel',category_txt);
% set(gca,'YTick',1:length(celllist),'YTickLabel',cellnamelist(sort_indx));
% xlabel('Stimuli');
% ylabel('Cells');
% xlim([0 size(norm_fr,1)]+0.5);
% ylim([0 size(norm_fr,2)]+0.5);
% title([TankName ', Finger Printing (n=' num2str(length(celllist)) ')'],'Interpreter','none','FontSize',14);
% colormap('hot');
% colorbar;
% 
% % print(gcf,[directory.save 'heatmap_population' TankName '.eps'],'-depsc2');
%      
     