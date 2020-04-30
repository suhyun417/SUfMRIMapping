% xcorr_continuous - calculate cross-correlogram from session data
% 
% written by Kenji Koyano, 12/13/2015
% last modified -12/14/2015


clear all;
close all;

%% Parameters
directory.cells     = '/procdata/koyanok/physiology/cells/';   % data to load
directory.save      = '/procdata/koyanok/physiology/xcorr/';

% list.tank  = {'Matcha151007';'Matcha151008';'Matcha151015';'Matcha151023';'Matcha151030'};
list.tank  = {'Dango170208'};
% list.block = {'Resting1'};
list.block = {'PView1'};

params.duration_xcorr    = [-100 100];   % cross-correlogram range, ms
params.binw = 1;                    % resolution, ms
params.shift     = 1000;                 % time-shift for shift predictor, ms
params.duration_baseline = [-100 100];   % baeline in correlogram, ms
params.duration_peak = [-10 10];         % range for peak detection, ms
params.threshold = 0.001;                % peak height threshold, percentile
params.threshold_sd = [norminv(params.threshold), norminv(1-params.threshold)]; % peak height threshold in SD
indx_baseline = (1+round((params.duration_baseline(1)-params.duration_xcorr(1))./params.binw)):...
                (diff(params.duration_xcorr)./params.binw+1-round((params.duration_xcorr(2)-params.duration_baseline(2))./params.binw));
indx_peak     = (1+round((params.duration_peak(1)-params.duration_xcorr(1))./params.binw)):...
                (diff(params.duration_xcorr)./params.binw+1-round((params.duration_xcorr(2)-params.duration_peak(2))./params.binw));

if exist('list','var')
    if isfield(list,'block')
        block_specified = true;
    else
        block_specified = false;
    end
end

if exist('list','var')
    if ~isfield(list,'tank')
        tmp_folderlist = dir(directory.cells);
        isd = [tmp_folderlist(:).isdir]; % only use folders
        list.tank = {tmp_folderlist(isd).name}';
        list.tank(ismember(list.tank,{'.','..'})) = []; % remove the . and .. from the list
        clear tmp_folderlist isd
    end
end

%% Loop for tank
for tankno = 1:length(list.tank)
    tankname = list.tank{tankno};
    disp(tankname);
    if ~block_specified
        tmp_folderlist = dir([directory.cells tankname filesep 'Resting*']);
        isd = [tmp_folderlist(:).isdir]; % only use folders
        list.block = {tmp_folderlist(isd).name}';
        list.block(ismember(list.block,{'.','..'})) = []; % remove the . and .. from the list
        clear tmp_folderlist isd
    end
    
    %% Loop for block
    for blockno = 1:length(list.block)
        blockname = list.block{blockno};
        disp(blockname);
        tmp_matlist = dir([directory.cells tankname filesep blockname]);
        isd = [tmp_matlist(:).isdir];
        list.cell = {tmp_matlist(~isd).name}';
        list.cell(ismember(list.cell,{'tankdat.mat'})) = []; % remove the . and .. from the list
        clear tmp_matlist isd
        
        %% Load timestamp data
        disp('loading data...')
        for cellno = 1:length(list.cell)
            disp(['  ' list.cell{cellno}])
            load([directory.cells tankname filesep blockname filesep list.cell{cellno}]);
            ts_all{cellno}      = celldata.ts;
            ch_all(cellno)      = celldata.info.ch;
            cluster_all(cellno) = celldata.info.cluster;
        end
        %% Calculate correlogram
        %% Loop for 1st cell
        disp('calculating cross correlogram...');
        for cellno1 = 1:length(list.cell)
            disp(['cell ' list.cell{cellno1} ', ' num2str(cellno1) '/' num2str(length(list.cell))]);
            for cellno2 = 1:length(list.cell)
                %% raw correlogram
                data1 = ts_all{cellno1};
                data2 = ts_all{cellno2};
                mat1 = repmat(data1,[1 length(data2)]);
                mat2 = repmat(data2',[length(data1) 1]);
                mat_subtracted = mat2-mat1;
                if cellno1==cellno2
                    mat_subtracted(mat_subtracted==0) = nan;
                end
                indx = mat_subtracted>(params.duration_xcorr(1)-0.5) & mat_subtracted<(params.duration_xcorr(2)+0.5);
                xcorr_dat.raw{cellno1,cellno2} = mat_subtracted(indx);
                xcorr_tmp = histc(mat_subtracted(indx),(params.duration_xcorr(1)-0.5):params.binw:(params.duration_xcorr(2)+0.5));
                if isempty(xcorr_tmp)
                    xcorr_tmp = zeros(diff(params.duration_xcorr)+1);
                end
                xcorr_dat.hist{cellno1,cellno2} = xcorr_tmp(1:end-1);
                xcorr_dat.n_spike(cellno1,cellno2)     = length(data1)+length(data2);
                xcorr_dat.n_spike_min(cellno1,cellno2) = nanmin([length(data1) length(data2)]);
                
                %% shift predictor
                shift_data2 = data2-params.shift;
                shift_data2(shift_data2<0) = data2(shift_data2<0)+max(shift_data2);
                shift_mat2 = repmat(shift_data2',[length(data1) 1]);
                mat_subtracted = shift_mat2-mat1;
                indx = mat_subtracted>(params.duration_xcorr(1)-0.5) & mat_subtracted<(params.duration_xcorr(2)+0.5);
                xcorr_dat.shift{cellno1,cellno2} = mat_subtracted(indx);
                xcorr_tmp = histc(mat_subtracted(indx),(params.duration_xcorr(1)-0.5):params.binw:(params.duration_xcorr(2)+0.5));
                if isempty(xcorr_tmp)
                    xcorr_tmp = zeros(diff(params.duration_xcorr)+1);
                end
                xcorr_dat.shift_hist{cellno1,cellno2} = xcorr_tmp(1:end-1);
                xcorr_dat.shift_baseline_mean(cellno1,cellno2) = nanmean(xcorr_tmp(indx_baseline));
                xcorr_dat.shift_baseline_sd(cellno1,cellno2)   = nanstd (xcorr_tmp(indx_baseline));
                
                %% subtract
                xcorr_dat.subtracted_hist{cellno1,cellno2} = xcorr_dat.hist{cellno1,cellno2} - xcorr_dat.shift_baseline_mean(cellno1,cellno2);
                xcorr_dat.zscore_hist{cellno1,cellno2} = xcorr_dat.subtracted_hist{cellno1,cellno2}./xcorr_dat.shift_baseline_sd(cellno1,cellno2);
                
                %% peak detection, statistics
                xcorr_tmp = xcorr_dat.zscore_hist{cellno1,cellno2};
                xcorr_tmp = xcorr_tmp(indx_peak);
                xcorr_left   = xcorr_tmp(1:-params.duration_peak(1));
                xcorr_center = xcorr_tmp(1-params.duration_peak(1));
                xcorr_right  = xcorr_tmp((2-params.duration_peak(1)):end);
                positive_peak = [max(xcorr_left), xcorr_center, max(xcorr_right)];
                negative_peak = [min(xcorr_left), xcorr_center, min(xcorr_right)];
                xcorr_dat.peak_count_positive(cellno1,cellno2,:) = positive_peak;
                xcorr_dat.peak_count_negative(cellno1,cellno2,:) = negative_peak;
                if any(isnan(positive_peak))
                    xcorr_dat.peak_time_positive(cellno1,cellno2,:)  = [nan nan nan];
                    xcorr_dat.peak_time_negative(cellno1,cellno2,:)  = [nan nan nan];
                else
                    xcorr_dat.peak_time_positive(cellno1,cellno2,:)  = [find(xcorr_left==positive_peak(1),1,'last')-length(xcorr_left)-1, 0, find(xcorr_right==positive_peak(3),1,'first')].*params.binw;
                    xcorr_dat.peak_time_negative(cellno1,cellno2,:)  = [find(xcorr_left==negative_peak(1),1,'last')-length(xcorr_left)-1, 0, find(xcorr_right==negative_peak(3),1,'first')].*params.binw;
                end
                xcorr_dat.peak_positive(cellno1,cellno2,:) = positive_peak>params.threshold_sd(2);
                xcorr_dat.peak_negative(cellno1,cellno2,:) = negative_peak<params.threshold_sd(1);
                rec_length = nanmax([data1; data2]);
                if isempty(rec_length)
                    rec_length = 60000;
                end
                xcorr_dat.ncc_positive(cellno1,cellno2,:) = positive_peak*sqrt(params.binw/rec_length); % neural correlation coefficient (Eggermont, 1992, JNP)
                xcorr_dat.ncc_negative(cellno1,cellno2,:) = negative_peak*sqrt(params.binw/rec_length);
                sum_right = sum(abs(xcorr_right(abs(xcorr_right)>params.threshold_sd(2))));
                sum_left  = sum(abs(xcorr_left (abs(xcorr_left) >params.threshold_sd(2))));
                sum_center = abs(xcorr_center(abs(xcorr_center) >params.threshold_sd(2)));
                if isempty(sum_right);  sum_right = 0; end;
                if isempty(sum_left);   sum_left = 0; end;
                if isempty(sum_center); sum_center = 0; end;
                if (sum_right+sum_left) ~= 0
                    xcorr_dat.ai(cellno1,cellno2) = (sum_right-sum_left)/(sum_right+sum_left);  % assymmetry index (Alonso et al., 1998, Nat. Neurosci)
                else
                    xcorr_dat.ai(cellno1,cellno2) = 0;
                end
                xcorr_dat.cs(cellno1,cellno2) = sum_right + sum_left + sum_center; % coupling strength (Hirabayashi et al., 2013, Science)
            end
        end
        xcorr_dat.list_cell    = list.cell;
        xcorr_dat.list_ch      = ch_all;
        xcorr_dat.list_cluster = cluster_all;
        %% Save .mat
        if ~exist([directory.save tankname],'dir')
            mkdir(directory.save, tankname);
        end
        if ~exist([directory.save tankname filesep blockname],'dir')
            mkdir([directory.save tankname], blockname);
        end
        save([directory.save tankname filesep blockname filesep 'xcorr_dat.mat'],'params','xcorr_dat');
        clear xcorr_dat list_cell ts_all ch_all cluster_all;
    end
end


