function [fsi] = calc_fr_from_multidays_for_fsi(multiday, cond_face, cond_obj)

% 2021/05/03 SHP
%       - modify to incorporate indices for face & object conditions 
%       - fixed the final check property: face selectivity criterion 
%           (before: cur_fsi < -0.33 : face selective neuron (negative))


% %% load data
% clear all
% % load /procdata/parksh/_macaque/Dav/_orgData/Davida180723/FPrint/6_037_2.mat
% load /procdata/parksh/_macaque/Dav/_orgData/Davida180723/FPrint/22_109_1.mat

%% parameters
response_time_win = [50 250];    % [0 200]; %[50 250];    % calculation window, ms
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
        fr(dayno,stimno)  = fr_w_the_stim;
    end
end

fr = nanmean(fr,1);
baseline_fr = nanmean(multiday.fr.base(:));


%% face-selective index (fsi)
% stimulus index for face & object category
ind_face = [];
for icond_face = 1:length(cond_face)
     ind_face = cat(2, ind_face, (cond_face(icond_face)-1)*10+1:cond_face(icond_face)*10);
end
ind_nonface = [];
for icond_obj = 1:length(cond_obj)
     ind_nonface = cat(2, ind_nonface, (cond_obj(icond_obj)-1)*10+1:cond_obj(icond_obj)*10);
end

fsi_fun = @(x,y) (x-y)/(abs(x)+abs(y))*sign(double(sign(x)>0|sign(y)>0)-0.5);
response_face    = nanmean(fr(ind_face)) -baseline_fr;  % face category
response_nonface = nanmean(fr(ind_nonface))-baseline_fr;  % non-face category 
fsi = fsi_fun(response_face,response_nonface);

disp(['FSI: ' num2str(fsi)]);


%% prepare firing rate and stim vector for statistics
time_raster = multiday.rasterwin(1):multiday.rasterwin(2);
indx_time_response   = time_raster>=response_time_win(1) & time_raster<=response_time_win(2);
indx_time_baseline   = time_raster>=baseline_time_win(1) & time_raster<=baseline_time_win(2);

fr_alltrial_response =[];
fr_alltrial_baseline =[];
stim_alltrial        =[];

for dayno =1:length(multiday.rasterdata)
    raster = multiday.rasterdata{dayno};
    stim   = multiday.stim{dayno};
    raster_during_response = raster(indx_time_response,:);
    raster_during_baseline = raster(indx_time_baseline,:);
    spikes_during_response = sum(raster_during_response,1);    % number of spikes
    spikes_during_baseline = sum(raster_during_baseline,1);    % number of spikes
    fr_during_response     = spikes_during_response/diff(response_time_win)*1000; % Hz
    fr_during_baseline     = spikes_during_baseline/diff(baseline_time_win)*1000; % Hz
    fr_alltrial_response   = [fr_alltrial_response; fr_during_response'];
    fr_alltrial_baseline   = [fr_alltrial_baseline; fr_during_baseline'];
    stim_alltrial          = [stim_alltrial; stim];
end

%% paired t-test (response vs baseline, each stimulus)
for stimno =1:60
    indx_stim = stim_alltrial==stimno;
    [~,p] = ttest(fr_alltrial_response(indx_stim),fr_alltrial_baseline(indx_stim));
    pval_ttest(stimno) = p;   
end

pval_ttest = pval_ttest./length(pval_ttest); % Bonferroni correction
min_p_ttest = min(pval_ttest);
best_stim = find(pval_ttest==min_p_ttest);

disp(['t-test: p<' num2str(min_p_ttest,'%.3f') ' at stim ' num2str(best_stim)]);


%% 1-way anova (across stimulus)
[p,tbl] = anovan(fr_alltrial_response, stim_alltrial,'display','off');
disp(['anova: p<' num2str(p,'%.3f')]);


%% firing rate
disp(['mean baseline firing rate is ' num2str(baseline_fr,'%.2f') ' Hz']);
disp(['mean response firing rate is ' num2str(mean(fr),   '%.2f') ' Hz, this is ' num2str(mean(fr)-baseline_fr,'%.2f') ' Hz higher than baseline']);
disp(['max  response firing rate is ' num2str(max(fr),    '%.2f') ' Hz, this is ' num2str(max(fr)- baseline_fr,'%.2f') ' Hz higher than baseline']);
disp(['min  response firing rate is ' num2str(min(fr),    '%.2f') ' Hz, this is ' num2str(min(fr)- baseline_fr,'%.2f') ' Hz different from baseline']);


%% check cell properties
disp(' ')
if fsi>0.33
    disp('face-selective response');
elseif fsi<-0.33
    disp('object-selective response'); %disp('face-selective response (negative)');
end

if min_p_ttest<0.05
    disp('significantly respond to at least one stimulus');
else
    disp('there is no significant response');
end

if p<0.05
    disp('selective response to the stimuli');
else
    disp('no selectivity for the stimuli')
end
    
if abs(mean(fr)-baseline_fr) > min_mean_fr
    disp('mean response is high enough');
else
    disp('mean respone is too weak');
end

if (abs(max(fr)- baseline_fr) > min_max_fr) || (abs(min(fr)- baseline_fr) > min_max_fr)
    disp('max response is high enough');
else
    disp('max respone is too weak');
end



