%% response calculation from rasterdata of multiday

function [matFR_fp nDate] = response_calc_SHP(filename)

%  dirFPrint = fullfile(dirDataNeural, '_orgData', 'Davida180515_other', 'FPrint');
%  aa = dir(fullfile(dirDataNeural, '_orgData', 'Davida180515_other', 'FPrint', '*.mat'));
load(filename)

window_response = [50 150]; % response calculation window. unit=ms. from stimulus onset

raster_allday =  horzcat(multiday.rasterdata{:});   % concatenate raster data across days
stim_allday   =  vertcat(multiday.stim{:});         % concatenate stimulus data across days

window_rasterdata = window_response - multiday.rasterwin(1);                    % response window in rasterdata matrix
raster_allday     = raster_allday(window_rasterdata(1):window_rasterdata(2),:);   % raster within the response window
response_alltrial = sum(raster_allday,1) ./ diff(window_response) .* 1000;      % response of each trial, unit=Hz

fr_response =[]; % vector for response

for ii = 1:length(multiday.stimlist)    % loop for each stimulus
    indx_stim = stim_allday==ii;        % index for the stimulus
    fr_response(ii) = nanmean(response_alltrial(indx_stim));    % average response of the stimulus across trials
end
    
matFR_fp = reshape(fr_response, 10, 6);
nDate = length(multiday.date);

end