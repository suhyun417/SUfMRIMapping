function generate_sort_report(data_path,sorted_file,save_path)

if strcmp(sorted_file(end-3:end),'.mat')
    sorted_file = sorted_file(1:end-4);
end

if nargin < 3
   save_path = [data_path filesep sorted_file '_sort_report']; 
end

try
    load([data_path sorted_file],'DSP*','WAV*')
catch
    load([data_path sorted_file '.mat'],'DSP*','WAV*')
end

%get channel list
var_list = whos('WAV*');
chan_N = nan(size(var_list));
for ii = 1:length(var_list)
   curr_var_name = char(var_list(ii).name);
   chan_N(ii) = str2double(curr_var_name(4:6));
end
chan_N = unique(chan_N);

for ii = 1:length(chan_N)
   curr_ch = chan_N(ii);
   if curr_ch < 10
       wave_name = ['WAV00' num2str(curr_ch)];
   elseif curr_ch < 100
       wave_name = ['WAV0' num2str(curr_ch)];
   else
       wave_name = num2str(curr_ch);
   end
   
   append_list = ['i' 'a' 'b' 'c' 'd'];
   
   clear spike_cells
   for jj = 1:length(append_list)
       curr_waves = [wave_name append_list(jj)];
       if exist(curr_waves,'var')
           eval(sprintf('spike_cells(jj,1) = {%s.Data};',curr_waves))
       end
   end  
   
   dspname = wave_name;
   dspname(1:3) = 'DSP';
   file_name = data_path;
   if strcmp(file_name(end),filesep)
       file_name(end) = [];
   end
   file_name = regexp(file_name,filesep,'split');
   file_name = char(file_name(end));
   
   get_spike_pca(spike_cells,[file_name ' ' dspname])
   if ~isempty(save_path)
       if ~isdir(save_path)
          mkdir(save_path) 
       end
%        saveas(gcf,[save_path dspname '.fig'],'fig')
       saveas(gcf,[save_path filesep dspname '.pdf'],'pdf')
   end    
   close(gcf)
end

