%% parameters

data_directory = '/procdata/parksh/Mat/_orgData/';

for iM = 2:6
    clear data
    
filename =  sprintf('MatMov1to6_Movie%02dData', iM); % 'AvaRight01_Movie06Data';

plot_range = [-20 20];
plot_resolution = 0.5;
trialno = 6;

%% load data
disp('loading data...')
load([data_directory filename '.mat']);
data = eval(filename);
clear(filename)
disp('done.')
fprintf('\n')

%% loop for cell combinations

ncells = length(data);
vec_nspikes = nan([ncells,1]);
mat_ratio_centerpeak1 = nan(ncells);
mat_ratio_centerpeak2 = nan(ncells);
mat_ccg = nan([ncells, ncells, diff(plot_range)/plot_resolution+1]);

for ii=1:ncells
    disp(['cell ' num2str(data(ii).chan(1)) '_' num2str(data(ii).chan(2))]);
    spike1  = data(ii).mov_spikes{trialno,2};
    nspike1 = length(spike1);
    if nspike1==0
        continue
    end
    vec_nspikes(ii) = nspike1;
    for jj=1:ncells
        disp(['  cell ' num2str(data(jj).chan(1)) '_' num2str(data(jj).chan(2))]);
        spike2 = data(jj).mov_spikes{trialno,2};
        nspike2 = length(spike2);
        if nspike2==0
            continue
        end
        
        spikemat1 = repmat(spike1, [1, nspike2]);
        spikemat2 = repmat(spike2',[nspike1, 1]);
        
        mat_subtracted = spikemat1 - spikemat2;
        mat_subtracted(mat_subtracted<plot_range(1) | mat_subtracted>plot_range(2)) = nan;
        
        ccg = hist(mat_subtracted(:),plot_range(1):plot_resolution:plot_range(2));
        center_bin = -plot_range(1)/plot_resolution;
        centerpeak_height = max(ccg(center_bin:center_bin+2));
        
        ratio_centerpeak1 = centerpeak_height/nspike1;
        ratio_centerpeak2 = centerpeak_height/nspike2;
        
        mat_ccg(ii,jj,:) = ccg;
        mat_ratio_centerpeak1(ii,jj) = ratio_centerpeak1;
        mat_ratio_centerpeak2(ii,jj) = ratio_centerpeak2;
        
    end
end

figure
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
imagesc(mat_ratio_centerpeak1)
set(gca, 'CLim', [0 0.3])
colorbar
title(filename)
dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';
print(gcf, fullfile(dirFig, ['xcorr_', filename]), '-r150', '-dtiff')

end


% cell selection
for iC = 1:length(mat_ratio_centerpeak1)
tempCount(iC,1) = sum(mat_ratio_centerpeak1(iC, :)>0.3);
tempCount(iC,2) = sum(mat_ratio_centerpeak1(:, iC)>0.3);
end
indCellValid = find(sum(tempCount,2)<3);
validChanID= cat(1, data(indCellValid).chan);
% in addition to this, some manual addition / exclusion is done


