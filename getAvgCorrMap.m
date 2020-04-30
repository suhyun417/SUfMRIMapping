for iChan=1:size(S,1)
mapR_catMov{iChan} = cat(4, S(iChan,:).mapR);
meanR_mov{iChan} = mean(mapR_catMov{iChan}, 4);
end

mapR_catAll = cat(4, meanR_mov{:});
meanR_all = mean(mapR_catAll, 4);

DSP.proc.fracvarmap_3d = meanR_all; %S(iChan).mapR.*-1; %S(6).mapR.*-1;


for iChan=1:size(S,1)
mapR_catMov_rect{iChan} = abs(cat(4, S(iChan,:).mapR));
meanR_mov_rect{iChan} = mean(mapR_catMov_rect{iChan}, 4);
end

mapR_catAll_rect = cat(4, meanR_mov_rect{:});
meanR_all_rect = mean(mapR_catAll_rect, 4);