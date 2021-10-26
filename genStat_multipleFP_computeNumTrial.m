% genStat_multipleFP_computeNumTrial.m
%
% 2021/02/17 SHP
% This code computes number of trials (viewings) for each subject for each
% of the movie 1-3

setSubj = {'Tor', 'Rho', 'Sig', 'Spi', 'Dav', 'Mat', 'Was', 'Dan', 'Moc'};
matNumTrial = NaN(length(setSubj), 3);
for iSubj = 1:length(setSubj)
    load(sprintf('/procdata/parksh/_macaque/%s/%s_movieTS_SU_indMov.mat', setSubj{iSubj}, setSubj{iSubj}));
    for iMovie = 1:3
        iCell = 1;
        if isempty(S(iCell,1).matFR)
            iCell = 2;
        end
        matNumTrial(iSubj, iMovie) = size(S(iCell,iMovie).matFR{1}, 2);
    end
end

paramStat.setSubj = setSubj;
resultsStat.matNumTrial = matNumTrial;
resultsStat.matNumTrial_description = 'subjects (rows) by movies 1-3 (columns)';
resultsStat.range = [min(matNumTrial(:)) max(matNumTrial(:))];
resultsStat.median = median(matNumTrial(:));
resultsStat.mean = mean(matNumTrial(:));

save('/procdata/parksh/_macaque/multipleFP_SU_numTrials.mat', 'paramStat', 'resultsStat')