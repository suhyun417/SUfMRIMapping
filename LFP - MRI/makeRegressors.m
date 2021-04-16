function makeRegressors
% This function makes regressors of the average power in each frequency
% band to use in fMRI analyses
% USAGE:    makeRegressorsaw
% INPUT:    BLP/CSP is a structure containing, in BLP.blp_dat, the 4-dimensional 
%           BLP/CSP data, such that the first dimension represents channel, the 
%           second dimension represents time, the third dimension represents
%           epoch and the fourth dimension represents bandwidth.
% OUTPUT:   vectors for each channel and bandwidth, representing the mean
%           power (BLP or CSP) per MRI volume, e.g. 'ch1_delta'

% user input
date = '11-08-08';   date2 = '110808';   sess = '_0002';   monkey = 'Varia';
code = 'OP1_27';     % code of monkey _ scan nr
shot = 1;            % number of shots in multishot EPI (1 for single shot)
pad = 0;             % number of -1000's added at the end of regressor in case BV stopped earlier than scanner; default 0
indiv = true;        % indiv = true means all channels get their own regressors
csp = false;         % csp = true means these are CSP regressors we're computing
groups = [3 7;11 15;19 23];         % groups of channels to compute the mean regressor of 
gr_names = {'upper','middle','lower'};
bw_names = {'delta_new';'theta_new';'alpha_new';'beta_new';'gamma_new';'MUA'};

% get to right channel index in matrix
% if csp
%    groups = (groups-1)/2;     % to account for missing 2 channels in CSP analysis
% else
%     groups = (groups+1)/2;    % to account for half channels we have
% end

% load input
cd(['/einstein0/USRlab/projects/scholvinckm/data/' monkey '/inside scanner/' date '/Matlab']);
load([date2 sess '_BLP_new']);
mkdir(['regressors' sess]);
cd(['/einstein0/USRlab/projects/scholvinckm/data/' monkey '/inside scanner/' date '/Matlab/regressors' sess]);

numchan = size(BLP.blp_dat,1);
numbw = size(BLP.blp_dat,4);

regressors  = squeeze(mean(BLP.blp_dat,2));
av_regr = squeeze(mean(regressors,2));

% % deal with multishot EPI: take mean BLP/CSP of n subsequent epochs (shots)
% regressNew=zeros(numchan,size(regressors,2)/shot,numbw);
% for r = 1:size(regressors,2)/shot
%     regressNew(:,r,:) = mean(regressors(:,(shot*r-(shot-1)):shot*r,:),2);
% end
% regressors = regressNew;
% % decrease 'remove' vector in multishot EPI
% removeNew=zeros(1,length(remove)/shot);
% for r = 1:length(remove)/shot
%     removeNew(r) = (remove(shot*r))/shot;
% end
% remove = removeNew;

% % 'covariation' of alpha and gamma
% covar = zeros(numchan,size(regressors,2));
% for chan = 1:numchan
%     max_alpha = max(regressors(chan,:,3),[],2);        % max alpha power per channel
%     norm_alpha = regressors(chan,:,3)/max_alpha;       % normalised alpha power per channel
%     max_gamma = max(regressors(chan,:,5),[],2);        % max gamma
%     norm_gamma = regressors(chan,:,5)/max_gamma;       % normalised gamma power
%     covar(chan,:) = 1 - (abs(norm_alpha-norm_gamma));  % 1 - abs diff btwn normalised alpha and gamma
%     tot_length = length(covar(chan,:))+length(remove); % add -1000's
%     regr = -1000*ones(1,tot_length);
%     regr(find(~ismember(1:tot_length,remove))) = covar(chan,:);
%     eval(['save ' code '_ch' num2str(2*chan-1) '_covary.rgr regr -ascii']);
% end

% make regressors per channel
if indiv
    for chan = 1:numchan
        for bw = 1:numbw
            regr1 = regressors(chan,:,bw);
            tot_length = length(regr1)+length(remove);
            regr = -1000*ones(1,tot_length);
            regr(find(~ismember(1:tot_length,remove))) = regr1;
            if pad ~= 0, regr(end+1:end+pad) = -1000; end
            if csp
                eval(['save ' code '_ch' num2str(chan) 'cspfilter_' bw_names{bw} '.rgr regr -ascii']);
                                        %num2str(2*chan+1) 
            else
                eval(['save ' code '_ch' num2str(2*chan) '_' bw_names{bw} '.rgr regr -ascii']);
                                        %num2str(2*chan-1)
            end
        end
    end
end

% % make regressors per group of channels
% for gr = 1:size(groups,1)
%     for bw = 1:numbw
%         for chan = groups(gr,1):groups(gr,2)
%             regr1 = regressors(chan,:,bw);
%             tot_length = length(regr1)+length(remove);
%             regr = -1000*ones(1,tot_length);
%             regr(find(~ismember(1:tot_length,remove))) = regr1;
%             gr_regr(chan,:) = regr;
%         end
%         gr_regr = squeeze(mean(gr_regr,1));
%         if csp
%             eval(['save ' code '_' gr_names{gr} 'cspfilter_' bw_names{bw} '.rgr regr -ascii']);
%         else
%             eval(['save ' code '_' gr_names{gr} '_' bw_names{bw} '.rgr regr -ascii']);
%         end
%         
%     end
% end  

% % plot average BLP/CSP per channel for each bandwidth
% for bw = 1:numbw
%     figure; bar(av_regr(:,bw));
%     title([code,' ',bw_names{bw}],'Interpreter','none'); xlabel('channel'); ylabel('average BLP');
%     set(gca,'XTick',1:16); set(gca,'XTicklabel',[1:2:31]);
%     saveas(gcf,[code,'_',bw_names{bw},'.fig']);
% end

copyfile(pwd, '/einstein0/data/lfpmri/rgr'); % copy regressors to blockana directory