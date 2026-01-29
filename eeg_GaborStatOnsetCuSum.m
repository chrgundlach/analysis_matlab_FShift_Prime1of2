function [varargout] = eeg_GaborStatOnsetCuSum(data, baselinedata, varargin)
%EEG_GABORSTATONSET Function to estimate the onset points of Gabor timecourses based on CuSum method
%   function to test all the Bayes Factor approaches around
%
%   input:
%       - data in 2D (one condition: time X participants) or 3D (condition X time X participants)
%       - baseline either 2D or 3D as above

%% defaultparameters
datadims = {'time','cons','RDK','subs'}; % default structure of data
pl.critical = 2.5;
pl.critical = 1.96;
pl.concols = num2cell(hex2rgb(["#F03F3F","#5CACEE"])',1);
% pl.concols = num2cell([255 133 4; 41 60 74; 25 138 131]'./255,1);


%% check for input
if ~(round(numel(varargin)/2) == numel(varargin)/2)
    error('Odd number of input arguments??')
end

for i = 1:2:length(varargin)
    Param = varargin{i};
    Value = varargin{i+1};
    if ~isstr(Param)
        error('Flag arguments must be strings')
    end
    Param = lower(Param);
    switch Param
        case 'datadims' % check for data structure and alter if difference
            if ~all(strcmp(lower(Value),lower(datadims)))
                t.permidx = cellfun(@(x) find(strcmp(lower(Value),x)),lower(datadims));
                data = permute(data,t.permidx);
                baselinedata = permute(baselinedata,t.permidx);
            end
        case 'conlabels' % define time window
            if iscellstr(Value)
                cons2disp = unique(Value);
                conRDKidx = nan(size(Value));
                for i_con = 1:numel(cons2disp)
                    conRDKidx(strcmp(Value,cons2disp(i_con)))=i_con;
                end
            end
        case 'rdkfreqs' %
            rdkfreqs = Value;
            rdkfreqsunq = unique(Value);
            rdkfreqidx = nan(size(Value));
            for i_freq = 1:numel(rdkfreqsunq)
                rdkfreqidx(Value==rdkfreqsunq(i_freq))=i_freq;
            end
        case 'timevec'
            timevec = Value;
    end
end

%% first z-transform the real data
% separately for the RDKs of the same frequency
% extract mean and sd for the baseline first which serve as conceptual null distribution
% having more than one SSVEP with different frequencies, this means that SSVEPs may differ per se due to 1/f distribution of
% amplitudes, hence M and SD need to be extracted for each RDK frequency
% preallocate memory
M_base = nan(size(data,[3 4]));
SD_base = M_base;
% loop across RDK frequencies
for i_freq = 1:numel(rdkfreqs)
    t.dat = [];
    % as frequencies for RDKs may differ across participants, do this seperately for each participant
    for i_sub = 1:size(data,4)
        % concatenate the relevant amplitudes
        t.dat = cat(2,t.dat,reshape(baselinedata(:,:,rdkfreqidx(i_sub,:)==i_freq,i_sub),1,[]));
    end
    % figure; histogram( t.dat,linspace(0,10,100)); title('baseline amplitudes'); xlabel('amplitude in \muV/m²')
    % figure; histogram( log(t.dat),linspace(0,log(10),100)); title('baseline amplitudes'); xlabel('amplitude in \muV/m²')

    % extract mean and sd for same frequency SSVEP across all participants
    t.m = mean(t.dat);
    t.sd = std(t.dat);
    % % extract mean and sd for same frequency SSVEP across all participants | log transformed data
    % t.m = mean(log(t.dat));
    % t.sd = std(log(t.dat));
    % now write everything back to have the original structure required for the subsequent calculations
    for i_sub = 1:size(data,4)
        M_base(rdkfreqidx(i_sub,:)==i_freq,i_sub)=t.m;
        SD_base(rdkfreqidx(i_sub,:)==i_freq,i_sub) = t.sd;
    end
end

% now z-transform the actual time resolved data
% expand mean and sd data to have the size of the post-cue data to allow for element by element calculations
M_base_exp = permute(repmat(M_base,[1 1 size(data,(1:2))]),[3,4,1,2]);
SD_base_exp = permute(repmat(SD_base,[1 1 size(data,(1:2))]),[3,4,1,2]);

% ztransform the data
zdata = (data-M_base_exp)./SD_base_exp;
% ztransform the data | log transform
% zdata = (log(data)-M_base_exp)./SD_base_exp;
% calculate cumulative sum
zdatacusum = cumsum(zdata,1);


%% now same as above but for the bootsamples + null distribution for condition shuffled data
nbootsamples = 5000;
% create index of participants to be drawn
bootsampleidx = reshape(randsample(size(data,4),size(data,4)*nbootsamples,true),size(data,4),[]);
% loop across bootsamples
[bs.zdatacusum bs.zdata bs.zdatacusum_shuffle bs.zdata_shuffle]=deal(nan([size(data) nbootsamples]));

h_wait = waitbar(0, '...extracting bootatrap samples...');  
for i_bs = 1:nbootsamples
    t.progress = (i_bs / nbootsamples);
    waitbar( t.progress, h_wait, sprintf('...extracting bootatrap samples...Progress: %d%%', round( t.progress * 100)));
    
    % pre-allocated memory
    bs.M_base = nan(size(data,[3 4]));
    bs.SD_base = bs.M_base;
    for i_freq = 1:size(rdkfreqs,2)
        t.dat = [];
        for i_sub = 1:size(data,4)
            t.dat = cat(2,t.dat,reshape(baselinedata(:,:,rdkfreqidx(bootsampleidx(i_sub,i_bs),:)==i_freq,bootsampleidx(i_sub,i_bs)),1,[]));
        end
        t.m = mean(t.dat);
        t.sd = std(t.dat);
        
        % % log transform
        % t.m = mean(log(t.dat));
        % t.sd = std(log(t.dat));
        % now write everything back
        for i_sub = 1:size(data,4)
            bs.M_base(rdkfreqidx(bootsampleidx(i_sub,i_bs),:)==i_freq,i_sub)=t.m;
            bs.SD_base(rdkfreqidx(bootsampleidx(i_sub,i_bs),:)==i_freq,i_sub) = t.sd;
        end
    end
    % now z-transform the actual time resolved data
    % expand mean and sd data
    bs.M_base_exp = permute(repmat(bs.M_base,[1 1 size(data,(1:2))]),[3,4,1,2]);
    bs.SD_base_exp = permute(repmat(bs.SD_base,[1 1 size(data,(1:2))]),[3,4,1,2]);

    % ztransform the data
    bs.zdata = (data(:,:,:,bootsampleidx(:,i_bs))-bs.M_base_exp)./bs.SD_base_exp;
    bs.zdatacusum(:,:,:,:,i_bs) = cumsum(bs.zdata,1);
    % % ztransform the data | log
    % bs.zdata = (log(data(:,:,:,bootsampleidx(:,i_bs)))-bs.M_base_exp)./bs.SD_base_exp;
    % bs.zdatacusum(:,:,:,:,i_bs) = cumsum(bs.zdata,1);

    % create empirical null distribution by shuffling condition labels across participants
    data_shuffle = nan(size(data));
    for i_sub = 1:size(data_shuffle,4)
        data_shuffle(:,:,:,i_sub)=data(:,randperm(size(data_shuffle,2)),:,i_sub);
    end
    % ztransform the empirical null data
    bs.zdata_shuffle = (data_shuffle(:,:,:,bootsampleidx(:,i_bs))-bs.M_base_exp)./bs.SD_base_exp;
    bs.zdatacusum_shuffle(:,:,:,:,i_bs) = cumsum(bs.zdata_shuffle,1);
    % % ztransform the empirical null data | logtransformed data
    % bs.zdata_shuffle = (log(data_shuffle(:,:,:,bootsampleidx(:,i_bs)))-bs.M_base_exp)./bs.SD_base_exp;
    % bs.zdatacusum_shuffle(:,:,:,:,i_bs) = cumsum(bs.zdata_shuffle,1);
end
close(h_wait);

%% plot data | cumsum + CI of bootstrap samples
pl.data = nan([size(zdatacusum,[1 3]) numel(cons2disp)]);
% loop across RDKs and conditions
for i_rdk = 1:size(zdatacusum,3)
    for i_con=1:numel(cons2disp)
        % sum all values? looks strange
        % pl.data(:,i_rdk,i_con)=sum(zdatacusum(:,conRDKidx(:,i_rdk)==i_con,i_rdk,:),[2 3 4]);
        % rather average z-cusum values across condition indices and participants
        pl.data(:,i_rdk,i_con)=mean(zdatacusum(:,conRDKidx(:,i_rdk)==i_con,i_rdk,:),[2 3 4]);
    end
end
% pl.data = squeeze(sum(pl.data,2));
pl.data = squeeze(mean(pl.data,2)); % average across rdks
% do the same for the bs data
pl.bsdata = nan([size(bs.zdatacusum,[1 3]) numel(cons2disp) size(bs.zdatacusum,[5])]);
for i_rdk = 1:size(bs.zdatacusum,3)
    for i_con=1:numel(cons2disp)
        % pl.data(:,i_rdk,i_con)=sum(zdatacusum(:,conRDKidx(:,i_rdk)==i_con,i_rdk,:),[2 3 4]);
        pl.bsdata(:,i_rdk,i_con,:)=mean(bs.zdatacusum(:,conRDKidx(:,i_rdk)==i_con,i_rdk,:,:),[2 3 4]);
    end
end
% pl.bsdata = squeeze(sum(pl.bsdata,2));
pl.bsdata = squeeze(mean(pl.bsdata,2));
% pl.bsdata: [time x condition x bootstrap_samples]
pl.bsci_low  = prctile(pl.bsdata, 2.5, 3);   % lower 2.5% across bootstrap samples
pl.bsci_high = prctile(pl.bsdata, 97.5, 3);  % upper 97.5%

figure('Position',[100 100 500 300]);
h.pl = {}; h.plsem=[];  h.plm = []; h.pls = []; h.plst = [];
for i_con = 1:size(pl.data,2)
    % plot SEM as boundary
    % create data
    pl.xconf = [timevec timevec(end:-1:1)] ;
    pl.yconf = [pl.bsci_high(:,i_con)' pl.bsci_low(end:-1:1,i_con)'];
    % plot
    h.plsem{i_con} = fill(pl.xconf,pl.yconf,pl.concols{i_con}','EdgeColor','none','FaceAlpha',0.3);
    hold on
    % plot mean lines
    h.plm{i_con}=plot(timevec, pl.data(:,i_con),'Color',pl.concols{i_con},'LineWidth',1);
end
plot(timevec,pl.data)
title('cusum z-values | real data')
xlabel('time in ms')
ylabel('standardized cumulative sum')
legend([h.plm{:}],cons2disp)
hline([-1 1]*pl.critical)


% permuted data with shuffled condition labels
% do the same for the bs data with the shuffled condition labels
pl.bsdata_shuffle = nan([size(bs.zdatacusum_shuffle,[1 3]) numel(cons2disp) size(bs.zdatacusum_shuffle,[5])]);
for i_rdk = 1:size(bs.zdatacusum_shuffle,3)
    for i_con=1:numel(cons2disp)
        % pl.data(:,i_rdk,i_con)=sum(zdatacusum(:,conRDKidx(:,i_rdk)==i_con,i_rdk,:),[2 3 4]);
        pl.bsdata_shuffle(:,i_rdk,i_con,:)=mean(bs.zdatacusum_shuffle(:,conRDKidx(:,i_rdk)==i_con,i_rdk,:,:),[2 3 4]);
    end
end
% pl.bsdata = squeeze(sum(pl.bsdata,2));
pl.bsdata_shuffle = squeeze(mean(pl.bsdata_shuffle,2));
% pl.bsdata: [time x condition x bootstrap_samples]
pl.bsci_low_shuffle  = prctile(pl.bsdata_shuffle, 2.5, 3);   % lower 2.5% across bootstrap samples
pl.bsci_high_shuffle = prctile(pl.bsdata_shuffle, 97.5, 3);  % upper 97.5%

figure('Position',[100 100 500 300]);
h.pl = {}; h.plsem=[];  h.plm = []; h.pls = []; h.plst = [];
for i_con = 1:size(pl.data,2)
    % plot SEM as boundary
    % create data
    pl.xconf = [timevec timevec(end:-1:1)] ;
    pl.yconf = [pl.bsci_high_shuffle(:,i_con)' pl.bsci_low_shuffle(end:-1:1,i_con)'];
    % plot
    h.plsem{i_con} = fill(pl.xconf,pl.yconf,pl.concols{i_con}','EdgeColor','none','FaceAlpha',0.3);
    hold on
    % plot mean lines
    % h.plm{i_con}=plot(timevec, pl.data(:,i_con),'Color',pl.concols{i_con},'LineWidth',1);
end
plot(timevec,pl.data)
title('cusum z-values | condition shuffled data (early)')
xlabel('time in ms')
ylabel('standardized cumulative sum')
legend([h.plsem{:}],cons2disp)
hline([-1 1]*pl.critical)

%% now do the BayesFactor benchmarking
% version 1
% do everything from scratch
% first calculate difference in timing between the two conditions "cons2disp"
% for real data, early and late shuffled data serving as conceptual null distributions
%
% this does treat the effects as paired differences

% create real effect data: timing differences between bootstrapped data
[bs.onsets bs.onsets_earlyshuffle bs.onsets_lateshuffle] = deal(nan(numel(cons2disp),nbootsamples));
% for real data
for i_con = 1:size(pl.bsdata,2)
    for i_bs = 1:size(pl.bsdata,3)
        % index first |value| >pl.critical
        t.idx = find(abs(pl.bsdata(:,i_con,i_bs))>pl.critical,1,"first");
        % corresponds to timevec?
        bs.onsets(i_con,i_bs) = timevec(t.idx);
    end
end
% for early shuffled data, based on condition shuffling during bootstrap cumsum calculation
for i_con = 1:size(pl.bsdata,2)
    for i_bs = 1:size(pl.bsdata,3)
        % index first |value| >pl.critical
        t.idx = find(abs(pl.bsdata_shuffle(:,i_con,i_bs))>pl.critical,1,"first");
        % corresponds to timevec?
        if ~isempty(t.idx)
            bs.onsets_earlyshuffle(i_con,i_bs) = timevec(t.idx);
        else
            bs.onsets_earlyshuffle(i_con,i_bs) = nan;
        end
    end
end
% for late shuffled data, shffling conditions randomly acroos bootstrap samples
for i_bs = 1:size(pl.bsdata_shuffle,3)
    % randomly swap condition labels
    t.cidx = randperm(size(pl.bsdata_shuffle,2));
    for i_con = 1:size(pl.bsdata_shuffle,2)
        % index first |value| >pl.critical
        t.idx = find(abs(pl.bsdata_shuffle(:, t.cidx(i_con),i_bs))>pl.critical,1,"first");
        % corresponds to timevec?
        if ~isempty(t.idx)
            bs.onsets_lateshuffle(i_con,i_bs) = timevec(t.idx);
        else
            bs.onsets_lateshuffle(i_con,i_bs) = nan;
        end
    end
end

% concat diff distributions
pl.diffdists = [diff(bs.onsets,1,1);  diff( bs.onsets_earlyshuffle,1,1);  diff( bs.onsets_lateshuffle,1,1)];

% display distributions
figure('Position',[100 100 500 500]);
pl.diffcols = num2cell([255 133 4; 41 60 74; 25 138 131]'./255,1);
pl.difflabels = {'original diff';'shuffled early';'shuffled late'};
for i_diff = 1:size(pl.diffdists,1)
    t.edges = linspace( ...
        min(pl.diffdists,[],"all")-abs(diff([min(pl.diffdists,[],"all") max(pl.diffdists,[],"all")]))*0.1, ...
        max(pl.diffdists,[],"all")+abs(diff([min(pl.diffdists,[],"all") max(pl.diffdists,[],"all")]))*0.1, ...
        100);
    histogram(pl.diffdists(i_diff,:),t.edges,'FaceColor',pl.diffcols{i_diff})
    hold on
end
title(sprintf('distributions of onset differences | %s - %s',cons2disp{end:-1:1}))
legend(pl.difflabels,'Location','SouthOutside','Orientation','horizontal')
legend('boxoff')

% calculate posterior of effect model
prob_H1_effect_model = sum(pl.diffdists(1,:) > 0) / length(pl.diffdists); % probability of difference larger > 0
% calculate posterior of null models
% probability of difference larger > 0
prob_H1_null_model_earlyshuffle   = sum(pl.diffdists(2,:) > 0) / length(pl.diffdists);
prob_H1_null_model_lateshuffle   = sum(pl.diffdists(3,:) > 0) / length(pl.diffdists);

% calculate posterior odds
% probability of difference larger > 0
odds_posterior_pos = prob_H1_effect_model / (1 - prob_H1_effect_model);
% probability of difference smaller < 0
odds_posterior_neg = (1-prob_H1_effect_model) / prob_H1_effect_model;

% caculate odds for the null distribution that serve as priors (no effect should be expected)
% crucially: if there is no effect expected, the prior odds should approach 1, i.e. 0.5/0.5
% probability of difference larger > 0
odds_prior_earlyshuffle_pos = prob_H1_null_model_earlyshuffle / (1 - prob_H1_null_model_earlyshuffle);
odds_prior_lateshuffle_pos = prob_H1_null_model_lateshuffle / (1 - prob_H1_null_model_lateshuffle);
% probability of difference smaller < 0
odds_prior_earlyshuffle_neg = (1-prob_H1_null_model_earlyshuffle) / prob_H1_null_model_earlyshuffle;
odds_prior_lateshuffle_neg = (1-prob_H1_null_model_lateshuffle) / prob_H1_null_model_lateshuffle;

% now calculate all combinations of Bayes Factors
% essentially it is  posterior odds by the prior odds: how much more likely is the actual effect as compared to the expected
BF10_early_pos = odds_posterior_pos / odds_prior_earlyshuffle_pos;
BF10_late_pos = odds_posterior_pos / odds_prior_lateshuffle_pos;
BF10_early_neg = odds_posterior_neg / odds_prior_earlyshuffle_neg; % neg = 1/pos
BF10_late_neg = odds_posterior_neg / odds_prior_lateshuffle_neg;

fprintf('\nearly shuffle | %s > %s | Bayes Factor = %.3f', cons2disp{end:-1:1},BF10_early_pos);
fprintf('\nearly shuffle | %s < %s | Bayes Factor = %.3f', cons2disp{end:-1:1},BF10_early_neg);
fprintf('\nlate  shuffle | %s > %s | Bayes Factor = %.3f', cons2disp{end:-1:1},BF10_late_pos);
fprintf('\nlate  shuffle | %s < %s | Bayes Factor = %.3f\n', cons2disp{end:-1:1},BF10_late_neg);


%% plot distribution of onsets | based on own cumsum data | Bayes directly contrasting both distributions | partly based on chatgpt/gemini [might be wrong]
% version 2
% also treats the test as dependent, pairwise contrast
% tests the hypothesis of a systematic difference in onset times
bs.onsets = nan(numel(cons2disp),nbootsamples);
for i_con = 1:size(pl.bsdata,2)
    for i_bs = 1:size(pl.bsdata,3)
        % index first |value| >pl.critical
        t.idx = find(abs(pl.bsdata(:,i_con,i_bs))>pl.critical,1,"first");
        % corresponds to timevec?
        bs.onsets(i_con,i_bs) = timevec(t.idx);
    end
end

% calculate probability of attended < unattended: across all pairs, how likely is the expected direction to be found?
p_H1 = mean(bs.onsets(1,:) < bs.onsets(2,:));  
% calculate probability of alternative hypothesis attended >= unattended
p_H0 = 1 - p_H1;
% calculate posterior odds: how much more likely is the alternative hypothesis than the null hypothesis
posterior_odds = p_H1 / p_H0;
% what odds would one expect? if non-informative: it should be 1 as p_H1_null ~ 0.5, p_H0_null ~ 0.5 --> prior_odds = 0.5/0.5 = 1
prior_odds = 1;  % flat prior
BF10 = posterior_odds / prior_odds;
fprintf('Posterior P(%s < %s) = %.3f\n',cons2disp{:}, p_H1);
fprintf('Bayes Factor = %.3f\n', BF10);

% display true data
figure('Position',[100 100 500 500]);
for i_con = 1:size(bs.onsets,1)
    t.edges = linspace( ...
        min(bs.onsets,[],"all")-abs(diff([min(bs.onsets,[],"all") max(bs.onsets,[],"all")]))*0.1, ...
        max(bs.onsets,[],"all")+abs(diff([min(bs.onsets,[],"all") max(bs.onsets,[],"all")]))*0.1, ...
        100);
    histogram(bs.onsets(i_con,:),t.edges,'FaceColor',pl.concols{i_con})
    hold on
end
title(sprintf('distributions of estimated onsets for |cusum Z| > %1.2f\n%s < %s BF10 = %1.3f', ...
    pl.critical,cons2disp{:},BF10))
legend(cons2disp,'Location','SouthOutside','Orientation','horizontal')
legend('boxoff')

% display true difference
figure('Position',[100 100 500 500]);
histogram(bs.onsets(2,:)-bs.onsets(1,:),100,'FaceColor',[0.3 0.9 0.6])
title(sprintf('diffs ofonsets for |cusum Z| > %1.2f\n%s MINUS %s BF10 = %1.3f', ...
    pl.critical,cons2disp{[2 1]},BF10))


% display true difference | colored by hypothesis label
diffdata = bs.onsets(2,:)-bs.onsets(1,:);
diffdata_edges = linspace(min(diffdata),max(diffdata),100);
figure('Position',[100 100 500 500]);
histogram(diffdata(diffdata<=0),diffdata_edges)
hold on
histogram(diffdata(diffdata>0),diffdata_edges)
title(sprintf('diffs ofonsets for |cusum Z| > %1.2f\n%s MINUS %s BF10 = %1.3f', ...
    pl.critical,cons2disp{[2 1]},BF10))
legend({'p_H_0';'p_H_1'},'Location','SouthOutside','Orientation','horizontal')
legend('boxoff')



%% plot distribution of onsets with fitted normal distribution according to Andreas/Sebastian/Norman
% first extract onsets of cucsum data
bs.onsets = nan(numel(cons2disp),nbootsamples);
for i_con = 1:size(pl.bsdata,2)
    for i_bs = 1:size(pl.bsdata,3)
        % index first |value| >pl.critical
        t.idx = find(abs(pl.bsdata(:,i_con,i_bs))>pl.critical,1,"first");
        % corresponds to timevec?
        bs.onsets(i_con,i_bs) = timevec(t.idx);
    end
end

% there are conceptually some issues:
% if the data are z standardized jointly:
%   how does this behave if both distributions are qualitatively different?
%   in case of a difference: doesn't this push the prior distribution out?

[BFattFASTERunatt] = bootstrap2BF_z(bs.onsets(1,:),bs.onsets(2,:), 1)
[BFunattFASTERatt] = bootstrap2BF_z(bs.onsets(2,:),bs.onsets(1,:), 1)


[BFattFASTERunatt2]= bootstrap2BF_mc(bs.onsets(1,:),bs.onsets(2,:), 1, 'test_type', 'directional_positive')
[BFunattFASTERatt2] =bootstrap2BF_mc(bs.onsets(1,:),bs.onsets(2,:), 1, 'test_type', 'directional_negative')


%% plot distribution of onsets with fitted normal distribution according to Andreas/Sebastian/Norman | adapted
% adapted to meet the issues above:
% don't do joint transformation, but use null distribution as reference for z-transform?
% now the problem is that both distributions may look different testing in both directions may yield totally different
% findings
%
% first extract onsets of cucsum data
bs.onsets = nan(numel(cons2disp),nbootsamples);
for i_con = 1:size(pl.bsdata,2)
    for i_bs = 1:size(pl.bsdata,3)
        % index first |value| >pl.critical
        t.idx = find(abs(pl.bsdata(:,i_con,i_bs))>pl.critical,1,"first");
        % corresponds to timevec?
        bs.onsets(i_con,i_bs) = timevec(t.idx);
    end
end

[BFattFASTERunatt] = bootstrap2BF_z_adapted(bs.onsets(1,:),bs.onsets(2,:), 1)
[BFunattFASTERatt] = bootstrap2BF_z_adapted(bs.onsets(2,:),bs.onsets(1,:), 1)


%% alternative approach
% simply calculate monte carlo derived p-value
bs.onsets = nan(numel(cons2disp),nbootsamples);
for i_con = 1:size(pl.bsdata,2)
    for i_bs = 1:size(pl.bsdata,3)
        % index first |value| >pl.critical
        t.idx = find(abs(pl.bsdata(:,i_con,i_bs))>pl.critical,1,"first");
        % corresponds to timevec?
        bs.onsets(i_con,i_bs) = timevec(t.idx);
    end
end

% monte carlo-based p-values
bs.diff_p_mc_pos = sum((bs.onsets(1,:)-bs.onsets(2,:))<0)/(nbootsamples+1);
bs.diff_p_mc_neg = sum((bs.onsets(2,:)-bs.onsets(1,:))<0)/(nbootsamples+1);
bs.diff_p_mc = min([bs.diff_p_mc_pos bs.diff_p_mc_neg])*2;
bs.diff_CI = prctile(bs.onsets(1,:)-bs.onsets(2,:), [2.5 97.5]);

figure('Position',[100 100 500 500]);
histogram(bs.onsets(1,:)-bs.onsets(2,:),100,"FaceColor",[0.8 0.1 0.1])
title(sprintf('onset diffs %s-%s for |cusum Z| > %1.0f\nP_M_o_n_t_e_C_a_r_l_o = %1.4f | CI = [%1.3f %1.3f]ms', ...
    cons2disp{:},pl.critical,bs.diff_p_mc, bs.diff_CI))


%% jackknife approach | part 1: extracting jackknife samples
% use jackknifing instead of
jckknifeidx = true(size(data,4));
jckknifeidx(diag(diag(jckknifeidx)))=false;
% loop across jackknifesamples
jk.zdata=nan([size(data)+[0 0 0 -1] size(jckknifeidx,2)]);
for i_jk = 1:size(jckknifeidx,2)
    jk.M_base = nan(numel(cons2disp),size(jckknifeidx,1)-1);
    jk.SD_base = jk.M_base;
    t.jkidx = find(jckknifeidx(:,i_jk));
    for i_freq = 1:size(rdkfreqs,2)
        t.dat = [];
        for i_sub = 1:numel(t.jkidx)
            t.dat = cat(2,t.dat,reshape(baselinedata(:,:,rdkfreqidx(t.jkidx(i_sub),:)==i_freq,t.jkidx(i_sub)),1,[]));
        end
        t.m = mean(t.dat);
        t.sd = std(t.dat);
        % now write everything back
        for i_sub = 1:numel(t.jkidx)
            jk.M_base(rdkfreqidx(t.jkidx(i_sub),:)==i_freq,i_sub)=t.m;
            jk.SD_base(rdkfreqidx(t.jkidx(i_sub),:)==i_freq,i_sub) = t.sd;
        end
    end
    % now z-transform the actual time resolved data
    % expand mean and sd data
    jk.M_base_exp = permute(repmat(jk.M_base,[1 1 size(data,(1:2))]),[3,4,1,2]);
    jk.SD_base_exp = permute(repmat(jk.SD_base,[1 1 size(data,(1:2))]),[3,4,1,2]);

    % ztransform the data
    jk.zdata = (data(:,:,:,t.jkidx)-jk.M_base_exp)./jk.SD_base_exp;
    jk.zdatacusum(:,:,:,:,i_jk) = cumsum(jk.zdata,1);

end
% bsdata
%% jackknife approach | part 2: plot data | cumsum + CI of bootstrap samples
pl.data = nan([size(zdatacusum,[1 3]) numel(cons2disp)]);
% loop across RDKs and conditions
for i_rdk = 1:size(zdatacusum,3)
    for i_con=1:numel(cons2disp)
        % sum all values? looks strange
        % pl.data(:,i_rdk,i_con)=sum(zdatacusum(:,conRDKidx(:,i_rdk)==i_con,i_rdk,:),[2 3 4]);
        % rather average z-cusum values across condition indices and participants
        pl.data(:,i_rdk,i_con)=mean(zdatacusum(:,conRDKidx(:,i_rdk)==i_con,i_rdk,:),[2 3 4]);
    end
end
% pl.data = squeeze(sum(pl.data,2));
pl.data = squeeze(mean(pl.data,2));
% do the same for the jk data
pl.jkdata = nan([size(jk.zdatacusum,[1 3]) numel(cons2disp) size(jk.zdatacusum,[5])]);
for i_rdk = 1:size(jk.zdatacusum,3)
    for i_con=1:numel(cons2disp)
        % pl.data(:,i_rdk,i_con)=sum(zdatacusum(:,conRDKidx(:,i_rdk)==i_con,i_rdk,:),[2 3 4]);
        pl.jkdata(:,i_rdk,i_con,:)=mean(jk.zdatacusum(:,conRDKidx(:,i_rdk)==i_con,i_rdk,:,:),[2 3 4]);
    end
end
% pl.jkdata = squeeze(sum(pl.jkdata,2));
pl.jkdata = squeeze(mean(pl.jkdata,2));
% pl.jkdata: [time x condition x jackknife_samples]

pl.jkdata_smulder = nan(size(pl.jkdata));
for i_sub = 1:size(pl.jkdata_smulder,3)
    pl.jkdata_smulder(:,:,i_sub) = size(pl.jkdata,3)*mean(pl.jkdata,3)-(size(pl.jkdata,3)-1)*pl.jkdata(:,:,i_sub);
end

pl.jkci_low  = prctile(pl.jkdata_smulder, 2.5, 3);   % lower 2.5% across bootstrap samples
pl.jkci_high = prctile(pl.jkdata_smulder, 97.5, 3);  % upper 97.5%
% pl.jkci_low  = prctile(pl.jkdata, 2.5, 3);   % lower 2.5% across bootstrap samples
% pl.jkci_high = prctile(pl.jkdata, 97.5, 3);  % upper 97.5%



figure('Position',[100 100 500 500]);
h.pl = {}; h.plsem=[];  h.plm = []; h.pls = []; h.plst = [];
for i_con = 1:size(pl.data,2)
    % plot SEM as boundary
    % create data
    pl.xconf = [timevec timevec(end:-1:1)] ;
    pl.yconf = [pl.jkci_high(:,i_con)' pl.jkci_low(end:-1:1,i_con)'];
    % plot
    h.plsem{i_con} = fill(pl.xconf,pl.yconf,pl.concols{i_con}','EdgeColor','none','FaceAlpha',0.3);
    hold on
    % plot mean lines
    h.plm{i_con}=plot(timevec, pl.data(:,i_con),'Color',pl.concols{i_con},'LineWidth',1);
end
plot(timevec,pl.data)
title('cusum z-values | real data')
xlabel('time in ms')
ylabel('standardized cumulative sum')
legend([h.plm{:}],cons2disp)
hline([-1 1]*pl.critical)

%% jackknife approach | part 3:  now extract smuldered onset points
jk.onsets = nan(numel(cons2disp),size(jckknifeidx,1));
for i_con = 1:size(pl.jkdata,2)
    for i_jk = 1:size(pl.jkdata,3)
        % index first |value| >pl.critical
        t.idx = find(abs(pl.jkdata(:,i_con,i_jk))>pl.critical,1,"first");
        % corresponds to timevec?
        jk.onsets(i_con,i_jk) = timevec(t.idx);
    end
end

jk.onsets_smulder = nan(size(jk.onsets));
for i_sub = 1:size(jk.onsets_smulder,2)
    jk.onsets_smulder(:,i_sub) = size(jk.onsets,2)*mean(jk.onsets,2)-(size(jk.onsets,2)-1)*jk.onsets(:,i_sub);
end

figure; histogram(jk.onsets_smulder(1,:),20); hold on; histogram(jk.onsets_smulder(2,:),20)
figure; histogram(diff(jk.onsets_smulder,1,1),100)






