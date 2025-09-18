function [varargout] = SSVEP_FShiftPrime1of2_GaborStatOnsetCuSum(data, baselinedata, varargin)
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
    % extract mean and sd for same frequency SSVEP across all participants
    t.m = mean(t.dat);
    t.sd = std(t.dat);
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
% calculate cumulative sum
zdatacusum = cumsum(zdata,1);


%% now same as above but for the bootsamples + null distribution for condition shuffled data
nbootsamples = 5000;
% create index of participants to be drawn
bootsampleidx = reshape(randsample(size(data,4),size(data,4)*nbootsamples,true),size(data,4),[]);
% loop across bootsamples
[bs.zdatacusum bs.zdata bs.zdatacusum_shuffle bs.zdata_shuffle]=deal(nan([size(data) nbootsamples]));

% h_wait = waitbar(0, '...extracting bootatrap samples...');  
for i_bs = 1:nbootsamples
    % t.progress = (i_bs / nbootsamples);
    % waitbar( t.progress, h_wait, sprintf('...extracting bootatrap samples...Progress: %d%%', round( t.progress * 100)));
    
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

    % create empirical null distribution by shuffling condition labels across participants
    data_shuffle = nan(size(data));
    for i_sub = 1:size(data_shuffle,4)
        data_shuffle(:,:,:,i_sub)=data(:,randperm(size(data_shuffle,2)),:,i_sub);
    end
    % ztransform the empirical null data
    bs.zdata_shuffle = (data_shuffle(:,:,:,bootsampleidx(:,i_bs))-bs.M_base_exp)./bs.SD_base_exp;
    bs.zdatacusum_shuffle(:,:,:,:,i_bs) = cumsum(bs.zdata_shuffle,1);
end
close(h_wait);

%% plot data 
pl.data = nan([size(zdatacusum,[1 3]) numel(cons2disp)]);
for i_rdk = 1:size(zdatacusum,3)
    for i_con=1:numel(cons2disp)
        % pl.data(:,i_rdk,i_con)=sum(zdatacusum(:,conRDKidx(:,i_rdk)==i_con,i_rdk,:),[2 3 4]);
        pl.data(:,i_rdk,i_con)=mean(zdatacusum(:,conRDKidx(:,i_rdk)==i_con,i_rdk,:),[2 3 4]);
    end
end
% pl.data = squeeze(sum(pl.data,2));
pl.data = squeeze(mean(pl.data,2));
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

figure('Position',[100 100 500 500]);
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
% do the same for the bs data
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

figure('Position',[100 100 500 500]);
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
% plot(timevec,pl.data)
title('cusum z-values | condition shuffled data')
xlabel('time in ms')
ylabel('standardized cumulative sum')
legend([h.plsem{:}],cons2disp)
hline([-1 1]*pl.critical)

%% plot distribution of onsets | based on own cumsum data | Bayes directly contrasting both distributions [might be wrong]
bs.onsets = nan(numel(cons2disp),nbootsamples);
for i_con = 1:size(pl.bsdata,2)
    for i_bs = 1:size(pl.bsdata,3)
        % index first |value| >pl.critical
        t.idx = find(abs(pl.bsdata(:,i_con,i_bs))>pl.critical,1,"first");
        % corresponds to timevec?
        bs.onsets(i_con,i_bs) = timevec(t.idx);
    end
end

bs.onsets_z = (bs.onsets(1,:)-mean(bs.onsets,"all"))/std(bs.onsets,1,"all");
bs.onsets_z(2,:) = (bs.onsets(2,:)-mean(bs.onsets,"all"))/std(bs.onsets,1,"all");

p_H1 = mean(bs.onsets_z(1,:) < bs.onsets_z(2,:));   % posterior probability RSVP < RDK paired?
p_H0 = 1 - p_H1;
posterior_odds = p_H1 / p_H0;
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
title(sprintf('distributions of estimated onsets for |cusum Z| > %1.2f\nBF10 = %1.3f', ...
    pl.critical,BF10))
legend(cons2disp,'Location','SouthOutside','Orientation','horizontal')
legend('boxoff')

% display zscored data
figure('Position',[100 100 500 500]);
for i_con = 1:size(bs.onsets_z,1)
    t.edges = linspace( ...
        min(bs.onsets_z,[],"all")-abs(diff([min(bs.onsets_z,[],"all") max(bs.onsets_z,[],"all")]))*0.1, ...
        max(bs.onsets_z,[],"all")+abs(diff([min(bs.onsets_z,[],"all") max(bs.onsets_z,[],"all")]))*0.1, ...
        100);
    histogram(bs.onsets_z(i_con,:),t.edges,'FaceColor',pl.concols{i_con})
    hold on
end
title(sprintf('distributions of estimated onsets for | joined z-scored |cusum Z| > %1.2f\nBF10 = %1.3f', ...
    pl.critical,BF10))
legend(cons2disp,'Location','SouthOutside','Orientation','horizontal')
legend('boxoff')


figure('Position',[100 100 500 500]);
histogram(bs.onsets(1,:)-bs.onsets(2,:),100,"FaceColor",[0.8 0.1 0.1])
title(sprintf('onset diffs %s-%s for |cusum Z| > %1.0f\nP_M_o_n_t_e_C_a_r_l_o = %1.4f', ...
    cons2disp{:},pl.critical,1-(sum((bs.onsets(1,:)-bs.onsets(2,:))<0)/(nbootsamples+1))))

%% plot distribution of onsets | based on own cumsum data | normalized by one condition
bs.onsets = nan(numel(cons2disp),nbootsamples);
% bs.onsets(1,:) = effect; bs.onsets(2,:) = null;
for i_con = 1:size(pl.bsdata,2)
    for i_bs = 1:size(pl.bsdata,3)
        % index first |value| >pl.critical
        t.idx = find(abs(pl.bsdata(:,i_con,i_bs))>pl.critical,1,"first");
        % corresponds to timevec?
        bs.onsets(i_con,i_bs) = timevec(t.idx);
    end
end

% joint z-standardization
% bs.onsets_z = (bs.onsets(1,:)-mean(bs.onsets,"all"))/std(bs.onsets,1,"all");
% bs.onsets_z(2,:) = (bs.onsets(2,:)-mean(bs.onsets,"all"))/std(bs.onsets,1,"all");
% % z-standardization by null distribution only
bs.onsets_z = (bs.onsets(1,:)-mean(bs.onsets(2,:),"all"))/std(bs.onsets(2,:),1,"all");
bs.onsets_z(2,:) = (bs.onsets(2,:)-mean(bs.onsets(2,:),"all"))/std(bs.onsets(2,:),1,"all");

posteriorsignedlikelyhood_effect = mean(bs.onsets_z(1,:) < 0);
signedlikelyhood_null = mean(bs.onsets_z(2,:) < 0);

odds_posterior = posteriorsignedlikelyhood_effect ./(1-posteriorsignedlikelyhood_effect);
odds_prior = signedlikelyhood_null ./(1-signedlikelyhood_null);

BF10 = odds_posterior / odds_prior;
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
title(sprintf('distributions of estimated onsets for |cusum Z| > %1.2f\nBF10 = %1.3f', ...
    pl.critical,BF10))
legend(cons2disp,'Location','SouthOutside','Orientation','horizontal')
legend('boxoff')

% display zscored data
figure('Position',[100 100 500 500]);
for i_con = 1:size(bs.onsets_z,1)
    t.edges = linspace( ...
        min(bs.onsets_z,[],"all")-abs(diff([min(bs.onsets_z,[],"all") max(bs.onsets_z,[],"all")]))*0.1, ...
        max(bs.onsets_z,[],"all")+abs(diff([min(bs.onsets_z,[],"all") max(bs.onsets_z,[],"all")]))*0.1, ...
        100);
    histogram(bs.onsets_z(i_con,:),t.edges,'FaceColor',pl.concols{i_con})
    hold on
end
title(sprintf('distributions of estimated onsets for | joined z-scored |cusum Z| > %1.2f\nBF10 = %1.3f', ...
    pl.critical,BF10))
legend(cons2disp,'Location','SouthOutside','Orientation','horizontal')
legend('boxoff')


figure('Position',[100 100 500 500]);
histogram(bs.onsets(1,:)-bs.onsets(2,:),100,"FaceColor",[0.8 0.1 0.1])
title(sprintf('onset diffs %s-%s for |cusum Z| > %1.0f\nP_M_o_n_t_e_C_a_r_l_o = %1.4f', ...
    cons2disp{:},pl.critical,1-(sum((bs.onsets(1,:)-bs.onsets(2,:))<0)/(nbootsamples+1))))



%% plot distribution of onsets with fitted normal distribution according to Andreas/Sebastian | based on own cumsum data 
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

[BFattFASTERunatt] = log10(bootstrap2BF_z(bs.onsets(1,:),bs.onsets(2,:), 1))
[BFunattFASTERatt] = log10(bootstrap2BF_z(bs.onsets(2,:),bs.onsets(1,:), 1))


bootstrap2BF_mc(bs.onsets(1,:),bs.onsets(2,:), 1, 'test_type', 'directional_positive')
bootstrap2BF_mc(bs.onsets(1,:),bs.onsets(2,:), 1, 'test_type', 'directional_negative')


% to avoid issues with scaling, do a joint z-normalization of both
% bootstrapped distributions together, thus maintaining any differences in
% mean and spread
bs.onsets_znormjoined = zscore([bs.onsets(1,:)'; bs.onsets(2,:)']);
bs.onsets_znorm = reshape(bs.onsets_znormjoined,[],2)';


% fit the bootstrapped distributions as normal distributions
bs.PD1 = fitdist(bs.onsets_znorm(1,:)','normal');
bs.PD2 = fitdist(bs.onsets_znorm(2,:)','normal');

% apply the normal distributions
x_values = -5:.01:5;
bs.y1 = pdf(bs.PD1,x_values);
bs.y2 = pdf(bs.PD2,x_values);

bs.fit_posteriorsignedlikelyhood_effect = sum(bs.y1(x_values > 0))./100;
bs.fit_signedlikelyhood_null = sum(bs.y2(x_values > 0))./100;

bs.fit_odds_posterior = bs.fit_posteriorsignedlikelyhood_effect ./(1-bs.fit_posteriorsignedlikelyhood_effect);
bs.fit_odds_prior = bs.fit_signedlikelyhood_null ./(1-bs.fit_signedlikelyhood_null);
    

bs.fit_BF = log10(bs.fit_odds_posterior/bs.fit_odds_prior);


figure;
subplot(2,1,1), histogram(bs.onsets_znorm(1,:)', 100, 'Normalization','pdf')
title('empirical')
hold on
plot(x_values, bs.y1, 'LineWidth',3)
xline(0)
xlim([-6 6])

subplot(2,1,2), histogram(bs.onsets_znorm(2,:)', 100, 'Normalization','pdf')
title('random')
hold on
plot(x_values, bs.y2, 'LineWidth',3)
xline(0)
xlim([-6 6])

%% do everything from scratch (according to gemini and chatgpt)
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

% display diff distributions
pl.diffdists = [diff(bs.onsets,1,1);  diff( bs.onsets_earlyshuffle,1,1);  diff( bs.onsets_lateshuffle,1,1)];

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
legend(pl.difflabels,'Location','SouthOutside','Orientation','horizontal')
legend('boxoff')

% calculate posteriors
prob_H1_effect_model = sum(pl.diffdists(1,:) > 0) / length(pl.diffdists); % probability of difference larger > 0
prob_H1_null_model_earlyshuffle   = sum(pl.diffdists(2,:) > 0) / length(pl.diffdists);
prob_H1_null_model_lateshuffle   = sum(pl.diffdists(3,:) > 0) / length(pl.diffdists);

odds_posterior = prob_H1_effect_model / (1 - prob_H1_effect_model);
odds_prior_earlyshuffle = prob_H1_null_model_earlyshuffle / (1 - prob_H1_null_model_earlyshuffle);
odds_prior_lateshuffle = prob_H1_null_model_lateshuffle / (1 - prob_H1_null_model_lateshuffle);

% BF10 = odds_posterior / odds_prior_earlyshuffle;
BF10 = odds_posterior / odds_prior_lateshuffle;



%% jackknife approach
jckknifeidx = true(size(data,4));
jckknifeidx(diag(diag(jckknifeidx)))=false;
% loop across jackknifesamples
jk.zdata=nan([size(data)+[0 0 0 -1] size(jckknifeidx,2)]);
jk.zdatacusum=jk.zdata;
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
