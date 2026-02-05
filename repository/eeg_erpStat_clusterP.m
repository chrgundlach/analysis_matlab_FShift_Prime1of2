function [varargout] = eeg_erpStat_clusterP(data1, data2, NPermutation,  direction)
%EEG_ERPSTAT_CLUSTERP Summary of this function goes here
%   Based on originally proposed for MRI data by Bullmore et al.
%   (1999) and for EEG/MEG analysis by Maris & Oostenveld (2007).
%   clustercorrection of running t-test
%   - finds clusters and calculates t_sum for each cluster
%   - creates empiric disribution by shuffling experimental condition and
%   extracting largest t_sum value
%
%
%   input:
%       - data1         datavector (points x measurements)
%       - data2         datavector (points x measurements)
%       - NPermutation  number of permutations (default: 1000)
%       - direction     [-1, 1, 2] [data1<data2, data1>data2, data1<data2 | data1>data2]
%
%   output:
%       - cluster: contains samples, summed t-values, cluster significance
%       - timecourses: contains data, diff, t_vals, p_vals, ci
%
%   example:
%       - [cluster1, timecourse1]=eeg_erpStat_clusterP(data1,data2,10000,2);
%   testdata:
%       data1 = randn(400,20); data1(41:100,:)=randn(60,20)+bsxfun(@times, repmat(tukeywin(60,0.75),1,20),((rand(20,1)*2)+0.5)');
%       data2 = randn(400,20); data2(231:380,:)=randn(150,20)+bsxfun(@times, repmat(tukeywin(150,0.95),1,20),((rand(20,1)*1)+0.5)');
%       data1 = randn(400,25); data1(41:100,:)=randn(60,25)+bsxfun(@times, repmat(tukeywin(60,0.75),1,25),((rand(25,1)*2)+0.5)');
%       data2 = randn(400,25); data2(231:380,:)=randn(150,25)+bsxfun(@times, repmat(tukeywin(150,0.95),1,25),((rand(25,1)*1)+0.5)');
%       figure; plot(mean(data1,2)); hold on; plot(mean(data2,2));

% (c) 2015,2018,2022,2024 - C.Gundlach, S Wehle
%%
if nargin<3
    NPermutation = 1000;
end

%% check for data consistency
if any(size(data1)~=size(data2))
    fprintf(1,'\n###\ndata1 is not of same size das data2\n###\n')
    return
end

%% actual analysis
% 1: running t-tests
switch direction
    case -1
        [tt_h,tt_p,tt_ci,tt_stats] = ttest(data1,data2,0.05,'left',2);
    case 1
        [tt_h,tt_p,tt_ci,tt_stats] = ttest(data1,data2,0.05,'right',2);
    case 2
        [tt_h,tt_p,tt_ci,tt_stats] = ttest(data1,data2,0.05,'both',2);
end

returnflag = 0;
if ~any(tt_h)
    fprintf(1,'\n###\nno significant differences found\n###\n')
    returnflag = 1;
end

% create additionsal timecourses variables
timecourse.data1 = mean(data1,2);
timecourse.data2 = mean(data2,2);
timecourse.diff = mean(data1-data2,2);
timecourse.t_vals = tt_stats.tstat;
timecourse.p_vals_raw = tt_p;
timecourse.h_raw = tt_h;
timecourse.h_corr = zeros(size(tt_h));
timecourse.ci = tt_ci;

% plot data for inspection
% figure; subplot(3,1,1);
% plot(timecourse.data1,'r'); hold on; plot(timecourse.data2,'b')
% if any(timecourse.h_raw)
%     t.yval = min([timecourse.data1(timecourse.h_raw==1); timecourse.data2(timecourse.h_raw==1)])-(diff(get(gca, 'yLim'))*0.1);
%     t.hvals = nan(size(timecourse.h_raw)); t.hvals(find(timecourse.h_raw))=t.yval;
%     plot(t.hvals,'k','LineWidth',3)
% end
% subplot(3,1,2);
% plot(data1,'r'); hold on; plot(data2,'b');
% subplot(3,1,3);
% plot(data1-data2,'Color',[0.6 0.6 0.6]); hold on; plot(timecourse.diff,'r'); grid on



if returnflag == 0;
    % 2: extract cluster with cluster t-values, sorted by cluster size
    cluster = eeg_erpStat_findcluster(tt_p,tt_stats.tstat);
    
    % 3: create empiric distribution for comparison
    % how many unique ways of permutation?
    % randomly assign conditions to each subject
    num_allperm = 2^size(data1,2); % total number of possible combinations
    if num_allperm<NPermutation
        fprintf(1,'\n###\nrunning only %1.0f permutation\n###\n',num_allperm)
        NPermutation = num_allperm;
    end
    
    % create matrix with permutation indices
    if size(data1,2)<21 % due to memory issues
        comb_index = (permn([0 1],size(data1,2)))==1;
        perm_index = randperm(size(comb_index,1),NPermutation);
    else 
        comb_index = false(NPermutation,size(data1,2));
        perm_index = 1:NPermutation;
        for i_perm = 1:NPermutation
            comb_index(i_perm,:) = randi([0 1],1,size(data1,2))==1;
            % check for unique solutions
            while any(~any(bsxfun(@minus,comb_index(1:i_perm-1,:), comb_index(i_perm,:)),2))
                comb_index(i_perm,:) = randi([0 1],1,size(data1,2))==1;
            end
        end
    end
    
    % loop for each permutation to create t-value distribution
    t.absi = round(linspace(0,NPermutation,21));
    t.perci = linspace(0,100,21);
    fprintf(1,'\n###\npermutation running\nprogress in percent:')
    data3=[data1 data2];
    t_distribution_t_sum = []; % t_sum values
    t_distribution_samples = []; % cluster size
    for i_p = 1:NPermutation
        if any( t.absi == i_p)
            fprintf(1,'%4.0f',t.perci(t.absi==i_p))
        end
        idx = [comb_index(perm_index(i_p),:) ~comb_index(perm_index(i_p),:)];
        switch direction
            case -1
                [t.h,t.p,t.ci,t.stats] = ttest(data3(:,idx),data3(:,~idx),0.05,'left',2);
            case 1
                [t.h,t.p,t.ci,t.stats] = ttest(data3(:,idx),data3(:,~idx),0.05,'right',2);
            case 2
                [t.h,t.p,t.ci,t.stats] = ttest(data3(:,idx),data3(:,~idx),0.05,'both',2);
        end
        % figure; plot(mean(data3(:,idx),2)); hold on; plot((mean(data3(:,~idx),2)))
%         % plot data for inspection
%         figure; subplot(3,1,1);
%         plot(timecourse.data1,'r'); hold on; plot(timecourse.data2,'b')
%         plot(mean(data3(:,idx),2),'r:'); plot((mean(data3(:,~idx),2)),'b:')
%         if any(timecourse.h_raw)
%             t.yval = min([timecourse.data1(timecourse.h_raw==1); timecourse.data2(timecourse.h_raw==1)])-(diff(get(gca, 'yLim'))*0.1);
%             t.hvals = nan(size(timecourse.h_raw)); t.hvals(find(timecourse.h_raw))=t.yval;
%             plot(t.hvals,'k','LineWidth',3)
%         end
%         subplot(3,1,2);
%         plot(data3(:,idx),'r'); hold on; plot(data3(:,~idx),'b');
%         subplot(3,1,3);
%         plot(data3(:,idx)-data3(:,~idx),'Color',[0.6 0.6 0.6]); hold on; 
%         plot(data1-data2,'Color',[0.9 0.6 0.6]); 
%         plot(timecourse.diff,'r','LineWidth',2); grid on
%         plot(mean(data3(:,idx)-data3(:,~idx),2),'k','LineWidth',2);
%         hline(0,'c')
        
        t_cluster = eeg_erpStat_findcluster(t.p,t.stats.tstat); % find clusters for each iteration
        
        t.vals = [min(t_cluster.negative.t_sum) max(t_cluster.positive.t_sum)];
        t.clustersize = [max(cellfun(@(x) numel(x),t_cluster.negative.samples)) max(cellfun(@(x) numel(x),t_cluster.positive.samples))];
        
        if ~isempty(t.vals)
            [junk t.idx] = max(abs(t.vals));
            t_distribution_t_sum(i_p)=t.vals(t.idx);
            t_distribution_samples(i_p) = max(t.clustersize);
        else
            t_distribution_t_sum(i_p)=0;
            t_distribution_samples(i_p) =0;
        end
    end
    fprintf(1,'...done\n###\n')
    
    % figure; hist(t_distribution_t_sum,50)
    % figure; plot(sort(t_distribution_t_sum))
    % figure; hist(t_distribution_samples,50)
    % figure; plot(sort(t_distribution_samples))
    
    % 4: compare each cluster summed t-values to t-distribution
    % calculate Monte Carlo sginificance level + determine cluster significance
    if ~isempty(cluster.negative.t_sum)
        for i_c = 1:size(cluster.negative.t_sum,2)
            t.t_distribution = sort([t_distribution_t_sum cluster.negative.t_sum(i_c)]);
            [t.xr t.xc]=find(t.t_distribution>=cluster.negative.t_sum(i_c),1,'first');
            cluster.negative.p(i_c)=t.xc/(NPermutation+1);
            if direction == 2
                cluster.negative.h(i_c)=cluster.negative.p(i_c)<0.025;
            else
                cluster.negative.h(i_c)=cluster.negative.p(i_c)<0.05;
            end
            if cluster.negative.h(i_c)
                timecourse.h_corr(cluster.negative.samples{i_c})=ones;
            end
        end
    else
        cluster.negative.p = [];
        cluster.negative.h = [];
    end
    if ~isempty(cluster.positive.t_sum)
        for i_c = 1:size(cluster.positive.t_sum,2)
            t.t_distribution = sort([t_distribution_t_sum cluster.positive.t_sum(i_c)]);
            [t.xr t.xc]=find(t.t_distribution<=cluster.positive.t_sum(i_c),1,'last');
            cluster.positive.p(i_c)=(NPermutation+1-t.xc)/(NPermutation+1);
            if direction == 2
                cluster.positive.h(i_c)=cluster.positive.p(i_c)<0.025;
            else
                cluster.positive.h(i_c)=cluster.positive.p(i_c)<0.05;
            end
            if cluster.positive.h(i_c)
                timecourse.h_corr(cluster.positive.samples{i_c})=ones;
            end
        end
    else
        cluster.positive.p = [];
        cluster.positive.h = [];
    end
else
    cluster.positive = [];
    cluster.negative = [];
end

% plot final outcome
% figure;
% plot(timecourse.data1,'r'); hold on; plot(timecourse.data2,'b')
% if any(timecourse.h_raw)
%     t.yval = min([timecourse.data1(timecourse.h_raw==1); timecourse.data2(timecourse.h_raw==1)])-(diff(get(gca, 'yLim'))*0.1);
%     t.hvals = nan(size(timecourse.h_raw)); t.hvals(find(timecourse.h_raw))=t.yval;
%     t.hvals_corr = nan(size(timecourse.h_corr)); t.hvals_corr(find(timecourse.h_corr))=t.yval;
%     plot(t.hvals,'k','LineWidth',3)
%     plot(t.hvals_corr,'r','LineWidth',3)
% end

%% output
varargout{1}=cluster;
switch nargout
    case 1
    case 2
        varargout{2}=timecourse;
end
end

