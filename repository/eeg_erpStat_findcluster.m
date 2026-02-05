function [ cluster ] = eeg_erpStat_findcluster(p_values, t_values, criterion)
%find cluster of temporally adjacent p-values for subsequent permutation
%testing
%   Based on originally proposed for MRI data by Bullmore et al.
%   (1999) and for EEG/MEG analysis by Maris & Oostenveld (2007).
%   input:
%       - p_values      vector of p-values (from running t_test)
%       - t_values      vector of related t-values
%       - criterion     criterion for cluster (default: 0.05)


% (c) 2015 - C.Gundlach

%%
if nargin < 3
    criterion = 0.05;
end

%% indices
p_ind = p_values<criterion;
t_pos_ind = t_values>0;
t_neg_ind = t_values<0;
p_bord = diff(p_ind);

%% create samples
% index cluster of positive and negative t-values
pos_cluster_ind = bwlabel(p_ind&t_pos_ind);
neg_cluster_ind = bwlabel(p_ind&t_neg_ind);

cluster.positive.samples = {};
cluster.negative.samples = {};
cluster.positive.t_sum = [];
cluster.negative.t_sum = [];
% find samples and sum t-values
for i_cl1 = 1:max(unique(pos_cluster_ind))
    cluster.positive.samples{i_cl1} = find(pos_cluster_ind==i_cl1);
    cluster.positive.t_sum(i_cl1) = sum(t_values(pos_cluster_ind==i_cl1));
end
for i_cl2 = 1:max(unique(neg_cluster_ind))
    cluster.negative.samples{i_cl2} = find(neg_cluster_ind==i_cl2);
    cluster.negative.t_sum(i_cl2) = sum(t_values(neg_cluster_ind==i_cl2));
end

% sort values
if  ~isempty(cluster.positive.t_sum)
    [t.t t.ind]=sort(cluster.positive.t_sum,2,'descend');
    cluster.positive.t_sum=cluster.positive.t_sum(t.ind);
    cluster.positive.samples=cluster.positive.samples(t.ind);
end
if  ~isempty(cluster.negative.t_sum)
    [t.t t.ind]=sort(cluster.negative.t_sum);
    cluster.negative.t_sum=cluster.negative.t_sum(t.ind);
    cluster.negative.samples=cluster.negative.samples(t.ind);
end


end