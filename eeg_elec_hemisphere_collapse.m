function [elec_index_out] = eeg_elec_hemisphere_collapse(channames,idx_hemisphere)
%EEG_ELEC_HEMISPHERE_COLLAPSE electrode index to symmetrically collapse
%across hemispheres
%   works for symmetrical EEG electrode setup with standard labeling
%   requires:
%       -channames: cells with channel labels
%       -idx_hemisphere:
%           1: bias left
%           2: corresponding righ index

% 2020 C.Gundlach
%% initial check for input
% check for channel labels
if (sum(cellfun(@ischar,channames)) ~= numel(channames)) &&  (sum(cellfun(@isstr,channames)) ~= numel(channames))
    fprintf(1,'\n###\nInput of channel labels needs to be cell with strings or chars!\n###\n')
    return
end

% check for hemisphere index
if isnumeric(idx_hemisphere)
    if ~any(idx_hemisphere ~= 1 | idx_hemisphere ~= 2)
        fprintf(1,'\n###\nInput of hemisphere needs to be either [1],[2],[1,2]!\n###\n')
    end
else
    fprintf(1,'\n###\nInput of hemisphere needs to be numeric!\n###\n')
    return
end

%% do correspondance matching
% index biased to left hemisphere
elec_index = nan(2,numel(channames));
elec_index(1,1:end) = 1:numel(channames);

% look for corresponding symmetric channels in opposite hemisphere
for i_chan = 1:numel(channames)
    % extract number from channel label
    % is there a number?
    chan_num_idx = regexp(channames{i_chan},'\d');
    if isempty(chan_num_idx) % no
        % no number: could be artifact or central channel
        if strcmpi(channames{i_chan}(end),'z')
            % must be central electrode --> copy index
            elec_index(2,i_chan) = elec_index(1,i_chan);
        else
            % unrecognized channel
            fprintf(1,'\n###\nunrecognized channel label: %s of channel number %1.0f\n###\n', channames{i_chan},i_chan)
            return
        end
    else
        % extract string and number part
        chan_num = str2double(channames{i_chan}(chan_num_idx));
        chan_str = channames{i_chan}(1:chan_num_idx-1);
        % odd number or even number?
        if mod(chan_num,2)==1
            % odd! electrode must be on left side --> look for electrode on right side
            chan2search = sprintf('%s%1.0f',chan_str,chan_num+1);
            elec_index(2,i_chan) = find(strcmpi(channames,chan2search));
        elseif mod(chan_num,2)==0
            % even! electrode must be on right side --> look for electrode on left side
             chan2search = sprintf('%s%1.0f',chan_str,chan_num-1);
             elec_index(2,i_chan) = find(strcmpi(channames,chan2search));
        end
    end
    % output of relevant vector
    elec_index_out = elec_index(idx_hemisphere,:);
end

end

