%% parameters
clearvars
F.Pathlocal             = 'E:\work\data\SSVEP_FShift_Probabil\';
F.Pathlocal             = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2025_FShift_Prime1of2\';

F.PathInEEG             = fullfile(F.Pathlocal, 'eeg\epoch_erp\');
F.PathInBehavior        = fullfile(F.Pathlocal, 'behavior\');
F.PathInSCADS           = fullfile(F.Pathlocal, 'eeg\SCADS_erp\');
F.PathOut               = fullfile(F.Pathlocal, 'eeg\erp\'); 
F.subjects              = arrayfun(@(x) sprintf('%02.0f',x),1:60,'UniformOutput',false)';
% F.sub2use               = [1:13 15:24];%:53;
% F.sub2use               = [1 3 4 5 6 7 9 10 11 13 15 18 20 21 22 23 24 25];%:53; % for subject 12, 14: eeg and behavior data don't match
F.sub2use               = [2:14 16:28];%:53;

F.trigger               = {[111] [112] [113]}; % {'target primed';'target nonprimed';'distractor'}


F.EEGChans              = 64;
F.ERPepoch              = [-0.4 0.9];
F.ERPbase               = [-0.1 0]; % erp baseline in s
F.CSD_flag              = 1; % 0 = no; 1 = yes
F.ERP_FiltFreq          = [0 13];

%F.TFAfreqs              = [5:(1/6):40];
F.conname_within        = 'event_type';
F.conname_withinlabel   = {'target primed';'target nonprimed';'distractor'};

%% start processing
%% loop across subjects
for i_sub = 1:numel(F.sub2use)
    %% load files
    % EEG
    EEG = pop_loadset('filename',sprintf('VP%s_e.set',F.subjects{F.sub2use(i_sub)}),'filepath',F.PathInEEG);
    prep_input = load(fullfile(F.PathInSCADS,sprintf('VP%s_Preprocess_summary.mat',F.subjects{F.sub2use(i_sub)})));
    % pop_eegplot(EEG,1,1,1)
    % behavior (loads latest file)
    t.files = dir(fullfile(F.PathInBehavior,sprintf('VP%s_timing*.mat',F.subjects{F.sub2use(i_sub)})));
    [t.val t.idx ]=max([t.files.datenum]);
    behavior = load(fullfile(F.PathInBehavior,t.files(t.idx).name));


    temp.files = dir(fullfile(F.PathInBehavior,sprintf('VP%s_timing*.mat',F.subjects{F.sub2use(i_sub)})));
    behavior.resp.experiment = repmat({[nan]},1,17);
    behavior.button_presses.experiment = repmat({[nan]},1,17);
    for i_fi = 1:numel(temp.files)
        temp.data_in{i_fi}=load(sprintf('%s%s',F.PathInBehavior,temp.files(i_fi).name));
        % extract relevant data
        try behavior.conmat.experiment = temp.data_in{i_fi}.conmat.experiment;
        end
        if any(strcmp(fieldnames(temp.data_in{i_fi}.resp),'experiment'))
            temp.index1 = find(~cellfun(@isempty,temp.data_in{i_fi}.resp.experiment));
            temp.index2 = cell2mat(cellfun(@(x) ~isempty(cell2mat({x(:).trialnumber})), temp.data_in{i_fi}.resp.experiment(temp.index1),'UniformOutput',false));
            behavior.resp.experiment(temp.index1(temp.index2))=temp.data_in{i_fi}.resp.experiment(temp.index1(temp.index2));
            behavior.button_presses.experiment(temp.index1(temp.index2))=temp.data_in{i_fi}.button_presses.experiment(temp.index1(temp.index2));
        end
    end

    
         
    %% do csd transform
    if F.CSD_flag==1
        if  i_sub == 1 % calculate CSD matrix
            CSD.chanmat=ExtractMontage('C:\Dropboxdata\Dropbox\work\matlab\software\toolboxes\CSD\resource\10-5-System_Mastoids_EGI129.csd',{EEG.chanlocs.labels}');
            [CSD.G,CSD.H] = GetGH(CSD.chanmat);
        end
        fprintf(1,['\' ...
            'n###\ncalculating CSD transform\n###\n'])
        for i_tr = 1:EEG.trials
            % csd of raw data
            EEG.data(:,:,i_tr)= CSDTransform(EEG.data(:,:,i_tr), CSD.G, CSD.H);
        end
    end
    % pop_eegplot(EEG,1,1,1)

    

    %% additional calculation and condition allocation
    % raw (for FFT)
    % [EEG_ep, t.indices] = pop_epoch(EEG, num2cell(unique(cell2mat(F.trigger))), F.ERPepoch, 'epochinfo', 'yes');
    EEG_ep = pop_select(EEG, 'time', F.ERPepoch);

    EEG_ep = pop_rmbase(EEG_ep,F.ERPbase.*1000);
    % filtered for ERP
    EEG_f = pop_eegfiltnew(EEG, F.ERP_FiltFreq(1), F.ERP_FiltFreq(2), 8*EEG.srate, 0, [], 0);
    % [EEG_fep, t.indices] = pop_epoch(EEG_f, num2cell(unique(cell2mat(F.trigger))), F.ERPepoch, 'epochinfo', 'yes');
    EEG_fep = pop_select(EEG_f, 'time', F.ERPepoch);
    EEG_fep = pop_rmbase(EEG_fep,F.ERPbase.*1000);


    t.behavior = horzcat(behavior.resp.experiment{cellfun(@(x) isstruct(x), behavior.resp.experiment)});
    t.prep_idx = prep_input.PreProc.trial_nrtrialinexp( ...
        prep_input.PreProc.trial_blink & prep_input.PreProc.trial_eyemov & prep_input.PreProc.trial_SCADS);
    t.prep_idx2 = find( ...
        prep_input.PreProc.trial_blink & prep_input.PreProc.trial_eyemov & prep_input.PreProc.trial_SCADS);
    t.behavior = t.behavior(t.prep_idx);


    % add some information
    t.ur_epoch = num2cell(t.prep_idx);
    [t.behavior.urepoch] = deal(t.ur_epoch{:});
    
    % prune information to respective event (as more than one event could have been shown)
    for i_event = 1:numel(t.behavior)
        % index of event
        t.idx = prep_input.PreProc.trial_nreventintrial(t.prep_idx2(i_event));
        % event RDK
        t.behavior(i_event).eventRDK = t.behavior(i_event).eventRDK(t.idx);
        % event type
        t.behavior(i_event).eventtype2 = F.conname_withinlabel{t.behavior(i_event).eventtype2(t.idx)};
        % event color
        t.behavior(i_event).eventcolor = behavior.RDK.RDK( t.behavior(i_event).eventRDK).col_label;
        % event frequency
        t.behavior(i_event).eventfreq = behavior.RDK.RDK( t.behavior(i_event).eventRDK).freq;
        % event direction
        t.behavior(i_event).eventdirection = t.behavior(i_event).eventdirection(t.idx);
        % event onset frames
        t.behavior(i_event).event_onset_frames = t.behavior(i_event).event_onset_frames(t.idx);
        % event onset times
        t.behavior(i_event).event_onset_times = t.behavior(i_event).event_onset_times(t.idx);
        % event_response_type
        t.behavior(i_event).event_response_type = t.behavior(i_event).event_response_type{t.idx};
        % event_response_RT
        t.behavior(i_event).event_response_RT = t.behavior(i_event).event_response_RT(t.idx);

        % onset times after cue
        t.behavior(i_event).postcue_onset = ...
            (t.behavior(i_event).event_onset_times-t.behavior(i_event).pre_cue_times)*1000;
       
    end

    %% some bookkeeping
    EP.behavior = t.behavior;
    EP.EEG_ep = EEG_ep;
    EP.EEG_fep = EEG_fep;
    EP.parameters = F;
    EP.RDK = behavior.RDK.RDK;
    
    %% save
    EP.savetime=datestr(now);
    if ~exist(F.PathOut); mkdir(F.PathOut); end
    fprintf(1,'|| saving file ||  %s\\VP%s_events_erp.mat ||\n', ...
        F.PathOut,F.subjects{F.sub2use(i_sub)})
    save(fullfile(F.PathOut, sprintf('VP%s_events_erp.mat',F.subjects{F.sub2use(i_sub)})), 'EP')
    
end