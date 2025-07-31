% Script for analyzing EEG data from SSVEP_FShiftProbabil
%
%
%    2024 - C. Gundlach



%% General Definitions 
clearvars
p.path=             '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\';
p.bdf_path=         '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\eeg\raw\';
p.set_path=         [p.path 'eeg\set\'];
p.epoch_path=       [p.path 'eeg\epoch_erp\'];
p.scads_path=       [p.path 'eeg\SCADS_erp\'];
% p.chanlocs_path=    ['C:\Users\psy05cvd\Dropbox\work\matlab\Auswertungsskripte\Analyzer_G\ChanLocs\BioSemi64_1020.epf'];
p.chanlocs_path=    ['C:\Users\EEG\Documents\MATLAB\lab_library\BS_Chanlocs\BioSemi64_1020.epf'];
p.mean_path=        [p.path 'eeg\mean\'];
p.subs=             cellfun(@(x) sprintf('%02.0f',x),num2cell(1:60),'UniformOutput', false)';
% p.subs2use=         [1 3:6 7 9 10 11 12];%  % participant 2,8 measurement cancelled due to bad behavior
% p.subs2use=         [1:13 15:24];%
p.subs2use=         [49:52];%
p.part=             {'1';'2';'3'};
p.cue =          {[10 11 12 16 17 18 19]; ... %RDK1 attended; RDK1 and RDK2 colors in periphery peri attended + unattended
                    [20 21 22 26 27 28 29]; ... %RDK2 attended; RDK1 and RDK2 colors in periphery peri attended + unattended
                    [30 31 32 36 37 38 39]; ... %RDK1 attended; RDK1 and RDK3 colors in periphery peri attended + irrelevant
                    [40 41 42 46 47 48 49]; ... %RDK2 attended; RDK2 and RDK3 colors in periphery peri attended + irrelevant
                    [50 51 52 56 57 58 59]; ... %RDK1 attended; RDK2 and RDK3 colors in periphery peri unattended + irrelevant
                    [60 61 62 66 67 68 69]};    %RDK2 attended; RDK1 and RDK3 colors in periphery peri unattended + irrelevant
p.events=           {[201] [202]}; % trigger for events in trial (up to two possible)
p.con1name =        'type';
p.con1label =       {'target';'distractor'};
p.cue_epoch=        [-1 3.5];
p.event_epoch=      [-1 1.5];
% p.epoch2an=         [-1 3.25]; % not so conservative
p.epoch2an=         [-0.5 1]; % not so conservative
p.resample=         256;

p.AnaScriptName=    'SSVEP_FShiftProbabil_preprocessing';

% flags
ImportFlag=         0; % set to 1 if files have ti be importet from raw .bdf files
EpochFlag=          1; % set to 1 if data has to be epoched
MeanFlag=           0; % set to 1 to create mean files with trials averaged for each subject
TopoFlag=           0;
if ImportFlag==0 && EpochFlag==0 && MeanFlag==0 && TopoFlag==0, help(AnaScriptName),return, end

%% Main Script 
% loop for subjects
for i_sub=1:numel(p.subs2use)
    FileName=sprintf('VP%s',p.subs{p.subs2use(i_sub)});
    if ~MeanFlag || ((MeanFlag && EpochFlag) || (ImportFlag && EpochFlag && MeanFlag) || (MeanFlag && ~TopoFlag)),end
    
    %% import
    if ImportFlag    %% Import Files and merge %%
        % load file
%         EEG=pop_biosig([p.bdf_path FileName '.bdf']);
        % load data
        temp.files = dir(sprintf('%sVP%s*.bdf',p.bdf_path,p.subs{p.subs2use(i_sub)}));
        
        for i_fi = 1:numel(temp.files)
            EEG(i_fi)=pop_readbdf(sprintf('%s%s',p.bdf_path,temp.files(i_fi).name),[],73,[]);
        end
        if numel(EEG)>1,EEG=pop_mergeset(EEG,1:numel(EEG),0); end
        %pop_eegplot(EEG,1,1,1)
        
        % replace boundary trigger by 999999
        if sum(strcmp({EEG.event.type},'boundary'))~=0
            EEG.event(strcmp({EEG.event.type},'boundary')).type = '9999999';
            for i_ev = 1:numel(EEG.event)
                EEG.event(i_ev).type = str2num(EEG.event(i_ev).type);
            end
        end
        
        % troubleshooting
        if any(unique([EEG.event.type])>16128 & unique([EEG.event.type])~=9999999)
            t.triggernew = num2cell([EEG.event.type]-16128);
            [EEG.event.type]=t.triggernew{:};
            [EEG.urevent.type]=t.triggernew{:};
        end
        
        % added median filter
        %EEG.data([65 66 67 68],:) = medfilt1(EEG.data([65 66 67 68],:),ceil(EEG.srate/40),length(EEG.data),2);
        
        % resample
        EEG=pop_resample(EEG,p.resample);
        
        %
%         figure; pwelch(EEG.data(29,:),EEG.srate*1,EEG.srate/2,256,EEG.srate)
        
        % check trigger
        trigg.all = unique(cell2mat({EEG.event.type}));
        trigg.freq = arrayfun(@(x) sum(cell2mat({EEG.event.type})==x), trigg.all);
        trigg.sum = [trigg.all; trigg.freq];
        trigg.events = cellfun(@(x) sum(ismember(cell2mat({EEG.event.type}), x)), p.events);
        
        
        % save in new format
        if ~exist(p.set_path); mkdir(p.set_path);end
        EEG = pop_saveset(EEG,'filename',[FileName '.set'], 'filepath', p.set_path);
        clear ('EEG');
    end
    
    %% epoch Bipolarize, Remove Baseline, Detrend, Blinks + EyeMovements Threshold, SCADS%%
    if EpochFlag
        % load
        EEG = pop_loadset('filename',[FileName '.set'], 'filepath', p.set_path);
        % pop_eegplot(EEG,1,1,1)

        % check trigger
        t.idx_t = find([EEG.event.type] == 201); % target
        t.idx_d = find([EEG.event.type] == 202); % target
        t.freq = {'target' numel(t.idx_t) sum([EEG.event(t.idx_t+1).type]==2048); ...
            'distractor' numel(t.idx_t) sum([EEG.event(t.idx_d+1).type]==2048)};
        
        % epoch data around cue
        EEG_cue = pop_epoch( EEG, num2cell(unique(cell2mat(p.cue))), p.cue_epoch, 'epochinfo', 'yes');

        % epoch data around event
        EEG_event = pop_epoch( EEG_cue, num2cell(unique(cell2mat(p.events))), p.event_epoch, 'epochinfo', 'yes');
        
        % some bookkeeping to keep track of discarded trials
        PreProc.trial_nr = 1:numel(EEG_event.epoch);
        PreProc.trial_nrtrialinexp = nan(1,numel(EEG_event.epoch));
        PreProc.trial_con= nan(1,numel(EEG_event.epoch));
        PreProc.trial_blink=true(1,numel(EEG_event.epoch));
        PreProc.trial_eyemov=true(1,numel(EEG_event.epoch));
        PreProc.trial_SCADS=true(1,numel(EEG_event.epoch));
        for i_ep = 1:numel(EEG_event.epoch)
            % which event trigger?
            t.t=cell2mat(EEG_event.epoch(i_ep).eventtype);
            t.idx = cell2mat(EEG_event.epoch(i_ep).eventlatency)==0;
            PreProc.trial_con(i_ep)=t.t(t.idx);
            % which trial in experiment?
            PreProc.trial_nrtrialinexp(i_ep) = ...
                EEG_cue.event([EEG_cue.event.urevent] == EEG_event.epoch(i_ep).eventurevent{t.idx}).epoch;           
        end
        % first or second event
        PreProc.trial_nreventintrial = ones(size(PreProc.trial_nrtrialinexp));
        PreProc.trial_nreventintrial([false diff(PreProc.trial_nrtrialinexp)==0]) = 2;
        
        % create single horizontal and single vertical eye channel by
        % subtracting channels (e.g. VEOG1-VEOG2) from each other to
        % increase SNR of blinks and eye movements
        EEG_event = eegF_Bipolarize(EEG_event);
        
        % remove linear drift and offset of EEG signals
        EEG_event = eegF_Detrend(EEG_event,p.epoch2an); %
        % pop_eegplot(EEG,1,1,1)
        
        % index trials with blinks and eye movements for later rejection
        [EEG_temp,Trials2Remove1] = eegF_RemoveBlinks(EEG_event,65,p.epoch2an,1); % letzte Parameter '1': Artefakte nur markiert
        PreProc.trial_blink(Trials2Remove1)=false;
        [EEG_temp,Trials2Remove2] = eegF_RemoveEyeMovements(EEG_event,[65 66],p.epoch2an,25,1);
        PreProc.trial_eyemov(Trials2Remove2)=false;
        EEG_event = pop_select( EEG_event,'trial',find(PreProc.trial_eyemov& PreProc.trial_blink));
        
        s
        
        % filter before SCADS (is this a good idea?)
%         EEG = pop_eegfiltnew(EEG, 0.7, 0, 16*EEG.srate, 0, [], 0);
        
        % run SCADS to statistically detect artifacts
        [EEG_event Trials2Del SumData] = SCADS_Pass1(EEG_event,p.epoch2an,[6 6 15],1);%,[6 5 15]);  % ansonsten SCADS_Pass1()
        t.index = find(PreProc.trial_blink & PreProc.trial_eyemov);
        PreProc.trial_SCADS(t.index(Trials2Del))=false;        
        
        % save trials marked as artifacts for potential later inspection
        if ~isempty(Trials2Del)
            EEG_SCADS = pop_select(EEG_event, 'trial', Trials2Del);
        else
            EEG_SCADS = pop_select(EEG_event, 'trial', 1);
        end
        % discard trials indexed as containing artifacts
        EEG_event = pop_rejepoch(EEG_event, EEG_event.reject.rejmanual, 0);
        % rereference to average reference
        EEG_event=pop_reref(EEG_event,[],'refstate',0);   %
        % pop_eegplot(EEG,1,1,1)
        
        if ~exist(p.scads_path); mkdir(p.scads_path);end
        if ~exist(p.epoch_path); mkdir(p.epoch_path);end
        save([p.scads_path FileName '_Preprocess_summary.mat'],'SumData','PreProc')
        pop_saveset(EEG_SCADS,[FileName '_SCADS.set'],p.scads_path);% Save Data
        pop_saveset(EEG_event, 'filename', [FileName '_e.set'], 'filepath', p.epoch_path);
        clear EEG EEG_event EEG_cue
    end
end

%% Mean Trials for every Subject%% TO DO: Artefakte checken und rauswerfen
if MeanFlag 
    for i_sub=1:size(Subject,1)
        disp(' '), disp(['Subject ' Subject(i_sub,:)])
        %Mitteln der Trials, in beliebiger Variable ablegen
        FileName=[exp_name '-' Subject(i_sub,:)];
        EEG = pop_loadset('filename',[FileName '_e.set'], 'filepath', epoch_path);
        if i_sub==1
            TrialsAveraged=zeros(size(Subject,1),size(Events,1));
            NbChannels=EEG.nbchan;
            NbSamples=EEG.pnts;
            AvgTrials=zeros(size(Events,1),NbChannels,NbSamples,size(Subject,1));
            Avg = zeros(size(Events,1),NbChannels,NbSamples,size(Subject,1));
        end
        for i_event=1:size(Events,1)
            % [EEG2,epind]=pop_selectevent(EEG,'type',Events{i_event,:},'deleteevents','off','deleteepochs','on');
            EEG2 = pop_selectevent( EEG, 'omittype',cell2mat(Events(cell2mat(Events)~=Events{i_event}))'...
                ,'deleteevents','off','deleteepochs','on','invertepochs','off');
            TrialsAveraged(i_sub,i_event)=size(EEG2.data,3);
            AvgTrials(i_event,:,:,i_sub)=mean(EEG2.data,3);
            
%             Avg(i_event,:,:,i_sub) = mean(EEG2.data(:,:,[EEG.event(1,epind).epoch]),3);
        end
    end
    %% Eine Datei pro Bedingung speichern! %%
    EEG=pop_selectevent(EEG,'epoch',1:size(Subject,1));
    Conditions=Events;
%     chkpath(mean_path,'linux',1);
    for i_cond=1:size(Conditions,1)
        EEG.data=squeeze(AvgTrials(i_cond,:,:,:));
%         EEG.data=squeeze(Avg(i_cond,:,:,:,i_freq));
        EEG.setname=[exp_name '-Mean-C' int2str(i_cond)];
        EEG=pop_saveset(EEG,'filename',[EEG.setname '.set'], 'filepath', mean_path);
    end
    TrialsAveraged; % Anzahl der übrigen Trials
    save([path 'EEG\Dokumentation\TrialsAveraged_conditionmeans.mat'],'TrialsAveraged')
    clear('EEG')
end