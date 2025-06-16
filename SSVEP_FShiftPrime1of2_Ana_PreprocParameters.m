%% script to analyse preprocessing parameters

clearvars

%% parameters
% general parameters
F.PathIn            = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\EEG\SCADS\'; % with FWHM 0.5
F.subjects          = cellfun(@(x) sprintf('%02.0f',x),num2cell(1:60),'UniformOutput', false)';
% F.subs2use          = [9 10 11 12];%
F.subs2use          = [1:13 15:52];
F.con2an            ={[10 ]; ... %RDK1 attended; RDK1 and RDK2 colors in periphery peri attended + unattended
                    [20 ]; ... %RDK2 attended; RDK1 and RDK2 colors in periphery peri attended + unattended
                    [30 ]; ... %RDK1 attended; RDK1 and RDK3 colors in periphery peri attended + irrelevant
                    [40 ]; ... %RDK2 attended; RDK2 and RDK3 colors in periphery peri attended + irrelevant
                    [50 ]; ... %RDK1 attended; RDK2 and RDK3 colors in periphery peri unattended + irrelevant
                    [60 ]};    %RDK2 attended; RDK1 and RDK3 colors in periphery peri unattended + irrelevant
F.con2an_label      = {'att RDK1 | RDK3 RDK4';
    'att RDK2 | RDK3 RDK4';
    'att RDK1 | RDK3 RDK5';
    'att RDK2 | RDK4 RDK5';
    'att RDK1 | RDK4 RDK5';
    'att RDK2 | RDK3 RDK5'};

%% loop across subjects to extract
for i_sub = 1:numel(F.subs2use)
    % load mat
    input = load(sprintf('%sVP%s_Preprocess_summary.mat',F.PathIn,F.subjects{F.subs2use(i_sub)}));
    % extract relevant data
    params.SCADS_InterpChansPerTrial(i_sub)=input.SumData.interpolated_channels_avgPERtrial;
    params.trialnr_rem_blinks(i_sub) = numel(find(input.PreProc.trial_blink(ismember(input.PreProc.trial_con,unique(cell2mat(F.con2an))))==0));
    params.trialnr_rem_eyemov(i_sub) = numel(find(input.PreProc.trial_eyemov(ismember(input.PreProc.trial_con,unique(cell2mat(F.con2an))))==0));
    params.trialnr_rem_SCADS(i_sub) = numel(find(input.PreProc.trial_SCADS(ismember(input.PreProc.trial_con,unique(cell2mat(F.con2an))))==0));
    params.trialnr_in_anaylsis(:,i_sub) = cellfun(...
        @(x) numel(find(input.PreProc.trial_blink & input.PreProc.trial_eyemov & input.PreProc.trial_SCADS & ...
        ismember(input.PreProc.trial_con, x))),...
        F.con2an)';
end

%% output
fprintf('interpolated channels per trial: M=%1.3f; Std=%1.3f\n',mean(params.SCADS_InterpChansPerTrial),std(params.SCADS_InterpChansPerTrial))
fprintf('removed trials with blinks: M=%1.3f; Std=%1.3f\n',mean(params.trialnr_rem_blinks),std(params.trialnr_rem_blinks))
fprintf('removed trials with eye movements: M=%1.3f; Std=%1.3f\n',mean(params.trialnr_rem_eyemov),std(params.trialnr_rem_eyemov))
fprintf('removed trials by SCADS: M=%1.3f; Std=%1.3f\n',mean(params.trialnr_rem_SCADS),std(params.trialnr_rem_SCADS))
for i_con = 1:numel(F.con2an_label)
    fprintf('number of trials that entered analysis for condition %s: M=%1.3f; Std=%1.3f\n',...
        F.con2an_label{i_con}, mean(params.trialnr_in_anaylsis(i_con,:),2), std(params.trialnr_in_anaylsis(i_con,:),1,2))
end
fprintf('average trial number per condition M = %1.3f; SD = %1.3f\n',...
     mean(mean(params.trialnr_in_anaylsis,2)), mean(std(params.trialnr_in_anaylsis,1,2)))
fprintf('overall rejection rate of trials M = %1.3f; SD = %1.3f\n',...
    mean(mean(100-params.trialnr_in_anaylsis./120.*100,1)),std(mean(100-params.trialnr_in_anaylsis./120.*100,1)))