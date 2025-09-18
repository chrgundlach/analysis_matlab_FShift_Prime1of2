%% script to test contrasts of SSVEP timecourse analysis
% load FShiftPrime1of2

clearvars
F.PathInEEG             = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2025_FShift_Prime1of2\EEG\TFA'; % with FWHM 0.5

F.Subs                  = arrayfun(@(x) sprintf('%02.0f',x),1:70,'UniformOutput',false)';
% F.Subs2use              = [1:13 15:21];
% changed experiment from participant 22 onwards (stimuli isoluminant to
% background and used other frequencies
% participant 42 has lower trial number
F.Subs2use              = [1:14 16:28]; % no sub 15
                        
F.TFA.baseline          = [-500 -250];

F.SSVEP_Freqs           = [17 20 23]; 
F.RDK_pos               = [0 0 0];
F.RDK_pos_label         = {'center';'center';'center'};



F.conlabel_att = {'att RDK1+2';'att RDK2+3'; 'att RDK3+1'};
F.conlabel_primedRDK = {'RDK1';'RDK2'; 'RDK3'};
F.conlabel_nonprimedRDK = {'RDK2';'RDK3'; 'RDK1'};
F.conRDKattended = logical([1 1 0; 0 1 1; 1 0 1]);
F.conRDKprimed = logical([1 0 0; 0 1 0; 0 0 1]);
F.conRDKnonprimed = logical([0 1 0; 0 0 1; 1 0 0]);
F.conRDKattended_label = repmat({'attended'},size(F.conRDKattended));
F.conRDKattended_label(F.conRDKattended==0) = {'not attended'};
F.conRDKprimed_label = F.conRDKattended_label;
F.conRDKprimed_label(F.conRDKprimed==1) = {'primed'};
F.conRDKprimed_label(F.conRDKnonprimed==1) = {'nonprimed'};


pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];
%% load data
for i_sub = 1:numel(F.Subs2use)
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%s_tfa.mat ||\n',...
        i_sub,numel(F.Subs2use),F.PathInEEG,F.Subs{F.Subs2use(i_sub)})
    
    temp.tfa = open(fullfile(F.PathInEEG,sprintf('VP%s_tfa.mat',F.Subs{F.Subs2use(i_sub)})));
    
    % convert to single
    temp.tfa.TFA.data_evo = single(temp.tfa.TFA.data_evo);
    temp.tfa.TFA.data_ind = single(temp.tfa.TFA.data_ind);
    temp.tfa.TFA.FFT.data_evo = single(temp.tfa.TFA.FFT.data_evo);
    temp.tfa.TFA.FFT.data_ind = single(temp.tfa.TFA.FFT.data_ind);
    
    
    % preallocate memory
    if i_sub == 1
        TFA.data_evo = single(nan([size(temp.tfa.TFA.data_evo),numel(F.Subs2use)]));
        TFA.data_ind = single(nan([size(temp.tfa.TFA.data_ind),numel(F.Subs2use)]));
        TFA.time = temp.tfa.TFA.t;
        TFA.frequency = temp.tfa.TFA.f;
        TFA.electrodes = temp.tfa.TFA.electrodes;
        TFA.con_trialnum = temp.tfa.TFA.con_trialnum;
        TFA.srate = temp.tfa.TFA.params.srate/2;
        TFA.fftdata_ind = nan([size(temp.tfa.TFA.FFT.data_ind),numel(F.Subs2use)]);
        TFA.fftdata_evo = nan([size(temp.tfa.TFA.FFT.data_evo),numel(F.Subs2use)]);
        TFA.ffttimewin = temp.tfa.TFA.FFT.timewin;
        TFA.fftfreqs = temp.tfa.TFA.FFT.freqs;

        TFA.Gabor_FWHM_freq = temp.tfa.TFA.params.gabor_FWHM_freq;
        TFA.Gabor_FWHM_time = 2*log(2)/(pi*TFA.Gabor_FWHM_freq)*500; % time in ms
    end
    
    % assign data
    TFA.data_evo(:,:,:,:,i_sub) = temp.tfa.TFA.data_evo; % evoked data
    TFA.data_ind(:,:,:,:,i_sub) = temp.tfa.TFA.data_ind; % induced data
    %     TFA(i_exp).data_bc(:,:,:,i_sub) = bsxfun(@minus, temp.tfa.TFA.data_induced, ...
    %         mean(temp.tfa.TFA.data(:,eeg_time2points(TFA.baseline(1),TFA.time):eeg_time2points(TFA.baseline(2),TFA.time),:,:,:),2));
%     TFA(i_exp).data_bc(:,:,:,:,i_sub) = 100*((bsxfun(@rdivide, temp.tfa.TFA.data, ...
%         mean(temp.tfa.TFA.data(:,eeg_time2points(F.TFA.baseline(1),TFA(i_exp).time):eeg_time2points(F.TFA.baseline(2),TFA(i_exp).time),:,:,:),2)))-1);
    TFA.fftdata_evo(:,:,:,:,i_sub) = temp.tfa.TFA.FFT.data_evo;
    TFA.fftdata_ind(:,:,:,:,i_sub) = temp.tfa.TFA.FFT.data_ind;
    TFA.RDK(i_sub) = temp.tfa.TFA.RDK;

      
    clear temp    
end

pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];


%% plot electrode head
% pl.elec2plot = {'C3';'CP3'};
pl.elec2plot = {'P5';'P7';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'P8';'P10';'PO4';'PO8';'O2';'I2'};
pl.elec2plot = {'P5';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'PO4';'PO8';'O2';'I2'};
% pl.elec2plot = {'P6';'P8';'P10';'PO4';'PO8';'O2';'I2';'POz';'Oz';'Iz';'O1'};
% pl.elec2plot = {'P5';'P7';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'O2'};
% pl.elec2plot = {'POz';'O1';'Oz';'I2';'Iz'};
pl.elec2plot = {'P8';'P10';'PO8';'PO4';'POz';'Oz';'O1'};
% pl.elec2plot = {'P7';'P9';'PO7';'PO3';'POz';'Oz';'O2'};


figure;
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
% topoplot(find(pl.elec2plot_i),TFA.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on');
topoplot([],TFA.electrodes(1:64),'whitebk','on','style','blank','electrodes', 'on', 'emarker2', {find(pl.elec2plot_i),'o','r',8});


%% calculate everything with running t-tests ans cluster correction | central
% plotting parameters
% pl.elec2plot = {'Oz';'Iz'};
% pl.elec2plot = {'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';% vis alpha II
pl.elec2plot = {'P5';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'PO4';'PO8';'O2';'I2'}; sav.chan_add = 'VisualLarge';
pl.elec2plot = {'P7';'P5';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P8';'P6';'PO4';'PO8';'P10';'O2';'I2'}; sav.chan_add = 'VisualLarge'; % suppression irrelevant study[to be used]


% cluster analysis
% pl.elec2plot = {TFA.electrodes(1:64).labels}';
% pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) startsWith({TFA(1).electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.freqrange=[-0.1 0.1];

pl.time_post = [0 1800];

pl.xlims=[-1000 1800]; % index time 2 plot
pl.xlims_i = dsearchn(TFA.time', pl.xlims');

pl.ylims = [-20 40];

pl.base = F.TFA.baseline;
% pl.base = [-500 -0];
pl.base_i = dsearchn(TFA.time', pl.base');

pl.sub2plot = 1:numel(F.Subs2use);
pl.sub2plot(ismember(F.Subs2use,5))=[]; % participant 5 has no ssveps


pl.data_ind = []; pl.data_evo = []; pl.data_ind_bc = []; pl.data_evo_bc = [];
pl.RDKlabel = {'RDK1';'RDK2';'RDK3'};
pl.RDKidx = [1 2 3];

pl.con2plot = {'primed';'nonprimed';'not attended'}; %F.conRDKprimed_label(1,:)
pl.concols = num2cell([255 133 4; 41 60 74; 25 138 131]'./255,1);

% extract data
for i_sub = 1:numel(pl.sub2plot)
    t.fidx = cell2mat(arrayfun(@(x) dsearchn(TFA.frequency', x+pl.freqrange'),[TFA.RDK(pl.sub2plot(i_sub)).RDK(pl.RDKidx).freq],'UniformOutput',false));
    pl.freqlabel = TFA.frequency(t.fidx(1,1):t.fidx(2,1))-TFA.RDK(pl.sub2plot(i_sub)).RDK(1).freq;
    for i_RDK = 1:size(t.fidx,2)
        for i_con = 1:numel(pl.con2plot)
            % which condition
            t.cidx = strcmp(F.conRDKprimed_label(:,i_RDK),pl.con2plot{i_con});
            % raw
            pl.data_ind(:,i_con,i_RDK,i_sub) = ...
                squeeze(mean(TFA.data_ind(t.fidx(1,i_RDK):t.fidx(2,i_RDK),:,pl.elec2plot_i,t.cidx ,pl.sub2plot(i_sub)),[3,1]));
            pl.data_evo(:,i_con,i_RDK,i_sub) = ...
                squeeze(mean(TFA.data_evo(t.fidx(1,i_RDK):t.fidx(2,i_RDK),:,pl.elec2plot_i,t.cidx ,pl.sub2plot(i_sub)),[3,1]));
            % baseline corrected
            pl.data_ind_bc(:,i_con,i_RDK,i_sub) = ...
                100*(...
                bsxfun(@rdivide, pl.data_ind(:,i_con,i_RDK,i_sub), ...
                squeeze(mean(TFA.data_ind(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.base_i(1):pl.base_i(2),pl.elec2plot_i,t.cidx ,pl.sub2plot(i_sub)),[3,1,2]))')...
                -1);
            pl.data_evo_bc(:,i_con,i_RDK,i_sub) = ...
                100*(...
                bsxfun(@rdivide, pl.data_evo(:,i_con,i_RDK,i_sub), ...
                squeeze(mean(TFA.data_evo(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.base_i(1):pl.base_i(2),pl.elec2plot_i,t.cidx ,pl.sub2plot(i_sub)),[3,1,2]))')...
                -1);
        end
    end
end


% extract onset values with cusum method | actual Benchmark
t.time_rt = pl.time_post;
t.time_rt_i = dsearchn(TFA.time', t.time_rt');
pl.testdata = pl.data_evo(t.time_rt_i(1):t.time_rt_i(2),:,:,:);
pl.basedata = pl.data_evo(pl.base_i(1):pl.base_i(2),:,:,:);
pl.rdkfreqs = cell2mat(arrayfun(@(x) [TFA.RDK(x).RDK(1:3).freq]',pl.sub2plot,'UniformOutput',false))';
% only two conditions primed vs nonprimed [1 2]
t.idx = [1 2];
% only two conditions primed vs unattended [1 3]
t.idx = [1 3];
% only two conditions nonprimed vs unattended [2 3]
t.idx = [2 3];

[onset] = eeg_GaborStatOnsetCuSum(pl.testdata(:,t.idx,:,:),pl.basedata(:,t.idx,:,:), ...
    'datadims',{'time','cons','RDK','subs'}, ...
    'conlabels',repmat(pl.con2plot(t.idx),1,size(t.fidx,2)), ...
    'rdkfreqs',pl.rdkfreqs, ...
    'timevec',TFA.time(t.time_rt_i(1):t.time_rt_i(2)));





% running ttests
t.time_rt = pl.time_post;
t.time_rt_i = dsearchn(TFA.time', t.time_rt');

t.permut_n = 5000;
clear cluster_runt timecourse_runt

% run cluster correction for tests against zero
for i_con = 1:size(pl.data,2)
    t.data = squeeze(pl.data(t.time_rt_i(1):t.time_rt_i(2),i_con,:));
    t.nulldata = repmat(squeeze(mean(pl.data(pl.base_i(1):pl.base_i(2),i_con,:),1))',[numel(t.time_rt_i(1):t.time_rt_i(2)),1]);
    [cluster_runt{i_con}, timecourse_runt{i_con}]=eeg_erpStat_clusterP(t.data,t.nulldata,t.permut_n,2);
end

% run cluster correction for diffs against zero
t.diffs = [1 2];
for i_diff = 1:size(t.diffs,1)
    t.data = squeeze(diff(squeeze(pl.data(t.time_rt_i(1):t.time_rt_i(2),t.diffs(i_diff,:),:)),1,2));
    t.nulldata = repmat(diff(squeeze(mean(pl.data(pl.base_i(1):pl.base_i(2),t.diffs(i_diff,:),:),1))),[numel(t.time_rt_i(1):t.time_rt_i(2)),1]);
    [cluster_runt{size(pl.data,2)+i_diff}, timecourse_runt{size(pl.data,2)+i_diff}]=eeg_erpStat_clusterP(t.data,t.nulldata,t.permut_n,2);
end

pl.mdata = mean(pl.data,3); % mean data
pl.semdata = std(pl.data,1,3)./sqrt(numel(pl.sub2plot));

pl.conlabel = pl.conunique;
pl.col = pl.concols;
pl.col2 = [0.6 0.6 0.6];
pl.line = {'-';'-'};
figure;
subplot(7,1,[1:4])
h.pl = {}; h.plsem=[];  h.plm = []; h.pls = []; h.plst = [];
for i_con = 1:numel(pl.conlabel)
    % data index
    pl.idx = pl.xlims_i(1):pl.xlims_i(2);
    
    % plot SEM as boundary
    % create data
    pl.xconf = [TFA.time(pl.idx) TFA.time(pl.idx(end:-1:1))] ;
    pl.yconf = [pl.mdata(pl.idx,i_con,:)'+pl.semdata(pl.idx,i_con,:)' ...
        pl.mdata(pl.idx(end:-1:1),i_con,:)'-pl.semdata(pl.idx(end:-1:1),i_con,:)'];
    % plot
    h.plsem{i_con} = fill(pl.xconf,pl.yconf,pl.col{i_con}','EdgeColor','none','FaceAlpha',0.3);
    hold on
    
    % plot mean lines    
    h.plm{i_con}=plot(TFA.time(pl.idx), pl.mdata(pl.idx,i_con,:),'Color',pl.col{i_con},'LineStyle',pl.line{i_con},'LineWidth',2);
    
end
grid on
set(gca,'XTickLabel',[])
xlim(pl.xlims)
if ~isempty(pl.ylims)
    ylim(pl.ylims)
end
ylabel('modulation in %')

% plot lines for significant effects
pl.sign_y = 1:size(timecourse_runt,2);

subplot(7,1,[5:7])

for i_con = 1:numel(pl.conlabel)
        
    % uncorrected
    pl.sigdata = timecourse_runt{ i_con}.h_raw.*pl.sign_y(i_con);
    pl.sigdata(timecourse_runt{ i_con}.h_raw==0)=nan;
    
    h.pls{i_con}=plot(TFA.time(t.time_rt_i(1):t.time_rt_i(2)), pl.sigdata,...
        'Color',pl.col2,'LineWidth',6);
    hold on
    
    % corrected
    pl.sigdata = timecourse_runt{ i_con}.h_corr.*pl.sign_y(i_con);
    pl.sigdata(timecourse_runt{ i_con}.h_corr==0)=nan;
    
    h.pls{i_con}=plot(TFA.time(t.time_rt_i(1):t.time_rt_i(2)), pl.sigdata,...
        'Color',pl.col{i_con},'LineWidth',6);
    
end

for i_diff = 1:size(t.diffs,1)
    t.idx = numel(pl.conlabel)+ i_diff;

    % uncorrected
    pl.sigdata = timecourse_runt{t.idx}.h_raw.*pl.sign_y(t.idx);
    pl.sigdata(timecourse_runt{ t.idx}.h_raw==0)=nan;
    
    h.pls{t.idx}=plot(TFA.time(t.time_rt_i(1):t.time_rt_i(2)), pl.sigdata,...
        'Color',pl.col2,'LineWidth',6);
    hold on
    
    % corrected
    pl.sigdata = timecourse_runt{t.idx}.h_corr.*pl.sign_y(t.idx);
    pl.sigdata(timecourse_runt{t.idx}.h_corr==0)=nan;
    
    h.pls{ t.idx}=plot(TFA.time(t.time_rt_i(1):t.time_rt_i(2)), pl.sigdata,...
        'Color',pl.col{t.diffs(i_diff,1)},'LineWidth',6);
    h.pls{t.idx}=plot(TFA.time(t.time_rt_i(1):t.time_rt_i(2)), pl.sigdata,...
        'Color',pl.col{t.diffs(i_diff,2)},'LineWidth',3);
    
end
xlabel('time in ms')

legend([h.pls{1:numel(pl.conlabel)}],pl.conlabel,'Location','SouthOutside','Orientation','horizontal')
legend('boxoff')

ylim(pl.sign_y([1 end])+[-1 1])
xlim(pl.xlims)

