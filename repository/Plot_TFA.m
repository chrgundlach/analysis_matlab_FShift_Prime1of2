%% plot TFA images
clearvars
F.PathInEEG             = '\eeg_data'; % with FWHM 0.5

F.Subs                  = arrayfun(@(x) sprintf('%02.0f',x),1:70,'UniformOutput',false)';
F.Subs2use              = [1:14 16:28]; % no sub 15
                        
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
F.conRDKprimed_label(F.conRDKprimed==1) = {'stim_cued'};
F.conRDKprimed_label(F.conRDKnonprimed==1) = {'top_down_cued'};


pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];
%% load data
for i_sub = 1:numel(F.Subs2use)
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%s_tfa.mat ||\n',...
        i_sub,numel(F.Subs2use),F.PathInEEG,F.Subs{F.Subs2use(i_sub)})
    
    temp.tfa = open(fullfile(F.PathInEEG,sprintf('VP%s_tfa.mat',F.Subs{F.Subs2use(i_sub)})));
    
    % convert to single
    temp.tfa.TFA.data_evo = single(temp.tfa.TFA.data_evo);
     
    
    % preallocate memory
    if i_sub == 1
        TFA.data_evo = single(nan([size(temp.tfa.TFA.data_evo),numel(F.Subs2use)]));
        TFA.time = temp.tfa.TFA.t;
        TFA.frequency = temp.tfa.TFA.f;
        TFA.electrodes = temp.tfa.TFA.electrodes;
        TFA.con_trialnum = temp.tfa.TFA.con_trialnum;
        TFA.srate = temp.tfa.TFA.params.srate/2;
        TFA.Gabor_FWHM_freq = temp.tfa.TFA.params.gabor_FWHM_freq;
        TFA.Gabor_FWHM_time = 2*log(2)/(pi*TFA.Gabor_FWHM_freq)*500; % time in ms
    end
    
    % assign data
    TFA.data_evo(:,:,:,:,i_sub) = temp.tfa.TFA.data_evo; % evoked data
    TFA.RDK(i_sub) = temp.tfa.TFA.RDK;

      
    clear temp    
end

pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];

%% plot grand mean topographies for all SSVEPs | RDKs
pl.elec2plot = {{'P7';'P5';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P8';'P6';'PO4';'PO8';'P10';'O2';'I2'} 'center'};

pl.elec2plot_i=cellfun(@(y) ...
    logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), y, 'UniformOutput',false)),1)),...
    pl.elec2plot(:,1), 'UniformOutput', false);

pl.freq2plot = F.SSVEP_Freqs;
pl.fidx = cell2mat(arrayfun(@(x) dsearchn(TFA.frequency',(x)'),pl.freq2plot,'UniformOutput',false));

pl.time2plot = [-500 -TFA.Gabor_FWHM_time];
pl.tidx = dsearchn(TFA.time', pl.time2plot');pl.tidx = pl.tidx(1):pl.tidx(2);

pl.sub2plot = 1:numel(F.Subs2use);
pl.sub2plot(ismember(F.Subs2use,5))=[]; % participant 5 has no ssveps


pl.data_evo = squeeze(mean(TFA.data_evo(pl.fidx,pl.tidx,:,:,pl.sub2plot),[1,2,4,5]));

pl.clim = [0 max(pl.data_evo)];

% actual plotting
figure('Position',[100 100 250 250]);

topoplot( pl.data_evo, TFA.electrodes(1:64), ...
   'numcontour', 0, 'maplimits',pl.clim,'conv','on','colormap',fake_parula,'whitebk','on', ...
   'emarker2',{find(any(cell2mat(pl.elec2plot_i),1)),'o','r',2,1})
% colormap(fake_parula)
% colormap(parula)
% colormap(jet)
% colormap(turbo)
% colormap(viridis)
% colormap("nebula")
% colormap("copper")
% colormap(cbrewer2("PuBuGn"))
colormap(cbrewer2("YlGnBu"))
% colormap(flipud(cbrewer2("YlGnBu")))
% colormap(cbrewer2("BuPu"))
% colormap(cbrewer2("Greys"))
colorbar

% exportgraphics(gcf,'figures/SSVEP_GrandMeanTopoBaseline.pdf','ContentType','vector')


%% plot grand mean Gabor data | spectra | for distinct frequencies (lookup of respective electrode cluster)

% large center as in tango | periphery: central and lateral 
pl.elec2plot = {{'P7';'P5';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P8';'P6';'PO4';'PO8';'P10';'O2';'I2'} 'center'};

pl.elec2plot_i=cellfun(@(y) ...
    logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), y, 'UniformOutput',false)),1)),...
    pl.elec2plot(:,1), 'UniformOutput', false);

pl.base = [-500 -TFA.Gabor_FWHM_time];
pl.base_i = dsearchn(TFA.time', pl.base');


pl.sub2plot = 1:numel(F.Subs2use);
pl.sub2plot(ismember(F.Subs2use,5))=[]; % participant 5 has no ssveps


% extract data
pl.data_evo = squeeze(mean(TFA.data_evo(:,pl.base_i(1):pl.base_i(2),pl.elec2plot_i{1},:,pl.sub2plot),[2,3,4]));


% plotting
figure;
set(gcf,'Position',[100 100 600 300],'PaperPositionMode','auto')
plot(TFA.frequency,pl.data_evo,'Color',[0.5 0.5 0.5],'LineWidth',1)
hold on;
plot(TFA.frequency,mean(pl.data_evo,2),'Color','k','LineWidth',2)

xlim([0 35])
xlabel('frequency in Hz')
ylabel('amplitude in \muV/m²')
title(sprintf('evoked GrandMean spectra for Gabor Transforms| N = %1.0f | FOI = %1.1f %1.1f %1.1f Hz', ...
    numel(pl.sub2plot), F.SSVEP_Freqs),'Interpreter','none')
vline(F.SSVEP_Freqs,'k:')
box on


% exportgraphics(gcf,'figures/SSVEP_GrandMeanGaborSpectraBaseline.pdf','ContentType','vector')


%% calculate everything with running t-tests and cluster correction | central
% plotting parameters
pl.elec2plot = {'P7';'P5';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P8';'P6';'PO4';'PO8';'P10';'O2';'I2'}; sav.chan_add = 'VisualLarge'; % suppression irrelevant study[to be used]
% cluster analysis
% pl.elec2plot = {TFA.electrodes(1:64).labels}';
% pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) startsWith({TFA(1).electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.freqrange=[-0.05 0.05];

pl.xlims=[-1000 1700]; % index time 2 plot
pl.xlims_i = dsearchn(TFA.time', pl.xlims');

pl.time_post = [0 1700];

pl.base = [-500 -TFA.Gabor_FWHM_time];
pl.base_i = dsearchn(TFA.time', pl.base');

pl.sub2plot = 1:numel(F.Subs2use);
pl.sub2plot(ismember(F.Subs2use,5))=[]; % participant 5 has no ssveps


pl.data_evo = []; pl.data_evo_bc = [];
pl.RDKlabel = {'RDK1';'RDK2';'RDK3'};
pl.RDKidx = [1 2 3];

pl.con2plot = {'stim_cued';'top_down_cued';'not attended'}; 
pl.concols = num2cell([255 133 4; 25 138 131; 41 60 74;]'./255,1);

% extract data
for i_sub = 1:numel(pl.sub2plot)
    t.fidx = cell2mat(arrayfun(@(x) dsearchn(TFA.frequency', x+pl.freqrange'),[TFA.RDK(pl.sub2plot(i_sub)).RDK(pl.RDKidx).freq],'UniformOutput',false));
    pl.freqlabel = TFA.frequency(t.fidx(1,1):t.fidx(2,1))-TFA.RDK(pl.sub2plot(i_sub)).RDK(1).freq;
    for i_RDK = 1:size(t.fidx,2)
        for i_con = 1:numel(pl.con2plot)
            % which condition
            t.cidx = strcmp(F.conRDKprimed_label(:,i_RDK),pl.con2plot{i_con});
            % raw
            pl.data_evo(:,i_con,i_RDK,i_sub) = ...
                squeeze(mean(TFA.data_evo(t.fidx(1,i_RDK):t.fidx(2,i_RDK),:,pl.elec2plot_i,t.cidx ,pl.sub2plot(i_sub)),[3,1]));
            % baseline corrected
            pl.data_evo_bc(:,i_con,i_RDK,i_sub) = ...
                100*(...
                bsxfun(@rdivide, pl.data_evo(:,i_con,i_RDK,i_sub), ...
                squeeze(mean(TFA.data_evo(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.base_i(1):pl.base_i(2),pl.elec2plot_i,t.cidx ,pl.sub2plot(i_sub)),[3,1,2]))')...
                -1);
        end
    end
end

% figure; plot(mean(pl.data_evo,[3,4]))

% actual plotting
pl.data = squeeze(mean(pl.data_evo_bc,[3]));

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
t.diffs = [1 2; 1 3; 2 3];
for i_diff = 1:size(t.diffs,1)
    t.data = squeeze(diff(squeeze(pl.data(t.time_rt_i(1):t.time_rt_i(2),t.diffs(i_diff,:),:)),1,2));
    t.nulldata = repmat(diff(squeeze(mean(pl.data(pl.base_i(1):pl.base_i(2),t.diffs(i_diff,:),:),1))),[numel(t.time_rt_i(1):t.time_rt_i(2)),1]);
    [cluster_runt{size(pl.data,2)+i_diff}, timecourse_runt{size(pl.data,2)+i_diff}]=eeg_erpStat_clusterP(t.data,t.nulldata,t.permut_n,2);
end

pl.mdata = mean(pl.data,3); % mean data
pl.semdata = std(pl.data,1,3)./sqrt(numel(pl.sub2plot));

pl.conlabel = pl.con2plot;
pl.col = pl.concols;
pl.col2 = [0.6 0.6 0.6];
pl.line = {'-';'-';'-'};
figure('Position',[100 100 800 500]);
subplot(8,1,[1:5])
h.pl = {}; h.plsem=[];  h.plm = []; h.pls = []; h.plst = []; h.plbf = [];
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
ylabel('modulation in %')

plot([pl.base;pl.base],get(gca,"YLim"),'k')

% plot lines for significant effects
pl.sign_y = 1:size(timecourse_runt,2);

subplot(8,1,[6:8])

for i_con = 1:numel(pl.conlabel)
        
    % % uncorrected
    % pl.sigdata = timecourse_runt{ i_con}.h_raw.*pl.sign_y(i_con);
    % pl.sigdata(timecourse_runt{ i_con}.h_raw==0)=nan;
    % 
    % h.pls{i_con}=plot(TFA.time(t.time_rt_i(1):t.time_rt_i(2)), pl.sigdata,...
    %     'Color',pl.col2,'LineWidth',6);
    % hold on
    
    % corrected
    pl.sigdata = timecourse_runt{ i_con}.h_corr.*pl.sign_y(i_con);
    pl.sigdata(timecourse_runt{ i_con}.h_corr==0)=nan;
    
    h.pls{i_con}=plot(TFA.time(t.time_rt_i(1):t.time_rt_i(2)), pl.sigdata,...
        'Color',pl.col{i_con},'LineWidth',6);
    hold on
    
    % add text
    text(pl.xlims(1)+diff(pl.xlims)*0.01,i_con,[pl.conlabel{i_con} ' vs 0'], ...
        'FontSize',8,'Interpreter','none')
    
end

for i_diff = 1:size(t.diffs,1)
    t.idx = numel(pl.conlabel)+ i_diff;

    % % uncorrected
    % pl.sigdata = timecourse_runt{t.idx}.h_raw.*pl.sign_y(t.idx);
    % pl.sigdata(timecourse_runt{ t.idx}.h_raw==0)=nan;
    % 
    % h.pls{t.idx}=plot(TFA.time(t.time_rt_i(1):t.time_rt_i(2)), pl.sigdata,...
    %     'Color',pl.col2,'LineWidth',6);
    % hold on
    
    % corrected
    pl.sigdata = timecourse_runt{t.idx}.h_corr.*pl.sign_y(t.idx);
    pl.sigdata(timecourse_runt{t.idx}.h_corr==0)=nan;
    
    h.pls{ t.idx}=plot(TFA.time(t.time_rt_i(1):t.time_rt_i(2)), pl.sigdata,...
        'Color',pl.col{t.diffs(i_diff,1)},'LineWidth',6);
    hold on
    h.pls{t.idx}=plot(TFA.time(t.time_rt_i(1):t.time_rt_i(2)), pl.sigdata,...
        'Color',pl.col{t.diffs(i_diff,2)},'LineWidth',3);

    % add text
    text(pl.xlims(1)+diff(pl.xlims)*0.01,t.idx,[pl.conlabel{t.diffs(i_diff,1)} ' vs ' pl.conlabel{t.diffs(i_diff,2)}], ...
        'FontSize',8,'Interpreter','none')
    
end
% 

legend([h.pls{1:numel(pl.conlabel)}],pl.conlabel,'Location','SouthOutside','Orientation','horizontal','Interpreter','none')
legend('boxoff')

ylim(pl.sign_y([1 end])+[-1 1])
xlim(pl.xlims)

xlabel('time in ms')



% save data for analysis in R
pl.idx = pl.xlims_i(1):pl.xlims_i(2);
R_Mat.time = repmat(TFA.time(pl.idx)',[1,numel(pl.conlabel),numel(pl.sub2plot)]);
R_Mat.con = repmat(pl.conlabel',[numel(pl.idx),1,numel(pl.sub2plot)]);
R_Mat.sub = permute(repmat(F.Subs2use(pl.sub2plot)',[1, numel(pl.idx),numel(pl.conlabel)]),[2 3 1]);
R_Mat.data = pl.data(pl.idx,:,:);
R_Mat.all = table(R_Mat.time(:),R_Mat.con(:),R_Mat.sub(:),R_Mat.data(:), ...
    'VariableNames',{'time';'condition';'participant';'SSVEP_mod'});

t.path = 'C:\Dropboxdata\Dropbox\work\R-statistics\experiments\ssvep_fshiftprime1of2\data_in';
% t.path = 'C:\Users\EEG\Documents\R\Christopher\analysis_R_ssvep_fshift_perirr\data_in';
t.datestr = datestr(now,'mm-dd-yyyy_HH-MM');
% write to textfile

% writetable(R_Mat.all,fullfile(t.path,sprintf('Gabor_TimeCourses_%s.csv',t.datestr)),'Delimiter',';')





%% running t-tests ans cluster correction | other contrasts
% plotting parameters
pl.elec2plot = {'P7';'P5';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P8';'P6';'PO4';'PO8';'P10';'O2';'I2'}; sav.chan_add = 'VisualLarge'; % suppression irrelevant study
% cluster analysis
% pl.elec2plot = {TFA.electrodes(1:64).labels}';
% pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) startsWith({TFA(1).electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

% summ
pl.contrasts = {
    {{'stim_cued','top_down_cued'};{'not attended'}},'Selectivity ((S+TD)-U)';
    {{'stim_cued','top_down_cued','not attended'};{[]}},'Total Activity (S+TD+U)';
    };


pl.freqrange=[-0.05 0.05];

pl.xlims=[-1000 1700]; % index time 2 plot
pl.xlims_i = dsearchn(TFA.time', pl.xlims');

pl.time_post = [0 1700];

pl.base = [-500 -TFA.Gabor_FWHM_time];
pl.base_i = dsearchn(TFA.time', pl.base');

pl.sub2plot = 1:numel(F.Subs2use);
pl.sub2plot(ismember(F.Subs2use,5))=[]; % participant 5 has no ssveps


pl.data_evo = [];  pl.data_evo_bc = [];
pl.RDKlabel = {'RDK1';'RDK2';'RDK3'};
pl.RDKidx = [1 2 3];

pl.con2plot = {'stim_cued';'top_down_cued';'not attended'}; 
% pl.concols = num2cell([255 133 4; 41 60 74; 25 138 131]'./255,1);
pl.concols = num2cell([255 40 40; 40 180 40; 180 40 180]'./255,1);


% extract data
for i_sub = 1:numel(pl.sub2plot)
    t.fidx = cell2mat(arrayfun(@(x) dsearchn(TFA.frequency', x+pl.freqrange'),[TFA.RDK(pl.sub2plot(i_sub)).RDK(pl.RDKidx).freq],'UniformOutput',false));
    pl.freqlabel = TFA.frequency(t.fidx(1,1):t.fidx(2,1))-TFA.RDK(pl.sub2plot(i_sub)).RDK(1).freq;
    for i_RDK = 1:size(t.fidx,2)
        for i_con = 1:numel(pl.con2plot)
            % which condition
            t.cidx = strcmp(F.conRDKprimed_label(:,i_RDK),pl.con2plot{i_con});
            % raw
            pl.data_evo(:,i_con,i_RDK,i_sub) = ...
                squeeze(mean(TFA.data_evo(t.fidx(1,i_RDK):t.fidx(2,i_RDK),:,pl.elec2plot_i,t.cidx ,pl.sub2plot(i_sub)),[3,1]));
            % baseline corrected
            pl.data_evo_bc(:,i_con,i_RDK,i_sub) = ...
                100*(...
                bsxfun(@rdivide, pl.data_evo(:,i_con,i_RDK,i_sub), ...
                squeeze(mean(TFA.data_evo(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.base_i(1):pl.base_i(2),pl.elec2plot_i,t.cidx ,pl.sub2plot(i_sub)),[3,1,2]))')...
                -1);
        end
    end
end

% figure; plot(mean(pl.data_evo,[3,4]))
pl.data = nan([size(pl.data_evo_bc,1), size(pl.contrasts,1) size(pl.data_evo_bc,4)]);

% actual calculation of contrasts with summation
for i_cont = 1:size(pl.contrasts,1)
    % calculate contrasts as specified
    % first data
    t.idx1 = ismember(pl.con2plot,pl.contrasts{i_cont,1}{1});
    pl.xdata = squeeze(sum(mean(pl.data_evo_bc(:,t.idx1,:,:),[3]),2));
    % index second data
    if ~isempty(pl.contrasts{i_cont,1}{2}{1})
        t.idx2 = ismember(pl.con2plot,pl.contrasts{i_cont,1}{2});
        pl.ydata = squeeze(sum(mean(pl.data_evo_bc(:,t.idx2,:,:),[3]),2));
        pl.data(:,i_cont,:) = pl.xdata - pl.ydata;
    else
        pl.data(:,i_cont,:) = pl.xdata;
    end
end

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



pl.mdata = mean(pl.data,3); % mean data
pl.semdata = std(pl.data,1,3)./sqrt(numel(pl.sub2plot));

pl.conlabel = pl.contrasts(:,2);
pl.col = pl.concols;
pl.col2 = [0.6 0.6 0.6];
pl.line = {'-';'-';'-'};
figure('Position',[100 100 800 500]);
subplot(7,1,[1:5])
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
ylabel('modulation in %')
plot([pl.base;pl.base],get(gca,"YLim"),'k')


% plot lines for significant effects
pl.sign_y = 1:size(timecourse_runt,2);

subplot(7,1,[6:7])

for i_con = 1:numel(pl.conlabel)
        
    % % uncorrected
    % pl.sigdata = timecourse_runt{ i_con}.h_raw.*pl.sign_y(i_con);
    % pl.sigdata(timecourse_runt{ i_con}.h_raw==0)=nan;
    % 
    % h.pls{i_con}=plot(TFA.time(t.time_rt_i(1):t.time_rt_i(2)), pl.sigdata,...
    %     'Color',pl.col2,'LineWidth',6);
    
    
    % corrected
    pl.sigdata = timecourse_runt{ i_con}.h_corr.*pl.sign_y(i_con);
    pl.sigdata(timecourse_runt{ i_con}.h_corr==0)=nan;
    
    h.pls{i_con}=plot(TFA.time(t.time_rt_i(1):t.time_rt_i(2)), pl.sigdata,...
        'Color',pl.col{i_con},'LineWidth',6);
    hold on
    % add text
    text(pl.xlims(1)+diff(pl.xlims)*0.01,i_con,[pl.conlabel{i_con} ' vs 0'],'FontSize',8)
    
end


grid on
xlabel('time in ms')

legend([h.pls{1:numel(pl.conlabel)}],pl.conlabel,'Location','SouthOutside','Orientation','horizontal')
legend('boxoff')

ylim(pl.sign_y([1 end])+[-1 1])
xlim(pl.xlims)




% exportgraphics(gcf,'figures/SSVEP_mod_timecourse_selectivity_summ.pdf','ContentType','vector')






