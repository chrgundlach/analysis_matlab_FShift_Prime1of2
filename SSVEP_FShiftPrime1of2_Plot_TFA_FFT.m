%% plot TFA images
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
% pl.elec2plot = {'P6';'P8';'P10';'PO4';'PO8';'O2';'I2';'POz';'Oz';'Iz';'O1'};
% pl.elec2plot = {'P5';'P7';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'O2'};
% pl.elec2plot = {'POz';'O1';'Oz';'O2';'Iz'};
% pl.elec2plot = {'P8';'P10';'PO8';'PO3';'POz';'Oz';'O1'};
% pl.elec2plot = {'P7';'P9';'PO7';'PO2';'POz';'Oz';'O2'};


figure;
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmp({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
% topoplot(find(pl.elec2plot_i),TFA.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on');
topoplot([],TFA.electrodes(1:64),'whitebk','on','style','blank','electrodes', 'on', 'emarker2', {find(pl.elec2plot_i),'o','r',8});

%% plot grand mean FFT data | spectra | for distinct frequencies (lookup of respective electrode cluster)

% large center as in tango | periphery: central and lateral 
pl.elec2plot = {{'P5';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'PO4';'PO8';'O2';'I2'}, 'center'};
pl.elec2plot = {{'P7';'P5';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P8';'P6';'PO4';'PO8';'P10';'O2';'I2'} 'center'};

pl.elec2plot_i=cellfun(@(y) ...
    logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), y, 'UniformOutput',false)),1)),...
    pl.elec2plot(:,1), 'UniformOutput', false);


% pl.time2plot = [1];
pl.time2plot = [1];
pl.sub2plot = 1:numel(F.Subs2use);
pl.sub2plot(ismember(F.Subs2use,5))=[]; % participant 5 has no ssveps


% extract data
pl.data_ind = squeeze(mean(TFA.fftdata_ind(:,pl.elec2plot_i{1},:,pl.time2plot,pl.sub2plot),[2,3,4]));
pl.data_evo = squeeze(mean(TFA.fftdata_evo(:,pl.elec2plot_i{1},:,pl.time2plot,pl.sub2plot),[2,3,4]));


% plotting
figure;
set(gcf,'Position',[100 100 600 600],'PaperPositionMode','auto')
subplot(2,1,1)
plot(TFA.fftfreqs,pl.data_ind,'Color',[0.5 0.5 0.5],'LineWidth',1)
hold on;
plot(TFA.fftfreqs,mean(pl.data_ind,2),'Color','k','LineWidth',2)

xlim([0 50])
xlabel('frequency in Hz')
ylabel('amplitude in \muV/m²')
title(sprintf('induced GrandMean FFT spectra | N = %1.0f | FOI = %1.1f %1.1f %1.1f Hz', ...
    numel(pl.sub2plot), F.SSVEP_Freqs),'Interpreter','none')
vline(F.SSVEP_Freqs,'k:')
% draw topography with electrode positions
h.a1 = axes('position',[0.72 0.75 0.15 0.15],'Visible','off');
topoplot(find(any(cell2mat(pl.elec2plot_i))),TFA(1).electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(any(cell2mat(pl.elec2plot_i),1)),'o','r',4,1});

subplot(2,1,2)
hold on;
plot(TFA.fftfreqs,pl.data_evo,'Color',[0.5 0.5 0.5],'LineWidth',1)
plot(TFA.fftfreqs,mean(pl.data_evo,2),'Color','k','LineWidth',2)
xlim([0 50])
xlabel('frequency in Hz')
ylabel('amplitude in \muV/m²')
title(sprintf('evoked GrandMean FFT spectra | N = %1.0f | FOI = %1.1f %1.1f %1.1f Hz', ...
    numel(pl.sub2plot), F.SSVEP_Freqs),'Interpreter','none')
vline(F.SSVEP_Freqs,'k:')
box on

% draw topography with electrode positions
h.a1 = axes('position',[0.72 0.28 0.15 0.15],'Visible','off');
topoplot(find(any(cell2mat(pl.elec2plot_i))),TFA(1).electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(any(cell2mat(pl.elec2plot_i),1)),'o','r',4,1});



%% plot Grand Mean FFT data | topoplot for different frequencies
pl.time2plot = [1:3];
pl.time2plot = [1];
pl.freq2plot = F.SSVEP_Freqs(3);
pl.freqrange=[-0.1 0.1];
pl.sub2sel = 1:numel(F.Subs2use);


t.idx = dsearchn(TFA.fftfreqs',(pl.freqrange+pl.freq2plot)');
pl.data_ind = squeeze(mean(TFA.fftdata_ind(t.idx(1):t.idx(2),:,:,pl.time2plot,pl.sub2sel),[1 3 4 ]));
pl.data_evo = squeeze(mean(TFA.fftdata_evo(t.idx(1):t.idx(2),:,:,pl.time2plot,pl.sub2sel),[1 3 4 ]));


figure;
set(gcf,'Position',[100 100 800 300],'PaperPositionMode','auto')

% induced
h.s(1) = subplot(1,2,1);
pl.mdata = mean(pl.data_ind,2,'omitnan');
pl.clim = [0 max(pl.mdata)];
topoplot(pl.mdata, TFA.electrodes(1:64), ...
    'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clim,'conv','on','colormap',fake_parula,...
    'whitebk','on');
title(sprintf('induced SSVEP amps for %1.1f +- [%1.1f  %1.1f] Hz | [%1.0f %1.0f]ms', ...
    pl.freq2plot, pl.freqrange, min([TFA.ffttimewin{pl.time2plot}])*1000, max([TFA.ffttimewin{pl.time2plot}])*1000))
colorbar

h.s(2) = subplot(1,2,2);
pl.mdata = mean(pl.data_evo,2,'omitnan');
pl.clim = [0 max(pl.mdata)];
topoplot(pl.mdata, TFA.electrodes(1:64), ...
    'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clim,'conv','on','colormap',fake_parula,...
    'whitebk','on');
title(sprintf('evoked SSVEP amps for %1.1f +- [%1.1f  %1.1f] Hz | [%1.0f %1.0f]ms', ...
    pl.freq2plot, pl.freqrange, min([TFA.ffttimewin{pl.time2plot}])*1000, max([TFA.ffttimewin{pl.time2plot}])*1000))
colorbar



%% plot Grand Mean FFT data | topoplot for all SSVEP frequencies
pl.time2plot = [1:3];
pl.time2plot = [1];
pl.freq2plot = F.SSVEP_Freqs;
pl.freqrange=[-0.1 0.1];
pl.sub2sel = 1:numel(F.Subs2use);


t.idx = arrayfun(@(x) dsearchn(TFA.fftfreqs', (pl.freqrange+x)'), pl.freq2plot, 'UniformOutput', false);

% extract data
pl.data_evo = [];
for i_freq = 1:numel(t.idx)
    pl.data_evo(:,i_freq) = squeeze(mean(TFA.fftdata_evo(t.idx{i_freq}(1):t.idx{i_freq}(2),:,:,pl.time2plot,pl.sub2sel),[1 3 4 5]))';
end
pl.data_evo(:,end+1)=mean(pl.data_evo,2);


figure;
set(gcf,'Position',[100 100 1100 300],'PaperPositionMode','auto')

h = [];
for i_freq = 1:size(pl.data_evo,2)
    h.s(i_freq) = subplot(1,size(pl.data_evo,2),i_freq);
    pl.clim = [0 max(pl.data_evo,[],"all")];
    pl.clim = [0 max(pl.data_evo(:,i_freq),[],"all")];
    topoplot(pl.data_evo(:,i_freq), TFA.electrodes(1:64), ...
        'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clim,'conv','on','colormap',fake_parula,...
        'whitebk','on');
    if i_freq < size(pl.data_evo,2) 
        title(sprintf('evo SSVEP %1.1f +- [%1.1f  %1.1f] Hz\n[%1.0f %1.0f]ms', ...
            pl.freq2plot(i_freq), pl.freqrange, min([TFA.ffttimewin{pl.time2plot}])*1000, max([TFA.ffttimewin{pl.time2plot}])*1000))
    else
        title(sprintf('evo SSVEPs | freq averaged\n[%1.0f %1.0f]ms', ...
             min([TFA.ffttimewin{pl.time2plot}])*1000, max([TFA.ffttimewin{pl.time2plot}])*1000))
    end
    colorbar
end


%% plot FFT data modulation | topoplot effects on central stimuli
pl.time_base = 1;
pl.time2plot = [3];
pl.pos2plot='center';
pl.freqrange=[-0.1 0.1];
pl.sub2sel = 1:numel(F.Subs2use);
pl.con2plot = {'primed';'nonprimed';'not attended'}; %F.conRDKprimed_label(1,:)

pl.data_ind = []; pl.data_evo = [];
t.data_ind = nan(numel(TFA.electrodes),numel(TFA.RDK(1).RDK),numel(pl.con2plot),numel(TFA.ffttimewin),numel(pl.sub2sel));
t.data_evo = t.data_ind;

for i_sub = 1:numel(pl.sub2sel)
    for i_rdk = 1:numel(TFA.RDK(1).RDK)
        % which SSVEP frequency?
        t.freq = TFA.RDK(pl.sub2sel(i_sub)).RDK(i_rdk).freq;
        t.fidx = dsearchn(TFA.fftfreqs',(pl.freqrange+t.freq)');
        for i_con = 1:numel(pl.con2plot)
            % which condition
            t.cidx = strcmp(F.conRDKprimed_label(i_rdk,:),pl.con2plot{i_con});
            % extract data
            t.data_ind(:,i_rdk,i_con,:,i_sub) = squeeze(mean(TFA.fftdata_ind(t.fidx(1):t.fidx(2),:, t.cidx,:,pl.sub2sel(i_sub)),[1,3]));
            t.data_evo(:,i_rdk,i_con,:,i_sub) = squeeze(mean(TFA.fftdata_evo(t.fidx(1):t.fidx(2),:, t.cidx,:,pl.sub2sel(i_sub)),[1,3]));
        end
    end
end

% collapse across RDKs
t.data_ind_coll = squeeze(mean(t.data_ind,2));
t.data_evo_coll = squeeze(mean(t.data_evo,2));

% baseline corrected data
t.data_ind_coll_bc = ((mean(t.data_ind(:,:,:,pl.time2plot,:),4) ./ t.data_ind(:,:,:,pl.time_base,:))-1)*100;
t.data_evo_coll_bc = ((mean(t.data_evo(:,:,:,pl.time2plot,:),4) ./ t.data_evo(:,:,:,pl.time_base,:))-1)*100;
% t.data_ind_coll_bc = mean(t.data_ind(:,:,:,pl.time2plot,:),4) - t.data_ind(:,:,:,pl.time_base,:);
% t.data_evo_coll_bc = mean(t.data_evo(:,:,:,pl.time2plot,:),4) - t.data_evo(:,:,:,pl.time_base,:);

% collapse across RDKs
t.data_ind_coll_bc = squeeze(mean(t.data_ind_coll_bc,2));
t.data_evo_coll_bc = squeeze(mean(t.data_evo_coll_bc,2));

t.data_ind_coll_bc_m = mean(t.data_ind_coll_bc,3);
t.data_evo_coll_bc_m = mean(t.data_evo_coll_bc,3);

% do ttests
[t.data_ind_coll_ttp t.data_evo_coll_ttp] = deal(nan(size(t.data_ind_coll_bc_m)));

for i_con = 1:size(t.data_ind_coll_ttp,2)
    [tt.h t.data_ind_coll_ttp(:,i_con) tt.ci tt.stats] = ttest(squeeze(t.data_ind_coll_bc(:,i_con,:))');
    [tt.h t.data_evo_coll_ttp(:,i_con) tt.ci tt.stats] = ttest(squeeze(t.data_evo_coll_bc(:,i_con,:))');
end

% evoked
figure;
set(gcf,'Position',[100 100 800 500],'PaperPositionMode','auto')
pl.conlabel = pl.con2plot;

% modulations
for i_con = 1:size(t.data_evo_coll_bc_m,2)
    h.s(i_con) = subplot(2,size(t.data_evo_coll_bc_m,2),i_con);
    pl.clim = [-1 1] *max(abs(t.data_evo_coll_bc_m),[],'all');
    % pl.clim = [-20 20];
    topoplot(t.data_evo_coll_bc_m(:,i_con), TFA.electrodes(1:64), ...
        'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clim,'conv','off','colormap',flipud(cbrewer2('RdBu')),...
        'whitebk','on');
    title(sprintf('evoked SSVEP mod | %s\n[%1.0f %1.0f]ms', ...
        pl.conlabel{i_con}, [min([TFA.ffttimewin{pl.time2plot}]) max([TFA.ffttimewin{pl.time2plot}])]*1000))
    colorbar
end

% uncorrected t-values
for i_con = 1:size(t.data_evo_coll_ttp,2)
    t.data = abs(log10(t.data_evo_coll_ttp(:,i_con)));
    t.clim = [0 max(t.data(:))];
    t.pcriterion = abs(log10(0.05));
    if max(t.data(:))<t.pcriterion
        % temp.colormap = repmat([0.5 0.5 0.5],100,1);
        t.colormap = repmat(linspace(1,0.3,1000)',1,3);
    else
        t.border = ceil((t.pcriterion / max(t.data(:)))*1000);
        % temp.colormap = [repmat([0.5 0.5 0.5],temp.border,1); [linspace(0,1,100-temp.border)' repmat(0,100-temp.border,1) linspace(1,0,100-temp.border)']];
        t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
    end

    h.s2(i_con)=subplot(2,size(t.data_evo_coll_ttp,2),i_con+size(t.data_ind_coll_ttp,2));
    topoplot( t.data, TFA.electrodes(1:64), ...
        'shading', 'interp', 'numcontour', 0, 'maplimits',t.clim,'conv','on','colormap',t.colormap,'whitebk','on');
    title(sprintf('respective t-values'))
    h.cb2 = colorbar;
    t.yticks = get(h.cb2,'YTick');
    set(h.cb2,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.clim(end),1,'last')), ...
        'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.clim(end),1,'last')))

end







%% extract amplitude values for FFT
% plotting parameters
pl.elec2plot = {{'P5';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'PO4';'PO8';'O2';'I2'}, 'center'};
pl.elec2plot = {{'P7';'P5';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P8';'P6';'PO4';'PO8';'P10';'O2';'I2'}, 'center'};

pl.elec2plot_i=cellfun(@(y) ...
    logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), y, 'UniformOutput',false)),1)),...
    pl.elec2plot(:,1), 'UniformOutput', false);


pl.time2plot = [1:3];
pl.time2baseline = [1];
pl.freqrange=[-0.1 0.1];
% pl.freqrange=[0 0];
pl.sub2plot = 1:numel(F.Subs2use);

pl.data_ind = []; pl.data_evo = [];
pl.RDKlabel = {'RDK1';'RDK2';'RDK3'};
pl.RDKposlabel2 = {'center';'center';'center'};
pl.timelabel = cellfun(@(x) vararg2str(x),TFA.ffttimewin(pl.time2plot),'UniformOutput',false);

% extract data
% ([RDK num], con, time, sub)
pl.RDK.data_ind = nan(numel(TFA.RDK(1).RDK),numel(F.conlabel_att),numel(pl.time2plot),numel(pl.sub2plot));
pl.RDK.data_evo = pl.RDK.data_ind;


% the ones that are fixed
pl.RDK.con = permute(repmat((1:numel(F.conlabel_att))',[1 size(pl.RDK.data_ind,[1 3 4])]),[2 1 3 4]);
pl.RDK.conlabel = permute(repmat((F.conlabel_att),[1 size(pl.RDK.data_ind,[1 3 4])]),[2 1 3 4]);
pl.RDK.RDK_id = repmat(pl.RDKlabel,[1 size(pl.RDK.data_ind, [2 3 4 ])]);
pl.RDK.RDK_isattended = permute(repmat(F.conRDKattended_label,[1 1 size(pl.RDK.data_ind,[3 4])]),[2 1 3 4]);
pl.RDK.RDK_isprimed = permute(repmat(F.conRDKprimed_label,[1 1 size(pl.RDK.data_ind,[3 4])]),[2 1 3 4]);
pl.RDK.timewin = permute(repmat(pl.timelabel,[1 size(pl.RDK.data_ind,[1 2 4])]), [2 3 1 4]);
pl.RDK.sub = permute(repmat(F.Subs2use(pl.sub2plot)',[1 size(pl.RDK.data_ind,[1 2 3])]),[2 3 4 1]);
pl.RDK.RDK_pos2 = repmat(pl.RDKposlabel2,[1 size(pl.RDK.data_ind, [2 3 4 ])]);

% the onses that change
pl.RDK.RDK_freq = pl.RDK.data_ind;
pl.RDK.RDK_color = repmat({''},size(pl.RDK.data_ind));
pl.RDK.RDK_electrodes = repmat({''},size(pl.RDK.data_ind));

for i_sub = 1:numel(pl.sub2plot)
    for i_rdk = 1:numel(TFA.RDK(1).RDK)
        t.freq = TFA.RDK(pl.sub2plot(i_sub)).RDK(i_rdk).freq; % which SSVEP frequency?
        t.fidx = dsearchn(TFA.fftfreqs',(pl.freqrange+t.freq)');
       
        % extract data
        pl.RDK.data_ind(i_rdk,:,:,i_sub) = squeeze(mean( ...
            TFA.fftdata_ind(t.fidx(1):t.fidx(2),pl.elec2plot_i{1},:,pl.time2plot,pl.sub2plot(i_sub)),[1,2] ... %induced
            ));
        pl.RDK.data_evo(i_rdk,:,:,i_sub) = squeeze(mean( ...
            TFA.fftdata_evo(t.fidx(1):t.fidx(2),pl.elec2plot_i{1},:,pl.time2plot,pl.sub2plot(i_sub)),[1,2] ... %evoked
            ));

        % write some data
        pl.RDK.RDK_freq(i_rdk,:,:,i_sub) = t.freq;
        pl.RDK.RDK_color(i_rdk,:,:,i_sub) = {TFA.RDK(pl.sub2plot(i_sub)).RDK(i_rdk).col_label};
        pl.RDK.RDK_electrodes(i_rdk,:,:,i_sub) = {vararg2str(pl.elec2plot(1,1))};
    end
end

% baseline corrected data
pl.RDK.data_ind_bc = 100*(bsxfun(@rdivide, pl.RDK.data_ind, pl.RDK.data_ind(:,:,1,:))-1);
pl.RDK.data_evo_bc = 100*(bsxfun(@rdivide, pl.RDK.data_evo, pl.RDK.data_evo(:,:,1,:))-1);
pl.RDK.data_ind_bc_sub = bsxfun(@minus,  pl.RDK.data_ind, pl.RDK.data_ind(:,:,1,:));
pl.RDK.data_evo_bc_sub = bsxfun(@minus, pl.RDK.data_evo, pl.RDK.data_evo(:,:,1,:));

% % trouble shooting do subtraction and modulation correspond?
% % pl.tdata = reshape(squeeze(pl.RDK.data_evo_bc(1:2,:,2,:)),[],size(pl.RDK.data_evo_bc,4));
% % pl.tdata(:,:,2) = reshape(squeeze(pl.RDK.data_evo_bc_sub(1:2,:,2,:)),[],size(pl.RDK.data_evo_bc,4)).*10;
% pl.tdata = reshape(squeeze(mean(pl.RDK.data_evo_bc(1:2,:,2,:),1)),[],size(pl.RDK.data_evo_bc,4));
% pl.tdata(:,:,2) = reshape(squeeze(mean(pl.RDK.data_evo_bc_sub(1:2,:,2,:),1)),[],size(pl.RDK.data_evo_bc,4)).*10;
% figure;
% for i_sub = 1:size(pl.tdata,2)
% %     plot([-0.25 +0.25]+i_sub, squeeze(pl.tdata(:,i_sub,:)),'Color',[0.3 0.3 0.3 0.5])
%     plot([-0.25 +0.25]+i_sub, sign(squeeze(pl.tdata(:,i_sub,:))),'Color',[0.3 0.3 0.3 0.5])
%     hold on
% end
% 
% pl.tdata = reshape(squeeze(diff(sign(pl.RDK.data_evo_bc(1:2,:,2,:)),1,1)),[],size(pl.RDK.data_evo_bc,4));
% pl.tdata(:,:,2) = reshape(squeeze(diff(sign(pl.RDK.data_evo_bc_sub(1:2,:,2,:)),1,1)),[],size(pl.RDK.data_evo_bc,4));
% figure;
% for i_sub = 1:size(pl.tdata,2)
%     plot([-0.25 +0.25]+i_sub, ...
%         squeeze(pl.tdata(:,i_sub,:))+repmat(((randn(size(pl.tdata,1),1)-0.5).*0.1),1,2), ...
%         'Color',[0.3 0.3 0.3 0.5])
% %     plot([-0.25 +0.25]+i_sub, sign(squeeze(pl.tdata(:,i_sub,:))),'Color',[0.3 0.3 0.3 0.5])
%     hold on
% end
% figure; plot(sign(pl.RDK.data_evo_bc_sub(:)), sign(pl.RDK.data_evo_bc(:)))
% figure; plot(sign(R_Mat.all_table.modulation_evoked), sign(R_Mat.all_table.subtraction_evoked))


% % do some interim plotting for checking everything
% t.tdata = cat(3, ...
%     [squeeze(mean(pl.RDK.data_evo(1,[1 3 5],2,:),2)) squeeze(mean(pl.RDK.data_evo(1,[2 4 6],2,:),2))], ...
%     [squeeze(mean(pl.RDK.data_evo(2,[2 4 6],2,:),2)) squeeze(mean(pl.RDK.data_evo(2,[1 3 5],2,:),2))]);
% figure; boxplot(mean(t.tdata,3))
% t.tidx = strcmp(pl.RDK.RDK_pos1,'center') & strcmp(pl.RDK.timewin,'[0.5 1.5] ') & pl.RDK.RDK_freq == 29 ;
% t.tdata_ev0 =  pl.RDK.data_evo_bc;  t.tdata_ev0(~t.tidx)= nan;
% t.tdata = cat(3, ...
%     [squeeze(mean(t.tdata_ev0(1,[1 3 5],2,:),2,'omitnan')) squeeze(mean(t.tdata_ev0(1,[2 4 6],2,:),2,'omitnan'))], ...
%     [squeeze(mean(t.tdata_ev0(2,[2 4 6],2,:),2,'omitnan')) squeeze(mean(t.tdata_ev0(2,[1 3 5],2,:),2,'omitnan'))]);
% figure; boxplot(mean(t.tdata,3,'omitnan'))
% figure; plot(mean(t.tdata,3,'omitnan')')
% t.tdata2 =diff(mean(t.tdata,3,'omitnan'),1,2);
% figure; boxplot(diff(mean(t.tdata,3,'omitnan'),1,2));
% figure; histfit(diff(mean(t.tdata,3,'omitnan'),1,2),20)
% [tt.h,tt.p,tt.ci,tt.stats] = ttest(diff(mean(t.tdata,3,'omitnan'),1,2));

R_Mat.all = [{'amplitude_induced','amplitude_evoked','modulation_induced','modulation_evoked','subtraction_induced','subtraction_evoked', ...
    'subjects', 'condition', 'time', ...
    'RDK_id', 'RDK_position', 'RDK_freq', 'RDK_color', 'RDK_isattended', 'RDK_isprimed', 'RDK_electrodes'}; ...
    num2cell([pl.RDK.data_ind(:) pl.RDK.data_evo(:) pl.RDK.data_ind_bc(:) pl.RDK.data_evo_bc(:) pl.RDK.data_ind_bc_sub(:) pl.RDK.data_evo_bc_sub(:) ...
    pl.RDK.sub(:) pl.RDK.con(:)]) ...
    pl.RDK.timewin(:) pl.RDK.RDK_id(:) pl.RDK.RDK_pos2(:) num2cell(pl.RDK.RDK_freq(:)) pl.RDK.RDK_color(:) ...
    pl.RDK.RDK_isattended(:) pl.RDK.RDK_isprimed(:) pl.RDK.RDK_electrodes(:)
    ];

R_Mat.all_table = cell2table(R_Mat.all(2:end,:), "VariableNames",R_Mat.all(1,:));


t.path = 'C:\Dropboxdata\Dropbox\work\R-statistics\experiments\ssvep_fshiftprime1of2\data_in';
% t.path = 'C:\Users\EEG\Documents\R\Christopher\analysis_R_ssvep_fshift_perirr\data_in';
t.datestr = datestr(now,'mm-dd-yyyy_HH-MM');
% write to textfile

% writetable(R_Mat.all_table,fullfile(t.path,sprintf('FFT_Amp_data_CenterLarge_%s.csv',t.datestr)),'Delimiter',';')


% % test analysis
% mean(R_Mat.all_table.modulation_evoked(strcmp(pl.RDK.RDK_isprimed,'primed')& strcmp(pl.RDK.timewin,'[0.5 1.5] ')))
% mean(R_Mat.all_table.modulation_evoked(strcmp(pl.RDK.RDK_isprimed,'nonprimed')& strcmp(pl.RDK.timewin,'[0.5 1.5] ')))
% mean(R_Mat.all_table.modulation_evoked(strcmp(pl.RDK.RDK_isprimed,'not attended')& strcmp(pl.RDK.timewin,'[0.5 1.5] ')))
% 
% mean(R_Mat.all_table.modulation_evoked(strcmp(pl.RDK.RDK_isprimed,'primed')& strcmp(pl.RDK.timewin,'[1 2] ')))
% mean(R_Mat.all_table.modulation_evoked(strcmp(pl.RDK.RDK_isprimed,'nonprimed')& strcmp(pl.RDK.timewin,'[1 2] ')))
% mean(R_Mat.all_table.modulation_evoked(strcmp(pl.RDK.RDK_isprimed,'not attended')& strcmp(pl.RDK.timewin,'[1 2] ')))

%% actual plotting data | TFA Grand Mean timecourse
% plotting parameters
% pl.elec2plot = {'Oz';'Iz'};
% pl.elec2plot = {'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';% vis alpha II
pl.elec2plot = {'P5';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'PO4';'PO8';'O2';'I2'}; sav.chan_add = 'VisualLarge'; % as in tango study [to be used]
% cluster analysis
% pl.elec2plot = {TFA.electrodes(1:64).labels}';
% pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) startsWith({TFA(1).electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.flims= TFA.frequency([1 end]);
pl.flims_i=dsearchn(TFA.frequency', pl.flims');

pl.xlims=[-500 1800]; % index time 2 plot
pl.xlims_i = dsearchn(TFA.time', pl.xlims');

pl.base = F.TFA.baseline;
% pl.base = [-500 TFA.Gabor_FWHM_time];
pl.base_i = dsearchn(TFA.time', pl.base');

pl.sub2sel = 1:numel(F.Subs2use);

pl.data_ind = []; pl.data_evo = []; pl.data_ind_bc = []; pl.data_evo_bc = [];
pl.conlabel = F.conlabel_att;
pl.RDKlabel = {'RDK1';'RDK2';'RDK3'};
% pl.timelabel = cellfun(@(x) vararg2str(x),TFA.ffttimewin(pl.time2plot),'UniformOutput',false);


pl.data_ind = squeeze(mean(mean(TFA.data_ind(pl.flims_i(1):pl.flims_i(2),pl.xlims_i(1):pl.xlims_i(2),pl.elec2plot_i,:,pl.sub2sel),3),4));
pl.data_evo = squeeze(mean(mean(TFA.data_evo(pl.flims_i(1):pl.flims_i(2),pl.xlims_i(1):pl.xlims_i(2),pl.elec2plot_i,:,pl.sub2sel),3),4));
pl.data_ind_bc = 100*(bsxfun(@rdivide, pl.data_ind, ...
    mean(squeeze(mean(mean(TFA.data_ind(pl.flims_i(1):pl.flims_i(2),pl.base_i(1):pl.base_i(2),pl.elec2plot_i,:,pl.sub2sel),3),4)),2))-1);
pl.data_evo_bc = 100*(bsxfun(@rdivide, pl.data_evo, ...
    mean(squeeze(mean(mean(TFA.data_evo(pl.flims_i(1):pl.flims_i(2),pl.base_i(1):pl.base_i(2),pl.elec2plot_i,:,pl.sub2sel),3),4)),2))-1);

figure;
set(gcf,'Position',[100 100 800 600],'PaperPositionMode','auto')
subplot(2,1,1)
pl.data = mean(pl.data_ind,3); pl.clims1=[0 max(pl.data(:))];
imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),TFA.frequency(pl.flims_i(1):pl.flims_i(2)),pl.data,pl.clims1)
colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
set(gca,'YDir','normal')
title(sprintf('Grand Mean TFA-amplitude | raw | induced\n for channel [%s]', vararg2str(pl.elec2plot)), ...
    'FontSize',8,'Interpreter','none')
xlabel('time in ms')
ylabel('frequency in Hz')
set(gca,'FontSize',8)
% set(gca,'ColorScale','log')
cb = colorbar();

subplot(2,1,2)
pl.data = mean(pl.data_evo,3); pl.clims1=[0 max(pl.data(:))];
imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),TFA.frequency(pl.flims_i(1):pl.flims_i(2)),pl.data,pl.clims1)
colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
set(gca,'YDir','normal')
title(sprintf('Grand Mean TFA-amplitude | raw | evoked\n for channel [%s]',vararg2str(pl.elec2plot)), ...
    'FontSize',8,'Interpreter','none')
xlabel('time in ms')
ylabel('frequency in Hz')
set(gca,'FontSize',8)
% set(gca,'ColorScale','log')
cb = colorbar();


figure;
set(gcf,'Position',[100 100 800 600],'PaperPositionMode','auto')
subplot(2,1,1)
pl.data = mean(pl.data_ind_bc,3); pl.clims1=[-1 1]* max(abs(pl.data(:)));
imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),TFA.frequency(pl.flims_i(1):pl.flims_i(2)),pl.data,pl.clims1)
colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
set(gca,'YDir','normal')
title(sprintf('Grand Mean TFA-amplitude | baseline corrected | induced\n for channel [%s]',vararg2str(pl.elec2plot)), ...
    'FontSize',8,'Interpreter','none')
xlabel('time in ms')
ylabel('frequency in Hz')
set(gca,'FontSize',8)
% set(gca,'ColorScale','log')
cb = colorbar();

subplot(2,1,2)
pl.data = mean(pl.data_evo_bc,3); pl.clims1=[-1 1]* max(abs(pl.data(:)));
imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),TFA.frequency(pl.flims_i(1):pl.flims_i(2)),pl.data,pl.clims1)
colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
set(gca,'YDir','normal')
title(sprintf('Grand Mean TFA-amplitude | raw | evoked\n for channel [%s]',vararg2str(pl.elec2plot)), ...
    'FontSize',8,'Interpreter','none')
xlabel('time in ms')
ylabel('frequency in Hz')
set(gca,'FontSize',8)
% set(gca,'ColorScale','log')
cb = colorbar();


%% actual plotting data | TFA timecourse | central stimuli
% plotting parameters
% pl.elec2plot = {'Oz';'Iz'};
% pl.elec2plot = {'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';% vis alpha II
pl.elec2plot = {'P5';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'PO4';'PO8';'O2';'I2'}; sav.chan_add = 'VisualLarge'; % as in tango study [to be used]
% pl.elec2plot = {'P7';'P5';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P8';'P6';'PO4';'PO8';'O2';'I2'}; sav.chan_add = 'VisualLarge'; % as in tango study [to be used]
% cluster analysis
% pl.elec2plot = {TFA.electrodes(1:64).labels}';
% pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) startsWith({TFA(1).electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.freqrange=[-5 5];

pl.xlims=[-500 2000]; % index time 2 plot
pl.xlims=[-500 1800]; % index time 2 plot
pl.xlims_i = dsearchn(TFA.time', pl.xlims');

pl.base = F.TFA.baseline;
% pl.base = [-1000 -250];
% pl.base = [-500 TFA.Gabor_FWHM_time];

pl.base_i = dsearchn(TFA.time', pl.base');

pl.sub2plot = 1:numel(F.Subs2use);
pl.sub2plot(ismember(F.Subs2use, [5])) = []; %

pl.data_ind = []; pl.data_evo = []; pl.data_ind_bc = []; pl.data_evo_bc = [];
pl.RDKlabel = {'RDK1';'RDK2';'RDK3'};
pl.RDKidx = [1 2 3];

pl.con2plot = {'primed';'nonprimed';'not attended'}; %F.conRDKprimed_label(1,:)
% pl.con2plot = unique(F.conRDKprimed_label);


% extract data
for i_sub = 1:numel(pl.sub2plot)
    t.fidx = cell2mat(arrayfun(@(x) dsearchn(TFA.frequency', x+pl.freqrange'),[TFA.RDK(pl.sub2plot(i_sub)).RDK(pl.RDKidx).freq],'UniformOutput',false));
    pl.freqlabel = TFA.frequency(t.fidx(1,1):t.fidx(2,1))-TFA.RDK(pl.sub2plot(i_sub)).RDK(1).freq;
    for i_RDK = 1:size(t.fidx,2)
        for i_con = 1:numel(pl.con2plot)
            % which condition
            t.cidx = strcmp(F.conRDKprimed_label(:,i_RDK),pl.con2plot{i_con});
            % raw
            pl.data_ind(:,:,i_con,i_RDK,i_sub) = ...
                squeeze(mean(TFA.data_ind(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.xlims_i(1):pl.xlims_i(2),pl.elec2plot_i, t.cidx ,pl.sub2plot(i_sub)),[3 4]));
            pl.data_evo(:,:,i_con,i_RDK,i_sub) = ...
                squeeze(mean(TFA.data_evo(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.xlims_i(1):pl.xlims_i(2),pl.elec2plot_i, t.cidx ,pl.sub2plot(i_sub)),[3 4]));
            % baseline corrected
            pl.data_ind_bc(:,:,i_con,i_RDK,i_sub) = ...
                100*(bsxfun(@rdivide, pl.data_ind(:,:,i_con,i_RDK,i_sub), ...
                mean(squeeze(mean(TFA.data_ind(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.base_i(1):pl.base_i(2),pl.elec2plot_i, t.cidx ,pl.sub2plot(i_sub)),[3 4])),2))-1);
            pl.data_evo_bc(:,:,i_con,i_RDK,i_sub) = ...
                100*(bsxfun(@rdivide, pl.data_evo(:,:,i_con,i_RDK,i_sub), ...
                mean(squeeze(mean(TFA.data_evo(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.base_i(1):pl.base_i(2),pl.elec2plot_i,t.cidx ,pl.sub2plot(i_sub)),[3 4])),2))-1);
        end
    end   
end


% raw
pl.data = squeeze(mean(pl.data_ind,[4,5]));

pl.clims1=[0 1].*repmat(max(pl.data,[],"all"),numel(pl.con2plot),1);
figure;
set(gcf,'Position',[100 100 1500 600],'PaperPositionMode','auto')
for i_con = 1:size(pl.data,3)
    subplot(size(pl.data,3),2,1+2*(i_con-1))
    imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),pl.freqlabel,pl.data(:,:,i_con),pl.clims1(i_con,:))
    colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    % colormap(gca, "turbo") % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    
    set(gca,'YDir','normal')
    if i_con == 1
        title(sprintf('TFA-amplitude | raw | induced | %s\n for channel [%s]',pl.con2plot{i_con}, vararg2str(pl.elec2plot)), ...
            'FontSize',8,'Interpreter','none')
    end
    
    if i_con ~= 3
%         set(gca,'XTickLabel',[])
    else
        xlabel('time in ms')
    end
    ylabel({pl.con2plot{i_con};'delta f_R_D_K in Hz'})
    xlim(pl.xlims)
    set(gca,'FontSize',8)
    % set(gca,'ColorScale','log')
    cb = colorbar();
    
end


pl.data = squeeze(mean(pl.data_evo,[4,5]));
pl.clims1=[0 1].*repmat(max(pl.data,[],"all"),numel(pl.con2plot),1);
for i_con = 1:size(pl.data,3)
    subplot(size(pl.data,3),2,2+2*(i_con-1))
    imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),pl.freqlabel,pl.data(:,:,i_con),pl.clims1(i_con,:))
    colormap(gca, fake_parula) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    % colormap(gca, "turbo") % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    % colormap(gca, plasma) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    

    set(gca,'YDir','normal')
    if i_con == 1
        title(sprintf('TFA-amplitude | raw | corrected | evoked \n for channel [%s]', vararg2str(pl.elec2plot)), ...
            'FontSize',8,'Interpreter','none')
    end
    
    if i_con ~= 3
%         set(gca,'XTickLabel',[])
    else
        xlabel('time in ms')
    end
    ylabel({pl.con2plot{i_con};'delta f_R_D_K in Hz'})
    xlim(pl.xlims)
    set(gca,'FontSize',8)
    % set(gca,'ColorScale','log')
    cb = colorbar();
    
end


% draw topography with electrode positions
h.a1 = axes('position',[0.425 0.75 0.15 0.15],'Visible','off');
topoplot(find(pl.elec2plot_i),TFA(1).electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',4,1});



% baseline corrected
pl.data = squeeze(mean(pl.data_ind_bc,[4,5]));

pl.clims1=[-1 1].*repmat(max(abs(pl.data),[],"all"),numel(pl.con2plot),1);

figure;
set(gcf,'Position',[100 100 1500 600],'PaperPositionMode','auto')
for i_con = 1:size(pl.data,3)
    subplot(size(pl.data,3),2,1+2*(i_con-1))
    imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),pl.freqlabel,pl.data(:,:,i_con),pl.clims1(i_con,:))
    colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    set(gca,'YDir','normal')
    if i_con == 1
        title(sprintf('TFA-amplitude | baseline corrected | induced | %s\n for channel [%s]',pl.con2plot{i_con}, vararg2str(pl.elec2plot)), ...
            'FontSize',8,'Interpreter','none')
    end
    
    if i_con ~= 3
%         set(gca,'XTickLabel',[])
    else
        xlabel('time in ms')
    end
    ylabel({pl.con2plot{i_con};'delta f_R_D_K in Hz'})
    xlim(pl.xlims)
    set(gca,'FontSize',8)
    % set(gca,'ColorScale','log')
    cb = colorbar();
    
end

pl.data = squeeze(mean(pl.data_evo_bc,[4,5]));
pl.clims1=[-1 1].*repmat(max(abs(pl.data),[],"all"),numel(pl.con2plot),1);
for i_con = 1:size(pl.data,3)
    subplot(size(pl.data,3),2,2+2*(i_con-1))
    imagesc(TFA.time(pl.xlims_i(1):pl.xlims_i(2)),pl.freqlabel,pl.data(:,:,i_con),pl.clims1(i_con,:))
    colormap(gca, flipud(cbrewer2('RdBu'))) % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
    set(gca,'YDir','normal')
    if i_con == 1
        title(sprintf('TFA-amplitude | baseline corrected | evoked \n for channel [%s]', vararg2str(pl.elec2plot)), ...
            'FontSize',8,'Interpreter','none')
    end
    
    if i_con ~= 3
%         set(gca,'XTickLabel',[])
    else
        xlabel('time in ms')
    end
    ylabel({pl.con2plot{i_con};'frequency in Hz'})
    xlim(pl.xlims)
    set(gca,'FontSize',8)
    % set(gca,'ColorScale','log')
    cb = colorbar();
    
end


% draw topography with electrode positions
h.a1 = axes('position',[0.425 0.75 0.15 0.15],'Visible','off');
topoplot(find(pl.elec2plot_i),TFA(1).electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',4,1});



%% actual plotting data | TFA timecourse lineplot | central stimuli
% plotting parameters
% pl.elec2plot = {'Oz';'Iz'};
% pl.elec2plot = {'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';% vis alpha II
pl.elec2plot = {'P5';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'PO4';'PO8';'O2';'I2'}; sav.chan_add = 'VisualLarge'; % as in tango study [to be used]
pl.elec2plot = {'P7';'P5';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P8';'P6';'PO4';'PO8';'O2';'I2'}; sav.chan_add = 'VisualLarge'; % as in tango study [to be used]
% cluster analysis
% pl.elec2plot = {TFA.electrodes(1:64).labels}';
% pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) startsWith({TFA(1).electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.freqrange=[-0.05 0.05];

pl.xlims=[-1000 1800]; % index time 2 plot
pl.xlims_i = dsearchn(TFA.time', pl.xlims');

pl.base = F.TFA.baseline;
% pl.base = [-500 -0];
pl.base_i = dsearchn(TFA.time', pl.base');

pl.sub2plot = 1:numel(F.Subs2use);
pl.sub2plot(ismember(F.Subs2use,5))=[]; % participant 5 has no ssveps

pl.data_ind = []; pl.data_evo = []; pl.data_ind_bc = []; pl.data_evo_bc = [];
pl.RDKlabel = {'RDK1';'RDK2';'RDK3'};
pl.RDKidx = [1 2 3];

pl.con2plot = {'primed';'nonprimed';'not attended'}; %F.conRDKprimed_label(1,:)

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
                squeeze(mean(TFA.data_ind(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.xlims_i(1):pl.xlims_i(2),pl.elec2plot_i,t.cidx ,pl.sub2plot(i_sub)),[3,1]));
            pl.data_evo(:,i_con,i_RDK,i_sub) = ...
                squeeze(mean(TFA.data_evo(t.fidx(1,i_RDK):t.fidx(2,i_RDK),pl.xlims_i(1):pl.xlims_i(2),pl.elec2plot_i,t.cidx ,pl.sub2plot(i_sub)),[3,1]));
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

% figure; plot(mean(pl.data_evo,[3,4]))



% actual plotting
pl.data = squeeze(mean(pl.data_ind,[3]));


pl.cons = permute(repmat(pl.con2plot',[size(pl.data,1),1,size(pl.data,3)]),[1,2,3]); pl.cons = pl.cons(:);
pl.times = repmat(TFA.time(pl.xlims_i(1):pl.xlims_i(2))',[1 numel(pl.con2plot) numel(pl.sub2plot)]); pl.times = pl.times(:);
pl.data = pl.data(:);

clear g
g(1,1)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(1,1).stat_summary('type','sem');
g(1,1).set_names('x','time in ms','y','amplitude in \muV','color','RDK');
g(1,1).set_title('amplitudes of raw collapsed SSVEPs (induced analysis)');
% g(1,1).axe_property('YLim',[4 6]);


pl.data = squeeze(mean(pl.data_evo,[3]));
pl.data = pl.data(:);
g(2,1)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(2,1).stat_summary('type','sem');
g(2,1).set_names('x','time in ms','y','amplitude in \muV','color','RDK');
g(2,1).set_title('amplitudes of raw collapsed SSVEPs (evoked analysis)');
% g(2,1).axe_property('YLim',[1 2.5]);

pl.data = squeeze(mean(pl.data_ind_bc,[3]));
pl.data = pl.data(:);
g(1,2)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(1,2).stat_summary('type','sem');
g(1,2).set_names('x','time in ms','y','modulation in %','color','RDK');
g(1,2).set_title('modulations of collapsed SSVEPs from baseline (induced analysis)');
% g(1,2).axe_property('YLim',[-5 5]);

pl.data = squeeze(mean(pl.data_evo_bc,[3]));
pl.data = pl.data(:);
g(2,2)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(2,2).stat_summary('type','sem');
g(2,2).set_names('x','time in ms','y','amplitude in \muV','color','RDK');
g(2,2).set_title('modulations of collapsed SSVEPs from baseline (evoked analysis)');
% g(2,2).axe_property('YLim',[-20 40]);

% g.stat_boxplot();
g.set_text_options('base_size',8,'Interpreter','Tex');
g.axe_property('YGrid','on','Box','on');
figure('Position',[100 100 1600 700]);
g.draw();


clear g
pl.data = squeeze(mean(pl.data_evo,[3]));
pl.data = pl.data(:);
g(1,1)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(1,1).stat_summary('type','sem');
g(1,1).set_names('x','time in ms','y','amplitude in \muV','color','RDK');
g(1,1).set_title('amplitudes of raw collapsed SSVEPs (evoked analysis)');
% g(1,1).axe_property('YLim',[1 2.5]);

pl.data = squeeze(mean(pl.data_evo_bc,[3]));
pl.data = pl.data(:);
g(2,1)=gramm('x',pl.times,'y',pl.data,'color',pl.cons);
g(2,1).stat_summary('type','sem');
g(2,1).set_names('x','time in ms','y','amplitude in \muV/m²','color','RDK');
g(2,1).set_title('modulations of collapsed SSVEPs from baseline (evoked analysis)');
% g(2,1).axe_property('YLim',[-20 40]);

% g.stat_boxplot();
g.set_text_options('base_size',8,'Interpreter','Tex');
g.axe_property('YGrid','on','Box','on');
figure('Position',[100 100 700 700]);
g.draw();




%% calculate everything with running t-tests and cluster correction | central
% plotting parameters
% pl.elec2plot = {'Oz';'Iz'};
% pl.elec2plot = {'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';% vis alpha II
% pl.elec2plot = {'P5';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'PO4';'PO8';'O2';'I2'}; sav.chan_add = 'VisualLarge'; % as in tango study 
pl.elec2plot = {'P7';'P5';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P8';'P6';'PO4';'PO8';'P10';'O2';'I2'}; sav.chan_add = 'VisualLarge'; % suppression irrelevant study[to be used]
% cluster analysis
% pl.elec2plot = {TFA.electrodes(1:64).labels}';
% pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) startsWith({TFA(1).electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.freqrange=[-0.05 0.05];

pl.xlims=[-1000 1800]; % index time 2 plot
pl.xlims_i = dsearchn(TFA.time', pl.xlims');

pl.time_post = [0 1800];

pl.base = F.TFA.baseline;
% pl.base = [-750 -250];
pl.base = [-500 -TFA.Gabor_FWHM_time];
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

% Bayesian running t-tests against zero
bayesResults.Zero = nan(numel(t.time_rt_i(1):t.time_rt_i(2)),size(pl.data,2));
for i_con = 1:size(pl.data,2)
    t.data = squeeze(pl.data(t.time_rt_i(1):t.time_rt_i(2),i_con,:));
    for i_time = 1:size(t.data,1)
        % Bayesian analysis for each time point
        bayesResults.Zero(i_time,i_con) = bf.ttest(t.data(i_time,:)); % default cauchy prior [sqrt(2)/2]
    end
end

% run cluster correction for diffs against zero
t.diffs = [1 2; 1 3; 2 3];
for i_diff = 1:size(t.diffs,1)
    t.data = squeeze(diff(squeeze(pl.data(t.time_rt_i(1):t.time_rt_i(2),t.diffs(i_diff,:),:)),1,2));
    t.nulldata = repmat(diff(squeeze(mean(pl.data(pl.base_i(1):pl.base_i(2),t.diffs(i_diff,:),:),1))),[numel(t.time_rt_i(1):t.time_rt_i(2)),1]);
    [cluster_runt{size(pl.data,2)+i_diff}, timecourse_runt{size(pl.data,2)+i_diff}]=eeg_erpStat_clusterP(t.data,t.nulldata,t.permut_n,2);
end

% Bayesian running t-tests for diffs
bayesResults.Diffs = nan(numel(t.time_rt_i(1):t.time_rt_i(2)),size(t.diffs,1));
for i_diff = 1:size(t.diffs,1)
    t.data = squeeze(diff(squeeze(pl.data(t.time_rt_i(1):t.time_rt_i(2),t.diffs(i_diff,:),:)),1,2));
    for i_time = 1:size(t.data,1)
        % Bayesian analysis for each time point
        bayesResults.Diffs(i_time,i_diff) = bf.ttest(t.data(i_time,:)); % default cauchy prior [sqrt(2)/2]
    end
end


pl.mdata = mean(pl.data,3); % mean data
pl.semdata = std(pl.data,1,3)./sqrt(numel(pl.sub2plot));

pl.conlabel = pl.con2plot;
pl.col = pl.concols;
pl.col2 = [0.6 0.6 0.6];
pl.line = {'-';'-';'-'};
figure('Position',[100 100 800 700]);
subplot(11,1,[1:5])
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


% plot lines for significant effects
pl.sign_y = 1:size(timecourse_runt,2);

subplot(11,1,[6:8])

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
    
    % add text
    text(pl.xlims(1)+diff(pl.xlims)*0.01,i_con,[pl.conlabel{i_con} ' vs 0'],'FontSize',8)
    
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

    % add text
    text(pl.xlims(1)+diff(pl.xlims)*0.01,t.idx,[pl.conlabel{t.diffs(i_diff,1)} ' vs ' pl.conlabel{t.diffs(i_diff,2)}],'FontSize',8)
    
end
% 

legend([h.pls{1:numel(pl.conlabel)}],pl.conlabel,'Location','SouthOutside','Orientation','horizontal')
legend('boxoff')

ylim(pl.sign_y([1 end])+[-1 1])
xlim(pl.xlims)

xlabel('time in ms')


% plot lines running Bayesian tests
% pl.BFrange=[-inf -10;-10 -3;-3 0;0 3;3 10;10 inf];
pl.BFrange=[0 1/10;1/10 1/3;1/3 1;1 3;3 10;10 inf];
t.cols = colormap(gca, flipud(cbrewer2('RdBu'))); % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
t.cols([1 end],:) = [0.2 0.2 0.8; 1 0.2 0.2];
pl.BFcols = t.cols(linspace(1,size(t.cols,1),size(pl.BFrange,1)),:);
pl.BFlabels = {'H0++';'H0+';'H0~';'H1~';'H1+';'H1++'};
pl.sign_y = 1:(size(bayesResults.Zero,2)+size(bayesResults.Diffs,2));
 
subplot(11,1,[9:11])
for i_con = 1:numel(pl.conlabel)
    for i_BFrange = 1:size(pl.BFrange,1)
        pl.bfdata = nan(size(bayesResults.Zero,1),1);
        pl.bfdata(bayesResults.Zero(:,i_con)>=pl.BFrange(i_BFrange,1)&bayesResults.Zero(:,i_con)<=pl.BFrange(i_BFrange,2))=pl.sign_y(i_con);
        % plot data
        h.plbf{i_BFrange}=plot(TFA.time(t.time_rt_i(1):t.time_rt_i(2)), pl.bfdata,...
            'Color',pl.BFcols(i_BFrange,:),'LineWidth',6);
        hold on
    end
            
    % add text
    text(pl.xlims(1)+diff(pl.xlims)*0.01,i_con,[pl.conlabel{i_con} ' vs 0'],'FontSize',8)
    
end

for i_diff = 1:size(t.diffs,1)
    t.idx = numel(pl.conlabel)+ i_diff;

    for i_BFrange = 1:size(pl.BFrange,1)
        pl.bfdata = nan(size(bayesResults.Diffs,1),1);
        pl.bfdata(bayesResults.Diffs(:,i_diff)>=pl.BFrange(i_BFrange,1)&bayesResults.Diffs(:,i_diff)<=pl.BFrange(i_BFrange,2))=pl.sign_y(t.idx);
        % plot data
        h.plbf{i_BFrange}=plot(TFA.time(t.time_rt_i(1):t.time_rt_i(2)), pl.bfdata,...
            'Color',pl.BFcols(i_BFrange,:),'LineWidth',6);
        hold on
    end

    % add text
    text(pl.xlims(1)+diff(pl.xlims)*0.01,t.idx,[pl.conlabel{t.diffs(i_diff,1)} ' vs ' pl.conlabel{t.diffs(i_diff,2)}],'FontSize',8)
    
end
% 

legend([h.plbf{1:end}],pl.BFlabels,'Location','SouthOutside','Orientation','horizontal')
legend('boxoff')

ylim(pl.sign_y([1 end])+[-1 1])
xlim(pl.xlims)

xlabel('time in ms')

%% running t-tests ans cluster correction | other contrasts
% plotting parameters
% pl.elec2plot = {'Oz';'Iz'};
% pl.elec2plot = {'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';
% pl.elec2plot = {'P9';'P10';'PO7';'PO8';'PO3';'PO4';'POz';'O1';'O2';'Oz';'I1';'I2';'Iz'}; sav.chan_add = 'VisualLarge';% vis alpha II
% pl.elec2plot = {'P5';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P6';'PO4';'PO8';'O2';'I2'}; sav.chan_add = 'VisualLarge'; % as in tango study [to be used]
pl.elec2plot = {'P7';'P5';'P9';'PO3';'PO7';'O1';'I1';'POz';'Oz';'Iz';'P8';'P6';'PO4';'PO8';'P10';'O2';'I2'}; sav.chan_add = 'VisualLarge'; % suppression irrelevant study
% cluster analysis
% pl.elec2plot = {TFA.electrodes(1:64).labels}';
% pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) startsWith({TFA(1).electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));
pl.elec2plot_i=logical(sum(cell2mat(cellfun(@(x) strcmpi({TFA.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1));

pl.contrasts = {
    {{'primed','nonprimed'};{'not attended'}},'Selectivity ((P+NP)-U)';
    {{'primed'};{'nonprimed'}},'Prime (P-NP)';
    {{'primed','nonprimed','not attended'};{[]}},'Total Activity (P+NP+U)';
    };

pl.freqrange=[-0.05 0.05];

pl.xlims=[-1000 1800]; % index time 2 plot
pl.xlims_i = dsearchn(TFA.time', pl.xlims');

pl.time_post = [0 1800];

% pl.base = F.TFA.baseline;
% pl.base = [-750 -250];
pl.base = [-500 -TFA.Gabor_FWHM_time];
pl.base_i = dsearchn(TFA.time', pl.base');

pl.sub2plot = 1:numel(F.Subs2use);
pl.sub2plot(ismember(F.Subs2use,5))=[]; % participant 5 has no ssveps


pl.data_ind = []; pl.data_evo = []; pl.data_ind_bc = []; pl.data_evo_bc = [];
pl.RDKlabel = {'RDK1';'RDK2';'RDK3'};
pl.RDKidx = [1 2 3];

pl.con2plot = {'primed';'nonprimed';'not attended'}; %F.conRDKprimed_label(1,:)
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

% figure; plot(mean(pl.data_evo,[3,4]))
pl.data = nan([size(pl.data_evo_bc,1), size(pl.contrasts,1) size(pl.data_evo_bc,4)]);

% actual calculation of contrasts
for i_cont = 1:size(pl.contrasts,1)
    % calculate contrasts as specified
    % first data
    t.idx1 = ismember(pl.con2plot,pl.contrasts{i_cont,1}{1});
    pl.xdata = squeeze(mean(pl.data_evo_bc(:,t.idx1,:,:),[2 3]));
    % index second data
    if ~isempty(pl.contrasts{i_cont,1}{2}{1})
        t.idx2 = ismember(pl.con2plot,pl.contrasts{i_cont,1}{2});
        pl.ydata = squeeze(mean(pl.data_evo_bc(:,t.idx2,:,:),[2 3]));
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

% Bayesian running t-tests against zero
bayesResults.Zero = nan(numel(t.time_rt_i(1):t.time_rt_i(2)),size(pl.data,2));
for i_con = 1:size(pl.data,2)
    t.data = squeeze(pl.data(t.time_rt_i(1):t.time_rt_i(2),i_con,:));
    for i_time = 1:size(t.data,1)
        % Bayesian analysis for each time point
        bayesResults.Zero(i_time,i_con) = bf.ttest(t.data(i_time,:)); % default cauchy prior [sqrt(2)/2]
    end
end

% % run cluster correction for diffs against zero
% t.diffs = [1 2; 1 3; 2 3];
% for i_diff = 1:size(t.diffs,1)
%     t.data = squeeze(diff(squeeze(pl.data(t.time_rt_i(1):t.time_rt_i(2),t.diffs(i_diff,:),:)),1,2));
%     t.nulldata = repmat(diff(squeeze(mean(pl.data(pl.base_i(1):pl.base_i(2),t.diffs(i_diff,:),:),1))),[numel(t.time_rt_i(1):t.time_rt_i(2)),1]);
%     [cluster_runt{size(pl.data,2)+i_diff}, timecourse_runt{size(pl.data,2)+i_diff}]=eeg_erpStat_clusterP(t.data,t.nulldata,t.permut_n,2);
% end

pl.mdata = mean(pl.data,3); % mean data
pl.semdata = std(pl.data,1,3)./sqrt(numel(pl.sub2plot));

pl.conlabel = pl.contrasts(:,2);
pl.col = pl.concols;
pl.col2 = [0.6 0.6 0.6];
pl.line = {'-';'-';'-'};
figure('Position',[100 100 800 700]);
subplot(9,1,[1:5])
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

% plot lines for significant effects
pl.sign_y = 1:size(timecourse_runt,2);

subplot(9,1,[6:7])

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
    
    % add text
    text(pl.xlims(1)+diff(pl.xlims)*0.01,i_con,[pl.conlabel{i_con} ' vs 0'],'FontSize',8)
    
end

% for i_diff = 1:size(t.diffs,1)
%     t.idx = numel(pl.conlabel)+ i_diff;
% 
%     % uncorrected
%     pl.sigdata = timecourse_runt{t.idx}.h_raw.*pl.sign_y(t.idx);
%     pl.sigdata(timecourse_runt{ t.idx}.h_raw==0)=nan;
% 
%     h.pls{t.idx}=plot(TFA.time(t.time_rt_i(1):t.time_rt_i(2)), pl.sigdata,...
%         'Color',pl.col2,'LineWidth',6);
%     hold on
% 
%     % corrected
%     pl.sigdata = timecourse_runt{t.idx}.h_corr.*pl.sign_y(t.idx);
%     pl.sigdata(timecourse_runt{t.idx}.h_corr==0)=nan;
% 
%     h.pls{ t.idx}=plot(TFA.time(t.time_rt_i(1):t.time_rt_i(2)), pl.sigdata,...
%         'Color',pl.col{t.diffs(i_diff,1)},'LineWidth',6);
%     h.pls{t.idx}=plot(TFA.time(t.time_rt_i(1):t.time_rt_i(2)), pl.sigdata,...
%         'Color',pl.col{t.diffs(i_diff,2)},'LineWidth',3);
% 
%     % add text
%     text(pl.xlims(1)+diff(pl.xlims)*0.01,t.idx,[pl.conlabel{t.diffs(i_diff,1)} ' vs ' pl.conlabel{t.diffs(i_diff,2)}],'FontSize',8)
% 
% end
xlabel('time in ms')

legend([h.pls{1:numel(pl.conlabel)}],pl.conlabel,'Location','SouthOutside','Orientation','horizontal')
legend('boxoff')

ylim(pl.sign_y([1 end])+[-1 1])
xlim(pl.xlims)


% plot lines running Bayesian tests
% pl.BFrange=[-inf -10;-10 -3;-3 0;0 3;3 10;10 inf];
pl.BFrange=[0 1/10;1/10 1/3;1/3 1;1 3;3 10;10 inf];
t.cols = colormap(gca, flipud(cbrewer2('RdBu'))); % magma, viridis, plasma, parula, fake_parula, jet, inferno, cbrewer2('RdBu'),flipud(cbrewer2('RdBu'))
t.cols([1 end],:) = [0.2 0.2 0.8; 1 0.2 0.2];
pl.BFcols = t.cols(linspace(1,size(t.cols,1),size(pl.BFrange,1)),:);
pl.BFlabels = {'H0++';'H0+';'H0~';'H1~';'H1+';'H1++'};
pl.sign_y = 1:(size(bayesResults.Zero,2));
 
subplot(9,1,[8:9])
for i_con = 1:numel(pl.conlabel)
    for i_BFrange = 1:size(pl.BFrange,1)
        pl.bfdata = nan(size(bayesResults.Zero,1),1);
        pl.bfdata(bayesResults.Zero(:,i_con)>=pl.BFrange(i_BFrange,1)&bayesResults.Zero(:,i_con)<=pl.BFrange(i_BFrange,2))=pl.sign_y(i_con);
        % plot data
        h.plbf{i_BFrange}=plot(TFA.time(t.time_rt_i(1):t.time_rt_i(2)), pl.bfdata,...
            'Color',pl.BFcols(i_BFrange,:),'LineWidth',6);
        hold on
    end
            
    % add text
    text(pl.xlims(1)+diff(pl.xlims)*0.01,i_con,[pl.conlabel{i_con} ' vs 0'],'FontSize',8)
    
end

% for i_diff = 1:size(t.diffs,1)
%     t.idx = numel(pl.conlabel)+ i_diff;
% 
%     for i_BFrange = 1:size(pl.BFrange,1)
%         pl.bfdata = nan(size(bayesResults.Diffs,1),1);
%         pl.bfdata(bayesResults.Diffs(:,i_diff)>=pl.BFrange(i_BFrange,1)&bayesResults.Diffs(:,i_diff)<=pl.BFrange(i_BFrange,2))=pl.sign_y(t.idx);
%         % plot data
%         h.plbf{i_BFrange}=plot(TFA.time(t.time_rt_i(1):t.time_rt_i(2)), pl.bfdata,...
%             'Color',pl.BFcols(i_BFrange,:),'LineWidth',6);
%         hold on
%     end
% 
%     % add text
%     text(pl.xlims(1)+diff(pl.xlims)*0.01,t.idx,[pl.conlabel{t.diffs(i_diff,1)} ' vs ' pl.conlabel{t.diffs(i_diff,2)}],'FontSize',8)
% 
% end
% % 

legend([h.plbf{1:end}],pl.BFlabels,'Location','SouthOutside','Orientation','horizontal')
legend('boxoff')

ylim(pl.sign_y([1 end])+[-1 1])
xlim(pl.xlims)

xlabel('time in ms')