%% plot TFA images
clearvars
F.PathInEEG             = '\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\eeg\erp'; 
% F.PathInEEG             = 'N:\AllgPsy\experimental_data\2024_FShiftPerIrr\eeg\erp'; 

F.Subs                  = arrayfun(@(x) sprintf('%02.0f',x),1:60,'UniformOutput',false)';
F.Subs2use              = [1:13 15:41 43:52]; 
                        % 1 to 22
                        % for subject 12, 14, 39: eeg and behavior data don't match



F.conds_within=         {'target';'distractor'};
F.conds_between=        [repmat({'bckgrd_off'},1,numel([1:13 15:21])) repmat({'bckgrd_iso'},1,numel([22:24]))];

%% load data
for i_sub = 1:numel(F.Subs2use)
    fprintf(1,'|| file %1.0f out of %1.0f || %s\\VP%s_events_erp.mat ||\n',...
        i_sub,numel(F.Subs2use),F.PathInEEG,F.Subs{F.Subs2use(i_sub)})
    
    temp.erp = open(fullfile(F.PathInEEG,sprintf('VP%s_events_erp.mat',F.Subs{F.Subs2use(i_sub)})));
    
    % preallocate memory
    if i_sub == 1
        EP.filtdata = repmat({[]},1,numel(F.Subs2use));
        EP.behavior = repmat({[]},1,numel(F.Subs2use));
        EP.RDK = repmat({[]},1,numel(F.Subs2use));
        EP.params = temp.erp.EP.parameters;
        EP.electrodes = temp.erp.EP.EEG_fep.chanlocs;
        EP.time = temp.erp.EP.EEG_fep.times;
        EP.srate = temp.erp.EP.EEG_fep.srate;
    end
    
    % assign data
    EP.filtdata{i_sub} = temp.erp.EP.EEG_fep;
    EP.behavior{i_sub} = temp.erp.EP.behavior;
    EP.RDK{i_sub} = temp.erp.EP.RDK;
    
    clear temp
    
end

pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];



%% plot ERPs exploratively in topo array
    
% loop through conditions defined in contrasts
pl.con_contrast = {... % contrasts by 1st dim; averaged across second dim
    %'trial_timing_type', {{'regular'}}; ...
    'eventtype', {{'target'};{'distractor'}}};
pl.sub2plot = 1:numel(F.Subs2use);


% preallocate data
pl.dat2plot = nan([size(EP.electrodes,2), ...           % 1st dim channels
    size(EP.time,2), ...                                % 2nd dim time
    cellfun(@(x) size(x,1), pl.con_contrast(:,2))', ... % 3rd to n dim = levels defind in pl.con_contrast
    numel(pl.sub2plot)]);                               % n+1 dim participants

% extract data across contrasts
% define how to loop through contrasts
t.cont = cellfun(@(x) 1:size(x,1), pl.con_contrast(:,2), 'UniformOutput', false);
t.contidx = combvec(t.cont{:});
% loop across participants
for i_sub = 1:numel(pl.sub2plot)
    t.behavior = EP.behavior{pl.sub2plot(i_sub)};
    t.epdata = EP.filtdata{pl.sub2plot(i_sub)}.data;

    % loop through contrasts
    for i_cont = 1:size(t.contidx,2)
        % define evaluation syntax
                
        % preallocate logical idx for trials
        t.idx = false(size(t.contidx,1),numel(t.behavior ));
        for i_contstep = 1:size(t.contidx,1)
            % index step
            t.idx(i_contstep,:) = ...
                any( ... % if it is more than one level, combine
                cell2mat( ...
                cellfun(@(x) strcmp({t.behavior.(pl.con_contrast{i_contstep,1})},x), ...
                pl.con_contrast{i_contstep,2}{t.contidx(i_contstep,i_cont)}, 'UniformOutput', false)' ...
                ) ...
                ,1);
        end
        % combine all contrast conditions
        t.idx_all = all(t.idx,1);
        
        % average across trials
        t.text = sprintf("%1.0f,",t.contidx(:,i_cont));
        t.evaltext = sprintf("pl.dat2plot(:,:,%si_sub) = mean(t.epdata(:,:,t.idx_all),3);",t.text);
        eval(t.evaltext)
    end
end

% reshape pl.data
pl.dat2plot_rs = mean(squeeze(pl.dat2plot),4);

pl.con_label = [pl.con_contrast{1,2}{:}];

plot_data_topoarray(EP.electrodes, pl.dat2plot_rs,'ERP','times',EP.time,'conds',pl.con_label)


%% plot ERPs exploratively for specified electrodes | event type

% loop through conditions defined in contrasts
pl.con_contrast = {... % contrasts by 1st dim; averaged across second dim
    %'trial_timing_type', {{'regular'}}; ...
%     'event_response_type', {{'hit'}}; ...
%     'evnt_type_label', {{'chroma+'}}; ...
%     'event_response_type', {{'hit','FA','error','miss'}}; ...
    'eventtype', {{'target'};{'distractor'}}};
% pl.sub2plot = 1:numel(F.Subs2use);

pl.sub2plot=[]; pl.sublabel=[];
pl.sub2plot{1} = find(F.Subs2use<22); pl.sublabel{1} = 'luminance offset';% luminance offset
pl.sub2plot{2} = find(F.Subs2use>21); pl.sublabel{2} = 'isoluminant';% isoluminant to background


% pl.elec2plot = {'P8';'PO8';'P10';'P7';'PO7';'P9'}; % for P1 component lateral !
% pl.elec2plot = {'P7';'PO7';'P9';'O1';'I1';'Oz'; 'Iz';'O2';'I2';'P8';'PO8';'P10';}; % for N2 component posterior!
% pl.elec2plot = {'PO7';'O1';'I1';'Oz'; 'Iz';'O2';'I2';'PO8';'PO3';'POz';'PO4'}; % for N2 component posterior! more central
% pl.elec2plot = {'PO7';'O1';'I1';'Oz'; 'Iz';'O2';'I2';'PO8'}; % for N2 component posterior! more central
pl.elec2plot = {'P3';'P1';'Pz';'P4';'P2';'POz';'PO3';'PO4'}; % for P300 component centro-parietal!
% pl.elec2plot = {'POz';'Oz';'O1';'O2';'Iz'}; % for N1 component centro-parietal
% pl.elec2plot = {'P6';'P8';'PO8';'P10';'P5';'P7';'PO7';'P9'}; % for N2 component lateral
% pl.elec2plot = {'POz'}; % early N2 SN?
% pl.elec2plot = {'CPz';'Cz'}; % early N2 SN?
% pl.elec2plot = {'FCz';'Fz'}; % early frontal?

pl.ylim = [-6 25];

% pl.con_label = {'target'; 'distractor'};

pl.elec2plot_i = ...
    any(cell2mat(cellfun(@(x) strcmp({EP.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1);

pl.time2plot = [-100 500]; % time in ms
pl.time2plot = [-100 700]; % time in ms
pl.time2plot = [-100 900]; % time in ms
pl.time2plot_i = dsearchn(EP.time', pl.time2plot');

pl.concols = num2cell([25 138 130; 41 60 74; 241 131 26]'./255,1);
pl.concols = num2cell(hex2rgb({'#005b8a';'#006651';'#943826'})',1);


for i_subg = 1:2

% preallocate data
pl.dat2plot = nan([sum(pl.elec2plot_i), ...             % 1st dim channels
    size(EP.time,2), ...                                % 2nd dim time
    cellfun(@(x) size(x,1), pl.con_contrast(:,2))', ... % 3rd to n dim = levels defind in pl.con_contrast
    numel(pl.sub2plot{i_subg})]);                               % n+1 dim participants

% extract data across contrasts
% define how to loop through contrasts
t.cont = cellfun(@(x) 1:size(x,1), pl.con_contrast(:,2), 'UniformOutput', false);
t.contidx = combvec(t.cont{:});
% loop across participants
for i_sub = 1:numel(pl.sub2plot{i_subg})
    t.behavior = EP.behavior{pl.sub2plot{i_subg}(i_sub)};
    t.epdata = EP.filtdata{pl.sub2plot{i_subg}(i_sub)}.data;

    % loop through contrasts
    for i_cont = 1:size(t.contidx,2)
        % define evaluation syntax
                
        % preallocate logical idx for trials
        t.idx = false(size(t.contidx,1),numel(t.behavior ));
        for i_contstep = 1:size(t.contidx,1)
            % index step
            t.idx(i_contstep,:) = ...
                any( ... % if it is more than one level, combine
                cell2mat( ...
                cellfun(@(x) strcmp({t.behavior.(pl.con_contrast{i_contstep,1})},x), ...
                pl.con_contrast{i_contstep,2}{t.contidx(i_contstep,i_cont)}, 'UniformOutput', false)' ...
                ) ...
                ,1);
        end
        % combine all contrast conditions
        t.idx_all = all(t.idx,1);
        
        % average across trials
        t.text = sprintf("%1.0f,",t.contidx(:,i_cont));
        t.evaltext = sprintf("pl.dat2plot(:,:,%si_sub) = mean(t.epdata(pl.elec2plot_i,:,t.idx_all),3);",t.text);
        eval(t.evaltext)
    end
end

% reshape pl.data
pl.dat2plot_rs = squeeze(mean(pl.dat2plot,1));
pl.dat2plot_rs(:,end+1,:)=pl.dat2plot_rs(:,1,:)-pl.dat2plot_rs(:,2,:);

pl.data_m = mean(pl.dat2plot_rs,3);
pl.data_sem = std(pl.dat2plot_rs,1,3)./sqrt(numel(pl.sub2plot{i_subg}));



figure;
set(gcf,'Position',[100 100 600 300],'PaperPositionMode','auto')
pl.idx = pl.time2plot_i(1):pl.time2plot_i(2);
h.plsem=[];  h.plm = []; h.pls = []; h.plst = [];

for i_con = 1:size(pl.data_m,2)
    % plot SEM as boundary
    % create data
    pl.xconf = [EP.time(pl.idx) EP.time(pl.idx(end:-1:1))] ;
    pl.yconf = [pl.data_m(pl.idx,i_con)' + pl.data_sem(pl.idx,i_con)' ...
        pl.data_m(pl.idx(end:-1:1),i_con)' - pl.data_sem(pl.idx(end:-1:1),i_con)'];
    % plot
    h.plsem{i_con} = fill(pl.xconf,pl.yconf,pl.concols{i_con}','EdgeColor','none','FaceAlpha',0.3);
    hold on
    
    % plot mean lines
    h.plm{i_con}=...
        plot(EP.time(pl.idx), pl.data_m(pl.idx,i_con),'Color',pl.concols{i_con},'LineWidth',2);
end
xlim(pl.time2plot)
if ~isempty(pl.ylim)
    ylim(pl.ylim)
end
xlabel('time in ms')
ylabel('amplitude in \muV/cm²')
title(sprintf('N=%1.0f | %s',numel(pl.sub2plot{i_subg}), pl.sublabel{i_subg}))
legend([h.plm{:}],[pl.con_contrast{2}{:} 'diff'],'Location','EastOutside','box','off')
grid on
set(gca, 'Ydir','reverse')


h.ax2 = axes('Position', [0.725, .63, .25, .25],'Visible','off');
topoplot(find(pl.elec2plot_i),EP.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',3,1});

set(gcf, 'Color', [1 1 1]);
end

%% plot ERPs exploratively for specified electrodes | event type  by between subject factor stimulus luminance 
% adapted to run with between subject contrasts

% loop through conditions defined in contrasts
pl.con_contrast = {... % contrasts by 1st dim; averaged across second dim
    'eventtype', {{'target'};{'distractor'}}; ...
%     'evnt_type_label', {{'chroma+'}}; ...
    'stim_luminance', {{'offset_to_bckgrd'};{'isolum__to_bckgrd'}}};
pl.sub2plot = 1:numel(F.Subs2use);
% pl.elec2plot = {'P8';'PO8';'P10';'P7';'PO7';'P9'}; % for P1 component lateral !
% pl.elec2plot = {'P7';'PO7';'P9';'O1';'I1';'Oz'; 'Iz';'O2';'I2';'P8';'PO8';'P10';}; % for N2 component posterior!
% pl.elec2plot = {'P3';'P1';'Pz';'P4';'P2';'POz';'PO3';'PO4'}; % for P300 component centro-parietal!
% pl.elec2plot = {'POz';'Oz';'O1';'O2';'Iz'}; % for N1 component centro-parietal
% pl.elec2plot = {'P6';'P8';'PO8';'P10';'P5';'P7';'PO7';'P9'}; % for N2 component lateral

% pl.elec2plot = {'POz'}; % early N2 SN?
% pl.elec2plot = {'CPz';'Cz'}; % early N2 SN?
% pl.elec2plot = {'FCz';'Fz'}; % early frontal?


pl.elec2plot_i = ...
    any(cell2mat(cellfun(@(x) strcmp({EP.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1);

pl.time2plot = [-100 500]; % time in ms
pl.time2plot = [-100 700]; % time in ms
pl.time2plot_i = dsearchn(EP.time', pl.time2plot');

pl.concols = num2cell([1 0.1 0.1; 0.3 0.3 1],2);
pl.concols(:,2) = num2cell([1 0.6 0.6; 0.7 0.7 1],2);
pl.con_label = {'hit';'miss'};


% preallocate data
pl.dat2plot = nan([sum(pl.elec2plot_i), ...             % 1st dim channels
    size(EP.time,2), ...                                % 2nd dim time
    cellfun(@(x) size(x,1), pl.con_contrast(:,2))', ... % 3rd to n dim = levels defind in pl.con_contrast
    numel(pl.sub2plot)]);                               % n+1 dim participants

% extract data across contrasts
% define how to loop through contrasts
t.cont = cellfun(@(x) 1:size(x,1), pl.con_contrast(:,2), 'UniformOutput', false);
t.contidx = combvec(t.cont{:});
% loop across participants
for i_sub = 1:numel(pl.sub2plot)
    t.behavior = EP.behavior{pl.sub2plot(i_sub)};
    t.epdata = EP.filtdata{pl.sub2plot(i_sub)}.data;

    % loop through contrasts
    for i_cont = 1:size(t.contidx,2)
        % define evaluation syntax
                
        % preallocate logical idx for trials
        t.idx = false(size(t.contidx,1),numel(t.behavior ));
        for i_contstep = 1:size(t.contidx,1)
            % index step
            t.idx(i_contstep,:) = ...
                any( ... % if it is more than one level, combine
                cell2mat( ...
                cellfun(@(x) strcmp({t.behavior.(pl.con_contrast{i_contstep,1})},x), ...
                pl.con_contrast{i_contstep,2}{t.contidx(i_contstep,i_cont)}, 'UniformOutput', false)' ...
                ) ...
                ,1);
        end
        % combine all contrast conditions
        t.idx_all = all(t.idx,1);
        
        % average across trials
        t.text = sprintf("%1.0f,",t.contidx(:,i_cont));
        t.evaltext = sprintf("pl.dat2plot(:,:,%si_sub) = mean(t.epdata(pl.elec2plot_i,:,t.idx_all),3,'omitnan');",t.text);
        eval(t.evaltext)
    end
end

% reshape pl.data
pl.dat2plot_rs = squeeze(mean(pl.dat2plot,1));

pl.data_m = mean(pl.dat2plot_rs,4,'omitnan');
% between subject effects: varying number of subjects for contrast: account for that!
t.isnan = isnan(pl.dat2plot_rs);
t.isnan_sum = sum(t.isnan,4);

pl.data_sem = std(pl.dat2plot_rs,1,4,'omitnan')./sqrt(t.isnan_sum);



figure;
set(gcf,'Position',[100 100 600 300],'PaperPositionMode','auto')
pl.idx = pl.time2plot_i(1):pl.time2plot_i(2);
h.plsem=[];  h.plm = []; h.pls = []; h.plst = [];
pl.legend = {};

for i_con1 = 1:size(pl.data_m,2)
    for i_con2 = 1:size(pl.data_m,3)
        % plot SEM as boundary
        % create data
        pl.xconf = [EP.time(pl.idx) EP.time(pl.idx(end:-1:1))] ;
        pl.yconf = [pl.data_m(pl.idx,i_con1,i_con2)' + pl.data_sem(pl.idx,i_con1,i_con2)' ...
            pl.data_m(pl.idx(end:-1:1),i_con1,i_con2)' - pl.data_sem(pl.idx(end:-1:1),i_con1,i_con2)'];
        % plot
%         h.plsem{i_con1,i_con2} = fill(pl.xconf,pl.yconf,pl.concols{i_con1,i_con2},'EdgeColor','none','FaceAlpha',0.3);
        hold on

        % plot mean lines
        h.plm{i_con1,i_con2}=...
            plot(EP.time(pl.idx), pl.data_m(pl.idx,i_con1,i_con2),'Color',pl.concols{i_con1,i_con2},'LineWidth',2);
        
        pl.legend{i_con1,i_con2} = [strcat(pl.con_contrast{1,2}{i_con1},'+', pl.con_contrast{2,2}{i_con2})];
    end
end
xlim(pl.time2plot)
xlabel('time in ms')
ylabel('amplitude in \muV/cm²')
legend([h.plm{:}],[pl.legend{:}],'Location','EastOutside','box','off','Interpreter','none')
grid on
set(gca, 'Ydir','reverse')


h.ax2 = axes('Position', [0.725, .63, .25, .25],'Visible','off');
topoplot(find(pl.elec2plot_i),EP.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',3,1});

set(gcf, 'Color', [1 1 1]);



%% plot ERPs exploratively for specified electrodes | event type  by between subject factor stimulus luminance | running t-tests
% adapted to run with between subject contrasts 
% running ANOVA doesn't work!

% loop through conditions defined in contrasts
pl.con_contrast = {... % contrasts by 1st dim; averaged across second dim
    'eventtype', {{'target'};{'distractor'}}; ...
%     'evnt_type_label', {{'chroma+'}}; ...
    'stim_luminance', {{'offset_to_bckgrd'};{'isolum__to_bckgrd'}}};
pl.sub2plot = 1:numel(F.Subs2use);
% pl.elec2plot = {'P8';'PO8';'P10';'P7';'PO7';'P9'}; % for P1 component lateral !
% pl.elec2plot = {'P7';'PO7';'P9';'O1';'I1';'Oz'; 'Iz';'O2';'I2';'P8';'PO8';'P10';}; % for N2 component posterior!
% pl.elec2plot = {'P3';'P1';'Pz';'P4';'P2';'POz';'PO3';'PO4'}; % for P300 component centro-parietal!
% pl.elec2plot = {'POz';'Oz';'O1';'O2';'Iz'}; % for N1 component centro-parietal
pl.elec2plot = {'P6';'P8';'PO8';'P10';'P5';'P7';'PO7';'P9'}; % for N2 component lateral

% pl.elec2plot = {'POz'}; % early N2 SN?
% pl.elec2plot = {'CPz';'Cz'}; % early N2 SN?
% pl.elec2plot = {'FCz';'Fz'}; % early frontal?


pl.elec2plot_i = ...
    any(cell2mat(cellfun(@(x) strcmp({EP.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1);

pl.time2plot = [-100 500]; % time in ms
pl.time2plot = [-100 700]; % time in ms
pl.time2plot_i = dsearchn(EP.time', pl.time2plot');

pl.concols = num2cell([1 0.1 0.1; 0.3 0.3 1],2);
pl.concols(:,2) = num2cell([1 0.6 0.6; 0.7 0.7 1],2);
pl.con_label = {'hit';'miss'};


% preallocate data
pl.dat2plot = nan([sum(pl.elec2plot_i), ...             % 1st dim channels
    size(EP.time,2), ...                                % 2nd dim time
    cellfun(@(x) size(x,1), pl.con_contrast(:,2))', ... % 3rd to n dim = levels defind in pl.con_contrast
    numel(pl.sub2plot)]);                               % n+1 dim participants

% extract data across contrasts
% define how to loop through contrasts
t.cont = cellfun(@(x) 1:size(x,1), pl.con_contrast(:,2), 'UniformOutput', false);
t.contidx = combvec(t.cont{:});
% loop across participants
for i_sub = 1:numel(pl.sub2plot)
    t.behavior = EP.behavior{pl.sub2plot(i_sub)};
    t.epdata = EP.filtdata{pl.sub2plot(i_sub)}.data;

    % loop through contrasts
    for i_cont = 1:size(t.contidx,2)
        % define evaluation syntax
                
        % preallocate logical idx for trials
        t.idx = false(size(t.contidx,1),numel(t.behavior ));
        for i_contstep = 1:size(t.contidx,1)
            % index step
            t.idx(i_contstep,:) = ...
                any( ... % if it is more than one level, combine
                cell2mat( ...
                cellfun(@(x) strcmp({t.behavior.(pl.con_contrast{i_contstep,1})},x), ...
                pl.con_contrast{i_contstep,2}{t.contidx(i_contstep,i_cont)}, 'UniformOutput', false)' ...
                ) ...
                ,1);
        end
        % combine all contrast conditions
        t.idx_all = all(t.idx,1);
        
        % average across trials
        t.text = sprintf("%1.0f,",t.contidx(:,i_cont));
        t.evaltext = sprintf("pl.dat2plot(:,:,%si_sub) = mean(t.epdata(pl.elec2plot_i,:,t.idx_all),3,'omitnan');",t.text);
        eval(t.evaltext)
    end
end

% reshape pl.data
pl.dat2plot_rs = squeeze(mean(pl.dat2plot,1));

pl.data_m = mean(pl.dat2plot_rs,4,'omitnan');
% between subject effects: varying number of subjects for contrast: account for that!
t.isnan = isnan(pl.dat2plot_rs);
t.isnan_sum = sum(t.isnan,4);

pl.data_sem = std(pl.dat2plot_rs,1,4,'omitnan')./sqrt(t.isnan_sum);

% different contrasts
pl.contrast(1).contrast = sprintf('fixed: %s; %s vs %s',pl.con_contrast{2,2}{1}{:}, pl.con_contrast{1,2}{1}{:}, pl.con_contrast{1,2}{2}{:});
t.subidx = ~squeeze(isnan(pl.dat2plot_rs(1,1,1,:)));
pl.contrast(1).datax = squeeze(pl.dat2plot_rs(:,1,1,t.subidx));
pl.contrast(1).datay = squeeze(pl.dat2plot_rs(:,2,1,t.subidx));
[tt.h tt.p tt.ci tt.stats] = ttest(pl.contrast(1).datax,pl.contrast(1).datay,'Dim',2);
pl.contrast(1).tt_h = tt.h; pl.contrast(1).tt_p = tt.p; pl.contrast(1).tt_t = tt.stats.tstat;

pl.contrast(2).contrast = sprintf('fixed: %s; %s vs %s',pl.con_contrast{2,2}{2}{:}, pl.con_contrast{1,2}{1}{:}, pl.con_contrast{1,2}{2}{:});
t.subidx = ~squeeze(isnan(pl.dat2plot_rs(1,1,2,:)));
pl.contrast(2).datax = squeeze(pl.dat2plot_rs(:,1,2,t.subidx));
pl.contrast(2).datay = squeeze(pl.dat2plot_rs(:,2,2,t.subidx));
[tt.h tt.p tt.ci tt.stats] = ttest(pl.contrast(2).datax,pl.contrast(2).datay,'Dim',2);
pl.contrast(2).tt_h = tt.h; pl.contrast(2).tt_p = tt.p; pl.contrast(2).tt_t = tt.stats.tstat;

pl.contrast(3).contrast = sprintf('fixed: %s; %s vs %s',pl.con_contrast{1,2}{1}{:}, pl.con_contrast{2,2}{1}{:}, pl.con_contrast{2,2}{2}{:});
t.subidx = ~squeeze(isnan(pl.dat2plot_rs(1,1,1,:)));
pl.contrast(3).datax = squeeze(pl.dat2plot_rs(:,1,1,t.subidx));
pl.contrast(3).datay = squeeze(pl.dat2plot_rs(:,1,2,~t.subidx));
[tt.h tt.p tt.ci tt.stats] = ttest2(pl.contrast(3).datax,pl.contrast(3).datay,'Dim',2);
pl.contrast(3).tt_h = tt.h; pl.contrast(3).tt_p = tt.p; pl.contrast(3).tt_t = tt.stats.tstat;

pl.contrast(4).contrast = sprintf('fixed: %s; %s vs %s',pl.con_contrast{1,2}{2}{:}, pl.con_contrast{2,2}{1}{:}, pl.con_contrast{2,2}{2}{:});
t.subidx = ~squeeze(isnan(pl.dat2plot_rs(1,2,1,:)));
pl.contrast(4).datax = squeeze(pl.dat2plot_rs(:,2,1,t.subidx));
pl.contrast(4).datay = squeeze(pl.dat2plot_rs(:,2,2,~t.subidx));
[tt.h tt.p tt.ci tt.stats] = ttest2(pl.contrast(4).datax,pl.contrast(4).datay,'Dim',2);
pl.contrast(4).tt_h = tt.h; pl.contrast(4).tt_p = tt.p; pl.contrast(4).tt_t = tt.stats.tstat;

pl.contrast(5).contrast = sprintf('%s vs %s for %s MINUS %s diff', ...
    pl.con_contrast{2,2}{1}{:}, pl.con_contrast{2,2}{2}{:}, pl.con_contrast{1,2}{1}{:}, pl.con_contrast{1,2}{2}{:});
t.subidx = ~squeeze(isnan(pl.dat2plot_rs(1,1,1,:)));
pl.contrast(5).datax = squeeze(pl.dat2plot_rs(:,1,1,t.subidx))-squeeze(pl.dat2plot_rs(:,2,1,t.subidx));
pl.contrast(5).datay = squeeze(pl.dat2plot_rs(:,1,2,~t.subidx))-squeeze(pl.dat2plot_rs(:,2,2,~t.subidx));
[tt.h tt.p tt.ci tt.stats] = ttest2(pl.contrast(5).datax,pl.contrast(5).datay,'Dim',2);
pl.contrast(5).tt_h = tt.h; pl.contrast(5).tt_p = tt.p; pl.contrast(5).tt_t = tt.stats.tstat;

figure;
set(gcf,'Position',[100 100 600 500],'PaperPositionMode','auto')
pl.idx = pl.time2plot_i(1):pl.time2plot_i(2);
h.plsem=[];  h.plm = []; h.pls = []; h.plst = [];
pl.legend = {};

tiledlayout(4,1)
nexttile([3,1])

for i_con1 = 1:size(pl.data_m,2)
    for i_con2 = 1:size(pl.data_m,3)
        % plot SEM as boundary
        % create data
        pl.xconf = [EP.time(pl.idx) EP.time(pl.idx(end:-1:1))] ;
        pl.yconf = [pl.data_m(pl.idx,i_con1,i_con2)' + pl.data_sem(pl.idx,i_con1,i_con2)' ...
            pl.data_m(pl.idx(end:-1:1),i_con1,i_con2)' - pl.data_sem(pl.idx(end:-1:1),i_con1,i_con2)'];
        % plot
%         h.plsem{i_con1,i_con2} = fill(pl.xconf,pl.yconf,pl.concols{i_con1,i_con2},'EdgeColor','none','FaceAlpha',0.3);
        hold on

        % plot mean lines
        h.plm{i_con1,i_con2}=...
            plot(EP.time(pl.idx), pl.data_m(pl.idx,i_con1,i_con2),'Color',pl.concols{i_con1,i_con2},'LineWidth',2);
        
        pl.legend{i_con1,i_con2} = [strcat(pl.con_contrast{1,2}{i_con1},'+', pl.con_contrast{2,2}{i_con2})];
    end
end
xlim(pl.time2plot)
% xlabel('time in ms')
ylabel('amplitude in \muV/cm²')
legend([h.plm{:}],[pl.legend{:}],'Location','northoutside','box','off','Interpreter','none')
grid on
set(gca, 'Ydir','reverse')




nexttile()
for i_contr = 1:size(pl.contrast,2)
    t.pldata = nan(1,numel(pl.contrast(i_contr).tt_h));
    t.pldata(pl.contrast(i_contr).tt_h==1) = i_contr;
    plot(EP.time(pl.idx), t.pldata(pl.idx),'LineWidth',2)
    hold on
end
set(gca, 'Ydir','reverse')
xlim(pl.time2plot)
ylim([0 size(pl.contrast,2)+1])

legend({pl.contrast.contrast}','Location','southoutside','box','off','Interpreter','none')

h.ax2 = axes('Position', [0.725, .80, .20, .20],'Visible','off');
topoplot(find(pl.elec2plot_i),EP.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',3,1});

set(gcf, 'Color', [1 1 1]);



% with difference curves
figure;
set(gcf,'Position',[100 100 600 700],'PaperPositionMode','auto')
pl.idx = pl.time2plot_i(1):pl.time2plot_i(2);
h.plsem=[];  h.plm = []; h.plm2 = []; h.pls = []; h.plst = [];
pl.legend = {};

tiledlayout(7,1)
nexttile([3,1])

for i_con1 = 1:size(pl.data_m,2)
    for i_con2 = 1:size(pl.data_m,3)
        % plot SEM as boundary
        % create data
        pl.xconf = [EP.time(pl.idx) EP.time(pl.idx(end:-1:1))] ;
        pl.yconf = [pl.data_m(pl.idx,i_con1,i_con2)' + pl.data_sem(pl.idx,i_con1,i_con2)' ...
            pl.data_m(pl.idx(end:-1:1),i_con1,i_con2)' - pl.data_sem(pl.idx(end:-1:1),i_con1,i_con2)'];
        % plot
%         h.plsem{i_con1,i_con2} = fill(pl.xconf,pl.yconf,pl.concols{i_con1,i_con2},'EdgeColor','none','FaceAlpha',0.3);
        hold on

        % plot mean lines
        h.plm{i_con1,i_con2}=...
            plot(EP.time(pl.idx), pl.data_m(pl.idx,i_con1,i_con2),'Color',pl.concols{i_con1,i_con2},'LineWidth',2);
        
        pl.legend{i_con1,i_con2} = [strcat(pl.con_contrast{1,2}{i_con1},'+', pl.con_contrast{2,2}{i_con2})];
    end
end
xlim(pl.time2plot)
% xlabel('time in ms')
ylabel('amplitude in \muV/cm²')
legend([h.plm{:}],[pl.legend{:}],'Location','northoutside','box','off','Interpreter','none')
grid on
set(gca, 'Ydir','reverse')

nexttile([3,1])
for i_contr = 1:size(pl.contrast,2)
    t.mdata = mean(pl.contrast(i_contr).datax(pl.idx,:),2) - mean(pl.contrast(i_contr).datay(pl.idx,:),2);

        % plot mean lines
    h.plm2{i_contr }=...
        plot(EP.time(pl.idx), t.mdata,'LineWidth',2);
    hold on

end
xlim(pl.time2plot)
% xlabel('time in ms')
ylabel('amplitude in \muV/cm²')
% legend([h.plm2{:}],{pl.contrast.contrast}','Location','northoutside','box','off','Interpreter','none')
grid on
set(gca, 'Ydir','reverse')



nexttile()
for i_contr = 1:size(pl.contrast,2)
    t.pldata = nan(1,numel(pl.contrast(i_contr).tt_h));
    t.pldata(pl.contrast(i_contr).tt_h==1) = i_contr;
    plot(EP.time(pl.idx), t.pldata(pl.idx),'LineWidth',2)
    hold on
end
set(gca, 'Ydir','reverse')
xlim(pl.time2plot)
ylim([0 size(pl.contrast,2)+1])

legend({pl.contrast.contrast}','Location','southoutside','box','off','Interpreter','none')

h.ax2 = axes('Position', [0.725, .83, .17, .17],'Visible','off');
topoplot(find(pl.elec2plot_i),EP.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',3,1});

set(gcf, 'Color', [1 1 1]);

%% plot ERPs exploratively for specified electrodes | hit vs miss

% loop through conditions defined in contrasts
pl.con_contrast = {... % contrasts by 1st dim; averaged across second dim
    'eventtype', {{'target'}}; ...
%     'evnt_type_label', {{'chroma+'}}; ...
    'event_response_type', {{'hit'};{'miss'}}};
pl.sub2plot = 1:numel(F.Subs2use);
% pl.elec2plot = {'P8';'PO8';'P10';'P7';'PO7';'P9'}; % for P1 component lateral !
% pl.elec2plot = {'P7';'PO7';'P9';'O1';'I1';'Oz'; 'Iz';'O2';'I2';'P8';'PO8';'P10';}; % for N2 component posterior!
pl.elec2plot = {'P3';'P1';'Pz';'P4';'P2';'POz';'PO3';'PO4'}; % for P300 component centro-parietal!
% pl.elec2plot = {'POz';'Oz';'O1';'O2';'Iz'}; % for N1 component centro-parietal
% pl.elec2plot = {'P6';'P8';'PO8';'P10';'P5';'P7';'PO7';'P9'}; % for N2 component lateral

% pl.elec2plot = {'POz'}; % early N2 SN?
% pl.elec2plot = {'CPz';'Cz'}; % early N2 SN?
% pl.elec2plot = {'FCz';'Fz'}; % early frontal?


pl.elec2plot_i = ...
    any(cell2mat(cellfun(@(x) strcmp({EP.electrodes.labels},x), pl.elec2plot, 'UniformOutput',false)),1);

pl.time2plot = [-100 500]; % time in ms
pl.time2plot = [-100 700]; % time in ms
pl.time2plot_i = dsearchn(EP.time', pl.time2plot');

pl.concols = num2cell([63 63 240; 240 63 63]'./255,1);
pl.con_label = {'hit';'miss'};


% preallocate data
pl.dat2plot = nan([sum(pl.elec2plot_i), ...             % 1st dim channels
    size(EP.time,2), ...                                % 2nd dim time
    cellfun(@(x) size(x,1), pl.con_contrast(:,2))', ... % 3rd to n dim = levels defind in pl.con_contrast
    numel(pl.sub2plot)]);                               % n+1 dim participants

% extract data across contrasts
% define how to loop through contrasts
t.cont = cellfun(@(x) 1:size(x,1), pl.con_contrast(:,2), 'UniformOutput', false);
t.contidx = combvec(t.cont{:});
% loop across participants
for i_sub = 1:numel(pl.sub2plot)
    t.behavior = EP.behavior{pl.sub2plot(i_sub)};
    t.epdata = EP.filtdata{pl.sub2plot(i_sub)}.data;

    % loop through contrasts
    for i_cont = 1:size(t.contidx,2)
        % define evaluation syntax
                
        % preallocate logical idx for trials
        t.idx = false(size(t.contidx,1),numel(t.behavior ));
        for i_contstep = 1:size(t.contidx,1)
            % index step
            t.idx(i_contstep,:) = ...
                any( ... % if it is more than one level, combine
                cell2mat( ...
                cellfun(@(x) strcmp({t.behavior.(pl.con_contrast{i_contstep,1})},x), ...
                pl.con_contrast{i_contstep,2}{t.contidx(i_contstep,i_cont)}, 'UniformOutput', false)' ...
                ) ...
                ,1);
        end
        % combine all contrast conditions
        t.idx_all = all(t.idx,1);
        
        % average across trials
        t.text = sprintf("%1.0f,",t.contidx(:,i_cont));
        t.evaltext = sprintf("pl.dat2plot(:,:,%si_sub) = mean(t.epdata(pl.elec2plot_i,:,t.idx_all),3);",t.text);
        eval(t.evaltext)
    end
end

% reshape pl.data
pl.dat2plot_rs = squeeze(mean(pl.dat2plot,1));

pl.data_m = mean(pl.dat2plot_rs,3);
pl.data_sem = std(pl.dat2plot_rs,1,3)./sqrt(numel(pl.sub2plot));



figure;
set(gcf,'Position',[100 100 600 300],'PaperPositionMode','auto')
pl.idx = pl.time2plot_i(1):pl.time2plot_i(2);
h.plsem=[];  h.plm = []; h.pls = []; h.plst = [];

for i_con = 1:size(pl.data_m,2)
    % plot SEM as boundary
    % create data
    pl.xconf = [EP.time(pl.idx) EP.time(pl.idx(end:-1:1))] ;
    pl.yconf = [pl.data_m(pl.idx,i_con)' + pl.data_sem(pl.idx,i_con)' ...
        pl.data_m(pl.idx(end:-1:1),i_con)' - pl.data_sem(pl.idx(end:-1:1),i_con)'];
    % plot
    h.plsem{i_con} = fill(pl.xconf,pl.yconf,pl.concols{i_con}','EdgeColor','none','FaceAlpha',0.3);
    hold on
    
    % plot mean lines
    h.plm{i_con}=...
        plot(EP.time(pl.idx), pl.data_m(pl.idx,i_con),'Color',pl.concols{i_con},'LineWidth',2);
end
xlim(pl.time2plot)
xlabel('time in ms')
ylabel('amplitude in \muV/cm²')
legend([h.plm{:}],pl.con_label,'Location','EastOutside','box','off')
grid on
set(gca, 'Ydir','reverse')


h.ax2 = axes('Position', [0.725, .63, .25, .25],'Visible','off');
topoplot(find(pl.elec2plot_i),EP.electrodes(1:64),'style','blank','electrodes', 'on','whitebk','on',...
    'emarker2',{find(pl.elec2plot_i),'o','r',3,1});

set(gcf, 'Color', [1 1 1]);


%% plot topagraphies | cue validity


% loop through conditions defined in contrasts
pl.con_contrast = {... % contrasts by 1st dim; averaged across second dim
    'trial_timing_type', {{'regular'}}; ...
    'event_response_type', {{'hit'}}; ...
    'cue_validity_label', {{'valid'};{'neutral'};{'invalid'}}};
pl.sub2plot = 1:numel(F.Subs2use);

pl.time2plot = [160 190]; % time in ms
pl.time2plot = [160 250]; % time in ms
pl.time2plot = [100 150]; % time in ms

% pl.time2plot = [100 190]; % time in ms
pl.time2plot = [250 350]; % time in ms
% pl.time2plot = [400 550]; % time in ms
pl.time2plot_i = dsearchn(EP.time', pl.time2plot');


% preallocate data
pl.dat2plot = nan([size(EP.electrodes,2), ...           % 1st dim channels
    size(EP.time,2), ...                                % 2nd dim time
    cellfun(@(x) size(x,1), pl.con_contrast(:,2))', ... % 3rd to n dim = levels defind in pl.con_contrast
    numel(pl.sub2plot)]);                               % n+1 dim participants

% extract data across contrasts
% define how to loop through contrasts
t.cont = cellfun(@(x) 1:size(x,1), pl.con_contrast(:,2), 'UniformOutput', false);
t.contidx = combvec(t.cont{:});
% loop across participants
for i_sub = 1:numel(pl.sub2plot)
    t.behavior = EP.behavior{pl.sub2plot(i_sub)};
    t.epdata = EP.filtdata{pl.sub2plot(i_sub)}.data;

    % loop through contrasts
    for i_cont = 1:size(t.contidx,2)
        % define evaluation syntax
                
        % preallocate logical idx for trials
        t.idx = false(size(t.contidx,1),numel(t.behavior ));
        for i_contstep = 1:size(t.contidx,1)
            % index step
            t.idx(i_contstep,:) = ...
                any( ... % if it is more than one level, combine
                cell2mat( ...
                cellfun(@(x) strcmp({t.behavior.(pl.con_contrast{i_contstep,1})},x), ...
                pl.con_contrast{i_contstep,2}{t.contidx(i_contstep,i_cont)}, 'UniformOutput', false)' ...
                ) ...
                ,1);
        end
        % combine all contrast conditions
        t.idx_all = all(t.idx,1);
        
        % average across trials
        t.text = sprintf("%1.0f,",t.contidx(:,i_cont));
        t.evaltext = sprintf("pl.dat2plot(:,:,%si_sub) = mean(t.epdata(:,:,t.idx_all),3);",t.text);
        eval(t.evaltext)
    end
end

% reshape pl.data
pl.idx = pl.time2plot_i(1):pl.time2plot_i(2);
pl.dat2plot_rs = squeeze(mean(pl.dat2plot(:,pl.idx,:,:,:,:),2));

pl.data_m = mean(pl.dat2plot_rs,3);

pl.con_label = [pl.con_contrast{end,2}{:}];

pl.clims = [-1 1]*max(abs(pl.data_m),[],"all");

figure;
set(gcf,'Position',[100 100 1100 200],'PaperPositionMode','auto')

for i_con = 1:size(pl.data_m,2)
    subplot(1,size(pl.data_m,2)+1,i_con)
    
    topoplot(pl.data_m(:,i_con), EP.electrodes(1:64), ...
            'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clims,'conv','on','colormap',flipud(cbrewer2('RdBu')),...
            'whitebk','on'); % 'colormap',fake_parula; 'colormap',flipud(cbrewer2('RdBu'))
    
    title(sprintf("%s\n[%1.0f %1.0f]ms",pl.con_label{i_con},pl.time2plot))
    colorbar
end
% mean across conditions
subplot(1,size(pl.data_m,2)+1,i_con+1)

topoplot(mean(pl.data_m,2), EP.electrodes(1:64), ...
    'shading', 'interp', 'numcontour', 0, 'maplimits',pl.clims,'conv','on','colormap',flipud(cbrewer2('RdBu')),...
    'whitebk','on'); % 'colormap',fake_parula; 'colormap',flipud(cbrewer2('RdBu'))

title(sprintf("mean con\n[%1.0f %1.0f]ms",pl.time2plot))
colorbar

%% export single trial erp data for R analysis
pl.sub2plot = 1:numel(F.Subs2use);
pl.erpparams = { ...
    'PD130', [100 190], {'P8';'PO8';'P10';'P7';'PO7';'P9'}; ...
    'N2', [250 330], {'P7';'PO7';'P9';'O1';'I1';'Oz'; 'Iz';'O2';'I2';'P8';'PO8';'P10'}; ...
    'P300', [400 550], {'P3';'P1';'Pz';'P4';'P2';'POz';'PO3';'PO4'} ...
    };

% append behavioral data
out.data = horzcat(EP.behavior{pl.sub2plot});

% add participant number
t.subs = [];
for i_sub = 1:numel(pl.sub2plot)
    t.subs = [t.subs; ...
        repmat( ...
        num2strcell(F.Subs2use(pl.sub2plot(i_sub))), ...
        numel(EP.behavior{pl.sub2plot(i_sub)}), ...
        1)];
end
[out.data.participants] = deal(t.subs{:});


% extract data
for i_erp = 1:size(pl.erpparams,1)
    t.erpdata = [];
    for i_sub = 1:numel(pl.sub2plot)
        % index electrodes
        pl.elec2plot_i = ...
            any(cell2mat(cellfun(@(x) strcmp({EP.electrodes.labels},x), pl.erpparams{i_erp,3}, 'UniformOutput',false)),1);
        pl.time2plot_i = dsearchn(EP.time', pl.erpparams{i_erp,2}');
        % extract erpdata
        t.erpdata = [t.erpdata; ...
            squeeze(mean(EP.filtdata{pl.sub2plot(i_sub)}.data(pl.elec2plot_i,pl.time2plot_i(1):pl.time2plot_i(2),:),[1,2]))];
    end
    % append values
    t.erpdata_cell = num2cell(t.erpdata);
    [out.data.(pl.erpparams{i_erp,1})] = deal(t.erpdata_cell{:});
    % add electrode values
    t.string = sprintf('%s_electrodes',pl.erpparams{i_erp,1});
    [out.data.(t.string)] = deal(vararg2str(pl.erpparams{i_erp,3}));
    % add time values
    t.string = sprintf('%s_timewindow',pl.erpparams{i_erp,1});
    [out.data.(t.string)] = deal(sprintf('[%1.0f %1.0f]ms',pl.erpparams{i_erp,2}));
end

t.path = 'C:\Users\psy05cvd\Dropbox\work\R-statistics\experiments\ssvep_fshiftprobabil\data';
t.datestr = datestr(now,'mm-dd-yyyy_HH-MM');
% t.filename = 'FFT_SSVEP_Amp_data_withoutBehav_sepRDK_peripheral_2elecs';
t.filename = 'ERP_data';
out.data_t = struct2table(out.data);
% remove some data
out.data_t = removevars(out.data_t,{'color_attended','eventRDK_color','event_color','button_presses_fr','button_presses_t'});

% write to textfile
% writetable(out.data_t,fullfile(t.path,sprintf('%s.csv',t.filename)),'Delimiter',';')


%% plot some erp images across time [cue validity]
% loop through conditions defined in contrasts
pl.con_contrast = {... % contrasts by 1st dim; averaged across second dim
    'trial_timing_type', {{'regular'}}; ...
    'event_response_type', {{'hit'}}; ...
%     'evnt_type_label', {{'chroma+'}}; ...
%     'event_response_type', {{'hit','FA','error','miss'}}; ...
    'cue_validity_label', {{'valid'};{'neutral'};{'invalid'}}};
pl.sub2plot = 1:numel(F.Subs2use);

pl.p_time = [25 25 -100 800]; % width step min max
pl.posScale = 1.1;

pl.pl.p_time_i = dsearchn(EP.time', pl.p_time(3:4)');

pl.concols = num2cell([240 63 240; 100 63 100]'./255,1);
pl.con_label = {'valid';'neutral';'invalid'};

pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];


% preallocate data
pl.dat2plot = nan([size(EP.electrodes,2), ...           % 1st dim channels
    size(EP.time,2), ...                                % 2nd dim time
    cellfun(@(x) size(x,1), pl.con_contrast(:,2))', ... % 3rd to n dim = levels defind in pl.con_contrast
    numel(pl.sub2plot)]);                               % n+1 dim participants

% extract data across contrasts
% define how to loop through contrasts
t.cont = cellfun(@(x) 1:size(x,1), pl.con_contrast(:,2), 'UniformOutput', false);
t.contidx = combvec(t.cont{:});
% loop across participants
for i_sub = 1:numel(pl.sub2plot)
    t.behavior = EP.behavior{pl.sub2plot(i_sub)};
    t.epdata = EP.filtdata{pl.sub2plot(i_sub)}.data;

    % loop through contrasts
    for i_cont = 1:size(t.contidx,2)
        % define evaluation syntax
                
        % preallocate logical idx for trials
        t.idx = false(size(t.contidx,1),numel(t.behavior ));
        for i_contstep = 1:size(t.contidx,1)
            % index step
            t.idx(i_contstep,:) = ...
                any( ... % if it is more than one level, combine
                cell2mat( ...
                cellfun(@(x) strcmp({t.behavior.(pl.con_contrast{i_contstep,1})},x), ...
                pl.con_contrast{i_contstep,2}{t.contidx(i_contstep,i_cont)}, 'UniformOutput', false)' ...
                ) ...
                ,1);
        end
        % combine all contrast conditions
        t.idx_all = all(t.idx,1);
        
        % average across trials
        t.text = sprintf("%1.0f,",t.contidx(:,i_cont));
        t.evaltext = sprintf("pl.dat2plot(:,:,%si_sub) = mean(t.epdata(:,:,t.idx_all),3);",t.text);
        eval(t.evaltext)
    end
end

% squeeze data
pl.dat2plot = squeeze(pl.dat2plot);

% calcualte contrast differences
t.contridx = nchoosek(1:size(pl.dat2plot,3),2);
pl.typeidx = [ones(1,size(pl.dat2plot,3)) ones(1,size(t.contridx,1))+1 ones(1,size(t.contridx,1))+2];
pl.dat2plot2 = pl.dat2plot;
for i_contr = 1:size(t.contridx,1)
    pl.dat2plot2(:,:,end+1,:) = pl.dat2plot(:,:,t.contridx(i_contr,1),:) - pl.dat2plot(:,:,t.contridx(i_contr,2),:);
    pl.con_label{end+1} = sprintf('%s-%s',pl.con_label{t.contridx(i_contr,1)},pl.con_label{t.contridx(i_contr,2)});
end
pl.con_label(pl.typeidx==3) = pl.con_label(pl.typeidx==2); % append for ttests

% plot amplitudes of all values
t.time=[];
t.timedot=[];
for i_st = 1:floor((pl.p_time(4)-pl.p_time(1)-pl.p_time(3))/pl.p_time(2))+1
    t.time(i_st,:)=pl.p_time(3)+pl.p_time(2)*(i_st-1)+[0 pl.p_time(1)];
    t.timedot(i_st,1:2)=dsearchn(EP.time',t.time(end,1:2)');
end
% t.t = get(0,'MonitorPositions');
% t.row = round(sqrt((size(t.time,1)+2)/(1/(t.t(1,4)/t.t(1,3)))));
% t.col = ceil(sqrt((size(t.time,1)+2)/(t.t(1,4)/t.t(1,3))));
t.t = [1920 1080];
t.col = ceil(sqrt((size(t.time,1)+2)/(1/(t.t(1)/t.t(2)))));
t.row = ceil(sqrt((size(t.time,1)+2)/(t.t(1)/t.t(2))));
% create plotdata

plotdata=[];
for i_pl = 1:size(t.time,1)
    plotdata(:,:,i_pl)=squeeze(mean(pl.dat2plot2(:,t.timedot(i_pl,1):t.timedot(i_pl,2),:,:),[2,4]));
end

% add ttest data
t.idx = size(pl.dat2plot2,3); % index where data should be appended
t.idx2 = size(pl.dat2plot,3); % index where contrasts start
for i_pl = 1:size(t.time,1)
    for i_tt = 1:size(t.contridx,1)
        t.data = squeeze(mean(pl.dat2plot2(:,t.timedot(i_pl,1):t.timedot(i_pl,2),t.idx2+i_tt,:),[2]));
        [tt.h tt.p tt.ci tt.stats] = ttest(t.data,0,'Dim',2);
        plotdata(:,t.idx+i_tt,i_pl) = tt.p;
    end
end


h.fig = []; h.sp = [];

for i_fig = 1:size(plotdata,2)
    h.fig(i_fig)=figure;
    set(gcf,'Position',[100 100 1200 700],'PaperPositionMode','auto')
    if pl.typeidx(i_fig) == 1
        t.lims = [-1 1]*max(abs(plotdata(:,pl.typeidx==1,:)),[],'all');
    elseif pl.typeidx(i_fig) == 2
        t.lims = [-1 1]*max(abs(plotdata(:,pl.typeidx==2,:)),[],'all');
        % t.lims = [-1 1]*max(abs(plotdata(:,i_fig,:)),[],'all');
    elseif pl.typeidx(i_fig) == 3
        % t.lims = [0 max(abs(log10(plotdata(:,pl.typeidx==3,:))),[],'all')];
        t.lims = [0 max(abs(log10(plotdata(:,i_fig,:))),[],'all')];
        
        % define colormaps
        t.pcriterion = abs(log10(0.05));
        if t.lims(2)<t.pcriterion
            t.colormap = repmat(linspace(1,0.3,1000)',1,3);
        else
            t.border = ceil((t.pcriterion / t.lims(2))*1000);
            % temp.colormap = [repmat([0.5 0.5 0.5],temp.border,1); [linspace(0,1,100-temp.border)' repmat(0,100-temp.border,1) linspace(1,0,100-temp.border)']];
            t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
        end
            
    end
    for i_spl = 1:size(plotdata,3)
        h.sp(i_spl)=subplot(t.row,t.col,i_spl);
        if pl.typeidx(i_fig) == 1 || pl.typeidx(i_fig) == 2
%             topoplot( plotdata(:,i_fig,i_spl), EP.electrodes(1:64), ...
%                 'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',t.lims, ...
%                 'electrodes','off','colormap',fake_parula,'whitebk','on');
            topoplot( plotdata(:,i_fig,i_spl), EP.electrodes(1:64), ...
                'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',t.lims, ...
                'electrodes','off','colormap',flipud(cbrewer2('RdBu')),'whitebk','on');
        elseif pl.typeidx(i_fig) == 3
            topoplot( abs(log10(plotdata(:,i_fig,i_spl))), EP.electrodes(1:64), ...
                'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',t.lims, ...
                'electrodes','off','colormap',t.colormap,'whitebk','on');
        end
        title(sprintf('[%1.0f %1.0f]',t.time(i_spl,1),t.time(i_spl,2)),'FontSize',6)
        t.pos = get(h.sp(i_spl),'Position');
        set(h.sp(i_spl),'Position',[t.pos(1:2)-(t.pos(3:4).*((pl.posScale-1)/2)) t.pos(3:4).*pl.posScale])
    end
    h.sp(i_spl+1)=subplot(t.row,t.col,i_spl+2);
    topoplot( [], EP.electrodes(1:64),  ...
        'style','blank','whitebk','on');
    title(sprintf('%s',pl.con_label{i_fig}),'FontSize',8)
    t.pos = get(h.sp(i_spl+1),'Position');
    set(h.sp(i_spl+1),'Position',[t.pos(1:2)-(t.pos(3:4).*((pl.posScale-1)/2)) t.pos(3:4).*pl.posScale])
    
    t.pos2 = get(h.sp(i_spl+1),'Position');
    t.pos3 = get(h.sp(i_spl+1),'OuterPosition');
    h.a1 = axes('position',[t.pos3(1) t.pos2(2) t.pos3(3) t.pos2(4)],'Visible','off');
    if pl.typeidx(i_fig) == 3
        colormap(gca, t.colormap)
        clim(t.lims);
        h.c3 = colorbar();
        t.pos4 = get(h.c3,'Position');
        set(h.c3,'Position',[t.pos4(1)+0.065 t.pos4(2)+(t.pos4(4)*(1/6)) t.pos4(3) t.pos4(4)*(2/3)])
        t.yticks = get(h.c3,'YTick');
        set(h.c3,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.lims(end),1,'last')), ...
            'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.lims(end),1,'last')))
    else
        colormap(gca,flipud(cbrewer2('RdBu')))
        clim(t.lims);
        h.c3 = colorbar();
        t.pos4 = get(h.c3,'Position');
        set(h.c3,'Position',[t.pos4(1)+0.065 t.pos4(2)+(t.pos4(4)*(1/6)) t.pos4(3) t.pos4(4)*(2/3)])
    end
    axcopy(h.fig(i_fig))
end



%% plot some erp images across time [hit or no hit]
% loop through conditions defined in contrasts
pl.con_contrast = {... % contrasts by 1st dim; averaged across second dim
    'trial_timing_type', {{'regular'}}; ...
    'evnt_type_label', {{'chroma+'}}; ...
    'event_response_type', {{'hit'};{'FA','error','miss'}}};
pl.sub2plot = 1:numel(F.Subs2use);

pl.p_time = [25 25 -100 800]; % width step min max
pl.posScale = 1.1;

pl.pl.p_time_i = dsearchn(EP.time', pl.p_time(3:4)');

pl.concols = num2cell([240 63 240; 100 63 100]'./255,1);
pl.con_label = {'hit';'error+FA+miss'};

pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];


% preallocate data
pl.dat2plot = nan([size(EP.electrodes,2), ...           % 1st dim channels
    size(EP.time,2), ...                                % 2nd dim time
    cellfun(@(x) size(x,1), pl.con_contrast(:,2))', ... % 3rd to n dim = levels defind in pl.con_contrast
    numel(pl.sub2plot)]);                               % n+1 dim participants

% extract data across contrasts
% define how to loop through contrasts
t.cont = cellfun(@(x) 1:size(x,1), pl.con_contrast(:,2), 'UniformOutput', false);
t.contidx = combvec(t.cont{:});
% loop across participants
for i_sub = 1:numel(pl.sub2plot)
    t.behavior = EP.behavior{pl.sub2plot(i_sub)};
    t.epdata = EP.filtdata{pl.sub2plot(i_sub)}.data;

    % loop through contrasts
    for i_cont = 1:size(t.contidx,2)
        % define evaluation syntax
                
        % preallocate logical idx for trials
        t.idx = false(size(t.contidx,1),numel(t.behavior ));
        for i_contstep = 1:size(t.contidx,1)
            % index step
            t.idx(i_contstep,:) = ...
                any( ... % if it is more than one level, combine
                cell2mat( ...
                cellfun(@(x) strcmp({t.behavior.(pl.con_contrast{i_contstep,1})},x), ...
                pl.con_contrast{i_contstep,2}{t.contidx(i_contstep,i_cont)}, 'UniformOutput', false)' ...
                ) ...
                ,1);
        end
        % combine all contrast conditions
        t.idx_all = all(t.idx,1);
        
        % average across trials
        t.text = sprintf("%1.0f,",t.contidx(:,i_cont));
        t.evaltext = sprintf("pl.dat2plot(:,:,%si_sub) = mean(t.epdata(:,:,t.idx_all),3);",t.text);
        eval(t.evaltext)
    end
end

% squeeze data
pl.dat2plot = squeeze(pl.dat2plot);

% calcualte contrast differences
t.contridx = nchoosek(1:size(pl.dat2plot,3),2);
pl.typeidx = [ones(1,size(pl.dat2plot,3)) ones(1,size(t.contridx,1))+1 ones(1,size(t.contridx,1))+2];
pl.dat2plot2 = pl.dat2plot;
for i_contr = 1:size(t.contridx,1)
    pl.dat2plot2(:,:,end+1,:) = pl.dat2plot(:,:,t.contridx(i_contr,1),:) - pl.dat2plot(:,:,t.contridx(i_contr,2),:);
    pl.con_label{end+1} = sprintf('%s-%s',pl.con_label{t.contridx(i_contr,1)},pl.con_label{t.contridx(i_contr,2)});
end
pl.con_label(pl.typeidx==3) = pl.con_label(pl.typeidx==2); % append for ttests

% plot amplitudes of all values
t.time=[];
t.timedot=[];
for i_st = 1:floor((pl.p_time(4)-pl.p_time(1)-pl.p_time(3))/pl.p_time(2))+1
    t.time(i_st,:)=pl.p_time(3)+pl.p_time(2)*(i_st-1)+[0 pl.p_time(1)];
    t.timedot(i_st,1:2)=dsearchn(EP.time',t.time(end,1:2)');
end
% t.t = get(0,'MonitorPositions');
% t.row = round(sqrt((size(t.time,1)+2)/(1/(t.t(1,4)/t.t(1,3)))));
% t.col = ceil(sqrt((size(t.time,1)+2)/(t.t(1,4)/t.t(1,3))));
t.t = [1920 1080];
t.col = ceil(sqrt((size(t.time,1)+2)/(1/(t.t(1)/t.t(2)))));
t.row = ceil(sqrt((size(t.time,1)+2)/(t.t(1)/t.t(2))));
% create plotdata

plotdata=[];
for i_pl = 1:size(t.time,1)
    plotdata(:,:,i_pl)=squeeze(mean(pl.dat2plot2(:,t.timedot(i_pl,1):t.timedot(i_pl,2),:,:),[2,4]));
end

% add ttest data
t.idx = size(pl.dat2plot2,3); % index where data should be appended
t.idx2 = size(pl.dat2plot,3); % index where contrasts start
for i_pl = 1:size(t.time,1)
    for i_tt = 1:size(t.contridx,1)
        t.data = squeeze(mean(pl.dat2plot2(:,t.timedot(i_pl,1):t.timedot(i_pl,2),t.idx2+i_tt,:),[2]));
        [tt.h tt.p tt.ci tt.stats] = ttest(t.data,0,'Dim',2);
        plotdata(:,t.idx+i_tt,i_pl) = tt.p;
    end
end

h.fig = []; h.sp = [];

for i_fig = 1:size(plotdata,2)
    h.fig(i_fig)=figure;
    set(gcf,'Position',[100 100 1200 700],'PaperPositionMode','auto')
    if pl.typeidx(i_fig) == 1
        t.lims = [-1 1]*max(abs(plotdata(:,pl.typeidx==1,:)),[],'all');
    elseif pl.typeidx(i_fig) == 2
        t.lims = [-1 1]*max(abs(plotdata(:,pl.typeidx==2,:)),[],'all');
        % t.lims = [-1 1]*max(abs(plotdata(:,i_fig,:)),[],'all');
    elseif pl.typeidx(i_fig) == 3
        % t.lims = [0 max(abs(log10(plotdata(:,pl.typeidx==3,:))),[],'all')];
        t.lims = [0 max(abs(log10(plotdata(:,i_fig,:))),[],'all')];

        % define colormaps
        t.pcriterion = abs(log10(0.05));
        if t.lims(2)<t.pcriterion
            t.colormap = repmat(linspace(1,0.3,1000)',1,3);
        else
            t.border = ceil((t.pcriterion / t.lims(2))*1000);
            % temp.colormap = [repmat([0.5 0.5 0.5],temp.border,1); [linspace(0,1,100-temp.border)' repmat(0,100-temp.border,1) linspace(1,0,100-temp.border)']];
            t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
        end
    end
    for i_spl = 1:size(plotdata,3)
        h.sp(i_spl)=subplot(t.row,t.col,i_spl);
        if pl.typeidx(i_fig) == 1 || pl.typeidx(i_fig) == 2
            topoplot(plotdata(:,i_fig,i_spl), EP.electrodes(1:64), ...
                'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',t.lims, ...
                'electrodes','off','colormap',flipud(cbrewer2('RdBu')),'whitebk','on');
        elseif pl.typeidx(i_fig) == 3
            topoplot( abs(log10(plotdata(:,i_fig,i_spl))), EP.electrodes(1:64), ...
                'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',t.lims, ...
                'electrodes','off','colormap',t.colormap,'whitebk','on');
        end
       
        title(sprintf('[%1.0f %1.0f]',t.time(i_spl,1),t.time(i_spl,2)),'FontSize',6)
        t.pos = get(h.sp(i_spl),'Position');
        set(h.sp(i_spl),'Position',[t.pos(1:2)-(t.pos(3:4).*((pl.posScale-1)/2)) t.pos(3:4).*pl.posScale])
    end
    h.sp(i_spl+1)=subplot(t.row,t.col,i_spl+2);
    topoplot( [], EP.electrodes(1:64),  ...
        'style','blank','whitebk','on');
    title(sprintf('%s',pl.con_label{i_fig}),'FontSize',8)
    t.pos = get(h.sp(i_spl+1),'Position');
    set(h.sp(i_spl+1),'Position',[t.pos(1:2)-(t.pos(3:4).*((pl.posScale-1)/2)) t.pos(3:4).*pl.posScale])
    
    t.pos2 = get(h.sp(i_spl+1),'Position');
    t.pos3 = get(h.sp(i_spl+1),'OuterPosition');
    h.a1 = axes('position',[t.pos3(1) t.pos2(2) t.pos3(3) t.pos2(4)],'Visible','off');
    if pl.typeidx(i_fig) == 3
        colormap(gca, t.colormap)
        clim(t.lims);
        h.c3 = colorbar();
        t.pos4 = get(h.c3,'Position');
        set(h.c3,'Position',[t.pos4(1)+0.065 t.pos4(2)+(t.pos4(4)*(1/6)) t.pos4(3) t.pos4(4)*(2/3)])
        t.yticks = get(h.c3,'YTick');
        set(h.c3,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.lims(end),1,'last')), ...
            'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.lims(end),1,'last')))
    else
        colormap(gca,flipud(cbrewer2('RdBu')))
        clim(t.lims);
        h.c3 = colorbar();
        t.pos4 = get(h.c3,'Position');
        set(h.c3,'Position',[t.pos4(1)+0.065 t.pos4(2)+(t.pos4(4)*(1/6)) t.pos4(3) t.pos4(4)*(2/3)])
    end
    axcopy(h.fig(i_fig))
end


%% plot some erp images across time [chroma]
% loop through conditions defined in contrasts
pl.con_contrast = {... % contrasts by 1st dim; averaged across second dim
    'trial_timing_type', {{'regular'}}; ...
    'event_response_type', {{'hit'}}; ...
%     'event_response_type', {{'hit','FA','error','miss'}}; ...
    'evnt_type_label', {{'chroma+'};{'chroma-'}}};
pl.sub2plot = 1:numel(F.Subs2use);

pl.p_time = [25 25 -100 800]; % width step min max
pl.posScale = 1.1;

pl.pl.p_time_i = dsearchn(EP.time', pl.p_time(3:4)');

pl.concols = num2cell([240 63 240; 100 63 100]'./255,1);
pl.con_label = {'chroma+';'chroma-'};

pl.plegend = {'1';'0.5';'0.25';'0.1';'0.05';'0.01';'0.001';'0.0001';'0.00001'};
pl.pcorrect = [0 abs(log10([0.5 0.25 0.1 0.05 0.01 0.001 0.0001 0.00001]))];


% preallocate data
pl.dat2plot = nan([size(EP.electrodes,2), ...           % 1st dim channels
    size(EP.time,2), ...                                % 2nd dim time
    cellfun(@(x) size(x,1), pl.con_contrast(:,2))', ... % 3rd to n dim = levels defind in pl.con_contrast
    numel(pl.sub2plot)]);                               % n+1 dim participants

% extract data across contrasts
% define how to loop through contrasts
t.cont = cellfun(@(x) 1:size(x,1), pl.con_contrast(:,2), 'UniformOutput', false);
t.contidx = combvec(t.cont{:});
% loop across participants
for i_sub = 1:numel(pl.sub2plot)
    t.behavior = EP.behavior{pl.sub2plot(i_sub)};
    t.epdata = EP.filtdata{pl.sub2plot(i_sub)}.data;

    % loop through contrasts
    for i_cont = 1:size(t.contidx,2)
        % define evaluation syntax
                
        % preallocate logical idx for trials
        t.idx = false(size(t.contidx,1),numel(t.behavior ));
        for i_contstep = 1:size(t.contidx,1)
            % index step
            t.idx(i_contstep,:) = ...
                any( ... % if it is more than one level, combine
                cell2mat( ...
                cellfun(@(x) strcmp({t.behavior.(pl.con_contrast{i_contstep,1})},x), ...
                pl.con_contrast{i_contstep,2}{t.contidx(i_contstep,i_cont)}, 'UniformOutput', false)' ...
                ) ...
                ,1);
        end
        % combine all contrast conditions
        t.idx_all = all(t.idx,1);
        
        % average across trials
        t.text = sprintf("%1.0f,",t.contidx(:,i_cont));
        t.evaltext = sprintf("pl.dat2plot(:,:,%si_sub) = mean(t.epdata(:,:,t.idx_all),3);",t.text);
        eval(t.evaltext)
    end
end

% squeeze data
pl.dat2plot = squeeze(pl.dat2plot);

% calcualte contrast differences
t.contridx = nchoosek(1:size(pl.dat2plot,3),2);
pl.typeidx = [ones(1,size(pl.dat2plot,3)) ones(1,size(t.contridx,1))+1 ones(1,size(t.contridx,1))+2];
pl.dat2plot2 = pl.dat2plot;
for i_contr = 1:size(t.contridx,1)
    pl.dat2plot2(:,:,end+1,:) = pl.dat2plot(:,:,t.contridx(i_contr,1),:) - pl.dat2plot(:,:,t.contridx(i_contr,2),:);
    pl.con_label{end+1} = sprintf('%s-%s',pl.con_label{t.contridx(i_contr,1)},pl.con_label{t.contridx(i_contr,2)});
end
pl.con_label(pl.typeidx==3) = pl.con_label(pl.typeidx==2); % append for ttests

% plot amplitudes of all values
t.time=[];
t.timedot=[];
for i_st = 1:floor((pl.p_time(4)-pl.p_time(1)-pl.p_time(3))/pl.p_time(2))+1
    t.time(i_st,:)=pl.p_time(3)+pl.p_time(2)*(i_st-1)+[0 pl.p_time(1)];
    t.timedot(i_st,1:2)=dsearchn(EP.time',t.time(end,1:2)');
end
% t.t = get(0,'MonitorPositions');
% t.row = round(sqrt((size(t.time,1)+2)/(1/(t.t(1,4)/t.t(1,3)))));
% t.col = ceil(sqrt((size(t.time,1)+2)/(t.t(1,4)/t.t(1,3))));
t.t = [1920 1080];
t.col = ceil(sqrt((size(t.time,1)+2)/(1/(t.t(1)/t.t(2)))));
t.row = ceil(sqrt((size(t.time,1)+2)/(t.t(1)/t.t(2))));
% create plotdata

plotdata=[];
for i_pl = 1:size(t.time,1)
    plotdata(:,:,i_pl)=squeeze(mean(pl.dat2plot2(:,t.timedot(i_pl,1):t.timedot(i_pl,2),:,:),[2,4]));
end

% add ttest data
t.idx = size(pl.dat2plot2,3); % index where data should be appended
t.idx2 = size(pl.dat2plot,3); % index where contrasts start
for i_pl = 1:size(t.time,1)
    for i_tt = 1:size(t.contridx,1)
        t.data = squeeze(mean(pl.dat2plot2(:,t.timedot(i_pl,1):t.timedot(i_pl,2),t.idx2+i_tt,:),[2]));
        [tt.h tt.p tt.ci tt.stats] = ttest(t.data,0,'Dim',2);
        plotdata(:,t.idx+i_tt,i_pl) = tt.p;
    end
end


h.fig = []; h.sp = [];

for i_fig = 1:size(plotdata,2)
    h.fig(i_fig)=figure;
    set(gcf,'Position',[100 100 1200 700],'PaperPositionMode','auto')
    if pl.typeidx(i_fig) == 1
        t.lims = [-1 1]*max(abs(plotdata(:,pl.typeidx==1,:)),[],'all');
    elseif pl.typeidx(i_fig) == 2
        t.lims = [-1 1]*max(abs(plotdata(:,pl.typeidx==2,:)),[],'all');
        % t.lims = [-1 1]*max(abs(plotdata(:,i_fig,:)),[],'all');
    elseif pl.typeidx(i_fig) == 3
        % t.lims = [0 max(abs(log10(plotdata(:,pl.typeidx==3,:))),[],'all')];
        t.lims = [0 max(abs(log10(plotdata(:,i_fig,:))),[],'all')];

        % define colormaps
        t.pcriterion = abs(log10(0.05));
        if t.lims(2)<t.pcriterion
            t.colormap = repmat(linspace(1,0.3,1000)',1,3);
        else
            t.border = ceil((t.pcriterion / t.lims(2))*1000);
            % temp.colormap = [repmat([0.5 0.5 0.5],temp.border,1); [linspace(0,1,100-temp.border)' repmat(0,100-temp.border,1) linspace(1,0,100-temp.border)']];
            t.colormap = [repmat(linspace(1,0.3,t.border)',1,3); [linspace(0,1,1000-t.border)' zeros(1000-t.border,1) linspace(1,0,1000-t.border)']];
        end
    end
    for i_spl = 1:size(plotdata,3)
        h.sp(i_spl)=subplot(t.row,t.col,i_spl);
        if pl.typeidx(i_fig) == 1 || pl.typeidx(i_fig) == 2
            topoplot(plotdata(:,i_fig,i_spl), EP.electrodes(1:64), ...
                'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',t.lims, ...
                'electrodes','off','colormap',flipud(cbrewer2('RdBu')),'whitebk','on');
        elseif pl.typeidx(i_fig) == 3
            topoplot( abs(log10(plotdata(:,i_fig,i_spl))), EP.electrodes(1:64), ...
                'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',t.lims, ...
                'electrodes','off','colormap',t.colormap,'whitebk','on');
        end
       
        title(sprintf('[%1.0f %1.0f]',t.time(i_spl,1),t.time(i_spl,2)),'FontSize',6)
        t.pos = get(h.sp(i_spl),'Position');
        set(h.sp(i_spl),'Position',[t.pos(1:2)-(t.pos(3:4).*((pl.posScale-1)/2)) t.pos(3:4).*pl.posScale])
    end
    h.sp(i_spl+1)=subplot(t.row,t.col,i_spl+2);
    topoplot( [], EP.electrodes(1:64),  ...
        'style','blank','whitebk','on');
    title(sprintf('%s',pl.con_label{i_fig}),'FontSize',8)
    t.pos = get(h.sp(i_spl+1),'Position');
    set(h.sp(i_spl+1),'Position',[t.pos(1:2)-(t.pos(3:4).*((pl.posScale-1)/2)) t.pos(3:4).*pl.posScale])
    
    t.pos2 = get(h.sp(i_spl+1),'Position');
    t.pos3 = get(h.sp(i_spl+1),'OuterPosition');
    h.a1 = axes('position',[t.pos3(1) t.pos2(2) t.pos3(3) t.pos2(4)],'Visible','off');
    if pl.typeidx(i_fig) == 3
        colormap(gca, t.colormap)
        clim(t.lims);
        h.c3 = colorbar();
        t.pos4 = get(h.c3,'Position');
        set(h.c3,'Position',[t.pos4(1)+0.065 t.pos4(2)+(t.pos4(4)*(1/6)) t.pos4(3) t.pos4(4)*(2/3)])
        t.yticks = get(h.c3,'YTick');
        set(h.c3,'YTick',pl.pcorrect(1:find(pl.pcorrect<t.lims(end),1,'last')), ...
            'YTickLabel',pl.plegend(1:find(pl.pcorrect<t.lims(end),1,'last')))
    else
        colormap(gca,flipud(cbrewer2('RdBu')))
        clim(t.lims);
        h.c3 = colorbar();
        t.pos4 = get(h.c3,'Position');
        set(h.c3,'Position',[t.pos4(1)+0.065 t.pos4(2)+(t.pos4(4)*(1/6)) t.pos4(3) t.pos4(4)*(2/3)])
    end
    axcopy(h.fig(i_fig))
end

%% run tfce analysis | contrast valid - invalid

% loop through conditions defined in contrasts
pl.con_contrast = {... % contrasts by 1st dim; averaged across second dim
    'trial_timing_type', {{'regular'}}; ...
    'event_response_type', {{'hit'}}; ...
%     'evnt_type_label', {{'chroma+'}}; ...
%     'event_response_type', {{'hit','FA','error','miss'}}; ...
%     'cue_validity_label', {{'valid'};{'neutral'};{'invalid'}}};
    'cue_validity_label', {{'valid'};{'invalid'}}};
pl.sub2plot = 1:numel(F.Subs2use);

% pl.con_label = {'valid';'neutral';'invalid'};
pl.con_label = {'valid';'invalid'};

pl.time2plot = [0 700]; % time in ms
pl.time2plot_i = dsearchn(EP.time', pl.time2plot');

% pl.concols = num2cell([25 138 130; 41 60 74; 241 131 26]'./255,1);
pl.concols = num2cell([25 138 130; 241 131 26]'./255,1);

% tfce parameters
p.e_h               = [0.66 2]; % tfce parameter
p.Samp              = 128; % data sampling rate (TFA.srate)
p.nperm             = 10000; % number of permutations, if feasible use 100000
p.plevel            = 0.05;p.tfce_olddate = '22-Mar-2024';
startdate           = datetime("today");
savnam_tfce_erp     = sprintf('tfce_erp_%1.2fe_%1.2fh_%s.mat',p.e_h(1),p.e_h(2),startdate);
eloc                = EP.electrodes;


% preallocate data
pl.dat2plot = nan([numel(EP.electrodes), ...             % 1st dim channels
    size(EP.time,2), ...                                % 2nd dim time
    cellfun(@(x) size(x,1), pl.con_contrast(:,2))', ... % 3rd to n dim = levels defind in pl.con_contrast
    numel(pl.sub2plot)]);                               % n+1 dim participants

% extract data across contrasts
% define how to loop through contrasts
t.cont = cellfun(@(x) 1:size(x,1), pl.con_contrast(:,2), 'UniformOutput', false);
t.contidx = combvec(t.cont{:});
% loop across participants
for i_sub = 1:numel(pl.sub2plot)
    t.behavior = EP.behavior{pl.sub2plot(i_sub)};
    t.epdata = EP.filtdata{pl.sub2plot(i_sub)}.data;

    % loop through contrasts
    for i_cont = 1:size(t.contidx,2)
        % define evaluation syntax
                
        % preallocate logical idx for trials
        t.idx = false(size(t.contidx,1),numel(t.behavior ));
        for i_contstep = 1:size(t.contidx,1)
            % index step
            t.idx(i_contstep,:) = ...
                any( ... % if it is more than one level, combine
                cell2mat( ...
                cellfun(@(x) strcmp({t.behavior.(pl.con_contrast{i_contstep,1})},x), ...
                pl.con_contrast{i_contstep,2}{t.contidx(i_contstep,i_cont)}, 'UniformOutput', false)' ...
                ) ...
                ,1);
        end
        % combine all contrast conditions
        t.idx_all = all(t.idx,1);
        
        % average across trials
        t.text = sprintf("%1.0f,",t.contidx(:,i_cont));
        t.evaltext = sprintf("pl.dat2plot(:,:,%si_sub) = mean(t.epdata(:,:,t.idx_all),3);",t.text);
        eval(t.evaltext)
    end
end

% reshape pl.data
pl.dat2plot_rs = squeeze(pl.dat2plot(:,pl.time2plot_i(1):pl.time2plot_i(2),:,:,:,:));

if exist([savnam_tfce_erp(1:regexp(savnam_tfce_erp,[date '.mat'])-1) p.tfce_olddate '.mat'],'file') 
    try % try to load previous file
        load([savnam_tfce_erp(1:regexp(savnam_tfce_erp,[date '.mat'])-1) p.tfce_olddate '.mat'])
    catch
        error('error when loading file')
    end
    results_tfce_erp = Results;
else % do tfce
    % tfce on two dimensions (channels x time), subjects have to be in dim 1
    % channel neighborhood calculated within ept_TFCE function

    results_tfce_erp= ept_TFCE(permute(pl.dat2plot_rs(:,:,1,:),[4 1 2 3]),permute(pl.dat2plot_rs(:,:,2,:),[4 1 2 3]),...
        eloc,'rsample',p.Samp,...
        'type','d','plots',0,...
        'e_h',p.e_h,'nperm',p.nperm,...
        'flag_tfce',true,'flag_ft',false,... %channel testing
        'savename',savnam_tfce_erp, ...
        'flag_save',true);
end

ChN = ept_ChN2(eloc);
res.cluster_erp = ept_calculateClusters(results_tfce_erp,ChN,0.05);

% illustrate results

clear h
for i_cluster = 1:numel(res.cluster_erp)
    h.fig(i_cluster)=figure;
    set(gcf,'Position',[100 100 1100 300],'PaperPositionMode','auto')
    tiledlayout(1,6,'TileSpacing','compact')

    % show time traces
    nexttile([1 3])
    % extract data
    t.sigchans = any(res.cluster_erp(i_cluster).cluster_locations,2);
    t.sigtimes = any(res.cluster_erp(i_cluster).cluster_locations,1);
    pl.dat_raw = pl.dat2plot(t.sigchans,:,:,:,:,:);
    pl.dat_raw_m = squeeze(mean(pl.dat_raw,[1,6]));
    % plot traces
    for i_trace = 1:size(pl.dat_raw_m,2)
        h.pl_trace(i_trace) = plot(EP.time, pl.dat_raw_m(:,i_trace),'Color',pl.concols{i_trace})
        hold on
    end
    % plot difference
    pl.dat_diff = -diff(pl.dat_raw_m,1,2);
    h.pl_trace(i_trace+1) = plot(EP.time, pl.dat_diff,'Color',[0.4 0.4 0.4]);

%     % add significant cluster difference [looks weird]
%     t.idx = repmat(res.cluster_erp(i_cluster).cluster_locations,[1,1,size(pl.dat2plot_rs,[3,4])]);
%     t.data = pl.dat2plot_rs;
%     t.data(~t.idx) = nan;
%     pl.data_cluster = mean(t.data,[1,3,4],"omitnan");
%     h.pl_trace(i_trace+2) = plot(EP.time(pl.time2plot_i(1):pl.time2plot_i(2)), pl.data_cluster,'Color',[0 0 0],'LineWidth',2)

    pl.dat_diff_sigcluster = pl.dat_diff(pl.time2plot_i(1):pl.time2plot_i(2));
    pl.dat_diff_sigcluster(~t.sigtimes) = nan;
    h.pl_trace(i_trace+2) = plot(EP.time(pl.time2plot_i(1):pl.time2plot_i(2)), pl.dat_diff_sigcluster,'Color',[0 0 0],'LineWidth',2)

    legend(h.pl_trace(1:end-1),[pl.con_label; sprintf('%s - %s',pl.con_label{:})],'Location','northwest')
    xlim([-100 pl.time2plot(2)])
    xlabel('time in ms')
    ylabel('amplitude in \muV/m²')
    grid("on")
    title(sprintf('traces for cluster %1.0f',i_cluster))
    

    % show three difference topos
    t.tidx = [find(t.sigtimes,1,'first') find(t.sigtimes,1,'last')]+pl.time2plot_i(1)-1;
    t.tidx2 = round(linspace(t.tidx(1),t.tidx(2),4)); t.tidx3 = [t.tidx2(1:3);t.tidx2(2:4)];
    pl.dat_raw = squeeze(mean(pl.dat2plot,6));
    pl.dat2plot_topo = [];
    for i_topo = 1:3
        pl.dat2plot_topo(:,i_topo) = -diff(squeeze(mean(pl.dat_raw(:,t.tidx3(1,i_topo):t.tidx3(2,i_topo),:),2)),1,2);
    end
    pl.clim = [-1 1]*max(abs(pl.dat2plot_topo),[],"all");
    for i_topo = 1:3
        nexttile
        topoplot(pl.dat2plot_topo(:,i_topo), eloc, ...
                'shading', 'flat', 'numcontour', 0, 'conv','on','maplimits',pl.clim, ...
                'colormap',flipud(cbrewer2('RdBu')),'whitebk','on', ...
                'emarker2',{find(t.sigchans),'+','g',4,1});
        title(sprintf('diff [%1.0f %1.0f]ms',EP.time(t.tidx3(:,i_topo))))
    end
    colorbar
end









