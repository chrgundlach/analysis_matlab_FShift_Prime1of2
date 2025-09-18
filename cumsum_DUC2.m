%% fits for the grand mean
% compmodelfit
cd '/Users/andreaskeil/Dropbox (UFL)/ak_own/Ducmodel/DUC2'
clear 
clc

thresh = 1.5; 

taxis4plot = -23:1000/512:1900-4;

%close all
load('Exp2_ModelingDataCorrected.mat')

% first: Pleasant Pictures
RDKtempP = Exp2_data{1,1}{1,2}{2,3}.*10^13;
RSVPtempP = Exp2_data{2,1}{1,2}{2,3}.*10^13;

RDKtempP = bslcorr_percent(RDKtempP, 400:500);
RSVPtempP = bslcorr_percent(RSVPtempP, 400:500);

RDKtempP = RDKtempP(501:end);
RSVPtempP = RSVPtempP(501:end); 

figure (101), 
plot(taxis4plot, RDKtempP), hold on, plot(taxis4plot, RSVPtempP), hold off, title('pleasant')

cusumRDKP_mean = cumsum(RDKtempP)./std(abs(cumsum(RDKtempP)));
cusumRSVPP_mean = cumsum(RSVPtempP)./std(abs(cumsum(RSVPtempP)));

figure(102), plot(taxis4plot, cusumRDKP_mean), hold on, plot(taxis4plot, cusumRSVPP_mean), axis([-.02 1000 -4 4]),  yline([-1 1]),hold off, title('pleasant')

latencyRDKP_mean = min(find(cusumRDKP_mean< -thresh));
latencyRSVPP_mean = min(find(cusumRSVPP_mean> thresh));

% second: Unpleasant
RDKtempU = Exp2_data{1,1}{1,2}{2,4}.*10^13;
RSVPtempU = Exp2_data{2,1}{1,2}{2,4}.*10^13;

RDKtempU = bslcorr_percent(RDKtempU, 400:500);
RSVPtempU = bslcorr_percent(RSVPtempU, 400:500);

RDKtempU = RDKtempU(501:end);
RSVPtempU = RSVPtempU(501:end);

figure(201),
plot(taxis4plot, RDKtempU), hold on, plot(taxis4plot, RSVPtempU), hold off, title(' unpleasant')

cusumRDKU_mean = cumsum(RDKtempU)./std(abs(cumsum(RDKtempU)));
cusumRSVPU_mean = cumsum(RSVPtempU)./std(abs(cumsum(RSVPtempU)));

figure(202), plot(taxis4plot, cusumRDKU_mean), hold on, plot(taxis4plot, cusumRSVPU_mean), axis([-.02 1000 -4 4]),  yline([-thresh thresh]), hold off, title('unpleasant')

latencyRDKU_mean = min(find(cusumRDKU_mean < -thresh));
latencyRSVPU_mean = min(find(cusumRSVPU_mean > thresh));

%% fits for the single subjects, bootstrapped grand means
cd '/Users/andreaskeil/Dropbox (UFL)/ak_own/Ducmodel/DUC2'
clc

thresh = 1.5; 

 cusumRDKU = nan(5000, 983);
 cusumRSVPU = nan(5000, 983);
 cusumRDKP = nan(5000, 983);
 cusumRSVPP = nan(5000, 983);

%close all
load('Exp2_ModelingDataCorrected.mat')

for bootstrapInd = 1:5000
    % make a new grand mean for each stimulus
    subjectvec = randi(32, 32, 1);

    RSVP_unpleasants = Exp2_data{2,1}{2,2}{2,4};
    RDK_unpleasants = Exp2_data{1,1}{2,2}{2,4};

    RSVP_pleasants = Exp2_data{2,1}{2,2}{2,3};
    RDK_pleasants = Exp2_data{1,1}{2,2}{2,3};

    RSVP_neutrals = Exp2_data{2,1}{2,2}{2,2};
    RDK_neutrals = Exp2_data{1,1}{2,2}{2,2};


    RSVPtempU = mean(RSVP_unpleasants(subjectvec,:)*10^13);
    RDKtempU = mean(RDK_unpleasants(subjectvec, :).*10^13);

    RDKtempP = mean(RDK_pleasants(subjectvec,:)*10^13);
    RSVPtempP = mean(RSVP_pleasants(subjectvec,:)*10^13);

    % baseline correction of bootstrap means
    RDKtempU = bslcorr_percent(RDKtempU, 400:500);
    RDKtempU = RDKtempU(501:end); 

    RSVPtempU = bslcorr_percent(RSVPtempU, 400:500);
    RSVPtempU = RSVPtempU(501:end); 

    RDKtempP = bslcorr_percent(RDKtempP, 400:500);
    RDKtempP = RDKtempP(501:end); 

    RSVPtempP = bslcorr_percent(RSVPtempP, 400:500);
    RSVPtempP = RSVPtempP(501:end);

    cusumRDKU(bootstrapInd, :) = cumsum(RDKtempU)./std(abs(cumsum(RDKtempU)));
    cusumRSVPU(bootstrapInd, :) = cumsum(RSVPtempU)./std(abs(cumsum(RSVPtempU)));

    cusumRDKP(bootstrapInd, :) = cumsum(RDKtempP)./std(abs(cumsum(RDKtempP)));
    cusumRSVPP(bootstrapInd, :)  = cumsum(RSVPtempP)./std(abs(cumsum(RSVPtempP)));

     %figure(101), 
     % subplot(2,1,1), plot(abs(cusumRDKP)), hold on, plot(abs(cusumRSVPP)), yline(thresh), hold off
     % subplot(2,1,2), plot(abs(cusumRDKU)), hold on, plot(abs(cusumRSVPU)), yline(thresh), pause, hold off


   if find(abs(cusumRDKU(bootstrapInd, :)) >thresh), latencyRDKU(bootstrapInd) = find(abs(cusumRDKU(bootstrapInd, :)) >thresh, 1, "first"); end
   if find(abs(cusumRSVPU(bootstrapInd, :)) >thresh),  latencyRSVPU(bootstrapInd) = find(abs(cusumRSVPU(bootstrapInd, :)) >thresh, 1, "first"); end
   if find(abs(cusumRDKP(bootstrapInd, :)) >thresh), latencyRDKP(bootstrapInd) = find(abs(cusumRDKP(bootstrapInd, :)) >thresh, 1, "first"); end
   if find(abs(cusumRSVPP(bootstrapInd, :)) >thresh), latencyRSVPP(bootstrapInd) = find(abs(cusumRSVPP(bootstrapInd, :)) >thresh, 1, "first"); end

   if round(bootstrapInd./100) == bootstrapInd/100, disp(num2str(bootstrapInd)), end

end

%% BFS

[BFP] = log10(bootstrap2BF_z(latencyRDKP, latencyRSVPP, 1))
[BFU] = log10(bootstrap2BF_z(latencyRDKU, latencyRSVPU, 1))

%% figures
clc
orange = [1 0.5 0.0];
teal = [0.2 0.8 1];

taxis4plot = -23:1000/512:1900-4;

figure(901), 
subplot(2,1,2), title('Unpleasant', 'FontSize', 16 ), hold on,  xlabel('Time in ms', 'FontSize', 16 ), ylabel('Standardized Cumulative Sum', 'FontSize', 16 )
s1 = shadedErrorBar(taxis4plot, cusumRSVPU_mean, std(cusumRSVPU), 'transparent', 1,'patchSaturation',0.55, 'lineprops', {'-', 'Color', orange}); axis([-.02 1000 -4.1 4.5]),  yline([-thresh thresh]),
s1.mainLine.LineWidth = 3;
s1.patch.FaceColor = orange; hold on

s2 = shadedErrorBar(taxis4plot, cusumRDKU_mean, std(cusumRDKU), 'transparent', 1,'patchSaturation',0.55, 'lineprops', {'-', 'Color', teal}); axis([-.02 1000 -4.1 4.5]),  yline([-thresh thresh]),
s2.mainLine.LineWidth = 3;
s2.patch.FaceColor = teal;
xline(taxis4plot(latencyRSVPU_mean))
legend('RSVP', 'RDK', 'fontsize', 18, 'Location', "northwest")
set(gca, 'FontSize', 12)

subplot(2,1,1),  title('Pleasant', 'FontSize', 16 ), hold on,  xlabel('Time in ms', 'FontSize', 16 ), ylabel('Standardized Cumulative Sum', 'FontSize', 16 )
s1 = shadedErrorBar(taxis4plot, cusumRSVPP_mean, std(cusumRSVPP), 'transparent', 1,'patchSaturation',0.55, 'lineprops', {'-', 'Color', orange}); axis([-.02 1000 -4.1 4.5]),  yline([-thresh thresh]),
s1.mainLine.LineWidth = 3;
s1.patch.FaceColor = orange; hold on

s2 = shadedErrorBar(taxis4plot, cusumRDKP_mean', std(cusumRDKP), 'transparent', 1,'patchSaturation',0.55, 'lineprops', {'-', 'Color', teal}); axis([-.02 1000 -4.1 4.5]),  yline([-thresh thresh]),
s2.mainLine.LineWidth = 3;
s2.patch.FaceColor = teal;
xline(taxis4plot(latencyRSVPP_mean))
xline(taxis4plot(latencyRDKP_mean))
legend('RSVP', 'RDK', 'fontsize', 18, 'Location', "northwest")
set(gca, 'FontSize', 12)
