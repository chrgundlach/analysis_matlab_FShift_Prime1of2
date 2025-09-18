function [BF] = bootstrap2BF_z(dist1,dist2, plotflag)
% This function takes two bootstrapped distributions (>1000 draws needed)
% dist1 and dist2, which reflect an effect model (dist1) and a null model
% (dist2). It computes then the BF as posterior oddsover prior odds foe the
% model entering the bootstrap in dist1. in many cases dist2 will be a
% permutation distribution representing the null model.
% dist1 and dist2 are vectors

% to avoid issues with scaling, do a joint z-normalization of both
% bootstrapped distributions together, thus maintaining any differences in
% mean and spread

% note: copyright by Andreas Keil (reveived 01/24)

dist1 = column(dist1); % make sure they are column vectors
dist2 = column(dist2);


% temp = zscore([dist1; dist2]);

% dist1z = temp(1:length(dist1));
% dist2z = temp(length(dist1)+1:end);


% rather normalize by second distribution
dist1z = (dist1-mean(dist2))/std(dist2);
dist2z = (dist2-mean(dist2))/std(dist2);


% fit the bootstrapped distributions as normal distributions
PD1 = fitdist(column(dist1z),'normal');
PD2 = fitdist(column(dist2z),'normal');

% apply the normal distributions
% x_values = -10:.01:10;
x_values = -(abs(PD1.mu)+5*PD1.sigma):.01:(abs(PD1.mu)+5*PD1.sigma);
x_values = -(abs(PD1.mu)+5*PD1.sigma):.001:(abs(PD1.mu)+5*PD1.sigma);
y1 = pdf(PD1,x_values);
y2 = pdf(PD2,x_values);

% plot for control, then turn off plotflag
if plotflag == 1
    figure;
    subplot(2,1,1), histogram(dist1z, 100, 'Normalization','pdf')
    title('empirical')
    hold on
    plot(x_values, y1, 'LineWidth',3)
    xline(0)
    xlim([-1 1]*max(abs(x_values)))

    subplot(2,1,2), histogram(dist2z, 100, 'Normalization','pdf')
    title('random')
    hold on
    plot(x_values, y2, 'LineWidth',3)
    xline(0)
    xlim([-1 1]*max(abs(x_values)))
end

posteriorsignedlikelyhood_effect = sum(y1(x_values > 0))./round(1/((max(x_values)-min(x_values))/numel(x_values)));
signedlikelyhood_null = sum(y2(x_values > 0))./round(1/((max(x_values)-min(x_values))/numel(x_values)));

odds_posterior = posteriorsignedlikelyhood_effect./(1-posteriorsignedlikelyhood_effect);
odds_prior = signedlikelyhood_null./(1-signedlikelyhood_null);


BF = odds_posterior/odds_prior;


%% column function
function [datout] = column(datin)
datout = datin(:);

