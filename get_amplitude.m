%% go to DR data

load('dataOut_Gratings_POOLED.mat')
load('gaussFit_results_POOLED.mat')

%get cells
cells_chosen = [1:1:length(fit_oriBW)];
cells_chosen = cells_chosen(dataOut.stats.global.responsive_cells_p001_fdr_average_index' & fit_oriBW>4 & fit_oriBW<90);

%get amplitude (2 ways)
cells_chosen_amp = max(dataOut.stats.global.response_ACTUAL_avg_vals,[],2)';
cells_chosen_amp = cells_chosen_amp(cells_chosen);

for i = 1:length(gaussFitResults)
    cells_chosen_fitAmp(i) = gaussFitResults{i}.a1;
end
cells_chosen_fitAmp = cells_chosen_fitAmp(cells_chosen);

stim_oris = [0 0 0 0 0 15 15 15 15 15 30 30 30 30 30 45 45 45 45 45 60 60 60 60 60 75 75 75 75 75 ...
    90 90 90 90 90 105 105 105 105 105 120 120 120 120 120 ...
    135 135 135 135 135 150 150 150 150 150 165 165 165 165 165];
stim_sfs = repmat([.05 .1 .2 .3 .4],1,12);
for i = 1:length(gaussFitResults)
    ori = gaussFitResults{i}.b1;
    [c,v] = min(abs(stim_oris-ori));
    ori_est = stim_oris(v);
    stim = find(stim_oris==ori_est & stim_sfs==pref_sf_estimate(i));
    cells_chosen_closestFitAmp(i) = dataOut.stats.global.response_ACTUAL_avg_vals(i,stim);
end
cells_chosen_closestFitAmp = cells_chosen_closestFitAmp(cells_chosen);

%get BW for chosen cells
cells_chosen_BW = fit_oriBW(cells_chosen);

%get average
mean_cells_chosen = mean(cells_chosen_amp);
mean_cells_chosen_fit = mean(cells_chosen_fitAmp);

range_bottom = 4;
range_top =14;
mean_range = mean(cells_chosen_amp(cells_chosen_BW<range_top & cells_chosen_BW>range_bottom));
mean_range_fit = mean(cells_chosen_fitAmp(cells_chosen_BW<range_top & cells_chosen_BW>range_bottom));

save('DR_amplitude_data.mat','cells_chosen_fitAmp','cells_chosen_closestFitAmp','cells_chosen','cells_chosen_BW');

%% go to Control data

load('dataOut_Gratings_POOLED.mat')
load('gaussFit_results_POOLED.mat')

%get cells
cells_chosen = [1:1:length(fit_oriBW)];
cells_chosen = cells_chosen(dataOut.stats.global.responsive_cells_p001_fdr_average_index' & fit_oriBW>4 & fit_oriBW<90);

%get amplitude (3 ways)
cells_chosen_amp = max(dataOut.stats.global.response_ACTUAL_avg_vals,[],2)';
cells_chosen_amp = cells_chosen_amp(cells_chosen);

for i = 1:length(gaussFitResults)
    cells_chosen_fitAmp(i) = gaussFitResults{i}.a1;
end
cells_chosen_fitAmp = cells_chosen_fitAmp(cells_chosen);

stim_oris = [0 0 0 0 0 15 15 15 15 15 30 30 30 30 30 45 45 45 45 45 60 60 60 60 60 75 75 75 75 75 ...
    90 90 90 90 90 105 105 105 105 105 120 120 120 120 120 ...
    135 135 135 135 135 150 150 150 150 150 165 165 165 165 165];
stim_sfs = repmat([.05 .1 .2 .3 .4],1,12);
for i = 1:length(gaussFitResults)
    ori = gaussFitResults{i}.b1;
    [c,v] = min(abs(stim_oris-ori));
    ori_est = stim_oris(v);
    stim = find(stim_oris==ori_est & stim_sfs==pref_sf_estimate(i));
    cells_chosen_closestFitAmp(i) = dataOut.stats.global.response_ACTUAL_avg_vals(i,stim);
end
cells_chosen_closestFitAmp = cells_chosen_closestFitAmp(cells_chosen);
%get BW for chosen cells
cells_chosen_BW = fit_oriBW(cells_chosen);

%get average
mean_cells_chosen = mean(cells_chosen_amp);
mean_cells_chosen_fit = mean(cells_chosen_fitAmp);

range_bottom = 4;
range_top =14;
mean_range = mean(cells_chosen_amp(cells_chosen_BW<range_top & cells_chosen_BW>range_bottom));
mean_range_fit = mean(cells_chosen_fitAmp(cells_chosen_BW<range_top & cells_chosen_BW>range_bottom));

save('CT_amplitude_data.mat','cells_chosen_fitAmp','cells_chosen_closestFitAmp','cells_chosen','cells_chosen_BW');


%% compare

DR = load('DR_amplitude_data.mat');
CT = load('CT_amplitude_data.mat');

[h,p]=ttest2(DR.cells_chosen_closestFitAmp,CT.cells_chosen_closestFitAmp);

range_bottom = 4;
range_top =14;

vals_range_DR = DR.cells_chosen_closestFitAmp(DR.cells_chosen_BW<range_top & DR.cells_chosen_BW>range_bottom);
vals_range_fit_DR = DR.cells_chosen_fitAmp(DR.cells_chosen_BW<range_top & DR.cells_chosen_BW>range_bottom);

vals_range_CT = CT.cells_chosen_closestFitAmp(CT.cells_chosen_BW<range_top & CT.cells_chosen_BW>range_bottom);
vals_range_fit_CT = CT.cells_chosen_fitAmp(CT.cells_chosen_BW<range_top & CT.cells_chosen_BW>range_bottom);

[h,p]=ttest2(vals_range_DR,vals_range_CT);

mean(vals_range_DR)
mean(vals_range_CT)


