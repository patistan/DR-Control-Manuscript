%% Get grating - NS selectivity

N = load('dataOut_NatScenes_POOLED.mat');
G = load('dataOut_Gratings_POOLED.mat');
sig_resp_ind = NaN(N.dataOut.totalNumCells,2);

% NG_selec_all = [];
for c = 1:N.dataOut.totalNumCells
    %get max response for NS
    c_resps_N = N.dataOut.stats.global.response_ACTUAL_avg_vals(c,:);
    c_resps_N_sig = c_resps_N(N.dataOut.stats.global.response_average_pval_fdr(c,:)<0.01);
    if isempty(c_resps_N_sig)
        sig_resp_ind(c,1)=0;
    else
        sig_resp_ind(c,1)= max(c_resps_N_sig);
    end
    
    %get max response for grating
    c_resps_G = G.dataOut.stats.global.response_ACTUAL_avg_vals(c,1:60);
    c_resps_G_sig = c_resps_G(G.dataOut.stats.global.response_average_pval_fdr(c,1:60)<0.01);
    if isempty(c_resps_G_sig)
        sig_resp_ind(c,2)=0;
    else
        sig_resp_ind(c,2)= max(c_resps_G_sig);
    end

%     %get and store difference
%     NS_selec = (c_resps_N - c_resps_G)/(c_resp_G+c_resp_N);
%     NG_selec_all = [NG_selec_all; NS_selec];
end

NG_selec_all = -diff(sig_resp_ind,1,2)./sum(sig_resp_ind,2);

save('NG_selec_all.mat','NG_selec_all')

figure
h=histogram(NG_selec_all,[-1:.1:1]);
bin_values = h.Values;
figure
bar([-.95:.1:.95],[bin_values/sum(bin_values)])
title(sprintf('all responsive cells (n=%i)',sum(bin_values)))
xlabel({'Natural Scene-Grating Selectivity','(max_N - max_G)/(max_N + max_G)','1 = response only to NatScene'})
ylabel('# of Responsive Neurons')
ylabel('% of Responsive Neurons')
saveas(gcf,'NGSI_responsive_neurons_hist.fig')
saveas(gcf,'NGSI_responsive_neurons_hist.png')
saveas(gcf,'NGSI_responsive_neurons_bar.fig')
saveas(gcf,'NGSI_responsive_neurons_bar.png')
