
type = 'DR';
G = load('dataOut_Gratings_POOLED.mat');
N = load('dataOut_NatScenes_POOLED.mat');

G_ind = G.dataOut.stats.global.responsive_cells_p001_fdr_average_index;
G_cells = G.dataOut.stats.global.responsive_cells_p001_fdr_average';
SFs = repmat([1 2 3 4 5],1,12);

pref_sf=[];
resp_sfs_store=[];
for n = G_cells
    n_sig_resps = zeros(1,60);
    n_sig_resps(G.dataOut.stats.global.response_average_pval_fdr(n,1:60)<0.01) = 1;
    n_pref_resp = max(G.dataOut.stats.global.response_ACTUAL_avg_vals(n,n_sig_resps==1));
    if isempty(n_pref_resp)
        continue
    end
    n_pref_stim = find(G.dataOut.stats.global.response_ACTUAL_avg_vals(n,1:60)==n_pref_resp);
    pref_sf(1,n) = SFs(n_pref_stim);
    resp_sfs{1,n} = unique(SFs(n_sig_resps==1));
    resp_sfs_store = [resp_sfs_store unique(SFs(n_sig_resps==1))];
end

%%%plotting things
bins_pref = [];
bins_resp=[];
for SF = 1:5
    bins_pref(SF) = sum(pref_sf==SF);
    bins_resp(SF) = sum(resp_sfs_store==SF);
end

bins_pref_prop = bins_pref./sum(bins_pref);
bins_resp_prop = bins_resp./sum(bins_resp);

figure
bar(bins_pref_prop)
xlabel('preferred SF')
ylabel('proportion of neurons')
title(sprintf('%s preferred SF dist',type))
ylim([0 .6])
saveas(gcf,sprintf('pref_SF_dist_%s',type))
figure
bar(bins_resp_prop)
xlabel('SF')
ylabel('proportion of neurons')
ylim([0 .6])
title(sprintf('%s SF response dist',type))
saveas(gcf,sprintf('SF_dist_%s',type))

save(sprintf('get_SF_dist_data_%s.mat',type),'pref_sf','resp_sfs','resp_sfs_store','bins_pref','bins_resp')