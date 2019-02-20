%% above df/f threshold

home = pwd;
cd ..
G = load('dataOut_Gratings_POOLED.mat');
N = load('dataOut_NatScenes_POOLED.mat');
cd(home)

% G = load('dataOut_Gratings_POOLED_DR.mat');
% N = load('dataOut_NatScenes_POOLED_DR.mat');

G_ind = G.dataOut.stats.global.responsive_cells_p001_fdr_average_index;
G_cells = G.dataOut.stats.global.responsive_cells_p001_fdr_average;
N_ind = N.dataOut.stats.global.responsive_cells_p001_fdr_average_index;
N_cells = N.dataOut.stats.global.responsive_cells_p001_fdr_average;
NSO_ind = zeros(N.dataOut.totalNumCells,1);
NSO_ind(G_ind==0 & N_ind==1)= 1;
NSO_cells = find(NSO_ind==1);
G_NS_ind = zeros(N.dataOut.totalNumCells,1);
G_NS_ind(G_ind==1 & N_ind==1)= 1;
G_NS_cells = find(G_NS_ind==1);

%%%%%%%%%%%%
% for NS
thresh = .1;

stim_reliability = [];
stim_reliability_store=[];
for n=1:N.dataOut.totalNumCells
    n_resps = squeeze(N.dataOut.responseMatrix_1(:,n,:));
    n_resps(N.dataOut.isRemovedBlock(:,n,:)==1) = NaN;
    for stim = 1:N.dataOut.totalNumStimuli
        stim_reliability_store{n,stim} = find(n_resps(:,stim)>thresh);
        stim_reliability(n,stim) = length(find(n_resps(:,stim)>thresh))/length(find(~isnan(n_resps(:,stim))));
    end
    
end

stim_reliability(N.dataOut.stats.global.response_average_pval_fdr>0.01) = NaN;
cells_resps = N.dataOut.stats.global.response_ACTUAL_avg_vals;
cells_resps(isnan(stim_reliability)) = NaN;
[cells_maxs,cells_inds] = max(cells_resps,[],2);
all_cell_reliability_max = [];
for n =1:N.dataOut.totalNumCells
    all_cell_reliability_max(n,1) = stim_reliability(n,cells_inds(n));
end
% stim_reliability(N.dataOut.stats.global.response_average_pval_fdr>0.01) = NaN;
% cells_resps = N.dataOut.stats.global.response_ACTUAL_avg_vals;
% cells_resps(isnan(stim_reliability)) = NaN;
% [cells_maxs,cells_inds] = max(cells_resps,[],2);
% cells_inds(isnan(cells_maxs))=NaN;
% store_cellinds = cells_inds;
% store_cellinds(isnan(cells_inds)) = [];
% all_cell_reliability_max = stim_reliability(:,cells_inds);
% all_cell_reliability_max = diag(all_cell_reliability_max);
%selected_cells = N.dataOut.stats.global.responsive_cells_p001_fdr_average;
%selected_cells = G_NS_cells;
selected_cells = NSO_cells;
%C_diag = selected_cells_reliability;
%C_diag = get_diagonal(get_diagonal>0);
%N_diag = selected_cells_reliability;

selected_cells_reliability_all = stim_reliability(selected_cells,:);
selected_cells_reliability_all = selected_cells_reliability_all(:);
selected_cells_reliability_all(isnan(selected_cells_reliability_all))=[];
N_diag = selected_cells_reliability_all;
%[h,p]=kstest2(N_diag,DR_save)

figure
histogram(N_diag,[0:0.05:1])
xlabel('reliability')
ylabel('number of neurons')
title(sprintf('DR NS reliability df/f thresh %.2f',thresh))
saveas(gcf,sprintf('D_NS_reliability_dffThresh_%.2f.fig',thresh))
saveas(gcf,sprintf('D_NS_reliability_dffThresh_%.2f.png',thresh))


% NSO_reliability_all = N_diag;
% GR_NS_reliability_all = N_diag;
% [h,p] = kstest2(NSO_reliability_all,GR_NS_reliability_all);
% 
% DR_reliability_all = N_diag;
% C_reliability_all = N_diag;
% [h,p] = kstest2(DR_reliability_all,C_reliability_all);

save('reliability_set2_DR_NS_all_SIG_data.mat','stim_reliability','stim_reliability_store')

%%%% for Gratings
thresh = .1;

stim_reliability = [];
stim_reliability_store=[];
for n=1:G.dataOut.totalNumCells
    n_resps = squeeze(G.dataOut.responseMatrix_1(:,n,:));
    n_resps(G.dataOut.isRemovedBlock(:,n,:)==1) = NaN;
    for stim = 1:G.dataOut.totalNumStimuli
        stim_reliability_store{n,stim} = find(n_resps(:,stim)>thresh);
        stim_reliability(n,stim) = length(find(n_resps(:,stim)>thresh))/length(find(~isnan(n_resps(:,stim))));
    end
    
end

stim_reliability(G.dataOut.stats.global.response_average_pval_fdr>0.01) = NaN;
cells_resps = G.dataOut.stats.global.response_ACTUAL_avg_vals;
cells_resps(isnan(stim_reliability)) = NaN;
[cells_maxs,cells_inds] = max(cells_resps,[],2);
all_cell_reliability_max = [];
for n =1:G.dataOut.totalNumCells
    all_cell_reliability_max(n,1) = stim_reliability(n,cells_inds(n));
end
stim_reliability(G.dataOut.stats.global.response_average_pval_fdr>0.01) = NaN;
save('reliability_set2_DR_GR_all_SIG_data.mat','stim_reliability','stim_reliability_store')

%% Get above reliability for top response of each responsive cell and plot distribution
home = pwd;
cd ..
G = load('dataOut_Gratings_POOLED.mat');
N = load('dataOut_NatScenes_POOLED.mat');
cd(home)
GR = load('reliability_set2_DR_GR_all_SIG_data.mat');
NR = load('reliability_set2_DR_NS_all_SIG_data.mat');

reliability_all = vertcat(NR.stim_reliability', GR.stim_reliability')';
G_ind = G.dataOut.stats.global.responsive_cells_p001_fdr_average_index;
G_cells = G.dataOut.stats.global.responsive_cells_p001_fdr_average;
N_ind = N.dataOut.stats.global.responsive_cells_p001_fdr_average_index;
N_cells = N.dataOut.stats.global.responsive_cells_p001_fdr_average;
responsive_cells = find(G_ind==1 | N_ind ==1);


R_maxs = [];
for n = 1: G.dataOut.totalNumCells
    n_resps = [N.dataOut.stats.global.response_ACTUAL_avg_vals(n,:) G.dataOut.stats.global.response_ACTUAL_avg_vals(n,:)];
    n_pvals = [N.dataOut.stats.global.response_average_pval_fdr(n,:) G.dataOut.stats.global.response_average_pval_fdr(n,:)];
    n_resps_sig = n_resps;
    n_resps_sig(n_pvals>0.01) = NaN;
    [n_max,n_max_ind]=max(n_resps_sig);
    n_max_R = reliability_all(n,n_max_ind);
    R_maxs(n) = n_max_R;
end
R_maxs_responsiveCells = R_maxs(responsive_cells)';
save('reliability_best_stim_data.mat','R_maxs','R_maxs_responsiveCells','responsive_cells')

h=histogram(R_maxs_responsiveCells,[0:.05:1]);
bins = h.Values;
bins_prop = bins./sum(bins);
figure
bar([0.025:0.05:0.975],bins_prop)
ylim([0 0.15])
xlabel('reliability')
ylabel('proportion of cells')
title('reliability to preferred stim (NS or GR)')
saveas(gcf,'reliability_best_stim_set1.fig')
saveas(gcf,'reliability_best_stim_set1.png')
saveas(gcf,'reliability_best_stim_set1.svg')

% stim_reliability = [];
% for n=1:DR.dataOut.totalNumCells
%     n_resps = squeeze(DR.dataOut.responseMatrix_1(:,n,:));
%     n_resps(DR.dataOut.isRemovedBlock(:,n,:)==1) = NaN;
%     for stim = 1:DR.dataOut.totalNumStimuli
%         stim_reliability_store{n,stim} = find(n_resps(:,stim)>thresh);
%         stim_reliability(n,stim) = length(find(n_resps(:,stim)>thresh))/length(find(~isnan(n_resps(:,stim))));
%     end
% end
% 
% stim_reliability(DR.dataOut.stats.global.response_average_pval_fdr>0.01) = NaN;
% cells_resps = DR.dataOut.stats.global.response_ACTUAL_avg_vals;
% cells_resps(isnan(stim_reliability)) = NaN;
% [cells_maxs,cells_inds] = max(cells_resps,[],2);
% cells_inds(isnan(cells_maxs))=NaN;
% test = cells_inds;
% test(isnan(cells_inds)) = [];
% test3 = DR.dataOut.stats.global.responsive_cells_p001_fdr_average;
% test2 = stim_reliability(test3,test);
% get_diagonal = diag(test2);
% DR_diag = get_diagonal(get_diagonal>0);
% 
% figure
% histogram(DR_diag,[0:0.05:1])
% xlabel('reliability')
% ylabel('number of neurons')
% title(sprintf('DR reliability df/f thresh %.2f',thresh))
% saveas(gcf,sprintf('DR_reliability_dffThresh_%.2f.fig',thresh))
% saveas(gcf,sprintf('DR_reliability_dffThresh_%.2f.png',thresh))
% 
% save(sprintf('C_DR_reliability_dffThresh_%.2f_data.mat',thresh),'C_diag','DR_diag')
% [h,p]=kstest2(C_diag,DR_diag)

%% above std of global baseline

session_order = {'2210_NC_180629','2231_1R1L_180628','2253_1R_180621','2260_1R_180622','2260_2R_180517','2263_2R_180624'};
i=1;
NSdataOuts{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/Control/%s_suite2p_processed/processed_suite2p/%s_analysis/NatScenes/%s_NatScenes_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
NSdataOuts{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/Control/%s_suite2p_processed/processed_suite2p/%s_analysis/NatScenes/%s_NatScenes_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
NSdataOuts{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/Control/%s_suite2p_processed/processed_suite2p/%s_analysis/NatScenes/WITHOUTEYETRACKING/%s_NatScenes_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
NSdataOuts{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/Control/%s_suite2p_processed/processed_suite2p/%s_analysis/NatScenes/%s_NatScenes_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
NSdataOuts{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/Control/%s_suite2p_processed/processed_suite2p/%s_analysis/NatScenes/%s_NatScenes_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
NSdataOuts{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/Control/%s_suite2p_processed/processed_suite2p/%s_analysis/NatScenes/%s_NatScenes_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;

i=1;
GdataOut{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/Control/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/%s_Grating_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
GdataOut{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/Control/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/%s_Grating_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
GdataOut{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/Control/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/WITHOUTEYETRACKING/%s_Grating_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
GdataOut{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/Control/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/%s_Grating_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
GdataOut{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/Control/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/%s_Grating_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
GdataOut{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/Control/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/%s_Grating_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;


std_times = 2;

for s = 1:length(GdataOut)
    %load dataOuts and re-run portion to get values that went into baseline
    G=GdataOut{s};
    N=NSdataOuts{s};
    dataOut = klab_dataOut_addstats_PS(N.dataOut,N.cfg_stats);
    session_order{2,s} = dataOut;
%     baseline_std = cellfun(@std,dataOut.stats.globalBaselines_store);
%     baseline_means = cellfun(@mean,dataOut.stats.globalBaselines_store);
%     basline_std_times_thresh = (baseline_means+std_times*baseline_std)./baseline_means;
    session_stim_reliability = [];
    for n = 1:dataOut.totalNumCells
        %%%%%%%%%get stds threshold from gray screens
        %a trials for a cell
        trials = dataOut.stats.globalBaselines_store(:,n);
        %get mean of gray screens for each trial
        trials_mean=cellfun(@mean,trials);
        %subtract and divide baseline from response of cell
        trials_std = std_times*cellfun(@std,trials);
        response_dff_thresh = std_times*trials_std./trials_mean;
        
        %%%%get neuron reliability for each stim
        n_resps = squeeze(dataOut.responseMatrix_1(:,n,:));
        n_resps(dataOut.isRemovedBlock(:,n,:)==1) = NaN;
      
        n_stim_reliability = [];
        for stim = 1:dataOut.totalNumStimuli
            stim_resps = n_resps(:,stim);
            trials_good = sum(~isnan(stim_resps));
            trials_passed = sum(stim_resps>response_dff_thresh);
            n_stim_reliability(stim) = trials_passed/trials_good;
        end
        session_stim_reliability(n,:) = n_stim_reliability;
    end
    session_stim_reliability_sig = session_stim_reliability;
    session_stim_reliability_sig(dataOut.stats.global.response_average_pval_fdr>0.01) = NaN;
    A.session_stim_reliability = session_stim_reliability;
    A.session_stim_reliability_sig=session_stim_reliability_sig;
    
    %get specific groups of cells    
    G_cells = G.dataOut.stats.global.responsive_cells_p001_fdr_average;
    G_ind = false(1,G.dataOut.totalNumCells);
    G_ind(G_cells)=1;
    
    N_cells = N.dataOut.stats.global.responsive_cells_p001_fdr_average;
    N_ind = false(1,N.dataOut.totalNumCells);
    N_ind(N_cells) = 1;
    
    NSO_ind = zeros(N.dataOut.totalNumCells,1);
    NSO_ind(G_ind==0 & N_ind==1)= 1;
    NSO_cells = find(NSO_ind==1);
    G_NS_ind = zeros(N.dataOut.totalNumCells,1);
    G_NS_ind(G_ind==1 & N_ind==1)= 1;
    G_NS_cells = find(G_NS_ind==1);
    
    A.G_cells = G_cells;
    A.N_cells = N_cells;
    A.NSO_cells = NSO_cells;
    A.G_NS_cells = G_NS_cells;
    
    session_order{3,s} = A;
    %get distribution for all NS responsive cells
    %chosen_reliability = session_stim_reliability_sig(N_cells,:);    
    
end

%%%%%STOP HERE AND DO REST BY HAND

reliability_values_allSess = [];
for s = 1:length(GdataOut)
    cells_chosen = session_order{3,s}.G_NS_cells;
%     cells_chosen = session_order{3,s}.NSO_cells;
%     cells_chosen = session_order{3,s}.G_NS_cells;
    rel_chosen = session_order{3,s}.session_stim_reliability_sig(cells_chosen,:);
    rel_sig_vector = (rel_chosen(~isnan(rel_chosen)));
    reliability_values_allSess = [reliability_values_allSess rel_sig_vector'];
end

all_NS_reliability = reliability_values_allSess;
all_NSO_reliability = reliability_values_allSess;
all_G_NS_reliability = reliability_values_allSess;

session_order2 = session_order;
session_order2(2,:) = [];

save('all_sess_sdvThresh_data_control.mat','session_order2','std_times','all_NS_reliability','all_NSO_reliability','all_G_NS_reliability')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% now for DR
session_order = {'2268_NC_180523_DR','2269_NC_180521_DR','2270_1R_180519_DR','2271_1R_180520_DR','2280_1L_180618_DR','2297_NC_180702_DR','2298_1R1L_180704_DR'};
i=1;
NSdataOuts{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/NatScenes/%s_NatScenes_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
NSdataOuts{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/NatScenes/%s_NatScenes_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
NSdataOuts{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/NatScenes/%s_NatScenes_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
NSdataOuts{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/NatScenes/%s_NatScenes_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
NSdataOuts{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/NatScenes/%s_NatScenes_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
NSdataOuts{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/NatScenes/%s_NatScenes_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
NSdataOuts{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/NatScenes/%s_NatScenes_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;

i=1;
GdataOut{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/%s_Grating_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
GdataOut{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/%s_Grating_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
GdataOut{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/%s_Grating_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
GdataOut{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/%s_Grating_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
GdataOut{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/%s_Grating_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
GdataOut{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/%s_Grating_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
GdataOut{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/%s_Grating_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;


std_times = 2;

for s = 1:length(GdataOut)
    %load dataOuts and re-run portion to get values that went into baseline
    G=GdataOut{s};
    N=NSdataOuts{s};
    dataOut = klab_dataOut_addstats_PS(N.dataOut,N.cfg_stats);
    session_order{2,s} = dataOut;
%     baseline_std = cellfun(@std,dataOut.stats.globalBaselines_store);
%     baseline_means = cellfun(@mean,dataOut.stats.globalBaselines_store);
%     basline_std_times_thresh = (baseline_means+std_times*baseline_std)./baseline_means;
    session_stim_reliability = [];
    for n = 1:dataOut.totalNumCells
        %%%%%%%%%get stds threshold from gray screens
        %a trials for a cell
        trials = dataOut.stats.globalBaselines_store(:,n);
        %get mean of gray screens for each trial
        trials_mean=cellfun(@mean,trials);
        %subtract and divide baseline from response of cell
        trials_std = std_times*cellfun(@std,trials);
        response_dff_thresh = std_times*trials_std./trials_mean;
        
        %%%%get neuron reliability for each stim
        n_resps = squeeze(dataOut.responseMatrix_1(:,n,:));
        n_resps(dataOut.isRemovedBlock(:,n,:)==1) = NaN;
      
        n_stim_reliability = [];
        for stim = 1:dataOut.totalNumStimuli
            stim_resps = n_resps(:,stim);
            trials_good = sum(~isnan(stim_resps));
            trials_passed = sum(stim_resps>response_dff_thresh);
            n_stim_reliability(stim) = trials_passed/trials_good;
        end
        session_stim_reliability(n,:) = n_stim_reliability;
    end
    session_stim_reliability_sig = session_stim_reliability;
    session_stim_reliability_sig(dataOut.stats.global.response_average_pval_fdr>0.01) = NaN;
    A.session_stim_reliability = session_stim_reliability;
    A.session_stim_reliability_sig=session_stim_reliability_sig;
    
    %get specific groups of cells    
    G_cells = G.dataOut.stats.global.responsive_cells_p001_fdr_average;
    G_ind = false(1,G.dataOut.totalNumCells);
    G_ind(G_cells)=1;
    
    N_cells = N.dataOut.stats.global.responsive_cells_p001_fdr_average;
    N_ind = false(1,N.dataOut.totalNumCells);
    N_ind(N_cells) = 1;
    
    NSO_ind = zeros(N.dataOut.totalNumCells,1);
    NSO_ind(G_ind==0 & N_ind==1)= 1;
    NSO_cells = find(NSO_ind==1);
    G_NS_ind = zeros(N.dataOut.totalNumCells,1);
    G_NS_ind(G_ind==1 & N_ind==1)= 1;
    G_NS_cells = find(G_NS_ind==1);
    
    A.G_cells = G_cells;
    A.N_cells = N_cells;
    A.NSO_cells = NSO_cells;
    A.G_NS_cells = G_NS_cells;
    
    session_order{3,s} = A;
    %get distribution for all NS responsive cells
    %chosen_reliability = session_stim_reliability_sig(N_cells,:);    
    
end

%%%%%STOP HERE AND DO REST BY HAND

reliability_values_allSess = [];
for s = 1:length(GdataOut)
    cells_chosen = session_order{3,s}.N_cells;
%     cells_chosen = session_order{3,s}.NSO_cells;
%     cells_chosen = session_order{3,s}.G_NS_cells;
    rel_chosen = session_order{3,s}.session_stim_reliability_sig(cells_chosen,:);
    rel_sig_vector = (rel_chosen(~isnan(rel_chosen)));
    reliability_values_allSess = [reliability_values_allSess rel_sig_vector'];
end

all_NS_reliability = reliability_values_allSess;
all_NSO_reliability = reliability_values_allSess;
all_G_NS_reliability = reliability_values_allSess;

session_order2 = session_order;
session_order2(2,:) = [];

save('all_sess_sdvThresh_data_DR.mat','session_order2','std_times','all_NS_reliability','all_NSO_reliability','all_G_NS_reliability')
