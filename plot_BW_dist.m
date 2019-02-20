%loaded dataOut_Gratings_POOLED
%loaded gaussFit_results_POOLED
%got the dimensions of each session from pool_neurons.m

%or run this
%% Pool gauss fit and BW for all cells
clear all

session_order = {'2268_NC_180523_DR','2269_NC_180521_DR','2270_1R_180519_DR','2271_1R_180520_DR','2280_1L_180618_DR','2297_NC_180702_DR','2298_1R1L_180704_DR'};
i=1;
GdataOut{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/%s_Grating_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
GdataOut{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/%s_Grating_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
GdataOut{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/%s_Grating_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
GdataOut{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/%s_Grating_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
GdataOut{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/%s_Grating_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
GdataOut{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/%s_Grating_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
GdataOut{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/%s_Grating_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;

load('dataOut_Gratings_POOLED.mat');
load('gaussFit_results_POOLED.mat')



%% Get plot
start = 1;
for s = 1:size(dimensions,1)
    %get the BW for just this session
    sess_BW = fit_oriBW(start:start+dimensions(s,2)-1);    
    %get BW just for responsive neurons
    sess_BW_sig = sess_BW(GdataOut{1,s}.dataOut.stats.global.responsive_cells_p001_fdr_average);
    sess_BW_sig(sess_BW_sig<4)=[];
    sess_BW_sig(sess_BW_sig>90) = [];
    store_sess_BW_sig{s} = sess_BW_sig;
    start = start+dimensions(s,2);
end

bins = [0:5:90];
session_BW_count=zeros(length(bins)-1,size(dimensions,1));
session_BW_count_prop=zeros(length(bins)-1,size(dimensions,1));
for s = 1:size(dimensions,1)
    h=histogram(store_sess_BW_sig{s},bins);
    session_BW_count(:,s) = h.Values';
    session_BW_count_prop(:,s) = h.Values'./length(store_sess_BW_sig{s});
end

save('BW_dist_data.mat','session_BW_count','store_sess_BW_sig');

%plot raw count
b = bar(session_BW_count,'stacked');
bar_width= b.BarWidth;
bar_width_new = bar_width/2;
set(gca,'XTick',[bar_width_new*1.25:bar_width_new*2.5:20])
xticklabels([0:5:90])
xlim([0 19])
xlabel('Orientation BW (degrees)')
ylabel('Number of Neurons')
legend({'Mouse 1','Mouse 2','Mouse 3','Mouse 4','Mouse 5','Mouse 6'})
set(gca,'FontSize',16)
title(sprintf('Orientation Bandwidth for Preferred SF (%i neurons)',sum(sum(session_BW_count))))
saveas(gcf,'BW_dist_count.fig')
saveas(gcf,'BW_dist_count.png')

%plot total proportion
session_prop = sum(session_BW_count,2)/sum(sum(session_BW_count,2));

save('BW_dist_data.mat','session_prop','-append')

b = bar(session_prop,'stacked','k');
bar_width= b.BarWidth;
bar_width_new = bar_width/2;
set(gca,'XTick',[bar_width_new*1.25:bar_width_new*2.5:20])
xticklabels([0:5:90])
xlim([0 19])
xlabel('Orientation BW (degrees)')
ylabel('Proportion of neurons')
set(gca,'FontSize',16)
title(sprintf('Orientation Bandwidth for Preferred SF (%i neurons)',sum(sum(session_BW_count))))
saveas(gcf,'BW_dist_prop.fig')
saveas(gcf,'BW_dist_prop.png')

% %plot proportion
% b = bar(session_BW_count_prop,'stacked');
% bar_width= b.BarWidth;
% bar_width_new = bar_width/2;
% set(gca,'XTick',[bar_width_new*1.25:bar_width_new*2.5:20])
% xticklabels([0:5:90])
% xlim([0 19])
% xlabel('Orientation BW (degrees)')
% ylabel('Number of Neurons')
% legend({'Mouse 1','Mouse 2','Mouse 3','Mouse 4','Mouse 5','Mouse 6'})
% set(gca,'FontSize',16)
% title(sprintf('Orientation Bandwidth for Preferred SF (%i neurons)',sum(sum(session_BW_count))))
% saveas(gcf,'BW_dist_count.fig')
% saveas(gcf,'BW_dist_count.png')

