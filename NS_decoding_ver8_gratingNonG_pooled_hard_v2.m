%% Decode using specific groups of neurons
clear all
close all

expID = 'POOLED';
folds = 10;
bins = 3;
pairs = [1 2;2 3;3 4;4 5;5 7;7 8;8 9;9 10;11 12;12 13;13 14;14 15;15 17;17 18;18 19;19 20];

home = pwd;
cd ../..
NS = load('dataOut_NatScenes_POOLED.mat');
G = load('dataOut_Gratings_POOLED.mat');
cd(home)

grating_cells = G.dataOut.stats.global.responsive_cells_p001_fdr_average';
grat_ind = zeros(1,G.dataOut.totalNumCells);
grat_ind(grating_cells) = 1;

nat_cells = NS.dataOut.stats.global.responsive_cells_p001_fdr_average';
nat_ind = zeros(1,NS.dataOut.totalNumCells);
nat_ind(nat_cells) = 1;

natOnly_cells = find(grat_ind == 0 & nat_ind ==1);
nsOnly_groups_all{1,1} = natOnly_cells;

%get random grating groups
for r = 1:30
    group_ind = sort(randsample(1:length(grating_cells),length(natOnly_cells)));
    grating_group = grating_cells(group_ind);
    grating_groups_all{r,1} = grating_group;
end

%% get each group
types = ({'Grating','nonGresp'});

%decode with groups identified when decoding natural scenes
for g=1:length(types)    
    if g==1 %for broad
        type = types{g};
        for i = 1:size(grating_groups_all,1)         
            selected_cells = grating_groups_all{i,1};
            [Pair_data,selected_cells] = NatScene_decoding_ver8_for10_simple_pooled_hard_v2(expID,type,selected_cells,bins,folds,pairs);
            for k = 1:size(pairs,1)
                mean_acc = mean(Pair_data{k,2}{1,bins});
                grating_groups_all{i,1} = selected_cells;
                grating_groups_all{i,k+1} = mean_acc;
            end
            fprintf('rep %d finished\n',i);            
        end

    elseif g==2 %for nonGresp
        type = types{g};
        for i = 1:size(nsOnly_groups_all,1)            
            selected_cells = nsOnly_groups_all{i,1};
            [Pair_data,selected_cells] = NatScene_decoding_ver8_for10_simple_pooled_hard_v2(expID,type,selected_cells,bins,folds,pairs);
            for k = 1:size(pairs,1)
                mean_acc = mean(Pair_data{k,2}{1,bins});
                nsOnly_groups_all{i,1} = selected_cells;
                nsOnly_groups_all{i,k+1} = mean_acc;
            end
            fprintf('rep %d finished\n',i);
        end
    end  
end

save('group_data_grat_nsOnly.mat','nsOnly_groups_all','grating_groups_all','pairs')

%plot groups against eachother using averages
mean_grat = mean(cell2mat(grating_groups_all(:,2:(size(pairs,1)+1))));
plot_nonG = cell2mat(nsOnly_groups_all(:,2:(size(pairs,1)+1)));
%mean_nsOnly = mean(cell2mat(nsOnly_groups_all(:,2:11)));
figure('Position',[100 200 800 600])
hold on
% if size(grating_groups_all,1)>1
%     scatter(repmat([1],1,30),cell2mat(grating_groups_all(:,2)))
% end
% if size(nsOnly_groups_all,1)>1
%     scatter(repmat([2],1,30),cell2mat(nsOnly_groups_all(:,2)))
% end
% scatter(repmat([3],1,30),cell2mat(diverse_groups_all(:,2)))
scatter(1:size(pairs,1),mean_grat,80,'b','filled','LineWidth',2)
scatter(1:size(pairs,1),plot_nonG,80,'g','filled','LineWidth',2)
ylim([0 1])
xlim([0 size(pairs,1)])
ylabel('Accuracy')
xlabel('Pair')
legend({'Grating','Non-Grating'},'Location','best')
plot([0,size(pairs,1)],[.6 .6],'--k')
set(gca,'FontSize',16)
title(sprintf('Decoding Accuracy (%d neurons per group)',size(grating_groups_all{1,1},2)))
saveas(gcf,sprintf('%s_NBdecoding_%istim_n%i_%ibins_allgroup_avg.fig',expID,20,size(grating_groups_all{1,1},2),bins));
saveas(gcf,sprintf('%s_NBdecoding_%istim_n%i_%ibins_allgroup_avg.png',expID,20,size(grating_groups_all{1,1},2),bins));

%% Decoding %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('specific_groups','dir')
    mkdir('specific_groups')
end

types = {'grating','nsOnly'};
groups{1} = group_grating;
groups{2} = group_nsOnly; 
%groups{3} = group_diverse;

for g = 1:length(groups)
    %do decoding with specific group
    type = types{g};
    selected_cells = groups{g};

    [Pair_data,selected_cells] = NatScene_decoding_ver8_for10_simple_pooled_hard_v2(expID,type,selected_cells,bins,folds,pairs);
%     accuracies_bin{g} = bins_accuracy{1,bins};
%     
%     %get confusion matrix    
%     real_v_guessed = AllFold_AllBins{2,bins}(:,2:3);
%     real_v_guessed_sorted = sort(real_v_guessed,2);
%     total_stim = size(AllFold_AllBins{1,bins}(1).RespMatrix,3);
%     confusion_matrix = zeros(total_stim);
%     for i = 1:total_stim
%         stimnum = real_v_guessed(find(real_v_guessed(:,1)== i),:);
%         for n = 1:total_stim
%             stimguessed = find(stimnum(:,2)==n);
%             confused_perc = length(stimguessed)/length(stimnum);
%             confusion_matrix(n,i) = confused_perc;
%         end
%     end
%     figure('Position',[100 200 1000 800])
%     imagesc(confusion_matrix)
%     h=colorbar;
%     ylabel(h, 'P(PS|S)')
%     colormap hot
%     caxis([0 1])
%     xlabel('True Stim (S)')
%     ylabel('Predicted Stim (PS)')
%     set(gca,'FontSize',16)
%     title(sprintf('Confusion Matrix(%d neurons, %s)',length(selected_cells),type))
%     saveas(gcf,sprintf('confusionMatrix_%istim_n%i_bin%i_%s.fig',total_stim,length(selected_cells),bins,type));
%     saveas(gcf,sprintf('confusionMatrix_%istim_n%i_bin%i_%s.png',total_stim,length(selected_cells),bins,type));
%     save(sprintf('confusionMatrix_%istim_n%i_bin%i_%s',total_stim,length(selected_cells),bins,type),'confusion_matrix','bins');
%     close all
end

% %plot groups against eachother
% figure('Position',[100 200 800 600])
% colors = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880],[0.6350, 0.0780, 0.1840]};
% for g = 1:length(groups)
%     scatter(repmat(g,1,folds),accuracies_bin{g},60,'MarkerEdgeColor',[colors{g}])
%     hold on
%     %scatter(g,mean(accuracies_2bin{g}),70,'k','LineWidth',2)
%     scatter(g,mean(accuracies_bin{g}),70,'LineWidth',2,'MarkerEdgeColor',[colors{g}],'MarkerFaceColor',[colors{g}])
% end
% xlim([0 4])
% xticks([1:1:3])
% xticklabels(types)
% ylabel('Accuracy')
% set(gca,'FontSize',16)
% title(sprintf('Decoding Accuracy (%d neurons per group)',length(selected_cells)))
% saveas(gcf,sprintf('NBdecoding_%istim_n%i_%ibins_allgroup.fig',size(confusion_matrix,2),length(selected_cells),bins));
% saveas(gcf,sprintf('NBdecoding_%istim_n%i_%ibins_allgroup.png',size(confusion_matrix,2),length(selected_cells),bins));


