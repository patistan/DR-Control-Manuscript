%% load different dataouts FOR NATURAL SCENES
clear all

session_order = {'2268_NC_180523_DR','2269_NC_180521_DR','2270_1R_180519_DR','2271_1R_180520_DR','2280_1L_180618_DR','2297_NC_180702_DR','2298_1R1L_180704_DR'};
i=1;
NSdataOuts{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/NatScenes/%s_NatScenes_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
NSdataOuts{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/NatScenes/%s_NatScenes_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
NSdataOuts{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/NatScenes/%s_NatScenes_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
NSdataOuts{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/NatScenes/%s_NatScenes_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
NSdataOuts{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/NatScenes/%s_NatScenes_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
NSdataOuts{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/NatScenes/%s_NatScenes_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;
NSdataOuts{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/NatScenes/%s_NatScenes_dataOut.mat',session_order{i},session_order{i},session_order{i}));i=i+1;

%get dimensions of all response matrices
dimensions = zeros(length(NSdataOuts),3);
for s=1:length(NSdataOuts)
    dimensions(s,1) = size(NSdataOuts{s}.dataOut.responseMatrix_1,1); %reps
    dimensions(s,2) = size(NSdataOuts{s}.dataOut.responseMatrix_1,2); %cells
    dimensions(s,3) = size(NSdataOuts{s}.dataOut.responseMatrix_1,3); %stim
end

%create a large dataOut 
dataOut.responseMatrix_1 = nan(max(dimensions(:,1)),sum(dimensions(:,2)),dimensions(1,3));
dataOut.isRemovedBlock = true(max(dimensions(:,1)),sum(dimensions(:,2)),dimensions(1,3)); %anything that isn't included later on will be a bad block by making 1 the default
dataOut.hasLocomotion = true(max(dimensions(:,1)),sum(dimensions(:,2)),dimensions(1,3)); %anything that isn't included later on will be a bad block by making 1 the default
dataOut.totalNumStimuli = dimensions(1,3);
dataOut.totalNumCells = sum(dimensions(:,2));
dataOut.stimulus_ID = 'NatScenes';
dataOut.stats.global.response_average_pval_fdr = nan(sum(dimensions(:,2)),dimensions(1,3));
dataOut.stats.global.response_average_MED = nan(sum(dimensions(:,2)),dimensions(1,3));
dataOut.stats.global.responsive_cells_p001_fdr_average_index = false(sum(dimensions(:,2)),1);
dataOut.stats.global.response_ACTUAL_avg_vals = nan(sum(dimensions(:,2)),dimensions(1,3));
dataOut.stats.global.response_ACTUAL_avg_sem = nan(sum(dimensions(:,2)),dimensions(1,3));

%fill in large dataOut
count = 0;
count2 = 0;
FIRST_CELL = 1;
for s=1:length(NSdataOuts)
    %fill in response matrix
    current = NSdataOuts{s}.dataOut.responseMatrix_1;
    current_dimensions = [size(current,1) size(current,2) size(current,3)];
    dataOut.responseMatrix_1(1:current_dimensions(1),FIRST_CELL:FIRST_CELL+current_dimensions(2)-1,1:current_dimensions(3)) = current;
    
    %fill in removed block
    current = NSdataOuts{s}.dataOut.isRemovedBlock;
    current_dimensions = [size(current,1) size(current,2) size(current,3)];
    dataOut.isRemovedBlock(1:current_dimensions(1),FIRST_CELL:FIRST_CELL+current_dimensions(2)-1,1:current_dimensions(3)) = current;
    
    %fill in locomotion
    current = NSdataOuts{s}.dataOut.hasLocomotion;
    current_dimensions = [size(current,1) size(current,2) size(current,3)];
    dataOut.hasLocomotion(1:current_dimensions(1),FIRST_CELL:FIRST_CELL+current_dimensions(2)-1,1:current_dimensions(3)) = current;

    %fill in pvalue for median responses (using blcok average)
    current = NSdataOuts{s}.dataOut.stats.global.response_average_pval_fdr;
    current_dimensions = [size(current,1) size(current,2)];
    dataOut.stats.global.response_average_pval_fdr(FIRST_CELL:FIRST_CELL+current_dimensions(1)-1,1:current_dimensions(2))= current;
    
    %fill in value for median responses (using blcok average)
    current = NSdataOuts{s}.dataOut.stats.global.response_average_vals;
    current_dimensions = [size(current,1) size(current,2)];
    dataOut.stats.global.response_average_MED(FIRST_CELL:FIRST_CELL+current_dimensions(1)-1,1:current_dimensions(2))= current;
    
    %fill in index of if cell is significant (pval<0.01) or not
    current = NSdataOuts{s}.dataOut.stats.global.response_average_pval_fdr;
    for c = 1:NSdataOuts{s}.dataOut.totalNumCells
        count = count+1;
        c_resps_pval = current(c,:);
        if min(c_resps_pval) <0.01
            dataOut.stats.global.responsive_cells_p001_fdr_average_index(count) = 1;
        end            
    end
    
    %fill in average value across trials after removing bad blocks
    current = NSdataOuts{s}.dataOut.responseMatrix_1;
    current2 = NSdataOuts{s}.dataOut.isRemovedBlock;    
    for c = 1:NSdataOuts{s}.dataOut.totalNumCells
        count2=count2+1;
        c_resps = squeeze(current(:,c,:));
        c_bad = squeeze(current2(:,c,:));
        c_resps(c_bad==1) = nan;
        c_resps_avg = nanmean(c_resps,1);
        c_resps_sem = nanstd(c_resps,1)./sqrt(sum(~isnan(c_resps),1));
        dataOut.stats.global.response_ACTUAL_avg_vals(count2,:) = c_resps_avg;
        dataOut.stats.global.response_ACTUAL_avg_sem(count2,:) = c_resps_sem; 
    end 
    
    FIRST_CELL = FIRST_CELL + dimensions(s,2);
end
allcells = [1:sum(dimensions(:,2))]';
dataOut.stats.global.responsive_cells_p001_fdr_average = allcells(dataOut.stats.global.responsive_cells_p001_fdr_average_index);

original_cellID = [];
for s = 1:length(NSdataOuts)
    original_cellID = [original_cellID 1:dimensions(s,2)];
end

save('dataOut_NatScenes_POOLED.mat','dataOut','dimensions','original_cellID','session_order');


%% load different dataouts FOR GRATINGS
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

%get dimensions of all response matrices
dimensions = zeros(length(GdataOut),3);
for s=1:length(GdataOut)
    dimensions(s,1) = size(GdataOut{s}.dataOut.responseMatrix_1,1); %reps
    dimensions(s,2) = size(GdataOut{s}.dataOut.responseMatrix_1,2); %cells
    dimensions(s,3) = size(GdataOut{s}.dataOut.responseMatrix_1,3); %stim
end

%create a large dataOut 
dataOut.responseMatrix_1 = nan(max(dimensions(:,1)),sum(dimensions(:,2)),dimensions(1,3));
dataOut.isRemovedBlock = true(max(dimensions(:,1)),sum(dimensions(:,2)),dimensions(1,3)); %anything that isn't included later on will be a bad block by making 1 the default
dataOut.hasLocomotion = true(max(dimensions(:,1)),sum(dimensions(:,2)),dimensions(1,3)); %anything that isn't included later on will be a bad block by making 1 the default
dataOut.totalNumStimuli = dimensions(1,3);
dataOut.totalNumCells = sum(dimensions(:,2));
dataOut.stimulus_ID = 'Gratings';
dataOut.stats.global.response_average_pval_fdr = nan(sum(dimensions(:,2)),dimensions(1,3));
dataOut.stats.global.response_average_MED = nan(sum(dimensions(:,2)),dimensions(1,3));
dataOut.stats.global.responsive_cells_p001_fdr_average_index = false(sum(dimensions(:,2)),1);
dataOut.stats.global.response_ACTUAL_avg_vals = nan(sum(dimensions(:,2)),dimensions(1,3));
dataOut.stats.global.response_ACTUAL_avg_sem = nan(sum(dimensions(:,2)),dimensions(1,3));

%fill in large dataOut
count = 0;
count2 = 0;
FIRST_CELL = 1;
for s=1:length(GdataOut)
    %fill in response matrix
    current = GdataOut{s}.dataOut.responseMatrix_1;
    current_dimensions = [size(current,1) size(current,2) size(current,3)];
    dataOut.responseMatrix_1(1:current_dimensions(1),FIRST_CELL:FIRST_CELL+current_dimensions(2)-1,1:current_dimensions(3)) = current;
    
    %fill in removed block
    current = GdataOut{s}.dataOut.isRemovedBlock;
    current_dimensions = [size(current,1) size(current,2) size(current,3)];
    dataOut.isRemovedBlock(1:current_dimensions(1),FIRST_CELL:FIRST_CELL+current_dimensions(2)-1,1:current_dimensions(3)) = current;
    
    %fill in locomotion
    current = GdataOut{s}.dataOut.hasLocomotion;
    current_dimensions = [size(current,1) size(current,2) size(current,3)];
    dataOut.hasLocomotion(1:current_dimensions(1),FIRST_CELL:FIRST_CELL+current_dimensions(2)-1,1:current_dimensions(3)) = current;

    %fill in pvalue for median responses (using blcok average)
    current = GdataOut{s}.dataOut.stats.global.response_average_pval_fdr;
    current_dimensions = [size(current,1) size(current,2)];
    dataOut.stats.global.response_average_pval_fdr(FIRST_CELL:FIRST_CELL+current_dimensions(1)-1,1:current_dimensions(2))= current;
    
     %fill in value for median responses (using blcok average)
    current = GdataOut{s}.dataOut.stats.global.response_average_vals;
    current_dimensions = [size(current,1) size(current,2)];
    dataOut.stats.global.response_average_MED(FIRST_CELL:FIRST_CELL+current_dimensions(1)-1,1:current_dimensions(2))= current;

   
    %fill in index of if cell is significant (pval<0.01) or not
    current = GdataOut{s}.dataOut.stats.global.response_average_pval_fdr;
    for c = 1:GdataOut{s}.dataOut.totalNumCells
        count = count+1;
        c_resps_pval = current(c,:);
        if min(c_resps_pval) <0.01
            dataOut.stats.global.responsive_cells_p001_fdr_average_index(count) = 1;
        end            
    end
    
    %fill in average value across trials after removing bad blocks
    current = GdataOut{s}.dataOut.responseMatrix_1;
    current2 = GdataOut{s}.dataOut.isRemovedBlock;
    for c = 1:GdataOut{s}.dataOut.totalNumCells
        count2=count2+1;
        c_resps = squeeze(current(:,c,:));
        c_bad = squeeze(current2(:,c,:));
        c_resps(c_bad==1) = nan;
        c_resps_avg = nanmean(c_resps,1);
        c_resps_sem = nanstd(c_resps,1)./sqrt(sum(~isnan(c_resps),1));
        dataOut.stats.global.response_ACTUAL_avg_vals(count2,:) = c_resps_avg;
        dataOut.stats.global.response_ACTUAL_avg_sem(count2,:) = c_resps_sem;     
    end 
    
    FIRST_CELL = FIRST_CELL + dimensions(s,2);
end
allcells = [1:sum(dimensions(:,2))]';
dataOut.stats.global.responsive_cells_p001_fdr_average = allcells(dataOut.stats.global.responsive_cells_p001_fdr_average_index);

original_cellID = [];
for s = 1:length(GdataOut)
    original_cellID = [original_cellID 1:dimensions(s,2)];
end

save('dataOut_Gratings_POOLED.mat','dataOut','dimensions','original_cellID','session_order');

%% Pool gauss fit and BW for all cells
clear all

session_order = {'2268_NC_180523_DR','2269_NC_180521_DR','2270_1R_180519_DR','2271_1R_180520_DR','2280_1L_180618_DR','2297_NC_180702_DR','2298_1R1L_180704_DR'};

i=1;
gaussFits{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/Gauss_fit/gaussFit_results.mat',session_order{i},session_order{i}));i=i+1;
gaussFits{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/Gauss_fit/gaussFit_results.mat',session_order{i},session_order{i}));i=i+1;
gaussFits{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/Gauss_fit/gaussFit_results.mat',session_order{i},session_order{i}));i=i+1;
gaussFits{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/Gauss_fit/gaussFit_results.mat',session_order{i},session_order{i}));i=i+1;
gaussFits{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/Gauss_fit/gaussFit_results.mat',session_order{i},session_order{i}));i=i+1;
gaussFits{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/Gauss_fit/gaussFit_results.mat',session_order{i},session_order{i}));i=i+1;
gaussFits{i} = load(sprintf('H:/ProcessedDataArchive/Pati/DarkRearing/DarkReared/%s_suite2p_processed/processed_suite2p/%s_analysis/Gratings/Gauss_fit/gaussFit_results.mat',session_order{i},session_order{i}));i=i+1;

%piece all variables together 
%all_1minCV = [];
fit_oriBW = [];
fit_prefOri = [];
fit_rsquare = [];
pref_ori_estimate = [];
pref_sf_estimate = [];
gaussFitResults ={};
gaussgof ={};


for s = 1:length(gaussFits)
    %all_1minCV = [all_1minCV gaussFits{s}.all_1minCV']; %rotate to be like other variables
    fit_oriBW = [fit_oriBW gaussFits{s}.fit_oriBW];
    fit_prefOri = [fit_prefOri gaussFits{s}.fit_prefOri];
    fit_rsquare = [fit_rsquare gaussFits{s}.fit_rsquare];
    pref_ori_estimate = [pref_ori_estimate gaussFits{s}.pref_ori_estimate];
    pref_sf_estimate = [pref_sf_estimate gaussFits{s}.pref_sf_estimate];
    gaussFitResults = [gaussFitResults gaussFits{s}.gaussFitResults];
    gaussgof = [gaussgof gaussFits{s}.gaussgof];
end

clear i s gaussFits
save('gaussFit_results_POOLED.mat');

%% Pool RF data
clear all

i=1;
RFdataS{i} = load('H:/ProcessedDataArchive/Pati/NatScene2/2209_2L_180302_suite2p_processed/processed_suite2p/2209_2L_180302_analysis/RF/RFM_data.mat');i=i+1;
RFdataS{i} = load('H:/ProcessedDataArchive/Pati/NatScene2/2210_NC_180228_suite2p_processed/processed_suite2p/2210_NC_180228_analysis/RF/RFM_data.mat');i=i+1;
RFdataS{i} = load('H:/ProcessedDataArchive/Pati/NatScene2/2212_NC_180227_suite2p_processed/processed_suite2p/2212_NC_180227_analysis/RF/RFM_data.mat');i=i+1;
RFdataS{i} = load('H:/ProcessedDataArchive/Pati/NatScene2/2231_1R1L_180304_suite2p_processed/processed_suite2p/2231_1R1L_180304_analysis/RF/RFM_data.mat');i=i+1;
RFdataS{i} = load('H:/ProcessedDataArchive/Pati/NatScene2/2189_1R_180317_suite2p_processed/processed_suite2p/2189_1R_180317_analysis/RF/RFM_data.mat');i=i+1;
RFdataS{i} = load('H:/ProcessedDataArchive/Pati/NatScene2/2228_1R_180321_suite2p_processed/processed_suite2p/2228_1R_180321_analysis/RF/RFM_data.mat');i=i+1;

i=1;
RFdataOut{i} = load('H:/ProcessedDataArchive/Pati/NatScene2/2209_2L_180302_suite2p_processed/processed_suite2p/2209_2L_180302_analysis/RF/2209_2L_180302_checkersRFM_dataOut.mat');i=i+1;
RFdataOut{i} = load('H:/ProcessedDataArchive/Pati/NatScene2/2210_NC_180228_suite2p_processed/processed_suite2p/2210_NC_180228_analysis/RF/2210_NC_180228_checkersRFM_dataOut.mat');i=i+1;
RFdataOut{i} = load('H:/ProcessedDataArchive/Pati/NatScene2/2212_NC_180227_suite2p_processed/processed_suite2p/2212_NC_180227_analysis/RF/2212_NC_180227_checkersRFM_dataOut.mat');i=i+1;
RFdataOut{i} = load('H:/ProcessedDataArchive/Pati/NatScene2/2231_1R1L_180304_suite2p_processed/processed_suite2p/2231_1R1L_180304_analysis/RF/2231_1R1L_180304_checkersRFM_dataOut.mat');i=i+1;
RFdataOut{i} = load('H:/ProcessedDataArchive/Pati/NatScene2/2189_1R_180317_suite2p_processed/processed_suite2p/2189_1R_180317_analysis/RF/2189_1R_180317_checkersRFM_dataOut.mat');i=i+1;
RFdataOut{i} = load('H:/ProcessedDataArchive/Pati/NatScene2/2228_1R_180321_suite2p_processed/processed_suite2p/2228_1R_180321_analysis/RF/2228_1R_180321_checkersRFM_dataOut.mat');i=i+1;

session_order = {'2209_2L_180302','2210_NC_180228','2212_NC_180227','2231_1R1L_180304','2189_1R_180317','2228_1R_180321'};

%get dimensions of all response matrices
dimensions = zeros(length(RFdataOut),3);
for s=1:length(RFdataOut)
    dimensions(s,1) = size(RFdataOut{s}.dataOut.responseMatrix_1,1); %reps
    dimensions(s,2) = size(RFdataOut{s}.dataOut.responseMatrix_1,2); %cells
    dimensions(s,3) = size(RFdataOut{s}.dataOut.responseMatrix_1,3); %stim
end

%create a large dataOut 
dataOut.responseMatrix_1 = nan(max(dimensions(:,1)),sum(dimensions(:,2)),dimensions(1,3));
dataOut.isRemovedBlock = true(max(dimensions(:,1)),sum(dimensions(:,2)),dimensions(1,3)); %anything that isn't included later on will be a bad block by making 1 the default
dataOut.hasLocomotion = true(max(dimensions(:,1)),sum(dimensions(:,2)),dimensions(1,3)); %anything that isn't included later on will be a bad block by making 1 the default
dataOut.totalNumStimuli = dimensions(1,3);
dataOut.totalNumCells = sum(dimensions(:,2));
dataOut.stimulus_ID = 'RFM';
dataOut.stats.global.response_average_pval_fdr = nan(sum(dimensions(:,2)),dimensions(1,3));
dataOut.stats.global.response_average_MED = nan(sum(dimensions(:,2)),dimensions(1,3));
dataOut.stats.global.responsive_cells_p001_fdr_average_index = false(sum(dimensions(:,2)),1);
dataOut.stats.global.response_ACTUAL_avg_vals = nan(sum(dimensions(:,2)),dimensions(1,3));
dataOut.stats.global.response_ACTUAL_avg_sem = nan(sum(dimensions(:,2)),dimensions(1,3));

RFdata.RFM_cells_index = false(sum(dimensions(:,2)),1);
RFdata.RFM_cells = [];
RFdata.radius_all = nan(sum(dimensions(:,2)),1);
RFdata.CM_coordinates = nan(sum(dimensions(:,2)),2);

%fill in large dataOut
count = 0;
count2 = 0;
FIRST_CELL = 1;
for s=1:length(RFdataOut)
    %fill in response matrix
    current = RFdataOut{s}.dataOut.responseMatrix_1;
    current_dimensions = [size(current,1) size(current,2) size(current,3)];
    dataOut.responseMatrix_1(1:current_dimensions(1),FIRST_CELL:FIRST_CELL+current_dimensions(2)-1,1:current_dimensions(3)) = current;
    
    %fill in removed block
    current = RFdataOut{s}.dataOut.isRemovedBlock;
    current_dimensions = [size(current,1) size(current,2) size(current,3)];
    dataOut.isRemovedBlock(1:current_dimensions(1),FIRST_CELL:FIRST_CELL+current_dimensions(2)-1,1:current_dimensions(3)) = current;
    
    %fill in locomotion
    current = RFdataOut{s}.dataOut.hasLocomotion;
    current_dimensions = [size(current,1) size(current,2) size(current,3)];
    dataOut.hasLocomotion(1:current_dimensions(1),FIRST_CELL:FIRST_CELL+current_dimensions(2)-1,1:current_dimensions(3)) = current;

    %fill in pvalue for median responses (using blcok average)
    current = RFdataOut{s}.dataOut.stats.global.response_average_pval_fdr;
    current_dimensions = [size(current,1) size(current,2)];
    dataOut.stats.global.response_average_pval_fdr(FIRST_CELL:FIRST_CELL+current_dimensions(1)-1,1:current_dimensions(2))= current;
    
     %fill in value for median responses (using blcok average)
    current = RFdataOut{s}.dataOut.stats.global.response_average_vals;
    current_dimensions = [size(current,1) size(current,2)];
    dataOut.stats.global.response_average_MED(FIRST_CELL:FIRST_CELL+current_dimensions(1)-1,1:current_dimensions(2))= current;

   
    %fill in index of if cell is significant (pval<0.01) or not
    current = RFdataOut{s}.dataOut.stats.global.response_average_pval_fdr;
    for c = 1:RFdataOut{s}.dataOut.totalNumCells
        count = count+1;
        c_resps_pval = current(c,:);
        if min(c_resps_pval) <0.01
            dataOut.stats.global.responsive_cells_p001_fdr_average_index(count) = 1;
        end            
    end
    
    %fill in average value across trials after removing bad blocks
    current = RFdataOut{s}.dataOut.responseMatrix_1;
    current2 = RFdataOut{s}.dataOut.isRemovedBlock;
    for c = 1:RFdataOut{s}.dataOut.totalNumCells
        count2=count2+1;
        c_resps = squeeze(current(:,c,:));
        c_bad = squeeze(current2(:,c,:));
        c_resps(c_bad==1) = nan;
        c_resps_avg = nanmean(c_resps,1);
        c_resps_sem = nanstd(c_resps,1)./sqrt(sum(~isnan(c_resps),1));
        dataOut.stats.global.response_ACTUAL_avg_vals(count2,:) = c_resps_avg;
        dataOut.stats.global.response_ACTUAL_avg_sem(count2,:) = c_resps_sem;     
    end 
    
       
    %fill in index for passing RFM criteria and store cells
    current_cells = RFdataS{s}.RFM_cells+(FIRST_CELL-1);    
    RFdata.RFM_cells_index(current_cells) = 1;
    RFdata.RFM_cells = [RFdata.RFM_cells current_cells];
    %fill in radius info
    current = RFdataS{s}.radius_all;
    RFdata.radius_all(current_cells) = current;
    %fill in coordinate info
    current = RFdataS{s}.CM_coordinates;
    RFdata.CM_coordinates(current_cells,:) = current;
    
    
    
    
    FIRST_CELL = FIRST_CELL + dimensions(s,2);
end
allcells = [1:sum(dimensions(:,2))]';
dataOut.stats.global.responsive_cells_p001_fdr_average = allcells(dataOut.stats.global.responsive_cells_p001_fdr_average_index);

original_cellID = [];
for s = 1:length(RFdataOut)
    original_cellID = [original_cellID 1:dimensions(s,2)];
end

save('dataOut_RFM_POOLED.mat','dataOut','RFdata','dimensions','original_cellID','session_order');










