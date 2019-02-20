
i=1;
path{i} = 'H:\ProcessedDataArchive\Pati\DarkRearing\DarkReared\2268_NC_180523_DR_suite2p_processed\processed_suite2p';i=i+1;
path{i} = 'H:\ProcessedDataArchive\Pati\DarkRearing\DarkReared\2269_NC_180521_DR_suite2p_processed\processed_suite2p';i=i+1;
path{i} = 'H:\ProcessedDataArchive\Pati\DarkRearing\DarkReared\2270_1R_180519_DR_suite2p_processed\processed_suite2p';i=i+1;
path{i} = 'H:\ProcessedDataArchive\Pati\DarkRearing\DarkReared\2271_1R_180520_DR_suite2p_processed\processed_suite2p';i=i+1;
path{i} = 'H:\ProcessedDataArchive\Pati\DarkRearing\DarkReared\2280_1L_180618_DR_suite2p_processed\processed_suite2p';i=i+1;
path{i} = 'H:\ProcessedDataArchive\Pati\DarkRearing\DarkReared\2297_NC_180702_DR_suite2p_processed\processed_suite2p';i=i+1;
path{i} = 'H:\ProcessedDataArchive\Pati\DarkRearing\DarkReared\2298_1R1L_180704_DR_suite2p_processed\processed_suite2p';i=i+1;
session_order = {'2268_NC_180523_DR','2269_NC_180521_DR','2270_1R_180519_DR','2271_1R_180520_DR','2280_1L_180618_DR','2297_NC_180702_DR','2298_1R1L_180704_DR'};
home = pwd;

thresh = 6; %%%%%%%%%%%%%%%this is where set threshold

for s = 1:size(path,2)
    expID = session_order{s};
    cd(path{s})
    load(sprintf('F_%s_suite2p_plane1_proc.mat',expID))
    cd(home)
    
    iscellindex = logical([dat.stat(:).iscell]);%Extract only cells
    cellnum = find(iscellindex); % get the cell number to use
    numcells = length(cellnum);

    allsignals = zeros(numcells,dat.stat(1).blockstarts(end));
    for i=1:numcells
        allsignals(i,dat.stat(cellnum(i)).st) = dat.stat(cellnum(i)).c;
    end
    
    %allsignals = allsignals(:,1:20000);
    cell_index = allsignals>thresh; 
    cell_index2 = sum(cell_index,2);
    cell_index3 = find(cell_index2>0);
    
    session_cells{1,s} = expID;
    session_cells{2,s} = sprintf('out of %i',length(cell_index2));
    session_cells{3,s} = length(cell_index3);
    %min(cell_index2)
end

save(sprintf('total_num_cells_thresh%i',thresh),'session_cells');

