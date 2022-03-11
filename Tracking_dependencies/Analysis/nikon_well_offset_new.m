function [map] = nikon_well_offset_new(conditions,num_sites, frac_imaging)

%%% Create list of shots
wellindex = [];
for i = 1:size(conditions,1)
    rowMat = cell2mat(conditions(i,2));
    colMat = cell2mat(conditions(i,3));
    siteMat = cell2mat(conditions(i,4));
    for row = rowMat
        for col = colMat
            for site = siteMat
                try
                    wellindex = [wellindex; [row col site]];
                catch
                    keyboard;
                end
            end
        end
    end
end



totalshots = join([string(wellindex(:,1)) string(wellindex(:,2)) string(wellindex(:,3))],'_');
well_i_row = min(wellindex(:,1));
well_i_col = min(wellindex(wellindex(:,1) == well_i_row,2));

remaining_wells = wellindex;
positive_pass = true;

%%% initialize with first well
order = [];
current_well = [num2str(well_i_row) '_' num2str(well_i_col)];
order = [order; strcat(current_well,"_",string([1:num_sites]'))];
remove_ind = remaining_wells(:,1) == well_i_row & remaining_wells(:,2) == well_i_col;
remaining_wells(remove_ind,:) = [];


%%% loop through remaining wells
while ~isempty(remaining_wells)
    if positive_pass
        next_col = min(remaining_wells(remaining_wells(:,1) == well_i_row,2));
        if ~isempty(next_col)
            well_i_col = next_col;
        else
            positive_pass = false;
            well_i_row = min(remaining_wells(:,1));
            well_i_col= max(remaining_wells(remaining_wells(:,1) == well_i_row,2));
            
        end
    else
        next_col = max(remaining_wells(remaining_wells(:,1) == well_i_row,2));
        if ~isempty(next_col)
            well_i_col = next_col;
        else
            positive_pass = true;          
            well_i_row =  min(remaining_wells(:,1));
            well_i_col= min(remaining_wells(remaining_wells(:,1) == well_i_row,2));
            
        end
    end
    current_well = [num2str(well_i_row) '_' num2str(well_i_col)];
    order = [order; strcat(current_well,"_",string([1:num_sites]'))];
    remove_ind = remaining_wells(:,1) == well_i_row & remaining_wells(:,2) == well_i_col;
    remaining_wells(remove_ind,:) = [];
    
end

%%% Generate dictionary of shots to give offsets (in units of fraction of
%%% frame time). First frame will have highest offset
keys = order;
vec = (0:length(keys)-1) * frac_imaging / length(keys);
vals = num2cell(1 - vec);

map = containers.Map(keys, vals);