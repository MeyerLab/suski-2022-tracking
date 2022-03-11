function [map] = IXM_well_offset(conditions,num_sites, frac_imaging)

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
well_i_col = min(wellindex(:,2));
well_i_row = max(wellindex(wellindex(:,2) == well_i_col,1));

remaining_wells = wellindex;
positive_pass = false;

%%% initialize with first well
order = [];
current_well = [num2str(well_i_row) '_' num2str(well_i_col)];
order = [order; strcat(current_well,"_",string([1:num_sites]'))];
remove_ind = remaining_wells(:,1) == well_i_row & remaining_wells(:,2) == well_i_col;
remaining_wells(remove_ind,:) = [];


%%% loop through remaining wells
while ~isempty(remaining_wells)
    if positive_pass
        next_row = min(remaining_wells(remaining_wells(:,2) == well_i_col,1));
        if ~isempty(next_row)
            well_i_row = next_row;
        else
            positive_pass = false;
            well_i_col = min(remaining_wells(:,2));
            well_i_row= max(remaining_wells(remaining_wells(:,2) == well_i_col,1));
            
        end
    else
        next_row = max(remaining_wells(remaining_wells(:,2) == well_i_col,1));
        if ~isempty(next_row)
            well_i_row = next_row;
        else
            positive_pass = true;          
            well_i_col =  min(remaining_wells(:,2));
            well_i_row= min(remaining_wells(remaining_wells(:,2) == well_i_col,1));
            
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