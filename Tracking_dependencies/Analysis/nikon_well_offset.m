function [map] = nikon_well_offset(wellindex,num_sites)


totalshots = join([string(wellindex(:,1)) string(wellindex(:,2)) string(wellindex(:,3))],'_');
well_i_row = min(wellindex(:,1));
well_i_col = min(wellindex(wellindex(:,1) == well_i_row,2));

positive_pass = true;

%%% initialize with first well
order = [];


while well_i_row <= 8
    
    current_well = [num2str(well_i_row) '_' num2str(well_i_col)];
    if any(contains(totalshots, [current_well '_1']))
        order = [order; strcat(current_well,"_",string([1:num_sites]'))];
    end   
    if positive_pass       
        well_i_col = well_i_col + 1;
        if well_i_col > 12
            well_i_col = well_i_col - 1;
            positive_pass = false;
            well_i_row = well_i_row + 1;
        end
    else
        well_i_col = well_i_col -1;
        if well_i_col < 1
            well_i_col = well_i_col + 1;
            positive_pass = true;
                well_i_row = well_i_row + 1;
        end
    end
    
end

keys = order;
vals = num2cell((length(keys):-1:1)/length(keys));

map = containers.Map(keys, vals);