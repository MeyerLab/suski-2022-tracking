function generation = count_generation(current_cell,genealogy, cellIDs)
    curr_ind = find(cellIDs == current_cell);
    mother_cell = genealogy(curr_ind);
    mother_ind = find(cellIDs == mother_cell);
    if ~isempty(mother_ind)
        generation = count_generation(mother_cell,genealogy,cellIDs)+1;
    else
        generation = 1;
    end
end

