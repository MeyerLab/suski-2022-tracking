function wellName = rowColumnToWellName(row, col)
% Converts wellname to row_col format
rowletters='ABCDEFGH';
wellName = sprintf('%c%02d',rowletters(row),col);

end

