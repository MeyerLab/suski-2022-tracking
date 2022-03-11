function row_column = wellNameToRowColumn(str)
%Converts the IX_Micros naming system to a different naming system
%   IX_micro names the wells "A06", but this function changes that to 1_6.
%   The input is the well name (eg A06) and the output is the new name (eg
%   1_6)
str=char(str);
rowletters='ABCDEFGH';
row=find(rowletters==str(1));
column=str2num(str(2:3));
row_column=[num2str(row),'_',num2str(column)];
end

