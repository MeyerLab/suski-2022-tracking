function [ found ] = findFirstInMat( mat )
%findInMat Finds first nonzero element in each row of matrix

found = NaN*ones(size(mat,1),1);
for i = 1:size(mat, 1)
    ind = find(mat(i,:),1,'first');
    if ~isempty(ind)
        found(i) = ind;
    end
end



end

