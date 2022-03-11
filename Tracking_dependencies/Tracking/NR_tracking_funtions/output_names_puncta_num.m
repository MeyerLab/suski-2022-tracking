function [ names ] = output_names_puncta_num( sigs,localbg,ring,punctacalc, punctaThresh, varThresh, punctaCountThresh)
names={1,2,3,4;'x','y','nuclear area','nuclear mass'};
length=4;

for i=1:size(sigs,2)
    length=length+1;
    names(1:2,length)={length;[sigs{i} 'mean']};
end

for i=1:size(sigs,2)
    length=length+1;
    names(1:2,length)={length;[sigs{i} 'median']};
end

for i=1:size(sigs,2)
    if localbg(i)
        length=length+1;
        names(1:2,length)={length;[sigs{i} 'block bg']};
    end
end


for i=1:size(sigs,2)
    if ring(i)
        length=length+1;
        names(1:2,length)={length;[sigs{i} 'cyto ring']};
    end
end

if punctacalc
    length = length+1;
    names(1:2,length) = {length; ['PCNA mean']};
    
    length = length+1;
    names(1:2,length) = {length; ['Filter mean']};
    
    length = length+1;
    names(1:2,length) = {length; ['Variance mean']};
    
    length = length+1;
    names(1:2,length) = {length; ['Variance std']};
    
    for i=1:size(punctaThresh,2)
        length = length+1;
        names(1:2,length) = {length; ['Filter Masked mean ' num2str(punctaThresh(i))]};
    end
    
    for i=1:size(punctaThresh,2)
        length = length+1;
        names(1:2,length) = {length; ['Filter Masked area ' num2str(punctaThresh(i))]};
    end
    
    for i=1:size(varThresh,2)
        length = length+1;
        names(1:2,length) = {length; ['Variance Masked mean ' num2str(punctaThresh(i))]};
    end
    
    for i=1:size(varThresh,2)
        length = length+1;
        names(1:2,length) = {length; ['Variance Masked area ' num2str(punctaThresh(i))]};
    end
    
    for i=1:size(punctaCountThresh,2)
        length = length+1;
        names(1:2,length) = {length; ['Puncta Num ' num2str(punctaCountThresh(i))]};
    end
    
end

