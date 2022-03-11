function [ names ] = output_names_multiIF_puncta( sigs,localbg,ring, puncta, thresh )
names={1,2,3,4,5,6,7,8 ;'x','y','nuclear area','jitterX','jitterY','row','col','fixed site'};
ind=8;

for i=1:size(sigs,2)
    for j=1:size(sigs{i},2)
        ind=ind+1;
        names(1:2,ind)={ind;[num2str(i) '_' sigs{i}{j} '_mean']};
    end
end
for i=1:size(sigs,2)
    for j=1:size(sigs{i},2)
        ind=ind+1;
        names(1:2,ind)={ind;[num2str(i) '_' sigs{i}{j} '_median']};
    end
end
for i=1:size(sigs,2)
    for j=1:size(sigs{i},2)
        if localbg{i}(j)
            ind=ind+1;
            names(1:2,ind)={ind;[num2str(i) '_' sigs{i}{j} '_block bg']};
        end
    end
end
for i=1:size(sigs,2)
    for j=1:size(sigs{i},2)
        if ring{i}(j)
            ind=ind+1;
            names(1:2,ind)={ind;[num2str(i) '_' sigs{i}{j} '_cyto ring']};
        end
    end
end
for i=1:size(sigs,2)
    for j=1:size(sigs{i},2)
        if puncta{i}(j)
            for t=1:length(thresh{i}{j})
                ind=ind+1;
                names(1:2,ind)={ind;[num2str(i) '_' sigs{i}{j} '_puncta area_' num2str(thresh{i}{j}(t))]};
            end
        end
    end
end
for i=1:size(sigs,2)
    for j=1:size(sigs{i},2)
        if puncta{i}(j)
            for t=1:length(thresh{i}{j})
                ind=ind+1;
                names(1:2,ind)={ind;[num2str(i) '_' sigs{i}{j} '_puncta intensity_' num2str(thresh{i}{j}(t))]};
            end
        end
    end
end
for i=1:size(sigs,2)
    for j=1:size(sigs{i},2)
        if puncta{i}(j)
            for t=1:length(thresh{i}{j})
                ind=ind+1;
                names(1:2,ind)={ind;[num2str(i) '_' sigs{i}{j} '_puncta number_' num2str(thresh{i}{j}(t))]};
            end
        end
    end
end



