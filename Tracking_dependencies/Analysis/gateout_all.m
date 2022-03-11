function [ sensor_gated ] = gateout_all( sensor, ind )
    name = fieldnames(sensor);
    for i = 1:length(name)
        if size(sensor.(name{i}),1) == length(ind)
            sensor_gated.(name{i})=sensor.(name{i})(ind,:,:);
        else
            sensor_gated.(name{i})=sensor.(name{i});
        end
    end
    

end

