function [mode_of_trace] = nanmode(vector_of_data)
new_vector=vector_of_data(~isnan(vector_of_data));
new_vector_rounded=round(new_vector,1);
mode_of_trace=mode(new_vector_rounded);
end

