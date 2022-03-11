function [ trace1_smooth ] = nansmoothm(trace1,smooth_factor,method)
%Smooths single vector, ignoring NaNs
%Inputs: Trace1=single row vector of data you want smoothed, by a moving average filter;
%smooth_factor=span for moving average, default is 5
if nargin<2
    smooth_factor=5;
end

trace1_smooth=trace1;
nan_index=find(~isnan(trace1));
smooth_trace_temp=smooth(trace1(nan_index),smooth_factor,method);
trace1_smooth(nan_index)=smooth_trace_temp;
end

