function [ corrected_traces, mask_traces, bad_traces ] = maskBlankFrame( traces, shots, frames, thresh)
%maskBlankFrame finds blank frames from traces grouped by shot and
%creates a NaN mask for trace values and outputs the mask applied to the
%input traces

unique_shots = unique(shots);
mask_traces = ones(size(traces));
bad_traces = false(size(traces,1),1);


for i = 1:length(unique_shots)
    shot_inds = find(ismember(shots, unique_shots{i}));
    shot_frames = traces(shot_inds,frames);
    %plot(shot_frames');
    frame_mean = nanmean(shot_frames,1);
    %plot(frame_mean');
    frame_mean_diff = diff(frame_mean);
    rel_diff = [ 0 frame_mean_diff./frame_mean(1:end-1) ];
    curve = [diff(rel_diff) 0];
    %plot(rel_diff');
    %plot(curve');
    blank_frames = curve > thresh;
%     shot_frames_corrected = ones(size(shot_frames));
%     shot_frames_corrected(:,blank_frames) = NaN;
    mask_traces(shot_inds,blank_frames) = NaN;  
    if(any(blank_frames))
        bad_traces(shot_inds) = true;
    end
end

corrected_traces = mask_traces.*traces;


end

