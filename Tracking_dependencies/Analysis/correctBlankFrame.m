function [ corrected_traces, bad_traces ] = correctBlankFrame( traces, shots, frames, thresh)
%correctBlankFrame removes blank frames from traces grouped by shot and
%replaces it with NaN

unique_shots = unique(shots);
corrected_traces = traces;
bad_traces = false(size(traces,1),1);

for i = 1:length(unique_shots)
    shot_inds = find(ismember(shots, unique_shots{i}));
    shot_frames = traces(shot_inds,frames);
    %plot(shot_frames');
    frame_mean = nanmean(shot_frames,1);
    frame_mean_diff = [0 diff(frame_mean)];
    blank_frames = frame_mean_diff < thresh;
    shot_frames_corrected = shot_frames;
    shot_frames_corrected(:,blank_frames) = NaN;
    corrected_traces(shot_inds,frames) = shot_frames_corrected;  
    if(any(blank_frames))
        bad_traces(shot_inds) = true;
    end
end




end

