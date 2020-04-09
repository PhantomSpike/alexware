function [start_time_ms] = get_triggers(synch_ch,min_trig_length_s,min_inter_trig_length_s,fs)
%function [start_time_ms] = get_triggers(synch_ch,min_trig_length_s,fs)
%>> INPUT >>
%synch_ch = The extracted synch channel that contains the triggers
%min_trig_length_s = The minmum length of one trigger in seconds. This is used to exclude accicental mini triggers
%fs = The sampling rate in seconds
%<< OUTPUT <<
%start_time_ms = The starting time of every trigger in ms relative to the
%beginning of the recording 
trig_length_samples = min_trig_length_s*fs;
min_inter_trig_length_samples = min_inter_trig_length_s*fs;

diff_sig = diff(synch_ch); %Find the difference between every n+1 sample - n sample. This tells us the the beginning/end of each sweep

%Find the sample no for the beginning of each sweep
start_ix = find(diff_sig>0); %Find where the synch channel increased i.e. the beginning of a trigger
end_ix = find(diff_sig<0); %Find where the synch channel decreased i.e. the end if a trigger

if numel(end_ix) > numel(start_ix) %Remove a faulty end_ix if the triggers started from a high value
    end_ix(1) = [];
end

if numel(start_ix) > numel(end_ix) %Remove faulty start_ix if there more than end_ix
    start_ix = start_ix(1:numel(end_ix));
end

start_ix = start_ix(:); %Force a column vector
end_ix = end_ix(:); %Force a column vector
% Fuse triggers that have a inter-trigger interval too short
% Quentin
ITI_sample = end_ix(1:end-1) - start_ix(2:end);
rm_idx = (ITI_sample > -min_inter_trig_length_samples);
tmp = start_ix(1);
start_ix = start_ix(2:end);
start_ix = [tmp; start_ix(~rm_idx)];
tmp = end_ix(end);
end_ix = end_ix(1:end-1);
end_ix = [end_ix(~rm_idx); tmp];
% / end Fusion

% Check if trigger is too long
diff_ix = end_ix - start_ix; %Find the length of each sweep in samples
normal_ix = diff_ix >= trig_length_samples;

if ismember(0,normal_ix)
    ix_problem = find(0 == normal_ix);
    for problem_no = 1:numel(ix_problem)
        warning('There were presentations which were interupted!');
        fprintf('Check out presentation number %.0f\n',ix_problem(problem_no));
    end
end
% / end check.

start_ix = start_ix(diff_ix >= trig_length_samples); %Keep only the triggers which have length >= minimum triger length
start_time_ms = (start_ix/fs).*1000; %Convert the starting sample numbers to times in ms