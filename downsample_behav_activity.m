function [trial_time_ds, forward_response_trial_ds, pause_response_trial_ds, reverse_response_trial_ds]...
                = downsample_behav_activity(dt, trial_dur, num_trials, forward_behav, pause_behav, reverse_behav)

    trial_time = 0:dt:trial_dur-dt;
    trial_time_ds = downsample(trial_time, 5); 

    forward_response_trial = reshape(forward_behav,[trial_dur/dt,num_trials]);
    forward_response_trial_ds = smoothdata(downsample(forward_response_trial, 5),'movmedian', 5);
    forward_response_trial_ds = forward_response_trial_ds';

    pause_response_trial = reshape(pause_behav, [trial_dur/dt, num_trials]);
    pause_response_trial_ds = smoothdata(downsample(pause_response_trial, 5),'movmedian', 5);
    pause_response_trial_ds = pause_response_trial_ds';

    reverse_response_trial = reshape(reverse_behav,[trial_dur/dt,num_trials]);
    reverse_response_trial_ds = smoothdata(downsample(reverse_response_trial, 5),'movmedian', 5);
    reverse_response_trial_ds = reverse_response_trial_ds'; 

end

