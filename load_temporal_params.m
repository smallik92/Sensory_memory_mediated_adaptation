function [dt, t_init, t_on, trial_dur, num_trials] = load_temporal_params()

% This function contains information that specifies the temporal axis of
% the experiment
% Outputs:
%   dt: time interval chosen for simulation
%   t_init: initial period of stimulus absence at the beginning of a trial
%   t_on: the period for which stimulus is present during a trial
%   trial_dur: the duration of one trial
%   num_trials: number of times the same stimulus repeats

    dt = 0.1;
    t_init = 5;
    t_on = 20;
    trial_dur = 60;
    num_trials = 25;

end