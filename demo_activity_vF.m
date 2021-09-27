%% Description

%This code demonstrates neural and behavioral activity synthesized by the
%proposed model when the system is excited by the same stimulus for 24
%repeats.
% Author: Sruti Mallik
% Date: 09/26/2021

%% Clear command window and workspace and close all figures

clc; clear all; close all;

%% Set parameters along the time axis

[dt, t_init, t_on, trial_dur, num_trials] = load_temporal_params();

t_off = trial_dur - (t_init + t_on);
t_index = 0:dt:num_trials*trial_dur - dt; %temporal axis at which simulation values will be synthesized

%% Set dimension of the model 

latent_dim = 2; %dimension of the latent space (m in the paper)
neural_dim = 20; %dimension of the 'neural' space (n in the paper)
overlap = 5; %overlap between the tuning curves that maps the 'neural' space to latent space

if overlap >= neural_dim
    error('Check the parameter inputs... overlap should be stictly less than the neural dimension');
end

%% Set parameters of the model

load('sample_params_data_fit.mat');

% set parameters for the dynamical decoder in the model

% parameters for the decoder that creates instantaneous latent
% representations
drift_nu = (params_50s(1)+params_4s(1))/2;
drift_nu_gamma = (params_50s(2)+params_4s(2))/2; 
tau_nu = (params_50s(3)+params_4s(3))/2; 
b_matrix = weighting_matrix(latent_dim, neural_dim, overlap);

% parameters for the decoder that creates residual memory
drift_gamma = (params_50s(4)+params_4s(4))/2; 
tau_gamma = (params_50s(5)+params_4s(5))/2;
beta = (params_50s(6)+params_4s(6))/2;

%% Specification of the nominal representation to be tracked (intensity is embedded in the nominal rep.)  

c0 = 0.15;
conc = 0.75;

if conc<c0
    z_c = 0.1 + 0.9*exp(-8*(conc - c0)^2);
else
    z_c = 0.35 + 0.65*exp(-2*(conc - c0)^2);
end

z_target = zeros(latent_dim, trial_dur/dt);
z_target(1, t_init/dt + 1: floor((t_init+t_on)/dt)) = z_c;

z = repmat(z_target, 1, num_trials);

%% Set up the optimal control problem 

aux_eig_A = -1e-5*ones(1, latent_dim); %place very slow poles for numerical stability

latent_penalty = 10; %penalty for inaccurate latent representation
energy_penalty = 2; %penalty on energy expenditure
deriv_penalty = 0.1; %penalty on change in activity

[A_final, B_final, Q_final, R_final] = set_up_optimal_control(latent_dim, neural_dim, tau_nu, ...
                                            drift_nu_gamma, drift_nu, b_matrix, tau_gamma, drift_gamma,...
                                            beta, aux_eig_A, latent_penalty, energy_penalty, deriv_penalty);
                                        
%% Solve Riccati Equation forward in time 

K0 = zeros(size(A_final));% initial point for the differential equation
[T,K] = ode45(@(t,K)mRiccati_F(t,K,A_final,B_final,R_final,Q_final), t_index, K0); 

%% Set parameters for scaling response from a.u. to a positive range
    
baseline = 0.25;
low_asymp = 0; up_asymp = 1;
curve = 7.5;

%% Set parameters for the low pass filter that maps instantaneous neural activity to Calcium dynamics

low_pass = -0.145; 
gain = 1;

%% Set noise level of the activity

noise_level = 0.01;
SEED = 7;
 
%% Simulate neural activity
[Calcium_response, x_response, nu_response, gamma_response] = simulate_neural_activity_vF(dt, t_init, t_on, trial_dur, t_index,...
                                                              K, z, noise_level, SEED, B_final, R_final,...
                                                              tau_nu, drift_nu, drift_nu_gamma, b_matrix,...
                                                              tau_gamma, drift_gamma, beta, ...
                                                              latent_dim, neural_dim, ...
                                                              overlap, baseline, low_asymp, up_asymp, curve, low_pass, gain);
                                                          
%% Load parameters for neural behavioral transformation 

behav_states = 3; % we consider three behavioral states i.e., forward, pause, reverse

[forward_var, pause_var, reverse_var, ...
stim_on_mean, stim_off_mean, transition_mat, valence_mat, relative_weights] = load_behav_params(); %load parameters for simulating behavioral activity

[forward_behav, pause_behav, reverse_behav, ~, ~] = simulate_behavioral_activity_vF(t_index, behav_states,...
                                                                                    nu_response, gamma_response,...
                                                                                    forward_var, pause_var, reverse_var, ...
                                                                                    stim_on_mean, stim_off_mean, transition_mat, valence_mat, relative_weights,...
                                                                                    noise_level, SEED);
%% Downsample behavioral activity to match experimental data

[trial_time_ds, forward_response_trial_ds, pause_response_trial_ds, reverse_response_trial_ds] =...
    downsample_behav_activity(dt, trial_dur, num_trials, forward_behav, pause_behav, reverse_behav);

%% Visualize Calcium activity across trials

figure, 
plot(t_index/60, Calcium_response, 'LineWidth', 2);
xlabel('Time in mins')
ylabel('Calcium Activity')
xlim([0, 25])
set(gca, 'FontName', 'Arial', 'FontSize', 14);

%% Peak activity averaged across trials

calcium_response_trial = reshape(Calcium_response, [trial_dur/dt, num_trials]);
calcium_response_trial = calcium_response_trial'; 

peaks = max(calcium_response_trial,[],2);
relative_peaks = peaks./repmat(max(peaks), 25,1);

figure, plot(1:num_trials, relative_peaks(1:end),'-s', 'LineWidth', 2 );
xlim([0,num_trials]); ylim([0, 1.01]);
set(gca, 'FontName', 'Arial', 'FontSize', 14);
xlabel('Trial #');
ylabel('Relative peak activity');

%% Behavioral activity in trials

group_size = 3; % number of trials that are averaged

for k = 1:8
    subplot(2,8,k)
    area(trial_time_ds/60, mean(forward_response_trial_ds((k-1)*group_size+2:k*group_size+1, :), 1), 'FaceColor', [0.6, 0.6, 0.6]);
    ylim([0, 1])
    set(gca, 'FontSize',14,'FontName', 'Arial'); 
    ylabel('Fwd. prob.')    
    
    subplot(2,8,k+8)
    area(trial_time_ds/60, mean(reverse_response_trial_ds((k-1)*group_size+2:k*group_size+1, :), 1), 'FaceColor', [0.84,0.46,0.46]);
    ylim([0, 1])
    set(gca, 'FontSize',14,'FontName', 'Arial'); 
    ylabel('Rev. prob.')    
end




