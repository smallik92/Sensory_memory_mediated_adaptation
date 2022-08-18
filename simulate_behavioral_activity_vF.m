function [forward_behav, pause_behav, reverse_behav, prob_sensory_on, prob_sensory_off] = simulate_behavioral_activity_vF(t_index, behav_states,...
                                                                                            nu_response, gamma_response,...
                                                                                            forward_var, pause_var, reverse_var, ...
                                                                                            stim_on_mean, stim_off_mean, transition_mat, valence_mat, relative_weights,...
                                                                                            noise_level, SEED)
    
    % This function simulates the behavioral activity from latent space inputs
    
    %% Set the seed
    rng(SEED);
                                                                                                                                                                               
    %% Retrieve dimension of latent space 
    latent_dim = size(nu_response, 1);
    
    %% Initialize vectors
    
    cum_sensory_input = zeros(latent_dim, length(t_index));
    
    prob_sensory_off = zeros(behav_states, length(t_index));
    prob_sensory_on = zeros(behav_states, length(t_index));
    
    likelihood = zeros(behav_states, length(t_index));
    likelihood(:,1) = 1/behav_states;
    
    for i = 2:length(t_index) 
        cum_sensory_input(:,i) = relative_weights(1)*nu_response(:,i) + relative_weights*gamma_response(:,i)+mvnrnd([0,0],noise_level*eye(latent_dim))';
        
        prob_sensory_off(1, i) = mvnpdf(cum_sensory_input(:,i), stim_off_mean, forward_var*eye(latent_dim));
        prob_sensory_off(2, i) = mvnpdf(cum_sensory_input(:,i), stim_off_mean, pause_var*eye(latent_dim));
        prob_sensory_off(3, i) = mvnpdf(cum_sensory_input(:,i), stim_off_mean, reverse_var*eye(latent_dim));
        
        prob_sensory_on(1, i) = mvnpdf(cum_sensory_input(:,i), stim_on_mean, forward_var*eye(latent_dim));
        prob_sensory_on(2, i) = mvnpdf(cum_sensory_input(:,i), stim_on_mean, pause_var*eye(latent_dim));
        prob_sensory_on(3, i) = mvnpdf(cum_sensory_input(:,i), stim_on_mean, reverse_var*eye(latent_dim));
        
        likelihood(1,i) = ((valence_mat(1,1)*prob_sensory_off(1,i)+valence_mat(1,2)*prob_sensory_on(1,i))/sum(valence_mat(1,:)))*...
                           transition_mat(1,:)*likelihood(:,i-1);
    
        likelihood(2,i) = ((valence_mat(2,1)*prob_sensory_off(2,i)+valence_mat(2,2)*prob_sensory_on(2,i))/sum(valence_mat(2,:)))*...
                           transition_mat(2,:)*likelihood(:,i-1);
    
        likelihood(3,i) = ((valence_mat(3,1)*prob_sensory_off(3,i)+valence_mat(3,2)*prob_sensory_on(3,i))/sum(valence_mat(3,:)))*...
                          transition_mat(3,:)*likelihood(:,i-1);
    
        temp1 = likelihood(1,i)/sum(likelihood(:,i));
        temp2 = likelihood(2,i)/sum(likelihood(:,i));
        temp3 = likelihood(3,i)/sum(likelihood(:,i));
        
        likelihood(1,i) = temp1;
        likelihood(2,i) = temp2;
        likelihood(3,i) = temp3;
    end
    
    forward_behav = likelihood(1, :);
    reverse_behav = likelihood(2, :);
    pause_behav = likelihood(3, :);
%     
end
