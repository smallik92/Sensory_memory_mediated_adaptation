function [Calcium_response, x_response, nu_response, gamma_response] = simulate_neural_activity_vF(dt, t_init, t_on, trial_dur, t_index,...
                                                              K, z, noise_level, SEED,B_final, R_final,...
                                                              tau_nu, drift_nu, drift_nu_gamma, b_matrix,...
                                                              tau_gamma, drift_gamma, beta, ...
                                                              latent_dim, neural_dim, overlap, baseline, low_asymp, up_asymp,...
                                                              curve, low_pass, gain)
                                                          
   rng(SEED);
   
   stim_A_num = ceil((overlap + neural_dim)/ 2);
   pure_stim_A = stim_A_num - overlap;  
   
   total_dim = neural_dim + 3*latent_dim ;
   
   x_response = zeros(neural_dim, length(t_index));
   nu_response = zeros(latent_dim, length(t_index));
   gamma_response = zeros(latent_dim, length(t_index));
   
   x_calcium_raw = zeros(neural_dim, length(t_index));
   x_scaled = zeros(neural_dim, length(t_index));
   x_scaled(:,1) = low_asymp + (up_asymp - low_asymp)*sigmoid(x_scaled(:,1), baseline, curve);
   
   for i = 1:length(t_index)-1
        Kt = reshape(K(i+1,:),[total_dim, total_dim]);
        W = inv(R_final)*B_final'*Kt; % feedback matrix
        W_mat{i} = -W;
        
        W_1 = -W(:,1:latent_dim);     
        W_2 = -W(:,latent_dim+1:2*latent_dim);
        W_3 = -W(:,2*latent_dim+1:2*latent_dim+neural_dim);
        W_4 = -W(:,neural_dim+2*latent_dim+1:end);
        
        gamma_response(:,i+1) = gamma_response(:,i)+dt*((drift_gamma/tau_gamma)*gamma_response(:,i) + (1/tau_gamma)*beta*nu_response(:,i)+ noise_level*randn);
        nu_response(:,i+1) = nu_response(:,i)+dt*((drift_nu/tau_nu)*nu_response(:,i)+ (drift_nu_gamma/tau_nu)*gamma_response(:,i)+(1/tau_nu)*b_matrix*x_response(:,i)+ noise_level*randn);
        x_response(:,i+1) = x_response(:,i)+dt*(W_1*gamma_response(:,i)+W_2*nu_response(:,i)+W_3*x_response(:,i)+ W_4*z(:,i)+ noise_level*randn);

    %   Scale to firing rate
        x_scaled(:,i+1) = low_asymp + (up_asymp - low_asymp)*sigmoid(x_response(:,i+1), baseline, curve);
            
   end
    
   %% Use a low pass filter to obtain Calcium dynamics
 
    for i = 1:length(t_index)-1 
        x_calcium_raw(:, i+1) = x_calcium_raw(:, i)+dt*(low_pass*x_calcium_raw(:, i)+ gain*x_scaled(:, i));
    end
    intensity = 4;
    
    %% Baseline Correction    

    F0 = median(mean(x_calcium_raw(1:stim_A_num,floor((10*trial_dur+t_init+t_on)/dt)+1:(11*trial_dur)/dt)));
    x_calcium_trace = mean(x_calcium_raw(1:stim_A_num, :));
    x_calcium_trace(:,1:75) = F0;
    Calcium_response = intensity*(x_calcium_trace/F0 - 1);
                                                          
end                                                          