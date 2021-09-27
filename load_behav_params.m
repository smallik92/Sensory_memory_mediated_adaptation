function [forward_var, pause_var, reverse_var, ...
    stim_on_mean, stim_off_mean, transition_mat, valence_mat, relative_weights] = load_behav_params() 

    forward_var = 0.15; 
    pause_var = 0.2;
    reverse_var = 0.3; 
    
    stim_on_mean = [1; 0];
    stim_off_mean = [0; 0];
    
    transition_mat = [0.75, 0.15, 0.1; 
                      0.1, 0.7, 0.2; 
                      0.15, 0.15, 0.75];
                  
    valence_mat = [0.3, 0.6;
                   0.35, 0.2;
                   0.35, 0.2];
               
    relative_weights = [0.5, 0.5];    
end