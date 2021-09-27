function [A_final, B_final, Q_final, R_final] = set_up_optimal_control(latent_dim, neural_dim, tau_nu, ...
                                            drift_nu_gamma, drift_nu, b_matrix, tau_gamma, drift_gamma,...
                                            beta, aux_eig_A, latent_penalty, energy_penalty, deriv_penalty)
                                        
% This function creates the matrices to set up the optimal control problem.

A_modified = [(drift_gamma/tau_gamma)*eye(latent_dim), (1/tau_gamma)*beta*eye(latent_dim), zeros(latent_dim, neural_dim);...
               (drift_nu_gamma/tau_nu)*eye(latent_dim), (drift_nu/tau_nu)*eye(latent_dim), (1/tau_nu)*b_matrix;...
               zeros(neural_dim, latent_dim), zeros(neural_dim, latent_dim), zeros(neural_dim)];
aux_A_matrix = diag(aux_eig_A);
A_final = blkdiag(A_modified, aux_A_matrix); 

B_final = [zeros(2*latent_dim, neural_dim); eye(neural_dim); zeros(latent_dim, neural_dim)];

Q = latent_penalty*eye(latent_dim); S = energy_penalty*eye(neural_dim); R= deriv_penalty*eye(neural_dim);
Q_modified = [Q zeros(latent_dim, neural_dim) -Q; ...
                zeros(neural_dim, latent_dim), S, zeros(neural_dim, latent_dim);...
                -Q zeros(latent_dim, neural_dim) Q];
Q_final = blkdiag(1e-5*eye(latent_dim), Q_modified);
R_final = R;
         