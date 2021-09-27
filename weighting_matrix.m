function [ b_Matrix ] = weighting_matrix( m_dim, n_dim, overlap, spread, intensity )
 
if nargin<5
    intensity = 10;
end

if nargin<4
    spread = 10;
end

stim_A_num = ceil((overlap + n_dim)/ 2);

stim_A_mean =  ceil(stim_A_num/2);
std_A = spread;
stim_A_weights = pdf('Normal', 1:n_dim, stim_A_mean, std_A);

stim_B_mean =  n_dim - stim_A_mean + 1; 
std_B = spread;
stim_B_weights = pdf('Normal', 1:n_dim, stim_B_mean, std_B);

b_Matrix = intensity * [stim_A_weights; stim_B_weights];

end

