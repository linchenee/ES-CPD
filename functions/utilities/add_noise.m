function [H_ups,H_downs] = add_noise(H_up,H_down,SNR)
% Add the Gaussian noise to the uplink and downlink channels
% --------------------------------------------------------
% Input:
%  H_up:    noiseless uplink channel
%  H_down:  noiseless uplink channel
%  SNR:     signal-to-noise ratio
% --------------------------------------------------------
% Output:
%  H_ups:   noisy uplink channel
%  H_downs: noisy downlink channel
% ---------------------------------------------
% Written by Lin Chen (linchenee@sjtu.edu.cn)

Noise_tens = randn(size(H_up))+1j*randn(size(H_up)); % Noise tensor for the uplink channel
sigma = (10^(-SNR/10))*norm(H_up(:),'fro')^2/(size(H_up,1)*size(H_up,2)*size(H_up,3)); % noise variance
H_ups = H_up+sqrt(sigma)*Noise_tens./sqrt(2);

Noise_tens = randn(size(H_down))+1j*randn(size(H_down)); % Noise tensor for the downlink channel
sigma1 = (10^(-SNR/10))*norm(H_down(:),'fro')^2/(size(H_down,1)*size(H_down,2)*size(H_down,3));
H_downs = H_down+sqrt(sigma1)*Noise_tens./sqrt(2);

end
