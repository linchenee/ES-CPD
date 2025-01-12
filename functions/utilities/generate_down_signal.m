function [S,Y] = generate_down_signal(H,Ant_tot)
% Generate the transmitted pilot symbol and received signal for the downlink channel
% --------------------------------------------------------
% Input:
%  H:       downlink channel
%  Ant_tot: antenna number at the base station
% --------------------------------------------------------
% Output:
%  S:       transmitted pilot symbol
%  Y:       received signal 
% ---------------------------------------------
% Written by Lin Chen (linchenee@sjtu.edu.cn)
% -------------------------------------------------------

S = randn(Ant_tot,Ant_tot)+1j*randn(Ant_tot,Ant_tot); 
% S = randn(Ant_tot,Q)+1j*randn(Ant_tot,Q); 
S = orth(S);
for i = 1:4
    Y(:,:,i) = S.'*squeeze(H(:,:,i));
end

end
