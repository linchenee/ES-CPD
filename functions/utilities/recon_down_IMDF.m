function [NMSE] = recon_down_IMDF(angle1,angle2,delay,S,Y,H_down,Ant_Hor,Ant_Ver,lambda_up,lambda_down,f,fc_down,fs)
% Reconstruct the downlink channel 
% --------------------------------------------------------
% Input:
%  angle1:  estimated angle-1 parametes of multipaths through the uplink channel 
%  angle2:  estimated angle-2 parametes of multipaths through the uplink channel 
%  delay:   estimated delay parametes of multipaths through the uplink channel 
%  S:       transmitted downlink pilot symbol
%  Y:       received downlink signal
%  H_down:  ground-truth downlink channel
%  Ant_Hor,Ant_Ver,lambda_up,lambda_down,f,fc_down,fs: see "system parameter.m" for their definitions.
% --------------------------------------------------------
% Output:
%  NMSE:    NMSE for the downlink channel reconstruction
% -------------------------------------------------------
% Written by Lin Chen (linchenee@sjtu.edu.cn)
% -------------------------------------------------------

%% Reconstruct the matrices A_dl, B_dl, and C_dl according to Equation (28) for the downlink channel
A_dl = exp(1j*(0:Ant_Hor-1)'*angle1.*lambda_up./lambda_down);
B_dl = exp(1j*(0:Ant_Ver-1)'*angle2.*lambda_up./lambda_down);
C_dl = exp(1j*2*pi*(f+fc_down)'*(delay./(2*pi*fs)));  

%% Reconstruct the matrix D_dl according to Equation (29) for the downlink channel
temp = S.'*kr(A_dl,B_dl);
D_dl = pinv((C_dl'*C_dl).*(temp'*temp))*(kr(C_dl,temp)'*tens2mat(Y,2));

%% Reconstruct the downlink channel
H_down_recon = zeros(size(H_down));
for i = 1:4
 H_down_recon(:,:,i) = kr(A_dl,B_dl)*diag(D_dl(:,i))*C_dl.';
end

NMSE=(norm(H_down_recon(:)-H_down(:),'fro')/norm(H_down(:),'fro'))^2;

end