function [NMSE] = recon_down_CPD(A_up,B_up,C_up,S,Y,H_down,P,Ant_Hor,Ant_Ver,lambda_up,lambda_down,f,fc_down,fs)
% Reconstruct the downlink channel 
% --------------------------------------------------------
% Input:
%  A_up:    estimated matrix A in the CPD through the uplink channel 
%  B_up:    estimated matrix B in the CPD through the uplink channel 
%  C_up:    estimated matrix C in the CPD through the uplink channel 
%  S:       transmitted downlink pilot symbol
%  Y:       received downlink signal
%  H_down:  ground-truth downlink channel
%  P:       the estimated number of dominant paths
%  Ant_Hor,Ant_Ver,lambda_up,lambda_down,f,fc_down,fs: see "system parameter.m" for their definitions.
% --------------------------------------------------------
% Output:
%  NMSE:    NMSE for the downlink channel reconstruction
% -------------------------------------------------------
% Written by Lin Chen (linchenee@sjtu.edu.cn)
% -------------------------------------------------------

angle1 = zeros(1,P);
angle2 = angle1;
delay = angle1;
for i = 1:P
    angle1(1,i) = PUMA(A_up(:,i)); % estimate the angle1 parametes of multipaths          
    angle2(1,i) = PUMA(B_up(:,i)); % estimate the angle2 parametes of multipaths           
    delay(1,i) = PUMA(C_up(:,i));  % estimate the delay parametes of multipaths         
end

%% Reconstruct the matrices A_dl, B_dl, and C_dl according to Equation (28) for the downlink channel
A_dl = exp(1j*(0:Ant_Hor-1)'*angle1.*lambda_up./lambda_down);
B_dl = exp(1j*(0:Ant_Ver-1)'*angle2.*lambda_up./lambda_down);
C_dl = exp(1j*2*pi*(f+fc_down)'*(delay./(2*pi*fs)));  

%% Reconstruct the matrix D_dl according to Equation (29) for the downlink channel
temp = S.'*kr(A_dl,B_dl);
D_dl = pinv((C_dl'*C_dl).*(temp'*temp)) * (kr(C_dl,temp)'*tens2mat(Y,2));

%% Reconstruct the downlink channel
H_down_recon = zeros(size(H_down));
for i = 1:4
    H_down_recon(:,:,i) = kr(A_dl,B_dl)*diag(D_dl(:,i))*C_dl.'; 
end

NMSE = (norm(H_down_recon(:)-H_down(:),'fro')/norm(H_down(:),'fro'))^2;

end