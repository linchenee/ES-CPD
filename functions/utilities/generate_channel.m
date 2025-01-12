function [H_up,H_down] = generate_channel(L,Seed,varphi,phi,delay,lambda_up,lambda_down,Dh_up,Dv_up,Dh_down,Dv_down,Ant_Hor,Ant_Ver,f,fc_up,fc_down)
%  Generate the ground-truth uplink and downlink channels
% --------------------------------------------------------
% Input:
%  L:      number of multipaths
%  Seed:   random seed
%  varphi: elevation angle parameters of multipaths
%  phi:    azimuth angle parameters of multipaths
%  delay:  delay parameters of multipaths
%  lambda_up,lambda_down,Dh_up,Dv_up,Dh_down,Dv_down,Ant_Hor,Ant_Ver,f,fc_up,fc_down: see "system parameter.m" for their definitions.
% --------------------------------------------------------
% Output:
%  H_up:   uplink channel
%  H_down: downlink channel 
% ---------------------------------------------
% Written by Lin Chen (linchenee@sjtu.edu.cn)
% -------------------------------------------------------

temp1 = (sin(varphi).*cos(phi));
temp2 = cos(varphi);

%% uplink steering vectors
SteerHor_up = exp(1j*2*pi/lambda_up*Dh_up*(0:Ant_Hor-1)'*temp1); % antenna steering vector in the horizontal direction 
SteerVer_up = exp(1j*2*pi/lambda_up*Dv_up*(0:Ant_Ver-1)'*temp2); % antenna steering vector in the vertical direction 
Steer_up = kr(SteerHor_up,SteerVer_up);
Delay_up = exp(-1j*2*pi*(f+fc_up)'*delay); % delay steering vector
%% downlink steering vectors
SteerHor_down = exp(1j*2*pi/lambda_down*Dh_down*(0:Ant_Hor-1)'*temp1);
SteerVer_down = exp(1j*2*pi/lambda_down*Dv_down*(0:Ant_Ver-1)'*temp2);
Steer_down = kr(SteerHor_down,SteerVer_down);
Delay_down = exp(-1j*2*pi*(f+fc_down)'*delay);

%% generalized the path-losses for the dual-polarized uplink and downlink channels
[loss_up1,loss_up2,loss_up3,loss_up4,loss_down1, ...
  loss_down2,loss_down3,loss_down4] = generate_pathloss(L,Seed);

%% uplink channel
H_up(:,:,1) = Steer_up*diag(loss_up1)*Delay_up.'; % in the (1,1)th polarization mode
H_up(:,:,2) = Steer_up*diag(loss_up2)*Delay_up.'; % in the (1,2)th polarization mode
H_up(:,:,3) = Steer_up*diag(loss_up3)*Delay_up.'; % in the (2,1)th polarization mode
H_up(:,:,4) = Steer_up*diag(loss_up4)*Delay_up.'; % in the (2,2)th polarization mode
%% downlink channel
H_down(:,:,1) = Steer_down*diag(loss_down1)*Delay_down.'; 
H_down(:,:,2) = Steer_down*diag(loss_down2)*Delay_down.';  
H_down(:,:,3) = Steer_down*diag(loss_down3)*Delay_down.'; 
H_down(:,:,4) = Steer_down*diag(loss_down4)*Delay_down.'; 

end
