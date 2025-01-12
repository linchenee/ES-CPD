function [loss_up1,loss_up2,loss_up3,loss_up4,loss_down1,loss_down2,loss_down3,loss_down4] = generate_pathloss(L,Seed)
% Generate the path-losses for the dual-polarized uplink and downlink channels
% --------------------------------------------------------
% Input:
%  L:          number of multipaths
%  Seed:       random seed
% --------------------------------------------------------
% Output:
%  loss_up1:   uplink path-losses in the (1,1)th polarization mode
%  loss_up2:   uplink path-losses in the (1,2)th polarization mode
%  loss_up3:   uplink path-losses in the (2,1)th polarization mode
%  loss_up4:   uplink path-losses in the (2,2)th polarization mode
%  loss_down1: downlink path-losses in the (1,1)th polarization mode
%  loss_down2: downlink path-losses in the (1,2)th polarization mode
%  loss_down3: downlink path-losses in the (2,1)th polarization mode
%  loss_down4: downlink path-losses in the (2,2)th polarization mode
% --------------------------------------------------------
% Written by Lin Chen (linchenee@sjtu.edu.cn)
% --------------------------------------------------------

K_factordB = 20; % K_factor: the power ratio between the LOS path and NLOS paths is set to 20dB.
Kf = 10^(K_factordB/10);
PowerdB = -100;
sigma = sqrt(10^(PowerdB/10));

%% uplink
rng(Seed);
loss_up = sigma*(randn(L,1)+randn(L,1)*1j)/sqrt(2);
temp1 = exp(1j*pi*randn(L,4));
loss_up1 = [sqrt(Kf/(Kf+1))*loss_up(1,1)*temp1(1,1); sqrt(1/(Kf+1))*loss_up(2:end,1).*temp1(2:end,1)];
loss_up2 = [sqrt(Kf/(Kf+1))*loss_up(1,1)*temp1(1,2); sqrt(1/(Kf+1))*loss_up(2:end,1).*temp1(2:end,2)];
loss_up3 = [sqrt(Kf/(Kf+1))*loss_up(1,1)*temp1(1,3); sqrt(1/(Kf+1))*loss_up(2:end,1).*temp1(2:end,3)];
loss_up4 = [sqrt(Kf/(Kf+1))*loss_up(1,1)*temp1(1,4); sqrt(1/(Kf+1))*loss_up(2:end,1).*temp1(2:end,4)];
 
%% downlink
rng(Seed+1);
loss_down = sigma*(randn(L,1)+randn(L,1)*1j)/sqrt(2);
temp2 = exp(1j*pi*randn(L,4));
loss_down1 = [sqrt(Kf/(Kf+1))*loss_down(1,1)*temp2(1,1); sqrt(1/(Kf+1))*loss_down(2:end,1).*temp2(2:end,1)];
loss_down2 = [sqrt(Kf/(Kf+1))*loss_down(1,1)*temp2(1,2); sqrt(1/(Kf+1))*loss_down(2:end,1).*temp2(2:end,2)];
loss_down3 = [sqrt(Kf/(Kf+1))*loss_down(1,1)*temp2(1,3); sqrt(1/(Kf+1))*loss_down(2:end,1).*temp2(2:end,3)];
loss_down4 = [sqrt(Kf/(Kf+1))*loss_down(1,1)*temp2(1,4); sqrt(1/(Kf+1))*loss_down(2:end,1).*temp2(2:end,4)];

end
