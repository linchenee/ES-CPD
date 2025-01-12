%% Run this code to obtain the NMSE results shown in Figure 2(a) for the proposed ES-CPD method.
clc
clear all
close all
addpath(genpath('functions'));
addpath(genpath('functions\utilities'));

system_parameter;

trials = 1000;         % number of Monte-Carlo trials
SNRs = [0,5,10,15,20]; % signal-to-noise ratio

flag_IMDF = true;      % IMDF method
flag_C_CPD = true;     % C-CPD method
flag_ES_CPD = true;    % ES-CPD method   

NMSE_IMDF = zeros(length(SNRs),trials); 
NMSE_C_CPD = NMSE_IMDF;
NMSE_ES_CPD = NMSE_IMDF;

for i = 1:length(SNRs)
    SNR = SNRs(i);
    fprintf('SNR=%d dB\n',SNR);
    parfor j = 1:trials
    % for j = 1:trials
      rng(2*j);
      fprintf('The %dth trial\n',j);
      L = randi(50);           % number of multipaths
      P = L;                   % estimated number of dominant paths
      varphi = pi*rand(1,L)/2; % elevation angle parameters of multipaths
      phi  = pi*rand(1,L);     % azimuth angle parameters of multipaths
      tau = 3e-7*rand(1,L);    % delay parameters of multipaths
      
      %% Generate the ground-truth uplink and downlink channels
      [H_up, H_down] = generate_channel(L,j,varphi,phi,tau,lambda_up,lambda_down,Dh_up,Dv_up,Dh_down,Dv_down,Ant_Hor,Ant_Ver,f,fc_up,fc_down);
      %% Add Gaussian noise to channels
      [H_up_noisy,H_down_noisy] = add_noise(H_up,H_down,SNR);
      %% Reshape the uplink channel to form the respective inputs for C-CPD, ES-CPD, and IMDF methods
      [H_up_noisy_4D, H_up_noisy_IMDF] = reshape_channel_for_CPDandIMDF(H_up_noisy,Nf,Ant_Hor,Ant_Ver,Ant_tot);
      %% Generate the transmitted pilot symbol, S, and the received signal, Y, for the downlink channel
      [S,Y] = generate_down_signal(H_down_noisy,Ant_tot);
      
      if (flag_IMDF)
          left = Ant_Ver;  
          right = Ant_Hor; 
          [~, wr, wx, wy] = IMDF_multi(H_up_noisy_IMDF, Nf+1, left, right, P);
          NMSE_IMDF(i,j) = recon_down_IMDF(-wy,-wx,wr,S,Y,H_down,Ant_Hor,Ant_Ver,lambda_up,lambda_down,f,fc_down,fs);
      end

      if (flag_C_CPD)
          maxiter = 20;
          tol = 1e-4;
          [A0,B0,C0,~] = C_CPD(H_up_noisy_4D,P,maxiter,tol);
          NMSE_C_CPD(i,j) = recon_down_CPD(A0,B0,C0,S,Y,H_down,P,Ant_Hor,Ant_Ver,lambda_up,lambda_down,f,fc_down,fs);
      end

      if (flag_ES_CPD)
          maxiter1 = 20;
          maxiter2 = 1;
          tol = 1e-4;
          step_size = 0.02;          
          Ninit = 2;         
          A_init = exp(1j*(0:(size(H_up_noisy_4D,1)-1))'*(-wy));
          B_init = exp(1j*(0:(size(H_up_noisy_4D,2)-1))'*(-wx));
          C_init = exp(1j*(0:(size(H_up_noisy_4D,3)-1))'*wr);
          [A1, B1, C1, ~] = ES_CPD(H_up_noisy_4D,P,maxiter1,maxiter2,tol,step_size,Ninit,A_init,B_init,C_init);
          NMSE_ES_CPD(i,j) = recon_down_CPD(A1,B1,C1,S,Y,H_down,P,Ant_Hor,Ant_Ver,lambda_up,lambda_down,f,fc_down,fs);
      end 
   end
end

NMSE_IMDF_ave = mean(NMSE_IMDF,2);
NMSE_C_CPD_ave = mean(NMSE_C_CPD,2);
NMSE_ES_CPD_ave = mean(NMSE_ES_CPD,2);

figure(1);
semilogy(SNRs,NMSE_IMDF_ave,'g-^',SNRs,NMSE_C_CPD_ave,'b-x',SNRs,NMSE_ES_CPD_ave,'r-*','LineWidth',2); 
xlabel('SNR (dB)','fontsize',12,'FontName','Times new roman');
ylabel('NMSE','fontsize',12,'FontName','Times new roman');
xlim([0 20]);
ylim([0.01 0.2]);
set(gca,'FontName','Times New Roman','FontSize',12);
h=legend('IMDF','C-CPD','ES-CPD');
set(h,'Interpreter','latex','fontsize',10,'FontName','Times new roman');