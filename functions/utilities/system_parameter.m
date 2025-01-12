%% --- System parameter setting ---
c0 = 3e8;                            % speed of light (c)
fc_up = 28.0e9;                      % uplink carrier frequency (f)
fc_down = 28.1e9;                    % downlink carrier frequency (f^{dl})
fs = 645.16129e3/2;                  % subcarrier spacing (\Delta_f)
Nf = 29;                             % number of subcarriers is equal to Nf+1=30
f = linspace(-Nf/2*fs,Nf/2*fs,Nf+1); % frquency samples
lambda_up = c0/fc_up;                % uplink wavelength (\lambda)
lambda_down = c0/fc_down;            % downink wavelength 

Ant_Hor = 8;                         % number of horizontal antennas (M_h)
Ant_Ver = 8;                         % number of vertical antennas (M_v)
Ant_tot = Ant_Hor * Ant_Ver;         % total number of UPA antennas (M)
Dh_up = 0.5*lambda_up;               % antenna spacing in the horizontal direction for the uplink channel (d_h)
Dv_up = 0.5*lambda_up;               % antenna spacing in the vertical direction for the uplink channel (d_v)
Dh_down = Dh_up;                     % antenna spacing in the horizontal direction for the downlink channel
Dv_down = Dv_up;                     % antenna spacing in the horizontal direction for the downlink channel
