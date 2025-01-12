function [H_CPD,H_IMDF] = reshape_channel_for_CPDandIMDF(H,Nf,Ant_Hor,Ant_Ver,Ant_tot)
% Reshape the channel to form the respective inputs for C-CPD, ES-CPD, and IMDF methods
% --------------------------------------------------------
% Input:
%  H:                          channel
%  Nf,Ant_Hor,Ant_Ver,Ant_tot: see "system parameter.m" for their definitions.
% --------------------------------------------------------
% Output:
%  H_CPD:                      used as the input for both the C-CPD and ES-CPD methods
%  H_IMDF:                     used as the input for the IMDF method
% ---------------------------------------------
% Written by Lin Chen (linchenee@sjtu.edu.cn)
% -------------------------------------------------------

H_CPD = zeros(Ant_Hor,Ant_Ver,Nf+1,4);
for i = 1:Ant_Hor
    for j = 1:Ant_Ver
      H_CPD(i,j,:,:) = H((i-1)*Ant_Ver+j,:,:);
   end
end

H_IMDF = zeros(2*(Nf+1),2*Ant_tot);
temp = 0;
for i = 1:2
   for j = 1:2
        temp = temp+1;
        H_IMDF((i-1)*(Nf+1)+1:i*(Nf+1), (j-1)*Ant_tot+1:j*Ant_tot) = H(:,:,temp).';
   end
end

end