function [Y] = PGD(X,W11,W1,p,r,step_size,option,temp1,temp2)
% Perform one iteration of PGD to solve Problem (14)
% --------------------------------------------------------
% Input:
%  X:         matrix of size n1xn2
%  W11:       weight matrix W, as illustrated in Equation (16) 
%  W1:        weight matrix W, as illustrated in Equation (16) 
%  p:         next power of 2 for faster fft
%  r:         tensor frontal rank
%  step_size: step size in PGD
%  option:    opt the convergence tolerance in SVD 
%  temp1:     it is equal to the matrix E*E' in Equation (18)
%  temp2:     it is equal to the matrix Y1_T*E' in Equation (18)
% --------------------------------------------------------
% Output:
%  Y:         matrix of size n1xn2
% ---------------------------------------------
% Written by Lin Chen (linchenee@sjtu.edu.cn)
% -------------------------------------------------------

[n1,n2] = size(X);

Z = Hankel_transform( X-2*step_size*(W11.*(X*temp1-temp2)) ); % Equation (18)

% Perform the rank-r truncated tensor SVD in Equation (19) and 
% the inverse Hankel transform in Euqation (20) jointly
Y = zeros(n1,n2);
for k = 1:n2
    [u,s,v] = svds(Z(:,:,k),r,'largest',option);
    temp = ifft(fft(u,p,1).*fft(conj(v),p,1),[],1);
    Y(:,k) = sum(W1.*temp(1:n1,:)*s,2);
end

end