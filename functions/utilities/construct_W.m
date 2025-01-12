function [W,W_tilde] = construct_W(n1,n2)
% Construct the matrices W and \tilde{W}, as illustrated in Equation (16)
% --------------------------------------------------------
% Input:
%  n1:      row number
%  n2:      column number
% --------------------------------------------------------
% Output:
%  W:       the matrix W, as illustrated in Equation (16)
%  W_tilde: the matrix \tilde{W}, as illustrated in Equation (16)
% ---------------------------------------------
% Written by Lin Chen (linchenee@sjtu.edu.cn)
% -------------------------------------------------------

if mod(n1,2) == 0
 W_tilde = repmat([1:n1/2,n1/2:-1:1]',1,n2);
else
 W_tilde = repmat([1:(n1+1)/2,(n1-1)/2:-1:1]',1,n2);
end
W = ones(n1,n2)./W_tilde;

end
