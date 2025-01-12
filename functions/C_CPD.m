function [out1,out2,out3,out4] = C_CPD(H,R,K,tol)
% The classical CANDECOMP/PARAFAC decomposition (C-CPD) method to solve Equation (8)
% --------------------------------------------------------
% Input:
%  H:      4D channel tensor of size I1xI2xI3xI4
%  R:      number of rank-1 terms in the CPD
%  K:      maximum iteration number for ALS
%  tol:    termination tolerance for ALS
% --------------------------------------------------------
% Output:
%  out1:   matrix A of size I1xR
%  out2:   matrix B of size I2xR
%  out3:   matrix C of size I3xR
%  out4:   matrix D of size I4xR
% ---------------------------------------------
% Written by Lin Chen (linchenee@sjtu.edu.cn)
% -------------------------------------------------------

I1 = size(H,1);          
I2 = size(H,2);
I3 = size(H,3);
I4 = size(H,4);
H1 = reshape(permute(H,[4,3,2,1]),I2*I3*I4,I1); % mode-1 unfolding of H
H2 = reshape(permute(H,[1,4,3,2]),I3*I4*I1,I2); % mode-2 unfolding of H                
H3 = reshape(permute(H,[2,1,4,3]),I4*I1*I2,I3); % mode-3 unfolding of H 
H4 = reshape(permute(H,[3,2,1,4]),I1*I2*I3,I4); % mode-4 unfolding of H 
Ninit = 10; % number of starting points used
fit_error_previous = realmax;

%% Fit the CPD with several initializations
for ninit = 1:Ninit
 %% Random initialization
 A = randn(I1,R);%+1j*randn(I1,R);
 B = randn(I2,R);%+1j*randn(I2,R);
 C = randn(I3,R);%+1j*randn(I3,R);
 % if I1>=R;A=orth(A);else;A=orth(A')';end
 % if I2>=R;B=orth(B);else;B=orth(B')';end
 % if I3>=R;C=orth(C);else;C=orth(C')';end

 %% Initialization using the method and code provided in reference [20].
 %% [20] Z. Zhou, J. Fang, L. Yang, H. Li, Z. Chen, and R. S. Blum, "Low-rank tensor 
 %% decomposition-aided channel estimation for millimeter wave MIMO-OFDM systems," in 
 %% IEEE Journal on Selected Areas in Communications, vol. 35, no. 7, pp. 1524-1538, July 2017.
 if ninit > Ninit/2
  A = take_svd(H1.',R);
  B = take_svd(H2.',R);
  C = take_svd(H3.',R);
 end

 D = (pinv(kr(A,B,C))*H4).'; % Initialize D by the least-squares
 
 fit_error = norm(H4-kr(A,B,C)*D.','fro');
 %% Iteration of ALS
 for k = 1:K
    fit_error_old = fit_error;
    % A = (pinv(kr(B,C,D))*X1).';
    % B = (pinv(kr(C,D,A))*X2).';
    % C = (pinv(kr(D,A,B))*X3).';
    % D = (pinv(kr(A,B,C))*X4).';
    A = (inv((B'*B).*(C'*C).*(D'*D))* (kr(B,C,D)'*H1) ).';
    B = (inv((C'*C).*(D'*D).*(A'*A))* (kr(C,D,A)'*H2) ).'; 
    C = (inv((D'*D).*(A'*A).*(B'*B))* (kr(D,A,B)'*H3) ).';
    D = (inv((A'*A).*(B'*B).*(C'*C))* (kr(A,B,C)'*H4) ).';
    fit_error = norm(H4-kr(A,B,C)*D.','fro');
    if (abs((fit_error-fit_error_old)/fit_error_old) < tol)
       % fprintf('Iteration number of C-CPD is %d\n', k);
       break;
    end
 end   
 if fit_error < fit_error_previous
  out1 = A; 
  out2 = B; 
  out3 = C; 
  out4 = D;
  fit_error_previous = fit_error;
 end
end    
end

function U_matr = take_svd(Matr,rank)
    [U_matr,S,~] = svd(Matr,'econ');
    [m,n] = size(U_matr);
    if n >= rank
        U_matr = U_matr(:,1:rank);%*sqrt(S(1:rank,1:rank))
    else
%             rng(10,'v5uniform');
        U_matr = [U_matr,randn(m,rank-n)+1j*randn(m,rank-n)];
    end
end