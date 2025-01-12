function [out1,out2,out3,out4] = ES_CPD(H,R,K,T,tol,step_size,Ninit,A_init,B_init,C_init)
% Exponential structured CANDECOMP/PARAFAC decomposition (ES-CPD) method to solve Equation (12)
% --------------------------------------------------------
% Input:
%  H:         4D channel tensor of size I1xI2xI3xI4
%  R:         number of rank-1 terms in CPD
%  K:         maximum iteration number for ALS
%  T:         maximum iteration number for PGD
%  tol:       termination tolerance for ALS
%  step_size: step size in PGD
%  Ninit:     number of starting points used. The first (Ninit-1) initializations 
%             are random, and the last initialization uses A_init, B_init, C_init.
%  A_init:    starting point for the matrix A in the last initialization 
%  B_init:    starting point for the matrix B in the last initialization 
%  C_init:    starting point for the matrix C in the last initialization 
% --------------------------------------------------------
% Output:
%  out1:      matrix A of size I1xR
%  out2:      matrix B of size I2xR
%  out3:      matrix C of size I3xR
%  out4:      matrix D of size I4xR
% -------------------------------------------------------
% Written by Lin Chen (linchenee@sjtu.edu.cn)
% -------------------------------------------------------

Rank = 1;             % tensor frontal rank
option.tol = 1e-3;    % opt the convergence tolerance in SVD
[I1,I2,I3,I4] = size(H);
p1 = 2^nextpow2(I1);  % next power of 2 for faster fft
p2 = 2^nextpow2(I2);
p3 = 2^nextpow2(I3);
[W1,~] = construct_W(I1,Rank);
[W2,~] = construct_W(I2,Rank);
[W3,~] = construct_W(I3,Rank);   
[W11,~] = construct_W(I1,R);
[W22,~] = construct_W(I2,R);
[W33,~] = construct_W(I3,R); 

H1 = reshape(permute(H,[4,3,2,1]),I2*I3*I4,I1); % mode-1 unfolding of H
H2 = reshape(permute(H,[1,4,3,2]),I3*I4*I1,I2); % mode-2 unfolding of H                
H3 = reshape(permute(H,[2,1,4,3]),I4*I1*I2,I3); % mode-3 unfolding of H 
H4 = reshape(permute(H,[3,2,1,4]),I1*I2*I3,I4); % mode-4 unfolding of H
H1_T = H1.';
H2_T = H2.';                   
H3_T = H3.'; 
H4_T = H4.';

fit_error_previous = realmax;
%% Fit the CPD with several initializations
for ninit = 1:Ninit
    %% Initialization
    rng(ninit,'v5uniform');
    A = randn(I1,R)+1j*randn(I1,R);
    B = randn(I2,R)+1j*randn(I2,R);
    C = randn(I3,R)+1j*randn(I3,R);
    if I1>=R; A=orth(A); else A=orth(A')'; end
    if I2>=R; B=orth(B); else B=orth(B')'; end
    if I3>=R; C=orth(C); else C=orth(C')'; end
    %% In this paper, we use the result of the IMDF method for the last initialization, 
    %% which can be adopted to use other methods as needed.
    if ninit == Ninit
        A = A_init;
        B = B_init;
        C = C_init;
    end
    D = H4_T*pinv(kr(A,B,C).'); % Initialize D by the least-squares
    fit_error = norm(H4_T-D*kr(A,B,C).','fro');
   
    %% Iteration of ALS
    for k = 1:K
        fit_error_old = fit_error;

        %% Update A by solving Equation (14) via the PGD
        E1 = kr(B,C,D).';
        temp1 = E1*E1';
        temp2 = H1_T*E1';
        %% The first iteration in PGD, i.e., t=1
        A = PGD(H1_T*pinv(E1),W11,W1,p1,Rank,step_size,option,temp1,temp2);
        A1 = A;
        fitA_old = norm(H1_T-A1*E1,'fro');
        %% Form the second iteration to the Tth iteration in PGD
        for t = 2:T
            A1 = PGD(A1,W11,W1,p1,Rank,step_size,option,temp1,temp2);
            fitA = norm(H1_T-A1*E1,'fro');
            if fitA > fitA_old
                break;
            else
                A = A1;
                fitA_old = fitA;
            end
        end

        %% Update B by solving Equation (14) via the PGD
        E2 = kr(C,D,A).';
        temp1 = E2*E2';
        temp2 = H2_T*E2';
        B = PGD(H2_T*pinv(E2),W22,W2,p2,Rank,step_size,option,temp1,temp2);
        B1 = B;
        fitB_old = norm(H2_T-B1*E2,'fro');
        for t = 2:T
            B1 = PGD(B1,W22,W2,p2,Rank,step_size,option,temp1,temp2);
            fitB = norm(H2_T-B1*E2,'fro');
            if fitB > fitB_old
                break;
            else
                B = B1;
                fitB_old = fitB;
            end
        end

        %% Update C by solving Equation (14) via the PGD
        E3 = kr(D,A,B).'; 
        temp1 = E3*E3';
        temp2 = H3_T*E3';
        C = PGD(H3_T*pinv(E3),W33,W3,p3,Rank,step_size,option,temp1,temp2);
        C1 = C;
        fitC_old = norm(H3_T-C1*E3,'fro');
        for t = 2:T
            C1 = PGD(C1,W33,W3,p3,Rank,step_size,option,temp1,temp2);
            fitC = norm(H3_T-C1*E3,'fro');
            if fitC > fitC_old
                break;
            else
                C = C1;
                fitC_old = fitC;
            end
        end

        %% Update D by the least-squares
        E4 = kr(A,B,C).';
        D = H4_T*pinv(E4); 
        
        %% Stop if the termination condition is met
        fit_error = norm(H4_T-D*E4,'fro');
        if (abs((fit_error-fit_error_old)/fit_error_old)<tol)
         % fprintf('Iteration number of ES-CPD is %d\n', k);
        break;
       end
    end

    %% Select this initialization if it is better than previous one
    if fit_error < fit_error_previous
      out1 = A; 
      out2 = B; 
      out3 = C; 
      out4 = D;
      fit_error_previous = fit_error;
    end
end

end