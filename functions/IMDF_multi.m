%% If you find this code useful, please cite our paper:
% C. Qian, X. Fu, N. D. Sidiropoulos and Y. Yang, "Tensor-Based Channel 
% Estimation for Dual-Polarized Massive MIMO Systems," in IEEE Transactions 
% on Signal Processing. doi: 10.1109/TSP.2018.2873506

%% Thanks.

function [Hhat, wr, wx, wy] = IMDF_multi(H, Mr, left, right, P)

fun = @(a,b) exp(1j*[0:a-1]'*b(:)');
NN = left*right;

if Mr>1
    i3 = 0;
    for i1 = 1:2
        for i2 = 1:2
            i3 = i3 + 1;
            Hbreve(:,i3) = reshape(H((i1-1)*Mr+1:i1*Mr,(i2-1)*NN+1:i2*NN),Mr*NN,1); %H4
        end
    end
    
    X = reshape(Hbreve, [Mr,left,right,size(Hbreve,2)]);
    [N1, N2, N3, nSnap] = size(X);
    [K1,L1,K2,L2,K3,L3] = findPairs_3d(N1,N2,N3,nSnap);
    
    Ztot = [];
    for i1 = 1:nSnap
        X1 = squeeze(X(:,:,:,i1));
        X1s = IMDFS3(X1, N1, N2, N3, K1, K2, K3);
        [LL1,LL2]=size(X1s);
        x1s = X1s(:);
        Y1s = reshape(conj(x1s(end:-1:1)), LL1, LL2);
        Z = [X1s, Y1s];
        Ztot = [Ztot, Z];
    end
    
    [U, S, V] = svd(Ztot);
    us = U(:, 1:P);
    J1 = kron(kron(eye(K1-1,K1), eye(K2)), eye(K3));
    J2 = kron(kron([zeros(K1-1,1),eye(K1-1)], eye(K2)), eye(K3));
    
    U1 = J1*us;
    U2 = J2*us;
    [Tsp, ~] = eig((U1'*U1)\U1'*U2);
    Ahat = us*Tsp;
    
    A11 = J1*Ahat;
    A12 = J2*Ahat;
    
    for i = 1:P
        wr(i) = angle(A11(:,i)'*A12(:,i));
    end
    
    J1 = kron(kron(eye(K1), eye(K2-1,K2)), eye(K3));
    J2 = kron(kron(eye(K1), [zeros(K2-1,1),eye(K2-1)]), eye(K3));
    A21 = J1*Ahat.'';
    A22 = J2*Ahat.'';
    for i = 1:P
        wx(i) = angle(A21(:,i)'*A22(:,i));
    end
    
    J1 = kron(kron(eye(K1), eye(K2)), eye(K3-1,K3));
    J2 = kron(kron(eye(K1), eye(K2)), [zeros(K3-1,1),eye(K3-1)]);
    A31 = J1*Ahat.'';
    A32 = J2*Ahat.'';
    for i = 1:P
        wy(i) = angle(A31(:,i)'*A32(:,i));
    end
    
    Ar = fun(Mr,wr);
    Ax = fun(left,wx);
    Ay = fun(right,wy);
    
    A = krb(Ay.'',krb(Ax.'',Ar));
    B = (pinv(A)*Hbreve).';

    H2 = A*B.';
    i3 = 0;
    for i1 = 1:2
        for i2 = 1:2
            i3 = i3 + 1;
            Hhat((i1-1)*Mr+1:i1*Mr,(i2-1)*NN+1:i2*NN) = reshape(H2(:,i3),Mr,NN);
        end
    end
    
elseif Mr==1
    
    i3 = 0;
    for i1 = 1:2
        for i2 = 1:2
            i3 = i3 + 1;
            Hbreve(:,i3) = reshape(H(i1,(i2-1)*NN+1:i2*NN),NN,1);
        end
    end
    X = reshape(Hbreve, [left,right,4]);
    [N1, N2, nSnap] = size(X);
    [K1,L1,K2,L2] = findPairs_2d(N1, N2);
    
    Ztot = [];
    for i = 1:4
        X1 = squeeze(X(:,:,i1));
        X1s = IMDFS2(X1, N1, N2, K1, K2);
        [LL1,LL2]=size(X1s);
        x1s = X1s(:);
        Y1s = reshape(conj(x1s(end:-1:1)), LL1, LL2);
        Z = [X1s, Y1s];
        Ztot = [Ztot, Z];
    end
    
    [U, S, V] = svd(Ztot);
    us = U(:, 1:P);
    J1 = kron(eye(K1-1,K1), eye(K2));
    J2 = kron([zeros(K1-1,1),eye(K1-1)], eye(K2));
    
    U1 = J1*us;
    U2 = J2*us;
    [Tsp, ~] = eig((U1'*U1)\U1'*U2);
    Ahat = us*Tsp;
    
    A11 = J1*Ahat.'';
    A12 = J2*Ahat.'';
    
    for i = 1:P
        wx(i) = angle(A11(:,i)'*A12(:,i));
    end
    
    J1 = kron(eye(K1), eye(K2-1,K2));
    J2 = kron(eye(K1), [zeros(K2-1,1),eye(K2-1)]);
    A21 = J1*Ahat.'';
    A22 = J2*Ahat.'';
    for i = 1:P
        wy(i) = angle(A21(:,i)'*A22(:,i));
    end
    
    Ar = ones(1,P);
    Ax = exp(1j*[0:left-1]'*wx(:)');
    Ay = exp(1j*[0:right-1]'*wy(:)');
    At = krb(Ay.'',Ax.'');
    B = (pinv(At)*Hbreve).';
    
    H2 = At*B.';
    i3 = 0;
    for i1 = 1:2
        for i2 = 1:2
            i3 = i3 + 1;
            Hhat(i1,(i2-1)*NN+1:i2*NN) = reshape(H2(:,i3),1,NN);
        end
    end
    
end


end

function Z=IMDFS3(X, N1, N2, N3, K1, K2, K3)
L1=N1-K1+1;
L2=N2-K2+1;
L3=N3-K3+1;
for i1=1:L1
    X1(:,:,:,i1)=X(i1:N1-L1+i1,:,:);
end
for i2=1:L2
    X2(:,:,:,:,i2)=X1(:,i2:N2-L2+i2,:,:);
end
for i3=1:L3
    Y(:,:,:,:,:,i3)=X2(:,:,i3:N3-L3+i3,:,:);
end
Y2=permute(Y, [3,2,1,6,5,4]);
Z=reshape(Y2, K1*K2*K3, L1*L2*L3);
end

function Z=IMDFS2(X, N1, N2, K1, K2)
L1=N1-K1+1;
L2=N2-K2+1;
for i1=1:L1
    X1(:,:,i1)=X(i1:N1-L1+i1,:);
end
for i2=1:L2
    Y(:,:,:,i2)=X1(:,i2:N2-L2+i2,:);
end
Y2=permute(Y, [2,1,3,4]);
Z=reshape(Y2, K1*K2, L1*L2);
end



function [i1,i2,k1,k2,l1,l2] = findPairs_3d(I,M,N,T)
% T is number of snapshots
if nargin<4
    T = 4;
end
%% find the optimal subarray size for 3D-IMDF
Fold = 0;
for i1 = floor(I/2):I
    
    i2 = I + 1 - i1;
    
    for k1 = floor(M/2):M
        
        k2 = M + 1 - k1;
        
        for l1 = floor(N/2):N
            
            l2 = N + 1 - l1;
            
            F = 0;
            if M == 1
                
                while min(2*l2,l1-1)>=F || min(l1*2,l2-1)>=F
                    F = F + 1;
                end
                
            elseif M == 2
                
                while (8*k2*l2>=F && (i1-1)*k1*l1>=F )
                    F = F + 1;
                end
                
            else
                while (2*T*i2*k2*l2>=F && i1*k1*(l1-1)>=F )
                    F = F + 1;
                end
            end
            if F>Fold
                i1new = i1;
                l1new = l1;
                k1new = k1;
                Fold = F - 1;
            end
        end
    end
end

i1 = i1new;
l1 = l1new;
k1 = k1new;
i2 = I + 1 - i1;
k2 = M + 1 - k1;
l2 = N + 1 - l1;
end


function [k1,k2,l1,l2] = findPairs_2d(M,N)
%% find the optimal k1 and l1

Fold = 0;
for k1 = floor(M/2):M
    k2 = M + 1 - k1;
    for l1 = floor(N/2):N
        l2 = N + 1 - l1;
        F = 0;
        if M == 1
            while min(2*l2,l1-1)>=F || min(l1*2,l2-1)>=F
                F = F + 1;
            end
        elseif M == 2
            while (2*k2*l2>=F && (k1-1)*l1>=F )
                F = F + 1;
            end
        else
            while (8*k2*l2>=F && k1*(l1-1)>=F )
                F = F + 1;
            end
        end
        if F>Fold
            l1new = l1;
            k1new = k1;
            Fold = F - 1;
        end
    end
end
Fold
l1 = l1new;
k1 = k1new;
k2 = M + 1 - k1;
l2 = N + 1 - l1;
end
