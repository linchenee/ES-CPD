function [X_mat] = tens2mat(X,mode)
[I,J,K] = size(X);
if mode == 1
    X_mat = reshape(permute(X,[3 1 2]),I*K,J);
elseif mode == 2
    X_mat = reshape(X,J*I,K);
elseif mode == 3
    X_mat = reshape(permute(X,[2 3 1]),K*J,I);
else
    error('Input argument mode must be 1, 2 or 3');
end
end