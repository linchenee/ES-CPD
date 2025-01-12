function c = khatri(a,b)
K = size(a,2);
for i = 1:K
    c(:,i) = kron(a(:,i), b(:,i));
end
end
