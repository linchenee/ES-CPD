% C. Qian, L. Huang, H. C. So, N. D. Sidiropoulos and J. Xie, "Unitary PUMA 
% algorithm for estimating the frequency of a complex sinusoid," IEEE Transactions 
% on Signal Processing, vol. 63, no. 20, pp. 5358-5368, Oct. 2015.

function f = PUMA(u)

N = length(u);

J1 = eye(N-1,N);
J2 = [zeros(N-1,1),eye(N-1)];

u1 = J1*u;
u2 = J2*u;

f = pinv(u1)*u2;

% w = tmp(u,round(2*N/3),f);
% f = exp(1j*w);

it = 0;
while it < 3
    A = J2 - J1*f;
    W = A*A';
    f = (u1'/W*u2)/(u1'/W*u1);
    it = it + 1;   
end

doa = asin(angle(f)/pi)*180/pi;
f = angle(f);

end


function [w,DOA] = tmp(us,m,z2)
[M,K] = size(us);
z = angle( eig(us(1:m,:)\us(M-m+1:end,:)) );
w = zeros(K,1);
l = floor((M-m-1)/2);
dt = [];
for i = 1:K
    d = ((-l:l) * 2*pi + z(i))/(M-m);
    dt = [dt,d];
end
for i = 1:K
    [~,I] = min(abs(dt - z2(i)));
    w(i) = dt(I);
end
DOA = sort( asin( w(:)/pi ) * 180/pi );
end
