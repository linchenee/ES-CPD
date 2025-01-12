function C = kr(A1,A2,varargin)
% Kathri Rao product of N matrices
% INPUTS:
% - Matrices A1 (I1xR) and A2 (I2xR)
% OPTIONAL INPUTS:
% - varargin is a list of matrices A3,A4,A5,..., all of them having R columns
% OUTPUT:
% - Matrix C is the Khatri-product of A1 and A2 (or A1,A2,A3,A4,A5) 
N=size(varargin,2);
I1=size(A1,1);
I2=size(A2,1);
R1=size(A1,2);
R2=size(A2,2);
if R1~=R2
    error('Input matrices must have the same number of columns')
end
% perform Khathri-Rao product of A1 and A2
C=zeros(I1*I2,R1);
for r=1:R1
    C(:,r)=reshape(A2(:,r)*A1(:,r).',I1*I2,1);
end
% Do a recursive call to compute the product with the remaining matrices
for p=1:N
    C=kr(C,varargin{p});
end
end   
