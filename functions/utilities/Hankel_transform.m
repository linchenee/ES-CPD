function [Y] = Hankel_transform(X)
% Matrix Hankel transform
% --------------------------------------------------------
% Input:
%  X:     matrix of size n1xn2
% --------------------------------------------------------
% Output:
%  Y:     tensor of size n3xn4xn2
% ---------------------------------------------
% Written by Lin Chen (linchenee@sjtu.edu.cn)
% -------------------------------------------------------

[n1,n2] = size(X);
n3 = (n1+mod(n1,2))/2;
n4 = n1+1-n3;
Y = zeros(n3,n4,n2);
for i = 1:n2
    Y(:,:,i) = hankel(X(1:n3,i),X(n3:n1,i));
end

end

%% Equivalent result
% function [Y] = Hankel_transform(X)
% [n1,n2] = size(X);
% if mod(n1,2) == 0
%     n3 = n1/2;
% else
%     n3 = (n1+1)/2;
% end
% n4 = n1+1-n3;
% Y = zeros(n3,n4,n2);
% for k = 1:n2
%     for i = 0:n3-1
%         for j = 0:n4-1
%             Y(i+1,j+1,k) = X(i+j+1,k); 
%         end
%     end
% end
% end
