function [X, S, y, cx] = nALMSDP(At, b, c, mb, tao)
y = zeros(length(b),1);
Y =[];
p = 2;
[Y, S, y, cx] = ALMSDP(At, b, c, mb, p, 1e-4, 13, y, Y);
% [V,D,~] = svd(Y);
% r = 1;
% [vS, ~] = eig(S, 'vector');
% v = vS(:,1);
% e = diag(D);
% while r < p && e(r+1) > 1e-3*e(1) && e(r+1)/e(r) < 10
%     r = r + 1;
% end
% if r == p - 1
%    [q, ~] = eigs(Y'*Y, 1, 'smallestreal');
%    U = v*q';
% elseif r < p - 1
%     p = r + 1;
%     Y = V(:,1:p)*D(1:p,1:p);
%     [q, ~] = eigs(Y'*Y, 1, 'smallestreal');
%     U = v*q';
% else
%     U = [zeros(mb,p) v];
%     Y = [Y zeros(mb,1)]; 
%     p = p + 1;
% end
for k = 1:10
Y = Y + 0.1*rand(mb,p);
for i = 1:mb
   Y(i,:) = Y(i,:)/norm(Y(i,:));
end
% y = zeros(length(b),1);
[Y, S, y, cx] = ALMSDP(At, b, c, mb, p, 1e-6, 13, y, Y);
end
X = Y*Y';
end