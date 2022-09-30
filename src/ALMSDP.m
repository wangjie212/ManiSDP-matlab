function [X, S, y, cx] = ALMSDP(At, b, c, mb, p, tao, y, Y)
C = reshape(c, mb, mb);
sigma = 1;
gama = 10;
MaxIter = 50;
tolgrad = 1e-8;
U = [];
for iter = 1:MaxIter
    [Y, ~, ~] = SDP_ALM_subprog(At, b, c, C, mb, p, sigma, y, Y, U, tolgrad);
    X = Y*Y';
    x = X(:);
    cx = x'*c;
    Axb = At'*x - b;
    neta = norm(Axb);
%     yC = C - reshape(y'*At', mb, mb);
%     E = yC - diag(diag(yC*X));
%     yAt = reshape(Axb'*At', mb, mb);
%     F = sigma*(diag(diag(yAt*X)) - yAt);
    y = y - sigma*Axb;
    yA = reshape(y'*At', mb, mb);
    DfX = C - yA;
    lamda = diag(DfX*X);
    S = DfX - diag(lamda);
    % [v, mineigS] = eigs(S, 1, 'smallestreal');
    [vS, dS] = eig(S, 'vector');
    v = vS(:,1);
    mineigS = dS(1);
    r = 1;
%     by = b'*y + sum(lamda);
%     gap = abs(cx-by)/abs(cx+by);
    [V,D,~] = svd(Y);
    e = diag(D);
    while r < p && e(r+1) > 1e-2*e(1) && e(r+1)/e(r) < 10
        r = r + 1;
    end
%     if r == p - 1
%        [q, ~] = eigs(Y'*Y, 1, 'smallestreal');
%        U = v*q';
%     elseif r < p - 1
%         p = r + 1;
%         Y = V(:,1:p)*D(1:p,1:p);
%         [q, ~] = eigs(Y'*Y, 1, 'smallestreal');
%         U = v*q';
%     else
%          U = [zeros(mb,p) v];
%          Y = [Y zeros(mb,1)]; 
%          p = p + 1;
%    end
%    Y = Y + 0.1*U;
%     for i = 1:mb
%         Y(i,:) = Y(i,:)/norm(Y(i,:));
%     end
%     else
%         U = [];
%     end
   disp(['ALM iteration ' num2str(iter) ': fval = ' num2str(cx,10) ', rank X = ' num2str(r) ', mineigS = ' num2str(mineigS) ', eta = ' num2str(neta) ', p = ' num2str(p) ', sigma = ' num2str(sigma)]);
%    if max([eta,abs(mineigS)]) < 1e-6
    if iter == 1 || neta > 0.5*eta
%       break;
%    else
%        if sigma*eta < 0.1
      sigma = min(sigma*gama, 1e4);
    end
    eta = neta;
end
end