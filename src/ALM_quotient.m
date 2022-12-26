function [Y, S, y, fval, error] = ALM_quotient(At, b, c, n)
C = reshape(c, n, n);
A = At';
p = 1;
% sigma = 1e-1;
sigma = 1e1;
gama = 2;
MaxIter = 300;
tolgrad = 1e-8;
tao = 1e-8;
y = zeros(length(b),1);
normb = 1 + norm(b);
Y = [];
U = [];

timespend = tic;
for iter = 1:MaxIter
    Y = quotient(At, A, b, c, C, n, p, sigma, y, Y, U, tolgrad);
    X = Y*Y';
    x = X(:);
    fval = c'*x;
    Axb = A*x - b;
    neta = norm(Axb)/normb;
    y = y - sigma*Axb;
    S = C - reshape(y'*A, n, n);
    [vS, dS] = eig(S, 'vector');
    mS = abs(min(dS))/(1+dS(end));
    by = b'*y;
    gap = abs(fval-by)/(abs(by)+abs(fval)+1);
    [V, D, ~] = svd(Y);
    if size(D, 2) > 1
        e = diag(D);
    else
        e = D(1);
    end
%    r = sum(e > 1e-3*e(1));
    r = sum(e > 1e-2*e(1));
    if r <= p - 1         
        Y = V(:,1:r)*diag(e(1:r));
        p = r;
    end
    nne = min(sum(dS < 0), 8);
    U = [zeros(n, p) vS(:,1:nne)];
    p = p + nne;
    Y = [Y zeros(n, nne)];
    % Y = [Y 0.2*vS(:,1:nne)];
    fprintf('Iter:%d, fval:%0.8f, gap:%0.1e, mineigS:%0.1e, pinf:%0.1e, r:%d, p:%d, sigam:%0.3f, time:%0.2fs\n', ...
             iter,    fval,       gap,       mS,       neta,       r,    p,    sigma,   toc(timespend));
    error = max([neta, gap, mS]);
    if error < tao
        break;
    end
    if iter == 1 || neta > 0.7*eta
        if sigma < 1e3
              sigma = gama*sigma;
        else
%              sigma = 1e1;
              sigma = 1e2;
        end
    end
    eta = neta;
end
end