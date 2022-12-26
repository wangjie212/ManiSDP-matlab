function [Y, S, y, fval, error] = ALMSDP0(At, b, c, n)
C = reshape(c, n, n);
A = At';
p = 1;
sigma = 1e-1;
sigma_min = 1e0; 
sigma_max = 1e7;
gama = 2;
MaxIter = 300;
tolgrad = 1e0;
tao = 1e-8;
y = zeros(length(b),1);
normb = 1 + norm(b);
Y = [];
U = [];

timespend = tic;
for iter = 1:MaxIter
    [Y,dinf] = SDP_ALM_subprog0(At, A, b, c, C, n, p, sigma, y, Y, U, tolgrad);
    X = Y*Y';
    x = X(:);
    fval = c'*x;
    Axb = A*x - b;
    pinf = norm(Axb)/normb;
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
    r = sum(e > 1e-3*e(1));
    fprintf('Iter:%d, fval:%0.8f, gap:%0.1e, mineigS:%0.1e, pinf:%0.1e, dinf:%0.1e, r:%d, p:%d, sigam:%0.3f, time:%0.2fs\n', ...
             iter,    fval,       gap,       mS,       pinf,   dinf,    r,    p,    sigma,   toc(timespend));
    error = max([pinf, gap, mS]);
    if error < tao
        break;
    else
        if r <= p - 1         
            Y = V(:,1:r)*diag(e(1:r));
            p = r;
        end
        nne = min(sum(dS < 0), 6);
        U = [zeros(n, p) vS(:,1:nne)];
        p = p + nne;
        Y = [Y zeros(n, nne)];
        % Y = [Y 0.1*vS(:,1:nne)];
    end
%     neta = pinf;
%     if iter == 1 || neta > 0.5*eta
%         if sigma > 1e1
%             sigma = 1e-1;
%         else
%             sigma = sigma*gama;
%         end
%     end
%     eta = neta;
    
    if pinf < dinf/1e2
          sigma = max(sigma/gama, sigma_min);
    else
          sigma = min(sigma*gama, sigma_max);
    end
    tolgrad = 1e1*pinf;
end
end