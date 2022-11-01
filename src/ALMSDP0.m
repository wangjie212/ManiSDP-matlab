function [Y, S, y, fval] = ALMSDP0(At, b, c, mb)
C = reshape(c, mb, mb);
p = 2;
sigma = 1;
gama = 2;
MaxIter = 300;
tolgrad = 1e-6;
tao = 1e-6;
y = zeros(length(b),1);
normb = 1+norm(b);
Y = [];
U = [];
for iter = 1:MaxIter
    [Y, ~, ~] = SDP_ALM_subprog0(At, b, c, C, mb, p, sigma, y, Y, U, tolgrad);
    X = Y*Y';
    x = X(:);
    fval = x'*c;
    Axb = At'*x - b;
    neta = norm(Axb)/normb;
    y = y - sigma*Axb;
    yA = reshape(y'*At', mb, mb);
    S = C - yA;
    [vS, dS] = eig(S, 'vector');
    mS = min(dS)/(1+dS(end));
%     by = b'*y + sum(lamda);
%     gap = abs(cx-by)/abs(cx+by);
    [V,D,~] = svd(Y);
    if size(D, 2) > 1
        e = diag(D);
    else
        e = D(1);
    end
    r = sum(e > 1e-1*e(1));
    if r <= p - 1         
        Y = V(:,1:r)*diag(e(1:r));
        p = r;
    end
    nne = min(sum(dS < 0), 4);
    p = p + nne;
    Y = [Y 0.1*vS(:,1:nne)];
    disp(['ALM iter ' num2str(iter) ': fval = ' num2str(fval,10) ', rank X = ' num2str(r) ', mS = ' num2str(mS) ', eta = ' num2str(neta) ', p = ' num2str(p) ', sigma = ' num2str(sigma)]);
    if max(neta, abs(mS)) < tao
        break;
    end
    if iter == 1 || neta > 0.5*eta
        if sigma < 1
              sigma = gama*sigma;
        else
              sigma = 1;
        end
    else
    end
    eta = neta;
end
end