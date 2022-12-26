function [Y, S, y, fval] = ALMSDP(At, b, c, mb)
C = reshape(c, mb, mb);
p = 2;
sigma = 1e-3;
gama = 2;
MaxIter = 300;
tolgrad = 1e-6;
tao = 1e-6;
y = zeros(length(b),1);
normb = 1+norm(b);
Y = [];
U = [];
for iter = 1:MaxIter
    [Y, ~, ~] = SDP_ALM_subprog(At, b, c, C, mb, p, sigma, y, Y, U, tolgrad);
    X = Y*Y';
    x = X(:);
    fval = x'*c;
    Axb = At'*x - b;
    neta = norm(Axb)/normb;
    y = y - sigma*Axb;
    yA = reshape(y'*At', mb, mb);
    eS = C - yA;
    lamda = sum(reshape(x.*eS(:), mb, mb));
    S = eS - diag(lamda);
    [vS, dS] = eig(S, 'vector');
    mS = min(dS)/(1+norm(S));
    sy = norm(eS*Y);
%     by = b'*y + sum(lamda);
%     gap = abs(cx-by)/abs(cx+by);
    [V,D,~] = svd(Y);
    e = diag(D);
    r = sum(e > 1e-3*e(1));
    if r <= p - 1         
        Y = V(:,1:r)*diag(e(1:r));
        p = r;
    end
    nne = max(min(sum(dS < 0), 8), 1);
    p = p + nne;
    Y = [Y 0.1*vS(:,1:nne)];
    for i = 1:mb
        Y(i,:) = Y(i,:)/norm(Y(i,:));
    end
    disp(['ALM iter ' num2str(iter) ': fval = ' num2str(fval,10) ', rank X = ' num2str(r) ', mS = ' num2str(mS) ', eta = ' num2str(neta) ', p = ' num2str(p) ', sigma = ' num2str(sigma)]);
    if max(neta, abs(mS)) < tao
        break;
    end
   if iter == 1 || neta > 0.5*eta
%         if gama*sigma > 1
%             sigma = 1e-3;
%         else
          if sigma < 1 || sy < 2
              sigma = gama*sigma;
          elseif sy >= 2
              sigma = 1e-3;
          end
   end
    eta = neta;
end
end