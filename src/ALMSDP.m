function [Y, S, y, fval] = ALMSDP(At, b, c, mb)
C = reshape(c, mb, mb);
p = 2;
sigma = 1e-0;
gama = 1;
MaxIter = 300;
tolgrad = 1e-6;
tao = 1e-6;
y = zeros(length(b),1);
Y = [];
U = [];
for iter = 1:MaxIter
    [Y, ~, ~] = SDP_ALM_subprog(At, b, c, C, mb, p, sigma, y, Y, U, tolgrad);
    X = Y*Y';
    x = X(:);
    fval = x'*c;
    Axb = At'*x - b;
    neta = norm(Axb)/(1+norm(b));
    y = y - sigma*Axb;
    yA = reshape(y'*At', mb, mb);
    DfX = C - yA;
    lamda = diag(DfX*X);
    S = DfX - diag(lamda);
    [vS, dS] = eig(S, 'vector');
    v = vS(:,1);
    mS = norm(min(dS,0))/(1+norm(S));
    r = 1;
%     by = b'*y + sum(lamda);
%     gap = abs(cx-by)/abs(cx+by);
    [V,D,~] = svd(Y);
    e = diag(D);
    while r < p && e(r+1) > 1e-3*e(1)
        r = r + 1;
    end
    if r == p - 1
       [vY, ~] = eig(Y'*Y, 'vector');
       q = vY(:,1);
       U = v*q';
    elseif r < p - 1
        p = r + 1;
        Y = V(:,1:p)*D(1:p,1:p);
        [vY, ~] = eig(Y'*Y, 'vector');
        q = vY(:,1);
        U = v*q';
    else
         U = [zeros(mb,p) v];
         Y = [Y zeros(mb,1)]; 
         p = p + 1;
    end
    Y = Y + 0.1*U;
    for i = 1:mb
         Y(i,:) = Y(i,:)/norm(Y(i,:));
    end
    disp(['ALM iter ' num2str(iter) ': fval = ' num2str(fval,10) ', rank X = ' num2str(r) ', mS = ' num2str(mS) ', eta = ' num2str(neta) ', p = ' num2str(p) ', sigma = ' num2str(sigma)]);
    if max(neta, mS) < tao
        break;
    end
   if iter == 1 || neta > 0.8*eta
        % sigma = min(sigma*gama, 1);
        if sigma*gama > 1
            sigma = 1e-3;
        else
            sigma = sigma*gama;
        end
   end
    eta = neta;
end
end