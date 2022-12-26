function [Xopt, S, y, fval] = ALMSDPAbcn(At, b, c, n)
C = reshape(c, n, n);
A = At';
p = 2;
sigma = 1e-3;
gama = 2;
MaxIter = 300;
tolgrad = 1e-8;
tao = 1e-6;
y = zeros(length(b),1);
Y = [];
opts.verbosity = 0;     % Set to 0 for no output, 2 for normal output
opts.maxinner = 20;     % maximum Hessian calls per iteration
opts.mininner = 5;
opts.tolgradnorm = tolgrad; % tolerance on gradient norm
opts.maxiter = 4;
normb = 1+norm(b);
timespend = tic;
for iter = 1:MaxIter
    Y = SDP_ALM_subprog(Y);
    X = Y*Y';
    x = X(:);
    fval = c'*x;
    Axb = At'*x - b;
    neta = norm(Axb)/normb;
    y = y - sigma*Axb;
    yA = reshape(A'*y, n, n);
    DfX = C - yA;
    %lamda = diag(DfX*X); % 避免求其他乘积，从而加速
    lamda = sum(reshape(x.*DfX(:), n, n)); % lamda = diag(DfX*X);
    S = DfX - diag(lamda);
    %[vS, dS] = eig(S,'vector');
    [vS, dS] = mineig(S);
    v = vS(:,1);
    mineigS = dS(1);
    by = b'*y + sum(lamda);
    gap = abs(fval-by)/(abs(by)+abs(fval)+1);
    [V, D, UY] = svd(Y);
    e = diag(D);
    % r = 1;
    % while r < p && e(r+1) > 1e-3*e(1)
    %   r = r + 1;
    % end
    r = sum(e > 1e-3*e(1)); % r = rank(Y)
    if r == p - 1
        % [vY, ~] = eig(Y'*Y, 'vector');
        % q = vY(:,1);
        q = UY(:,end);
        U = v*q';
    elseif r < p - 1
        p = r + 1;
        Y = V(:,1:p)*diag(e(1:p));
        % [vY, ~] = eig(Y'*Y, 'vector');
        % q = vY(:,1);
        % U = v*q';
        U = [zeros(n,r) v];
    else
        U = [zeros(n,p) v];
        Y = [Y zeros(n,1)];
        p = p + 1;
    end
    Y = Y + 0.1*U;
    % for i = 1:n
    %     Y(i,:) = Y(i,:)/norm(Y(i,:));
    % end
    Y = Y';
    nrms = sqrt(sum(Y.^2, 1));
    Y = bsxfun(@times, Y, 1./nrms);
    Y = Y';
    %disp(['ALM iter ' num2str(iter) ': fval = ' num2str(fval,10) ', rank X = ' num2str(r) ', mS = ' num2str(mS) ', eta = ' num2str(neta) ', p = ' num2str(p) ', sigma = ' num2str(sigma)]);
    fprintf('Iter:%d, fval:%0.8f, gap:%0.1e, mineigS:%0.1e, pinf:%0.1e, r:%d, p:%d, sigam:%0.3f, time:%0.2fs\n', ...
             iter,    fval,        gap,       mineigS,       neta,       r,    p,    sigma,   toc(timespend));
    if max(neta, abs(mineigS)) < tao
        break;
    end
    if iter == 1 || neta > 0.7*eta
        % sigma = min(sigma*gama, 1);
        if sigma*gama > 1
            sigma = 1e-3;
        else
            sigma = sigma*gama;
        end
    end
    eta = neta;
end
Xopt = Y*Y';

    function  Y = SDP_ALM_subprog(Y0)
        problem.M = obliquefactory(p, n, true);
        problem.costgrad = @costgrad;
        problem.hess = @hess;
        Y = trustregions(problem, Y0, opts);

        function [f, G, store] = costgrad(Y, store)
            X0 = Y*Y';
            x0 = X0(:);
            Axb0 = At'*x0 - b - y/sigma;
            f = c'*x0 + 0.5*sigma*(Axb0'*Axb0);
            yA0 = reshape(A'*Axb0, n, n);
            S0 = C + sigma*yA0;
            store.S = S0;
            store.G = 2*S0*Y;
            G = problem.M.egrad2rgrad(Y, store.G);
        end

        % If you want to, you can specify the Riemannian Hessian as well.
        function [He, store] = hess(Y, Ydot, store)
            Xdot = Y*Ydot';
            xdot = Xdot(:);
            AxbdotA = A'*(At'*xdot);
            yAdot = reshape(AxbdotA, n, n);
            H = 2*store.S*Ydot + 4*sigma*(yAdot*Y);
            He = problem.M.ehess2rhess(Y, store.G, H, Ydot);
        end

    end

end
