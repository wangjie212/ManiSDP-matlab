function [Xopt, S, y, fval] = ALMSDPAbcK(A, b, c, K)
n = K.s;
C = reshape(c, n, n);
At = A';
p = 2;
sigma = 1e-3;
gama = 2;
MaxIter = 300;
tolgrad = 1e-8;
tao = 1e-6;
y = zeros(length(b),1);
Y = [];
opts.verbosity = 0;      % Set to 0 for no output, 2 for normal output
opts.maxinner = 20;     % maximum Hessian calls per iteration
opts.mininner = 5;
opts.tolgradnorm = tolgrad; % tolerance on gradient norm
opts.maxiter = 4;
normb = norm(b);
for iter = 1:MaxIter
    [Y, ~, ~] = SDP_ALM_subprog(Y);
    X = Y*Y';
    x = X(:);
    fval = x'*c;
    Axb = (x'*At)' - b;
    neta = norm(Axb)/(1+normb);
    y = y - sigma*Axb;
    yA = reshape(y'*A, n, n);
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
        U = [zeros(n,p) v];
        Y = [Y zeros(n,1)];
        p = p + 1;
    end
    Y = Y + 0.1*U;
    for i = 1:n
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
Xopt = Y*Y';

    function  [Y, fval, info] = SDP_ALM_subprog(Y0)
        problem.M = obliquefactory(p, n, true);
        problem.costgrad = @costgrad;
        function [f, G, store] = costgrad(Y, store)
            X0 = Y*Y';
            x0 = X0(:);
            cx = x0'*c;
            Axb0 = (x0'*At)' - b - y/sigma;
            f = cx + sigma/2*(Axb0'*Axb0);
            AxbA = Axb0'*A;
            yA0 = reshape(AxbA, n, n);
            S0 = C + sigma*yA0;
            store.S = S0;
            store.G = 2*S0*Y;
            G = problem.M.egrad2rgrad(Y, store.G);
        end

        % If you want to, you can specify the Riemannian Hessian as well.
        problem.hess = @hess;
        function [He, store] = hess(Y, Ydot, store)
            H = 2*store.S*Ydot;
            Xdot = Y*Ydot'+ Ydot*Y';
            xdot = Xdot(:);
            AxbdotA = (xdot'*At)*A;
            yAdot = reshape(AxbdotA, n, n);
            H = H + 2*sigma*(yAdot*Y);
            He = problem.M.ehess2rhess(Y, store.G, H, Ydot);
        end

        [Y, fval, info] = trustregions(problem, Y0, opts);
    end
end