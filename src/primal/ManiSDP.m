% This function solves general linear SDPs.
%  Min  <C, X>
%  s.t. A(X) = b
%       X >= 0.

function [X, obj, data] = ManiSDP(At, b, c, n, options)

if ~isfield(options,'p0'); options.p0 = 1; end
if ~isfield(options,'AL_maxiter'); options.AL_maxiter = 1000; end
if ~isfield(options,'gama'); options.gama = 2; end
if ~isfield(options,'sigma0'); options.sigma0 = 1e-2; end
if ~isfield(options,'sigma_min'); options.sigma_min = 1e-1; end
if ~isfield(options,'sigma_max'); options.sigma_max = 1e7; end
if ~isfield(options,'tol'); options.tol = 1e-8; end
if ~isfield(options,'theta'); options.theta = 1e-1; end
if ~isfield(options,'delta'); options.delta = 8; end
if ~isfield(options,'alpha'); options.alpha = 0.1; end
if ~isfield(options,'tolgradnorm'); options.tolgrad = 1e-8; end
if ~isfield(options,'TR_maxinner'); options.TR_maxinner = 30; end
if ~isfield(options,'TR_maxiter'); options.TR_maxiter = 4; end
if ~isfield(options,'tao'); options.tao = 0.25; end
if ~isfield(options,'line_search'); options.line_search = 0; end
if ~isfield(options,'solver'); options.solver = 0; end

fprintf('ManiSDP is starting...\n');
fprintf('SDP size: n = %i, m = %i\n', n, size(b,1));

A = At';
p = options.p0;
sigma = options.sigma0;
gama = options.gama;
y = zeros(length(b),1);
normb = 1 + norm(b);
if isfield(options, 'Y0') 
    Y = options.Y0; 
else
    Y = [];
end
U = [];
% fac_size = [];
% seta = [];
problem.cost = @cost;
problem.grad = @grad;
problem.hess = @hess;
opts.verbosity = 0;     % Set to 0 for no output, 2 for normal output
opts.maxinner = options.TR_maxinner;     % maximum Hessian calls per iteration
opts.maxiter = options.TR_maxiter;
opts.tolgradnorm = options.tolgrad;

data.status = 0;
timespend = tic;
for iter = 1:options.AL_maxiter
%     fac_size = [fac_size; p];
    problem.M = euclideanfactory(n, p);
    if ~isempty(U)
        Y = line_search(Y, U);
    end
    [Y, ~, info] = trustregions(problem, Y, opts);
    gradnorm = info(end).gradnorm;
    X = Y*Y';
    x = X(:);
    obj = c'*x;
    Axb = A*x - b;
    pinf = norm(Axb)/normb;
    y = y - sigma*Axb;
    S = reshape(c - At*y, n, n);
    S = 0.5*(S+S');
    [vS, dS] = eig(S, 'vector');
    dinf = max(0, -dS(1))/(1+dS(end));
    by = b'*y;
    gap = abs(obj-by)/(abs(by)+abs(obj)+1);
    [V, D, ~] = svd(Y);
    if size(D, 2) > 1
        e = diag(D);
    else
        e = D(1);
    end
    r = sum(e >= options.theta*e(1));
    fprintf('Iter %d, obj:%0.8f, gap:%0.1e, pinf:%0.1e, dinf:%0.1e, gradnorm:%0.1e, r:%d, p:%d, sigma:%0.3f, time:%0.2fs\n', ...
             iter,    obj,       gap,       pinf,       dinf,   gradnorm,    r,    p,    sigma,   toc(timespend));
    eta = max([pinf, gap, dinf]);
%     seta = [seta; eta];
    if eta < options.tol
        fprintf('Optimality is reached!\n');
        break;
    end
    if mod(iter, 10) == 0
        if iter > 20 && gap > gap0 && pinf > pinf0 && dinf > dinf0
            data.status = 2;
            fprintf('Slow progress!\n');
            break;
        else
            gap0 = gap;
            pinf0 = pinf;
            dinf0 = dinf;
        end
    end
    if r <= p - 1         
        Y = V(:,1:r)*diag(e(1:r));
        p = r;
    end
    nne = min(sum(dS < 0), options.delta);
    if options.line_search == 1
        U = [zeros(n, p) vS(:,1:nne)];
    end
    p = p + nne;
    if options.line_search == 1
        Y = [Y zeros(n, nne)];
    else
        Y = [Y options.alpha*vS(:,1:nne)];
    end
    if pinf < options.tao*gradnorm
          sigma = max(sigma/gama, options.sigma_min);
    else
          sigma = min(sigma*gama, options.sigma_max);
    end
%    tolgrad = pinf;
end
data.S = S;
data.y = y;
data.gap = gap;
data.pinf = pinf;
data.dinf = dinf;
data.gradnorm = gradnorm;
data.time = toc(timespend);
% data.fac_size = fac_size;
% data.seta = seta;
if data.status == 0 && eta > options.tol
    data.status = 1;
    fprintf('Iteration maximum is reached!\n');
end

fprintf('ManiSDP: optimum = %0.8f, time = %0.2fs\n', obj, toc(timespend));

    function nY = line_search(Y, U)
         alpha = 0.2;
         cost0 = co(Y);
         i = 1;
         nY = Y + alpha*U;
         while i <= 15 && co(nY) - cost0 > -1e-3
              alpha = 0.8*alpha;
              nY = Y + alpha*U;
              i = i + 1;
         end
    end

   function val = co(Y)
        X = Y*Y';
        x = X(:);
        Axb = A*x - b - y/sigma;
        val = c'*x + sigma/2*(Axb'*Axb);
   end

   function [f, store] = cost(Y, store)
        X = Y*Y';
        x = X(:);
        Axb = A*x - b - y/sigma;
        f = c'*x + sigma/2*(Axb'*Axb);
    end
    
    function [G, store] = grad(Y, store)
        S = reshape(c+sigma*At*Axb, n, n);
        G = 2*S*Y;
    end

    function [H, store] = hess(Y, U, store)
        YU = U*Y';
        AyU = reshape(A'*(At'*YU(:)), n, n);
        H = 2*S*U + 4*sigma*(AyU*Y);
    end
end