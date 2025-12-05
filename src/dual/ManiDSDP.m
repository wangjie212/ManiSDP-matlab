% This function solves SDPs using the dual approach:
%    sup  <C, X> + <c, w>
%    s.t. A(X) + B(w) = b
%         X in S_+^{n}
%         w in R^l

function [X, obj, data] = ManiDSDP(A, b, c, K, options)

% if ~isfield(options,'p0'); options.p0 = ceil(log(length(b))); end
if ~isfield(options,'p0'); options.p0 = 1; end
if ~isfield(options,'ADMM_maxiter'); options.ADMM_maxiter = 1000; end
if ~isfield(options,'gama'); options.gama = 2; end
if ~isfield(options,'sigma0'); options.sigma0 = 1e-1; end
if ~isfield(options,'sigma_min'); options.sigma_min = 1e-2; end
if ~isfield(options,'sigma_max'); options.sigma_max = 1e7; end
if ~isfield(options,'tol'); options.tol = 1e-8; end
if ~isfield(options,'theta'); options.theta = 1e-2; end
if ~isfield(options,'delta'); options.delta = 8; end
if ~isfield(options,'alpha'); options.alpha = 0.01; end
if ~isfield(options,'tolgradnorm'); options.tolgradnorm = 1e-8; end
if ~isfield(options,'TR_maxinner'); options.TR_maxinner = 20; end
if ~isfield(options,'TR_maxiter'); options.TR_maxiter = 4; end
if ~isfield(options,'tau1'); options.tau1 = 0.1; end
if ~isfield(options,'tau2'); options.tau2 = 1; end
if ~isfield(options,'line_search'); options.line_search = 1; end

n = K.s;
fprintf('ManiSDP is starting...\n');
fprintf('SDP size: n = %i, m = %i\n', n, size(b,1));
warning('off', 'manopt:trs_tCG_cached:memory');
normc = 1 + norm(c);
B = A(:, 1:K.f);
A = A(:, K.f+1:end);
cf = c(1:K.f);
c = c(K.f+1:end);

if ~isfield(options,'dAAt'); options.dAAt = diag(A*A'); end
iA = (sparse(1:size(b,1),1:size(b,1),options.dAAt)\A)'; % (inv(A*A')*A)'
bA = iA*b;
iAB = iA*B;

p = options.p0;
sigma = options.sigma0; 
gama = options.gama;
x = zeros(n^2, 1);
w = zeros(K.f, 1);
Y = [];
U = [];

problem.costgrad = @costgrad;
problem.hess = @hess;
opts.verbosity = 0;     % Set to 0 for no output, 2 for normal output
opts.maxinner  = options.TR_maxinner;    % maximum Hessian calls per iteration
opts.maxiter   = options.TR_maxiter;
opts.tolgradnorm = options.tolgradnorm;
% opts.trscache = false;
data.status = 0;
timespend = tic;
for iter = 1:options.ADMM_maxiter
    problem.M = euclideanfactory(n, p);
    if ~isempty(U)
        Y = line_search(Y, U);
    end
    [Y, ~, info] = trustregions(problem, Y, opts);
    gradnorm = info(end).gradnorm;
    S = Y*Y';
    sc = S(:) - c;
    y = iA'*sc;
    As = A'*y - sc;
    Af = B'*y - cf;
    pinf = (norm(As)+norm(Af))/normc;
    by = b'*y;
    x = x + sigma*(iAB*(Af - w/sigma) + A'*(iA'*(As - x/sigma)) - As);
    w = w - sigma*Af;
    X = reshape(x + bA, n, n);
    [vX, dX] = eig(X, 'vector');
    obj = c'*(x + bA) + cf'*w;
    dinf = max(0, -dX(1))/(1+abs(dX(end)));
    gap = abs(obj-by)/(1+abs(obj)+abs(by));
    [V, D, ~] = svd(Y);
    if size(D, 2) > 1
        e = diag(D);
    else
        e = D(1);
    end
    r = sum(e > options.theta*e(1));
    fprintf('Iter %d, obj:%0.8f, gap:%0.1e, pinf:%0.1e, dinf:%0.1e, gradnorm:%0.1e, r:%d, p:%d, sigma:%0.3f, time:%0.2fs\n', ...
             iter,    obj,       gap,       pinf,       dinf,       gradnorm,       r,    p,    sigma, toc(timespend));
    eta = max([gap, pinf, dinf]);
    if eta < options.tol
        fprintf('Optimality is reached!\n');
        break;
    end
    if mod(iter, 20) == 0
        if iter > 50 && gap > gap0 && pinf > pinf0 && dinf > dinf0
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
    nne = min(sum(dX < 0), options.delta);
    if options.line_search == 1
        U = [zeros(n, p) vX(:,1:nne)];
    end
    p = p + nne;
    if options.line_search == 1
        Y = [Y zeros(n, nne)];
    else
        Y = [Y options.alpha*vX(:,1:nne)];
    end
    if pinf < options.tau1*gradnorm
          sigma = max(sigma/gama, options.sigma_min);
    elseif pinf > options.tau2*gradnorm
          sigma = min(sigma*gama, options.sigma_max);
    end
end
data.X = X;
data.y = y;
data.S = S;
data.w = w;
data.gap = gap;
data.pinf = pinf;
data.dinf = dinf;
data.gradnorm = gradnorm;
data.time = toc(timespend);
if data.status == 0 && eta > options.tol
    data.status = 1;
    fprintf('Iteration maximum is reached!\n');
end

fprintf('ManiDSDP: optimum = %0.8f, time = %0.2fs\n', obj, toc(timespend));
        
    function val = co(Y)
        S = Y*Y';
        sc = S(:) - c;
        y = iA'*sc;
        As = A'*y - sc - x/sigma;
        Af = B'*y - cf - w/sigma;
        val = b'*y + 0.5*sigma*(As'*As + Af'*Af);
    end

    function nY = line_search(Y, U)
         alpha = 1;
         cost0 = co(Y);
         i = 1;
         nY = Y + alpha*U;
         while i <= 15 && co(nY) - cost0 > -1e-3
              alpha = 0.8*alpha;
              nY = Y + alpha*U;
              i = i + 1;
         end
    end

    function [f, G, store] = costgrad(Y, store)
        S = Y*Y';
        sc = S(:) - c;
        y = iA'*sc;
        As = A'*y - sc - x/sigma;
        Af = B'*y - cf - w/sigma;
        f = b'*y + 0.5*sigma*(As'*As + Af'*Af);
        X = reshape(bA + sigma*(iAB*Af + A'*(iA'*As) - As), n, n);
        G = 2*X*Y;
    end

    function [H, store] = hess(Y, U, store)      
        YU = U*Y';
        yAU = reshape((YU(:)'*iA)*A, n, n);
        H = 2*X*U + 2*sigma*(U*(Y'*Y) + Y*(U'*Y)) + 4*sigma*(reshape((YU(:)'*iAB)*iAB' + (yAU(:)'*iA)*A, n, n) - 2*yAU)*Y;
    end
end