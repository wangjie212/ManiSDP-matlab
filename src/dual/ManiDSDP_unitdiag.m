% This function solves SDPs with dual unit-diagonal (diag(S) = 1) using the dual approach:
%    sup  <C, X> + <c, w>
%    s.t. A(X) + B(w) = b
%         X in S_+^{n}
%         w in R^l
%         diag(S) = 1

function [X, obj, data] = ManiDSDP_unitdiag(A, b, c, K, options)

% if ~isfield(options,'p0'); options.p0 = 2; end
if ~isfield(options,'p0'); options.p0 = ceil(log(length(b))); end
if ~isfield(options,'ADMM_maxiter'); options.ADMM_maxiter = 300; end
if ~isfield(options,'gama'); options.gama = 2; end
if ~isfield(options,'sigma0'); options.sigma0 = 1e-3; end
if ~isfield(options,'sigma_min'); options.sigma_min = 1e-3; end
if ~isfield(options,'sigma_max'); options.sigma_max = 1e7; end
if ~isfield(options,'tol'); options.tol = 1e-8; end
if ~isfield(options,'theta'); options.theta = 1e-3; end
if ~isfield(options,'delta'); options.delta = 8; end
if ~isfield(options,'alpha'); options.alpha = 0.1; end
if ~isfield(options,'tolgradnorm'); options.tolgradnorm = 1e-8; end
if ~isfield(options,'TR_maxinner'); options.TR_maxinner = 20; end
if ~isfield(options,'TR_maxiter'); options.TR_maxiter = 4; end
if ~isfield(options,'tau1'); options.tau1 = 1e1; end
if ~isfield(options,'tau2'); options.tau2 = 1e2; end
if ~isfield(options,'line_search'); options.line_search = 0; end

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

p = options.p0;
sigma = options.sigma0; 
gama = options.gama;
x = zeros(n^2, 1);
w = zeros(K.f, 1);
eG = zeros(p, n);
Y = [];
U = [];
fac_size = [];
seta = [];

problem.cost = @cost;
problem.grad = @grad;
problem.hess = @hess;
opts.verbosity = 0;     % Set to 0 for no output, 2 for normal output
opts.maxinner  = options.TR_maxinner;    % maximum Hessian calls per iteration
opts.maxiter   = options.TR_maxiter;
opts.tolgradnorm = options.tolgradnorm;
% opts.trscache = false;
data.status = 0;
timespend = tic;
for iter = 1:options.ADMM_maxiter
    fac_size = [fac_size; p];
    problem.M = obliquefactoryNTrans(p, n);
    if ~isempty(U)
        Y = line_search(Y, U);
    end
    [Y, ~, info] = trustregions(problem, Y, opts);
    gradnorm = info(end).gradnorm;
    S = Y'*Y;
    sc = S(:) - c;
    y = iA'*sc;
    As = A'*y - sc;
    Af = B'*y - cf;
    pinf = (norm(As)+norm(Af))/normc;
    by = b'*y;
    x = x - sigma*As;
    w = w - sigma*Af;
    eX = reshape(x + bA, n, n);
    z = sum(S.*eX);
    X = eX - diag(z);
    [vX, dX] = eig(X, 'vector');
    % [vX, dX] = eigs(X, options.delta, 'smallestreal');
    % dX = diag(dX);
    obj = c'*eX(:) + cf'*w + sum(z);
    dinf = max(0, -dX(1))/(1+abs(dX(end)));
    gap = abs(obj-by)/(1+abs(obj)+abs(by));
    [~, D, V] = svd(Y);
    e = diag(D);
    r = sum(e > options.theta*e(1));
    fprintf('Iter %d, obj:%0.8f, gap:%0.1e, pinf:%0.1e, dinf:%0.1e, gradnorm:%0.1e, r:%d, p:%d, sigma:%0.3f, time:%0.2fs\n', ...
             iter,    obj,       gap,       pinf,       dinf,       gradnorm,       r,    p,    sigma,   toc(timespend));
    eta = max([gap, pinf, dinf]);
    seta = [seta; eta];
    if eta < options.tol
        fprintf('Optimality is reached!\n');
        break;
    end
    if mod(iter, 50) == 0
        if iter > 100 && gap > gap0 && pinf > pinf0 && dinf > dinf0
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
        Y = V(:,1:r)'.*e(1:r);
        p = r;
    end
    nne = max(min(sum(dX < 0), options.delta), 1);
    if options.line_search == 1
        U = [zeros(p, n); vX(:,1:nne)'];
    end
    p = p + nne;
    if options.line_search == 1
        Y = [Y; zeros(nne, n)];
    else
        Y = [Y; options.alpha*vX(:,1:nne)'];
        Y = Y./sqrt(sum(Y.^2));
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
data.fac_size = fac_size;
data.seta = seta;
if data.status == 0 && eta > options.tol
    data.status = 1;
    fprintf('Iteration maximum is reached!\n');
end

fprintf('ManiDSDP: optimum = %0.8f, time = %0.2fs\n', obj, toc(timespend));
        
    function val = co(Y)
        S = Y'*Y;
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
         nY = nY./sqrt(sum(nY.^2));
         while i <= 15 && co(nY) - cost0 > -1e-3
              alpha = 0.8*alpha;
              nY = Y + alpha*U;
              nY = nY./sqrt(sum(nY.^2));
              i = i + 1;
         end
    end

      function [f, store] = cost(Y, store)
        S = Y'*Y;
        sc = S(:) - c;
        y = iA'*sc;
        As = A'*y - sc - x/sigma;
        Af = B'*y - cf - w/sigma;
        f = b'*y + 0.5*sigma*(As'*As + Af'*Af);
      end

      function [G, store] = grad(Y, store)
        X = reshape(bA - sigma*As, n, n);
        eG = 2*Y*X;
        G = eG - Y.*sum(Y.*eG);
      end

    function [H, store] = hess(Y, U, store)
        YU = Y'*U;
        yAU = reshape(A'*(iA'*YU(:)), n, n);
        eH = 2*U*X - 4*sigma*(Y*yAU) + 2*sigma*((Y*U')*Y+(Y*Y')*U);
        H = eH - Y.*sum(Y.*eH) - U.*sum(Y.*eG);
    end

    function M = obliquefactoryNTrans(n, m)
        M.dim = @() (n-1)*m;
        M.inner = @(x, d1, d2) d1(:)'*d2(:);
        M.norm = @(x, d) norm(d);
        M.typicaldist = @() pi*sqrt(m);
        M.proj = @(X, U) U - X.*sum(X.*U);
        M.tangent = @(X, U) U - X.*sum(X.*U);

        M.retr = @retraction;
        function y = retraction(x, d)            
            xtd = x + d;
            y = xtd./sqrt(sum(xtd.^2));
        end

        M.rand = @() random(n, m);
        M.lincomb = @matrixlincomb;
        M.zerovec = @(x) zeros(n, m);
        M.transp = @(x1, x2, d) d - x2.*sum(x2.*d);
        
        function x = random(n, m)
            x = randn(n, m);
            x = x./sqrt(sum(x.^2));
        end
    end
end
