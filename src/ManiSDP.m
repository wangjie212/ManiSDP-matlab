% This function solves general linear SDPs.
% Min  <C, X>
% s.t. A(X) = b,
%      X >= 0.

function [X, obj, data] = ManiSDP(At, b, c, n, options)

if ~isfield(options,'p0'); options.p0 = 1; end
if ~isfield(options,'AL_maxiter'); options.AL_maxiter = 300; end
if ~isfield(options,'gama'); options.gama = 2; end
if ~isfield(options,'sigma0'); options.sigma0 = 1e-2; end
if ~isfield(options,'sigma_min'); options.sigma_min = 1e-1; end
if ~isfield(options,'sigma_max'); options.sigma_max = 1e7; end
if ~isfield(options,'tol'); options.tol = 1e-8; end
if ~isfield(options,'theta'); options.theta = 1e-1; end
if ~isfield(options,'delta'); options.delta = 8; end
if ~isfield(options,'alpha'); options.alpha = 0.2; end
if ~isfield(options,'tolgradnorm'); options.tolgrad = 1e-8; end
if ~isfield(options,'TR_maxinner'); options.TR_maxinner = 50; end
if ~isfield(options,'TR_maxiter'); options.TR_maxiter = 4; end
if ~isfield(options,'tao'); options.tao = 0.25; end
if ~isfield(options,'line_search'); options.line_search = 1; end
if ~isfield(options,'solver'); options.solver = 0; end

fprintf('ManiSDP is starting...\n');
fprintf('SDP size: n = %i, m = %i\n', n, size(b,1));

C = reshape(c, n, n);
A = At';
p = options.p0;
sigma = options.sigma0;
gama = options.gama;
y = zeros(length(b),1);
normb = 1 + norm(b);
Y = [];
U = [];
fac_size = [];
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
    fac_size = [fac_size; p];
    problem.M = euclideanfactory(n, p);
    if ~isempty(U)
        Y = line_search(Y, U);
    end
    if options.solver == 0
        [Y, ~, info] = trustregions(problem, Y, opts);
    elseif options.solver == 1
        [Y, ~, info] = arc(problem, Y, opts);
    elseif options.solver == 2
        [Y, ~, info] = steepestdescent(problem, Y, opts);
    elseif options.solver == 3
        [Y, ~, info] = conjugategradient(problem, Y, opts);
    elseif options.solver == 4
        [Y, ~, info] = barzilaiborwein(problem, Y, opts);
    elseif options.solver == 5
        [Y, ~, info] = rlbfgs(problem, Y, opts);
    else
        fprintf('Solver is not supported!\n');
        return;
    end
    gradnorm = info(end).gradnorm;
    X = Y*Y';
    x = X(:);
    obj = c'*x;
    Axb = A*x - b;
    pinf = norm(Axb)/normb;
    y = y - sigma*Axb;
    S = C - reshape(y'*A, n, n);
    [vS, dS] = eig(S, 'vector');
    dinf = abs(min(dS))/(1+dS(end));
    by = b'*y;
    gap = abs(obj-by)/(abs(by)+abs(obj)+1);
    [V, D, ~] = svd(Y);
    if size(D, 2) > 1
        e = diag(D);
    else
        e = D(1);
    end
    r = sum(e > options.theta*e(1));
    fprintf('Iter %d, obj:%0.8f, gap:%0.1e, pinf:%0.1e, dinf:%0.1e, gradnorm:%0.1e, r:%d, p:%d, sigma:%0.3f, time:%0.2fs\n', ...
             iter,    obj,       gap,       pinf,       dinf,   gradnorm,    r,    p,    sigma,   toc(timespend));
    eta = max([pinf, gap, dinf]);
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
data.fac_size = fac_size;
if data.status == 0 && eta > options.tol
    data.status = 1;
    fprintf('Iteration maximum is reached!\n');
end

fprintf('ManiSDP: optimum = %0.8f, time = %0.2fs\n', obj, toc(timespend));

%         function Y = line_search(Y, U, t)
%             X = Y*Y';
%             D = U*U';
%             YU = Y*U' + U*Y';
%             q0 = A*X(:) - b;
%             q1 = A*YU(:);
%             q2 = A*D(:);
%             aa = sigma/2*norm(q2)^2;
%             bb = sigma*q1'*q2;
%             cc = c'*D(:) - (y - sigma*q0)'*q2 + sigma/2*norm(q1)^2;
%             dd = c'*YU(:) - (y - sigma*q0)'*q1;
%             alpha_min = 0.02;
%             alpha_max = 0.5;
%             sol = vpasolve(4*aa*t^3 + 3*bb*t^2 + 2*cc*t + dd == 0, t, [alpha_min alpha_max]);
%             alpha = [alpha_min;eval(sol);alpha_max];
%             [~,I] = min(aa*alpha.^4 + bb*alpha.^3 + cc*alpha.^2 + dd*alpha);
%             Y = Y + alpha(I)*U;
%         end

%         function Y = line_search(Y, U)
%              alpha = [0.02;0.04;0.06;0.08;0.1;0.2];
%              val = zeros(length(alpha),1);
%              for i = 1:length(alpha)
%                 val(i) = co(Y + alpha(i)*U);
%              end
%              [~,I] = min(val);
%              Y = Y + alpha(I)*U;
%         end

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
        S = C + sigma*reshape(Axb'*A, n, n);
        G = 2*S*Y;
    end

    function [H, store] = hess(Y, U, store)
        YU = U*Y';
        AyU = reshape(YU(:)'*At*A, n, n);
        H = 2*S*U + 4*sigma*(AyU*Y);
    end
end