function [X, S, y, fval] = ALMSDPNT_ALLFUN(At, b, c, n)
% Solve standard linear SDP problem via Augmented Lagrangian Method based
% on Manifold Optimization
%
% Author: WangJie, HuLiangbing
% Data: 2022-07-29
%
% Change log:
%   HLB 2022-10-24:
%       problem using cost, grad and hess function, while ALMSDPNT_EIGV2
%       using costgrad and hess function.
%
%   HLB 2022-10-26:
%       change neta to sqrt(Axb'*Axb), not divided by normb = 1+norm(b).
%
%   HLB 2022-10-29:
%       1) lamda can be calculated by lamda = sum(Y0.*(Y0*DfX)); but suit for
%       small problems;
%       2) change "isempty(d)" to  "isempty(d) && neta<tao".
%       3) update sigma based on pinf and dinf(==trustregion.info(end).gradnorm).
%
%   HLB 2022-11-01:
%       put all trustregions call functions into one file.


C = reshape(c, n, n);
A = At';
p = 2;
sigma = 1e-3; sigma_min=1e-3; sigma_max = 1e3;
gamma = 2; eta1 = 1-10*eps; eta2 = 1+10*eps; h4 = 1;
it_pinf = 0; it_dinf = 0;
MaxIter = 30000;
tolgrad = 1e-8;
tao = 1e-6;
egrad = zeros(p,n);
y = zeros(length(b),1);
Axby = y;
normb = 1+norm(b);
lamda = zeros(1,n);
gap = 1e3;
r = p - 1;
Y = [];

problem.cost = @cost;
problem.grad = @grad;
problem.hess = @hess;
% opts.subproblemsolver = @trs_gep; % Call your favorite solver.
opts.verbosity = 0;     % Set to 0 for no output, 2 for normal output
opts.maxinner  = 20;    % maximum Hessian calls per iteration
opts.mininner  = 5;
opts.maxiter   = 4;
opts.tolgradnorm = tolgrad; % tolerance on gradient norm

timespend = tic;
for iter = 1:MaxIter
    problem.M = obliquefactoryNTrans(p, n);
    [Y0, Sd, dinf, fval, Axb, Axby, lamda]= trustregionsALM(problem, Y, opts);
    pinf = sqrt(Axb'*Axb)/normb;
    y = Axby;
    S = (Sd + Sd')/2;
    % [vS, dS] = eig(S);
    % dS = diag(dS);
    [vS, dS] = mineig(S);      % for large eigen decomposition problem
    %[vS, dS] = mineigimkl(S); % for small eigen decomposition problem
    % dinf = info(end).gradnorm;
    d  = dS(dS<0);
    if isempty(d) && pinf<tao  % S没有负的特征值，结束
        fprintf('Iter:%d, fval:%0.8f, gap:%0.1e, mineigS:%0.1e, pinf:%0.1e, dinf:%0.1e, r:%d, p:%d, sigam:%0.3f, time:%0.2fs\n', ...
            iter,    fval,       gap,       min(dS),       pinf,       dinf,       r,    p,    sigma,   toc(timespend));
        break;
    end
    rmiusS = min(length(d), 8); % 取小，防止Y增加太大
    v = vS(:,1:rmiusS)'; % 取主要的不大于8个负特征值对应的特征向量组
    mineigS = abs(d(1))/(1+dS(end));
    by = b'*y + sum(lamda);
    gap = abs(fval-by)/(abs(by)+abs(fval)+1);
    [~, D, V] = svd(Y0);
    e = diag(D);
    r = sum(e > 1e-3*e(1)); % r = rank(Y)
    if r <= p - 1
        % Y0 = diag(e(1:r))*V(:,1:r)';  % r个主分量
        Y0 = V(:,1:r)'.*e(1:r);  % r个主分量
        p = r;
    end
    p = p + rmiusS;
    Y = [Y0; 0.1*v]; % U = [zeros(p,n) ; v]; Y = Y + 0.1*U;
    Y = Y./sqrt(sum(Y.^2)); % nrms = sqrt(sum(Y.^2, 1)); Y = bsxfun(@times, Y, 1./nrms);

    fprintf('Iter:%d, fval:%0.8f, gap:%0.1e, mineigS:%0.1e, pinf:%0.1e, dinf:%0.1e, r:%d, p:%d, sigam:%0.3f, time:%0.2fs\n', ...
        iter,    fval,       gap,       mineigS,       pinf,       dinf,       r,    p,    sigma,   toc(timespend));
    if max(pinf, mineigS) < tao
        break;
    end
    %     if iter == 1 || neta > 0.7*eta
    %         % sigma = min(sigma*gama, 1);
    %         if sigma*gama > 1
    %            sigma = 1e-3;
    %         else
    %            sigma = sigma*gama;
    %         end
    %     end
    % eta = neta;

    if pinf/dinf < eta1
        it_pinf = it_pinf + 1; it_dinf = 0;
        if it_pinf >= h4
            sigma = max(sigma/gamma,sigma_min); it_pinf = 0;
        end
    elseif pinf/dinf > eta2
        it_dinf = it_dinf + 1; it_pinf = 0;
        if it_dinf >= h4
            sigma = min(sigma*gamma,sigma_max); it_dinf = 0;
        end
    end
    X = Y0'*Y0;
end

    function [f, store] = cost(Y, store)
        X0 = Y'*Y;
        x0 = X0(:);
        Axb = At'*x0 - b;
        Axby = y - sigma*Axb;
        fval = c'*x0;
        f = fval + 0.5/sigma*(Axby'*Axby);
    end

    function [G, store] = grad(Y, store)
        yA0 = reshape(A'*Axby, n, n);
        S = C - yA0;
        egrad = Y*S;
        % G = egrad - bsxfun(@times, Y, sum(Y.*egrad, 1));
        lamda = sum(Y.*egrad);
        G = 2*(egrad - Y.*lamda);
        store.Sdual = S - diag(lamda);
        store.fval = fval;
        store.Axb = Axb;
        store.Axby = Axby;
        store.lamda = lamda;
        % G = 2*Y*store.Sdual;
    end

    function [He, store] = hess(Y, Ydot, store)
        Xdot = Y'*Ydot;
        xdot = Xdot(:);
        AxbdotA = A'*(At'*xdot);
        yAdot = reshape(AxbdotA, n, n);
        H = 2*Ydot*S + 4*sigma*(Y*yAdot);
        % He = H - bsxfun(@times, Y, sum(Y.*H, 1)) - bsxfun(@times, Ydot, sum(Y.*egrad, 1)); % He = problem.M.ehess2rhess(Y, eG, H, Ydot);
        He = H - Y.*sum(Y.*H) - 2*Ydot.*lamda;
    end

    function M = obliquefactoryNTrans(n, m)
        % M.name = @() sprintf('Oblique manifold OB(%d, %d)', n, m);
        M.dim = @() (n-1)*m;
        M.inner = @(x, d1, d2) d1(:)'*d2(:);
        M.norm = @(x, d) norm(d(:));
        % M.dist = @(x, y) norm(real(2*asin(.5*sqrt(sum(trnsp(x - y).^2, 1)))));
        M.typicaldist = @() pi*sqrt(m);
        M.proj = @(X, U) U - X.*sum(X.*U);
        M.tangent = @(X, U) U - X.*sum(X.*U); %M.proj;

        M.retr = @retraction;
        % Retraction on the oblique manifold
        function y = retraction(x, d, t)
            xtd = x + d;
            y = xtd./sqrt(sum(xtd.^2)); % y = normalize_column
            % s(x + td);
        end

        M.rand = @() random(n, m);
        M.lincomb = @matrixlincomb;
        M.zerovec = @(x) zeros(n, m);
        M.transp = @(x1, x2, d) d - x2.*sum(x2.*d); %M.proj(x2, d);

        % Uniform random sampling on the sphere.
        function x = random(n, m)
            % x = normalize_columns(randn(n, m));
            x = randn(n, m);
            x = x./sqrt(sum(x.^2, 1));
        end
    end

    function [x, Sdual, norm_grad, fval, Axb, Axby, lamda] = trustregionsALM(problem, x, options)
        rho_prime = 0.1;
        rho_regularization = 1e3;
        options.storedepth = 2;
        options.subproblemsolver = @trs_tCG_cached;
        options.debug = false;
        M = problem.M;

        % If no initial point x is given by the user, generate one at random.
        if isempty(x)
            x = M.rand();
        end

        Delta_bar = M.typicaldist();
        Delta = Delta_bar / 8; % Initialize trust-region radius

        % Create a store database and get a key for the current x
        storedb = StoreDB(options.storedepth);
        key = storedb.getNewKey();

        k = 0;
        % accept tracks if the proposed step is accepted (true) or declined (false)
        accept = true;

        % Initialize solution and companion measures: f(x), fgrad(x)
        [fx, fgradx, Sdual, fval, Axb, Axby, lamda] = getCostGrad(problem, x, storedb, key);
        norm_grad = M.norm(x, fgradx);

        consecutive_TRplus = 0;
        consecutive_TRminus = 0;

        while true
            % Target gradient norm attained
            if  norm_grad < options.tolgradnorm || k >= options.maxiter
                break;
            end

            % Solve TR subproblem with solver specified by options.subproblemsolver
            trsinput = struct('x', x, 'fgradx', fgradx, 'Delta', Delta, ...
                'accept', accept);
            trsoutput = trs_tCG_cached(problem, trsinput, options, ...
                storedb, key);
            eta = trsoutput.eta;
            Heta = trsoutput.Heta;
            limitedbyTR = trsoutput.limitedbyTR;

            % Compute the tentative next iterate (the proposal)
            x_prop = M.retr(x, eta);
            key_prop = storedb.getNewKey();

            % Compute the function value of the proposal
            fx_prop = getCost(problem, x_prop, storedb, key_prop);

            % Will we accept the proposal or not?
            % Check the performance of the quadratic model against the actual cost.
            rhonum = fx - fx_prop;
            vecrho = M.lincomb(x, 1, fgradx, .5, Heta);
            rhoden = -M.inner(x, eta, vecrho);
            % rho_noreg = rhonum/rhoden;
            rho_reg_offset = max(1, abs(fx)) * eps * rho_regularization;
            rhonum = rhonum + rho_reg_offset;
            rhoden = rhoden + rho_reg_offset;
            model_decreased = (rhoden >= 0);
            rho = rhonum / rhoden;

            if rho < 1/4 || ~model_decreased || isnan(rho)
                %trstr = 'TR-';
                Delta = Delta/4;
                consecutive_TRplus = 0;
                consecutive_TRminus = consecutive_TRminus + 1;
            elseif rho > 3/4 && limitedbyTR
                % trstr = 'TR+';
                Delta = min(2*Delta, Delta_bar);
                consecutive_TRminus = 0;
                consecutive_TRplus = consecutive_TRplus + 1;
            else
                % Otherwise, keep the TR radius constant.
                consecutive_TRplus = 0;
                consecutive_TRminus = 0;
            end

            % Choose to accept or reject the proposed step based on the model
            % performance. Note the strict inequality.
            if model_decreased && rho > rho_prime
                accept = true;
                % accstr = 'acc';
                % We accept the step: no need to keep the old cache.
                storedb.removefirstifdifferent(key, key_prop);
                x = x_prop;
                key = key_prop;
                fx = fx_prop;
                [fgradx, Sdual, fval, Axb, Axby, lamda] = getGradient(problem, x, storedb, key);
                norm_grad = M.norm(x, fgradx);
            else
                % We reject the step: no need to keep cache related to the
                % tentative step.
                storedb.removefirstifdifferent(key_prop, key);
                accept = false;
                % accstr = 'REJ';
            end
            k = k + 1;
            % Make sure we don't use too much memory for the store database.
            storedb.purge();
        end  % of TR loop (counter: k)
        % Return the best cost reached
        cost = fx;
    end

    function [grad, Sdual, fval, Axb, Axby, lamda] = getGradient(problem, x, storedb, key)
        % Allow omission of the key, and even of storedb.
        if ~exist('key', 'var')
            if ~exist('storedb', 'var')
                storedb = StoreDB();
            end
            key = storedb.getNewKey();
        end

        % Contrary to most similar functions, here, we get the store by
        % default. This is for the caching functionality described below.
        store = storedb.getWithShared(key);
        store_is_stale = false;

        % If the gradient has been computed before at this point (and its
        % memory is still in storedb), then we just look up the value.
        force_grad_caching = true;
        if force_grad_caching && isfield(store, 'grad__')
            grad  = store.grad__;
            Sdual = store.Sdual__;
            fval =  store.fval__;
            Axb = store.Axb__;
            Axby = store.Axby__;
            lamda = store.lamda__;
            return;
        end

        if isfield(problem, 'grad')
            %% Compute the gradient using grad.

            % Check whether this function wants to deal with storedb or not.
            switch nargin(problem.grad)
                case 1
                    grad = problem.grad(x);
                case 2
                    [grad, store] = problem.grad(x, store);
                    Sdual = store.Sdual;
                    fval =  store.fval ;
                    Axb = store.Axb ;
                    Axby = store.Axby ;
                    lamda = store.lamda;
                case 3
                    % Pass along the whole storedb (by reference), with key.
                    grad = problem.grad(x, storedb, key);
                    % The store structure in storedb might have been modified
                    % (since it is passed by reference), so before caching
                    % we'll have to update (see below).
                    store_is_stale = true;
                otherwise
                    up = MException('manopt:getGradient:badgrad', ...
                        'grad should accept 1, 2 or 3 inputs.');
                    throw(up);
            end
        end

        % If we are not sure that the store structure is up to date, update.
        if store_is_stale
            store = storedb.getWithShared(key);
        end

        % Cache here.
        if force_grad_caching
            store.grad__ = grad;
            store.Sdual__ = Sdual;
            store.fval__ =  fval ;
            store.Axb__  =  Axb  ;
            store.Axby__ =  Axby ;
            store.lamda__ = lamda;
        end
        storedb.setWithShared(store, key);
    end

    function [cost, grad, Sdual, fval, Axb, Axby, lamda] = getCostGrad(problem, x, storedb, key)
        if ~exist('key', 'var')
            if ~exist('storedb', 'var')
                storedb = StoreDB();
            end
            key = storedb.getNewKey();
        end

        % Contrary to most similar functions, here, we get the store by
        % default. This is for the caching functionality described below.
        store = storedb.getWithShared(key);
        store_is_stale = false;

        % Check if the cost or gradient are readily available from the store.
        force_grad_caching = true;
        if isfield(store, 'cost__')
            cost = store.cost__;
            if force_grad_caching && isfield(store, 'grad__')
                grad = store.grad__;
                Sdual = store.Sdual__;
                fval =  store.fval__ ;
                Axb = store.Axb__ ;
                Axby = store.Axby__ ;
                lamda = store.lamda__ ;
                return;
            else
                grad = getGradient(problem, x, storedb, key); % caches grad
                return;
            end
        end
        % If we get here, the cost was not previously cached, but maybe the
        % gradient was?
        if force_grad_caching && isfield(store, 'grad__')
            grad = store.grad__;
            Sdual = store.Sdual__;
            fval =  store.fval__ ;
            Axb = store.Axb__ ;
            Axby = store.Axby__ ;
            lamda = store.lamda__ ;
            cost = getCost(problem, x, storedb, key); % this call caches cost
            return;
        end

        % Neither the cost nor the gradient were available: let's compute both.
        if ~isfield(problem, 'costgrad')
            cost = getCost(problem, x, storedb, key);
            [grad, Sdual, fval, Axb, Axby, lamda] = getGradient(problem, x, storedb, key);
            store_is_stale = true;
        end
        if store_is_stale
            store = storedb.getWithShared(key);
        end
        % Cache here.
        store.cost__ = cost;
        if force_grad_caching
            store.grad__ = grad;
            store.Sdual__ = Sdual;
            store.fval__ =  fval ;
            store.Axb__  =  Axb  ;
            store.Axby__ =  Axby ;
            store.lamda__ = lamda;
        end
        storedb.setWithShared(store, key);
    end
end