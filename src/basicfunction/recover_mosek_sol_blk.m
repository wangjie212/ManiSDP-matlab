function [Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res, K)
if ~isfield(K, 'f') && isscalar(K.s)
    n               = K.s;
    idxtril         = tril(true(n,n));
    Xopt            = zeros(n,n);
    Xopt(idxtril)   = res.sol.itr.barx;
    Xopt            = Xopt + Xopt';
    Xopt            = Xopt - 0.5*diag(diag(Xopt));
    Sopt            = zeros(n,n);
    Sopt(idxtril)   = res.sol.itr.bars;
    Sopt            = Sopt + Sopt';
    Sopt            = Sopt - 0.5*diag(diag(Sopt));
    
    Xopt            = {Xopt};
    Sopt            = {Sopt};
else
    if isfield(K, 'f')
        Xopt  = cell(length(K.s)+1, 1);
        Xopt{1} = res.sol.itr.xx;
    else
        Xopt  = cell(length(K.s), 1);
    end
    Sopt  = cell(length(K.s), 1);
    cid   = 0;
    for i = 1:length(K.s)
        n               = K.s(i);
        ndelta          = triangle_number(n);
        idxtril         = tril(true(n,n));
        Xopti           = zeros(n,n);
        Xopti(idxtril)  = res.sol.itr.barx(cid+1:cid+ndelta);
        Xopti           = Xopti + Xopti';
        Xopti           = Xopti - 0.5*diag(diag(Xopti));
        Sopti           = zeros(n,n);
        Sopti(idxtril)  = res.sol.itr.bars(cid+1:cid+ndelta);
        Sopti           = Sopti + Sopti';
        Sopti           = Sopti - 0.5*diag(diag(Sopti));        
        cid             = cid + ndelta;
        if isfield(K, 'f')
            Xopt{i+1}     = Xopti;
        else
            Xopt{i}     = Xopti;
        end
        Sopt{i}     = Sopti;
    end
end
yopt            = res.sol.itr.y;
obj             = [res.sol.itr.pobjval;res.sol.itr.dobjval];
end