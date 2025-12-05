function M = multiblockmanifold(pset, nset, nob)
    nelems = length(pset);
    M.dim = @() (pset(1:nob)-1)*nset(1:nob)' + pset(nob+1:end)*nset(nob+1:end)';
    
    M.inner = @inner;
    function val = inner(x, u, v)
           val = innerc(nelems, u, v);
    end

    M.norm = @(x, d) sqrt(M.inner(x, d, d));
    if nelems > nob
        M.typicaldist = @() sqrt(pi*sum(nset(1:nob)) + pset(nob+1:end)*nset(nob+1:end)');
    else
        M.typicaldist = @() sqrt(pi*sum(nset(1:nob)));
    end

    M.proj = @proj;
    M.tangent = @proj;
    function v = proj(x, u)
       v = projc(nelems, nob, x, u);
    end
   
    M.retr = @retr;
    function y = retr(x, u, t)
        y = retrc(nelems, nob, x, u);
    end
    
    M.lincomb = @lincomb;
    function v = lincomb(x, a1, u1, a2, u2)
       v = lincombc(nelems, a1, u1, a2, u2);
    end

    M.rand = @rand;
    function x = rand()
        x = randc(nelems, nob, pset, nset);
    end

    M.zerovec = @zerovec;
    function u = zerovec(x)
         u = zerovecc(nelems, x);
    end
end
