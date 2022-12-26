function [X, cx] = mani_nlp(px, C, par)
m = length(par.b);
p = 2;
A = par.At';
C = reshape(par.c, par.Nx, par.Nx);
options.maxtime = inf;
sigma = 1e-3;
gama = 13;
Y = msubs(par.v, par.x, px);
Y = [Y zeros(size(Y,1),1)];
yk = zeros(m,1);

MaxIter = 20;
for iter = 1:MaxIter
    [Y, fval, info] = SDP_ALM_subprog(A, par.At, par.b, C, par.c, par.Nx, p, sigma, yk, Y);
    X = Y*Y';
    z = X(:);
    cx = z'*par.c;
    Axb = (z'*par.At)' - par.b;
    if norm(Axb) < 1e-4
        break;
    else
        disp(['Iter ' num2str(iter) ': fval = ' num2str(cx,10)]);
        yk = yk + 2*Axb*sigma;
        sigma = min(sigma*gama, 1e4);
    end
end

end