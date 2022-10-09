function [lsol, gap] = DRS(sb, sB, M1, M2, mb, fval, sol, ogap, tol, flag)
maxiter = 2000;
minprog = 1e-8;
checkiter = 20;
vmb = mb*(mb+1)/2;
if isempty(sol)
    psd = zeros(mb, mb);
    sol = zeros(vmb+size(sB,2), 1);
else
    psd = Vec2Mat(sol(1:vmb), mb);
end
ssb = sb;
ssb(1) = ssb(1) - fval;
ssb = M2*ssb;
i = 1;
gap = 1;
if isempty(ogap)
    ogap = 0.2;
end
prog = 1;
while i <= maxiter && prog > minprog
    [V,D] = eig(psd, 'vector');
    psd = V*diag(max(0,D))*V';
    psol = sol;
    psol(1:vmb) = Mat2Vec(psd);
    lsol = 2*psol - sol;
    lsol = M1*lsol + ssb;
    sol = sol + 1*(lsol - psol);
    psd = Vec2Mat(sol(1:vmb), mb);
    if mod(i,checkiter) == 0
        minEig = min(eig(Vec2Mat(lsol(1:vmb), mb)));
        gap = - minEig*mb/abs(fval);
        prog = abs(gap-ogap);
        ogap = gap;
        % gap = - minEig*3/abs(fval);
        disp(['DRS iteration ' num2str(i) ': gap <= ' num2str(gap)]);
    end
    i = i + 1;
end
