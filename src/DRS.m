function [lsol, gap] = DRS(sb, sB, M1, M2, mb, fval, sol, ogap, tol, flag)
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
if isempty(ogap)
    ogap = 0.2;
end
i = 1;
gap = 1;
% while i <= niter
while (i <= 50 || (gap > tol && gap <= ogap) || (flag == 1 && gap > tol)) && i <= 1000
    [V,D] = eig(psd, 'vector');
    psd = V*diag(max(0,D))*V';
    psol = sol;
    psol(1:vmb) = Mat2Vec(psd);
    lsol = 2*psol - sol;
    lsol = M1*lsol + ssb;
    sol = sol + 1*(lsol - psol);
    psd = Vec2Mat(sol(1:vmb), mb);
    if mod(i,10) == 0
        minEig = min(eig(Vec2Mat(lsol(1:vmb), mb)));
        gap = - minEig*mb/abs(fval);
        % gap = - minEig*3/abs(fval);
        disp(['DRS iteration ' num2str(i) ': gap <= ' num2str(gap)]);
    end
    if mod(i,50) == 0
        ogap = ogap/2;
    end
    i = i + 1;
end
