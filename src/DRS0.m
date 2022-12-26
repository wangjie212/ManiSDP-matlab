function [sol, gap] = DRS0(sb, sB, M1, M2, mb, fval, sol, tol)
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
gap = 1;
ogap = 2;
i = 1;
while abs(gap - ogap) > 1e-6 && gap > tol
    [V,D] = eig(psd, 'vector');
    psd = V*diag(max(0,D))*V';
    psol = sol;
    psol(1:vmb) = Mat2Vec(psd);
    lsol = 2*psol - sol;
    lsol = M1*lsol + ssb;
    sol = sol + 1*(lsol - psol);
    psd = Vec2Mat(sol(1:vmb), mb);
    if mod(i,200) == 0
        ogap = gap;
        minEig = min(eig(Vec2Mat(lsol(1:vmb), mb)));
        gap = - minEig*mb/abs(fval);
    %    disp(['DRS iteration ' num2str(i) ': gap <= ' num2str(gap)]);
    end
    i = i + 1;
end
