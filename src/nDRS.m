% function [sol, fval, gap] = nDRS(sb, sB, M1, M2, mb, fval, lb, tol)
% ub = fval;
% i = 1;
% while (ub - lb)/abs(lb) > 2*tol
%     fval = (ub+lb)/2;
%     [sol, gap] = DRS0(sb, sB, M1, M2, mb, fval, [], tol);
%     if gap > 1.5*tol
%         ub = fval;
%     elseif gap > 0
%         lb = fval;
%     else
%         lb = fval*(1+gap);
%     end
%     disp(['DRS iteration ' num2str(i) ': gap <= ' num2str(gap) ', fval = ' num2str(fval) ', lb = ' num2str(lb) ', ub = ' num2str(ub)]);
%     i = i + 1;
% end
% end

function [sol, fval, gap] = nDRS(sb, sB, M1, M2, mb, fval, sol, lb, tol)
vmb = mb*(mb+1)/2;
if isempty(sol)
    psd = zeros(mb, mb);
    sol = zeros(vmb+size(sB,2), 1);
else
    psd = Vec2Mat(sol(1:vmb), mb);
end
ssb = sb;
ssb(1) = ssb(1) - fval;
i = 1;
niter = 500;
while (fval - lb)/abs(lb) > 2*tol
    [V,D] = eig(psd, 'vector');
    psd = V*diag(max(0,D))*V';
    psol = sol;
    psol(1:vmb) = Mat2Vec(psd);
    lsol = 2*psol - sol;
    lsol = M1*lsol + M2*ssb;
    sol = sol + 1*(lsol - psol);
    psd = Vec2Mat(sol(1:vmb), mb);
    if mod(i,niter) == 0
        minEig = min(eig(Vec2Mat(lsol(1:vmb), mb)));
        gap = - minEig*mb/abs(fval);
        disp(['DRS iteration ' num2str(i) ': gap <= ' num2str(gap) ', fval = ' num2str(fval) ', lb = ' num2str(lb)]);
        if gap < 1.5*tol
           if gap > 0
               lb = fval;
           else
               lb = fval*(1+gap);
           end
           nfval = lb*(1-1e-2);
        else
           lb = max(lb, fval*(1+gap));
           nfval = (lb+fval)/2;
        end
        ssb(1) = ssb(1) + fval - nfval;
        fval = nfval;
    end
    i = i + 1;
end
% disp(['DRS-SOS: fval = ' num2str(fval) ', lb = ' num2str(lb)]);
end
