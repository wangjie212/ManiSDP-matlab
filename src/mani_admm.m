function [X, fval] = mani_admm(SDP, At, b, c, mb, sb, sB, M1, M2, maniIter)
fprintf('ManiMom is starting...\n');
tic
options.printlevel = 0;
options.tol = 1e-6;
gap = 1;
i = 1;
bsol = [];
bgap = [];
while i <= maniIter && gap > 1e-2
y = zeros(length(b),1);
[X, fval] = ALMSDP(At, b, c, mb, 1e-4, y, []);
[sol, gap] = DRS(sb, sB, M1, M2, mb, fval, bsol, bgap);
if isempty(bsol) || fval < best
    best = fval;
    bX = X;
    bsol = sol;
    bgap = gap;
end
i = i + 1;
end
fprintf('ADMM+ is starting...\n');
S = Vec2Mat(bsol(1:mb*(mb+1)/2), mb);
[obj,aX,~,~,~] = admmplus(SDP.blk, SDP.At, SDP.C, SDP.b, [], [], [], [], [], options, {bX}, [], [], {S});
fprintf('\n');
% [V,D] = eig(aX{1});
% Y = [sqrt(D(end,end))*V(:,end) sqrt(D(end-1,end-1))*V(:,end-1)];
% y = zeros(length(b),1);
% [X, fval] = ALMSDP(At, b, c, mb, 1e-4, y, Y);
% if (obj(1)-fval)/max(1,abs(obj(1))) > 1e-6
%     aX = {X};
% end
% options.tol = 1e-6;
% [obj,aX,~,~,~] = admmplus(SDP.blk, SDP.At, SDP.C, SDP.b, [], [], [], [], [], options, aX, [], [], aS);
% fprintf('\n');
time = toc;
X = aX{1};
fval = obj(1);
disp(['ManiMom: optimum = ' num2str(fval) ', time = ' num2str(time) 's']);
end