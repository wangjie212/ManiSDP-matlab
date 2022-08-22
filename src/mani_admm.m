function [X, fval] = mani_admm(SDP, At, b, c, Nx)
tic
options.printlevel = 0;
options.tol = 1e-4;
while options.tol > 1e-7
if options.tol == 1e-4
    yk = zeros(length(b),1);
    Y = [];
    aS = [];
else
    [V,D] = eig(aX{1});
    Y = [sqrt(D(end,end))*V(:,end) sqrt(D(end-1,end-1))*V(:,end-1)];
end
[X, ~] = ALMSDP(At, b, c, Nx, yk, Y);
[obj,aX,~,yk,aS] = admmplus(SDP.blk, SDP.At, SDP.C, SDP.b, [], [], [], [], [], options, {X}, [], yk, aS);
options.tol = options.tol/100;
end
time = toc;
X = aX{1};
fval = obj(1);
fprintf('\n');
disp(['ManiMoment: optimum = ' num2str(fval) ', time = ' num2str(time) 's']);
end