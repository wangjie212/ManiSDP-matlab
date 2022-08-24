function [X, fval] = mani_admm(SDP, At, b, c, mb)
tic
options.printlevel = 3;
options.tol = 1e-4;
while options.tol > 1e-7
if options.tol == 1e-4
    Y = [];
    aS = [];
else
    [V,D] = eig(aX{1});
    Y = [sqrt(D(end,end))*V(:,end) sqrt(D(end-1,end-1))*V(:,end-1)];
end
y = zeros(length(b),1);
[X, y, fval] = ALMSDP(At, b, c, mb, 10*options.tol, y, Y);
if isempty(aS) || fval < obj(1) - 1e-3
    aX = {X};
end
[obj,aX,~,y,aS] = admmplus(SDP.blk, SDP.At, SDP.C, SDP.b, [], [], [], [], [], options, aX, [], [], aS);
options.tol = options.tol/100;
end
time = toc;
X = aX{1};
fval = obj(1);
fprintf('\n');
disp(['ManiMoment: optimum = ' num2str(fval) ', time = ' num2str(time) 's']);
end