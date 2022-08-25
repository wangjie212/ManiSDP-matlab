function [Sol, fval] = manisos(SDP, At, b, c, mb)
tic
options.printlevel = 0;
options.tol = 1e-4;
while options.tol > 1e-7
y = zeros(length(b),1);
if options.tol == 1e-4
    Y = [];
    aX = [];
    ay = [];
else
    [V,D] = eig(aS{1});
    Y = [sqrt(D(end,end))*V(:,end) sqrt(D(end-1,end-1))*V(:,end-1)];
end
[S{1}, ~] = ALMSDP(At, b, c, mb, y, Y);
S{2} = zeros(length(SDP.C{2}),1);
[obj,aX,~,ay,aS] = admmplus(SDP.blk, SDP.At, SDP.C, SDP.b, [], [], [], [], [], options, aX, [], ay, S);
options.tol = options.tol/100;
end
time = toc;
Sol = aS{1};
fval = - obj(1);
fprintf('\n');
disp(['ManiSOS: optimum = ' num2str(fval) ', time = ' num2str(time) 's']);
end