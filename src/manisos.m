function [M, fval] = manisos(SDP, At, b, c, mb, sb, sB, M1, M2, maniIter)
fprintf('ManiSOS is starting...\n');
tic
gap = 1;
i = 1;
bsol = [];
bgap = [];
while i <= maniIter && gap > 1e-2
y = zeros(length(b),1);
[M, fval] = ALMSDP(At, b, c, mb, 1e-4, y, []);
[sol, gap] = DRS(sb, sB, M1, M2, mb, fval, bsol, bgap);
if isempty(bsol) || fval < best
    best = fval;
    bM = M;
    bsol = sol;
    bgap = gap;
end
i = i + 1;
end
fprintf('spADMM is starting...\n');
S = [bsol;best];
vM = mattovec(M);
y = SDP.At{1}\vM;
sM = [vM; SDP.At{2}*y+SDP.C{2}];
[fval, ~, ~, M] = spadmm(SDP.blk, SDP.At, SDP.C, SDP.b, SDP.dA, S, sM);
time = toc;
disp(['ManiSOS: optimum = ' num2str(fval) ', time = ' num2str(time) 's']);
end