function fval = DRSPOP(At, b, c, mb, sb, sB, M1, M2)
fprintf('DRSPOP is starting...\n');
tic
y = zeros(length(b),1);
[~, fval] = ALMSDP(At, b, c, mb, 2, 1e-4, y, []);
[sol, gap] = DRS(sb, sB, M1, M2, mb, fval, [], [], 1e-2, 200, 0);
lb = fval*(1+gap);
[~, fval, gap] = nDRS(sb, sB, M1, M2, mb, fval, sol, lb, 1e-6);
time = toc;
disp(['DRSPOP: optimum = ' num2str(fval) ', gap <= ' num2str(gap) ', time = ' num2str(time) 's']);
end