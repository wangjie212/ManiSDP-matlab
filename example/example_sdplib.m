filename = "D:\project\manisdp-matlab\data\sdplib\gpp100.dat-s";
[At,b,c,K] = fromsdpa(filename);
mb = K.s;

rng(0);
tic
[~, ~, ~, fval, emani] = ALMSDPNT_EIGV3(At(:,1), b(1), full(c), mb);
tmani = toc;

fprintf('ManiSDP: optimum = %0.8f, eta = %0.1e, time = %0.2fs\n', fval, emani, tmani);
