clc
clear

filename = "D:\project\manisdp-matlab\data\gpp500-1.dat-s";
[At,b,c,K] = fromsdpa(filename);

rng(0);
tic
[X, S, y, fval] = ALMSDPNT_EIGV2(At, b, c, K.s);
tNTS_EIGV3 = toc
