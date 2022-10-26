clc
clear

filename = "D:\ADMM_for_SDPmaster\ADMM_for_SDP\problems\gpp500-1.dat-s";
[At,b,c,K] = fromsdpa(filename);

randn('state',1);
tic
[X, S, y, fval] = ALMSDPNT_EIGV3(At, b, c, K.s);
tNTS_EIGV3 = toc
