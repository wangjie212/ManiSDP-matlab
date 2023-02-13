% clear all
% clc

rng(1);
N  = 15; % code length
x  = msspoly('x',N); % symbolic decision variables using SPOTLESS
% e  = 1e-4*randn(N,1);
% f  = x'*e; % 偏置加速收敛,只取一个最优解
f = 0;
for k=1:N-2 % 副瓣平方和
    f = f + sum(x(1:N-k).*x(1+k:N))^2;
end
h  = x.^2 - 1; % equality constraints of binary code design

%% Relax BQP into an SDP
problem.vars            = x;
problem.objective       = f;
problem.equality        = h; 
kappa                   = 3; % relaxation order
[SDP,info]              = dense_sdp_relax_binary(problem,kappa);
At = SDP.sedumi.At;
b = SDP.sedumi.b;
c = SDP.sedumi.c;
K = SDP.sedumi.K;
mb = K.s;

rng(0);
tic
[X, fval, data] = ManiSDP_unitdiag(At, b, c, mb);
tmani = toc;