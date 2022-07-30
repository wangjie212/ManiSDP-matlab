clc
clear 

%randn('state',1);
fname =  'data/gpp250-1.dat-s';
%fname =  'data/control3.dat-s';
fname =  'data/hinf8.dat-s';
fname =  'data/truss2.dat-s';
%   truss8 
%   136.3546

[At,b,c,Nx,m]=fromsdpaSM(fname);

p = 10;

A = At';
C = reshape(c, Nx, Nx);


%% ALM方法参数设置
sigma = 8.5;
gama = 100;
Y = [];
yk = zeros(m,1);

%% 迭代循环
tic
MaxIter = 1000;
for iter = 1:MaxIter
    [Y, fval, info] = SDP_ALM_subprog(A, At, b, C, c, Nx, p, sigma, yk, Y);
    X = Y*Y';
    x = X(:);
    cx = x'*c;
    Axb = (x'*At)' - b ;
    if norm(Axb) < 1e-6
        break;
    else
        disp(['Iter:' num2str(iter) ' ,fval=' num2str(cx,10)])
        yk = yk + 2*Axb*sigma;
        sigma = sigma * gama;
        sigma = min(sigma,10000);
    end
end
toc

fvalend = cx

