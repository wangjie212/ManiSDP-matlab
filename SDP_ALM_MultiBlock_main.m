clc
clear 

%randn('state',1);
fname =  'data/gpp250-1.dat-s';
fname =  'data/control3.dat-s';
%fname =  'data/hinf11.dat-s';
fname =  'data/truss8.dat-s';
%   truss8 
%   136.3546

[At0,b,c0,K]=fromsdpa(fname);
Nx = K.l + sum(K.s);
Nz = length(K.s) + K.l;
m = size(b,1);
p = 10;

%% 表示为1个单个大阵
At = sparse(Nx*Nx, m);
xtt = zeros(Nx, Nz);
xtt((0:K.l-1)*Nx+(1:K.l)) = 1;
for ii = 1:length(K.s)
    offset = K.l+sum(K.s(1:ii-1));
    xtt(offset+1:offset+K.s(ii),K.l+ii) = 1;
end
xttt = repmat(xtt,p,1);
xttt = reshape(xttt,Nx, p*Nz);
IdXX = find((xttt*xttt')>0);
for ii = 1:m
    At(IdXX,ii) = At0(:,ii);
end
A = At';
c = sparse(Nx*Nx, 1);
c(IdXX) = c0;
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
    end
end
toc

fvalend = cx
s = K.s
l = K.l
