function [Xr,pobj] = local_mani(Xmat,Cmat,par,opt,roundonly)
if nargin < 5
    roundonly = false;
end

if nargin < 4
    opt = [1,2]; % default round the first two eigenvectors
end

    nhypo   = length(opt);
    x       = bqp_round(Xmat, opt, par.d);
    fopt    = zeros(nhypo,1);
    for i = 1:nhypo
        [Xi,fopti] = mani_nlp(x(:,i),Cmat,par);
        Xopt{i}     =  Xi;
        fopt(i)       = fopti;
    end
    [pobj,idxmin]    = min(fopt);
    Xr           = {Xopt{idxmin}};
end


function x = bqp_round(Xmat,opt,d)
nhypo       = length(opt);
Xmat        = Xmat{1};
X           = Xmat(1:d+1,1:d+1); % order 0 and 1 moment matrix
[V,D]       = eig(X);
[~,idx]     = sort(diag(D),'descend');
V           = V(:,idx);
xraw        = V(:,opt);
x           = zeros(d,nhypo);
for i = 1:nhypo
    xrawi   = xraw(:,i);
    xrawi   = xrawi/xrawi(1); % normalize first entry to be 1
    xi      = xrawi(2:end);
    
    xi      = sign(xi);
    xi(xi==0) = 1.0;
    x(:,i)  = xi;
end
end


