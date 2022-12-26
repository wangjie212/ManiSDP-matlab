% sdpt format data

function [blk, At, C, b, dA] = SOStoSDP_C(f, h, x, d)
pop = [f; h];
n = length(x);
m = length(h);
coe = cell(1, m+1);
supp = cell(1, m+1);
lt = zeros(1, m+1);
dg = zeros(1, m);
[~, fsupp, fcoe] = decomp(pop);
for k = 1:m+1
    [~, loc, coe{k}] = find(fcoe(k,:));
    lt(k) = length(loc);
    supp{k} = fsupp(loc,:)';
    if k > 1
        dg(k-1) = deg(pop(k), x);
    end
end

fbasis = get_basis(n, d);
flb = size(fbasis, 2);
hlb = zeros(m, 1);
hbasis = cell(1, m);
for k = 1:m
    hbasis{k} = get_basis(n, 2*d-dg(k));
    hlb(k) = size(hbasis{k}, 2);
end
sp = get_basis(n, 2*d);
lsp = size(sp, 2);

b = sparse(lsp, 1);
for i = 1:lt(1)
    locb = bfind(sp, lsp, supp{1}(:,i), n);
    b(locb) = coe{1}(i);
end
At{1,1} = sparse(flb*(flb+1)/2, lsp);
dA = zeros(lsp, 1);
for i = 1:flb
    for j = i:flb
        bi = fbasis(:,i) + fbasis(:,j);
        locb = bfind(sp, lsp, bi, n);
        if i == j
            At{1,1}(j*(j+1)/2, locb) = 1;
            dA(locb) = dA(locb) + 1;
        else
            At{1,1}(i+j*(j-1)/2, locb) = sqrt(2);
            dA(locb) = dA(locb) + 2;
        end
    end
end
At{2,1} = sparse(sum(hlb)+1,lsp);
loc = 0;
for k = 1:m
    for i = 1:hlb(k)
        for j = 1:lt(k+1)
            bi = hbasis{k}(:,i) + supp{k+1}(:,j);
            locb = bfind(sp, lsp, bi, n);
            At{2,1}(loc+i, locb) = coe{k+1}(j);
        end
    end
    loc = loc + hlb(k);
end
At{2,1}(end,1) = 1;
blk{1,1} = 's';
blk{1,2} = flb;
blk{2,1} = 'u';
blk{2,2} = size(At{2,1},1);
C{1,1} = sparse(flb,flb);
C{2,1} = sparse(size(At{2,1},1),1);
C{2,1}(end) = -1;
end
