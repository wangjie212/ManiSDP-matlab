% A: correspond to the PSD constraint
% B: correspond to the h equality constraints

function [A, B, b] = SOStoSDP(f, h, x, d)
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
flb = nchoosek(n+d, d);
hlb = zeros(m, 1);
hbasis = cell(1, m);
for k = 1:m
%    hlb(k) = nchoosek(n+d-ceil(dg(k)/2), d-ceil(dg(k)/2));
%    hbasis{k} = get_basis(n, d-ceil(dg(k)/2));
    hbasis{k} = get_basis(n, 2*d-dg(k));
    hlb(k) = size(hbasis{k}, 2);
end
% sA = flb + sum(hlb);

sp = get_basis(n, 2*d);
% sp = unique(sp, 'columns');
lsp = size(sp, 2);
A = cell(1, lsp);
b = sparse(lsp, 1);
for i = 1:lt(1)
    locb = bfind(sp, lsp, supp{1}(:,i), n);
    b(locb) = coe{1}(i);
end
for i = 1:lsp
    A{i} = sparse(flb, flb);
end
for i = 1:flb
    for j = i:flb
        bi = fbasis(:,i) + fbasis(:,j);
        locb = bfind(sp, lsp, bi, n);
        A{locb}(i,j) = 1;
    end
end
B = sparse(lsp, sum(hlb));
loc = 0;
for k = 1:m
    for i = 1:hlb(k)
        for j = 1:lt(k+1)
            bi = hbasis{k}(:,i) + supp{k+1}(:,j);
            locb = bfind(sp, lsp, bi, n);
            B(locb,loc+i) = coe{k+1}(j);
        end
    end
    loc = loc + hlb(k);
end
end
