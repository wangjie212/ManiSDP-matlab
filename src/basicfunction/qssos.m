% Generate the second-order SOS-SDP relaxation for the quartic sphere program:
% Min  coe'*[x]_4
% s.t. x^2 = 1.
% Output sedumi format data.
% dAAt: the diagnal of A*A'

function [A, b, c, K, dAAt] = qssos(n, coe)
sp2 = get_basis(n, 2);
sp4 = get_basis(n, 4);
mb = size(sp2, 2);
lsp = size(sp4, 2);
row = zeros(mb^2+mb*(n+1)+1, 1);
col = zeros(mb^2+mb*(n+1)+1, 1);
val = ones(mb^2+mb*(n+1)+1, 1);
dAAt = zeros(lsp, 1);
row(1) = 1;
col(1) = 1;
ind = 2;
for i = 1:mb
    for j = 1:n+1
        temp = zeros(n,1);
        if j < n + 1
            temp(j) = 2;
        else
            val(ind) = -1;
        end
        bi = sp2(:,i) + temp;
        locb = bfind(sp4, lsp, bi, n);
        row(ind) = locb;
        col(ind) = i + 1;
        ind = ind + 1;
    end
end
for i = 1:mb
    for j = i:mb
        bi = sp2(:,i) + sp2(:,j);
        locb = bfind(sp4, lsp, bi, n);
        row(ind) = locb;
        col(ind) = i*mb + j + 1;
        if j > i
            row(ind+1) = locb;
            col(ind+1) = j*mb + i + 1;
            dAAt(locb) = dAAt(locb) + 2;
            ind  = ind + 2;
        else
            dAAt(locb) = dAAt(locb) + 1;
            ind = ind + 1;
        end
    end
end
A = sparse(row,col,val,lsp,mb^2+mb+1);
b = coe;
K.f = mb + 1;
K.s = mb;
c = zeros(mb^2+mb+1, 1);
c(1) = 1;
end
