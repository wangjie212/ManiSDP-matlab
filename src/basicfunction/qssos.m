% Generate the second-order SOS-SDP relaxation for the quartic sphere program:
% Min  coe'*[x]_4
% s.t. x^2 = 1.
% Output sedumi format data.
% dAAt: the diagnal of A*A'

function [A, b, dAAt, mb] = qssos(n, coe)
sp = get_basis(n, 4);
mb = nchoosek(n+2, 2) - n;
lsp = size(sp, 2);
row = ones(mb^2, 1);
col = zeros(mb^2, 1);
val = ones(mb^2, 1);
dAAt = zeros(lsp, 1);
dAAt(1) = mb;
col(1:mb) = (0:(mb-1))*mb+(1:mb);
ind = mb + 1;
for i = 1:mb
    for j = i+1:mb
        bi = sp(:,i) + sp(:,j);
        locb = bfind(sp, lsp, bi, n);
        row(ind) = locb;
        row(ind+1) = locb;
        col(ind) = (i-1)*mb+j;
        col(ind+1) = (j-1)*mb+i;
        dAAt(locb) = dAAt(locb) + 2;
        ind = ind + 2;
    end
end
A = sparse(row,col,val,lsp,mb^2);
b = coe;
end
