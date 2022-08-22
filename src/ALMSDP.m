function [X, fval] = ALMSDP(At, b, c, Nx, yk, Y)
m = length(b);
p = 2;
A = At';
C = reshape(c, Nx, Nx);
sigma = 1e-3;
gama = 13;
MaxIter = 20;
for iter = 1:MaxIter
    [Y, ~, ~] = SDP_ALM_subprog(A, At, b, C, c, Nx, p, sigma, yk, Y);
    X = Y*Y';
    z = X(:);
    cx = z'*c;
    Axb = A*z - b;
    if norm(Axb) < 1e-4
        break;
    else
        % disp(['Iter ' num2str(iter) ': fval = ' num2str(cx,10)]);
        yk = yk + 2*Axb*sigma;
        sigma = min(sigma*gama, 1e4);
    end
end
fval = cx;
end