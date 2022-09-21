function [X, cx] = ALMSDP(At, b, c, mb, p, tao, y, Y)
C = reshape(c, mb, mb);
sigma = 1e-2;
gama = 13;
MaxIter = 30;
tolgrad = 1e-3;
for iter = 1:MaxIter
    [Y, ~, ~] = SDP_ALM_subprog(At, b, C, c, mb, p, sigma, y, Y, tolgrad);
    X = Y*Y';
    z = X(:);
    cx = z'*c;
    Axb = At'*z - b;
    y = y + 2*Axb*sigma;
   if norm(Axb) < tao
       break;
   else
        sigma = min(sigma*gama, 1e7);
        disp(['ALM iteration ' num2str(iter) ': fval = ' num2str(cx,10)]);
   end
end
end