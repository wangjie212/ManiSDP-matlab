%CQPSOS  Second-order complex SOS relaxation of a unit-modulus CQP.
%
%   [A,b,c,K,dAAt] = cqpsos(Q,e,n)
%
% generates SeDuMi data for
%
%       minimize    z'*Q*z + real(e'*z)
%       subject to  abs(z_i)^2 = 1,  i = 1,...,n,
%
% using the order-2 complex SOS relaxation and the holomorphic basis
%
%       v_2(z) = (z^alpha : alpha in N^n, |alpha| <= 2).
%
% For the problem without a linear term, use e = [] or e = zeros(n,1).
% The matrix Q must be Hermitian.  The linear-term convention real(e'*z)
% reduces to e'*z when e and z are real, as in bqpsos.m.
%
% The returned data are in the standard SeDuMi primal form
%
%       minimize    c'*x
%       subject to  A*x = b,  x in K.
%
% Here x = [gamma; vec(H)], H is an mb-by-mb Hermitian PSD Gram matrix,
% and minimizing c'*x = -gamma maximizes the SOS lower bound gamma.
% SeDuMi is told that H is Hermitian by K.scomplex = 1.  Rows 2:end are
% complex equalities (both their real and imaginary parts are imposed).
%
% Outputs
%   A,b,c,K : SeDuMi data; call [x,y,info] = sedumi(A,b,c,K).
%   dAAt    : diagonal of A*A'.  Different rows have disjoint supports.

%
% The certificate represented by a feasible x is
%
%       z'*Q*z + real(e'*z) - gamma = v_2(z)'*H*v_2(z)
%
% modulo the relations conjugate(z_i)*z_i = 1.

function [A,b,c,K,dAAt] = cqpsos(Q,e,n)

if nargin < 1
    error('cqpsos:NotEnoughInputs','Q is required.');
end
if nargin < 3 || isempty(n)
    n = size(Q,1);
end
if nargin < 2 || isempty(e)
    e = zeros(n,1);
end

if ~isscalar(n) || n < 1 || n ~= floor(n)
    error('cqpsos:InvalidDimension','n must be a positive integer.');
end
if ~isnumeric(Q) || ~isequal(size(Q),[n,n])
    error('cqpsos:InvalidQ','Q must be an n-by-n numeric matrix.');
end
if ~isnumeric(e) || numel(e) ~= n
    error('cqpsos:InvalidLinearTerm','e must have n entries.');
end
e = e(:);

hermerr = norm(Q-Q','fro');
qscale = max(1,norm(Q,'fro'));
if hermerr > 1e-12*qscale
    error('cqpsos:NonHermitianQ', ...
        'Q must be Hermitian so that the objective is real-valued.');
end
% Remove harmless roundoff asymmetry.
Q = (Q+Q')/2;

%% Holomorphic monomials of total degree at most two
mb = nchoosek(n+2,2);
sp = zeros(n,mb);
p = 1;                         % constant monomial
for i = 1:n                    % degree-one monomials
    p = p+1;
    sp(i,p) = 1;
end
for i = 1:n                    % degree-two monomials
    for j = i:n
        p = p+1;
        sp(i,p) = sp(i,p)+1;
        sp(j,p) = sp(j,p)+1;
    end
end
if p ~= mb
    error('cqpsos:InternalBasisError','Unexpected basis size.');
end

%% Quotient coefficient map H -> coefficients of v_2(z)'*H*v_2(z)
% On |z_i|=1, conjugate(z)^alpha*z^beta = z^(beta-alpha).
% Only one of delta and -delta is retained.  The retained representative
% is the one whose first nonzero entry is positive; Hermitian symmetry
% supplies the conjugate coefficient at -delta.
maxnnz = mb*(mb+1)/2;
row = zeros(maxnnz,1);
col = zeros(maxnnz,1);
val = ones(maxnnz,1);

dsp = zeros(1,n);              % delta = 0 is always row 1
dmap = containers.Map('KeyType','char','ValueType','double');
dmap(delta_key(zeros(n,1))) = 1;

p = 0;
for a = 1:mb
    for q = 1:mb
        delta = sp(:,q)-sp(:,a);
        first = find(delta~=0,1,'first');
        if isempty(first) || delta(first)>0
            key = delta_key(delta);
            if isKey(dmap,key)
                loc = dmap(key);
            else
                loc = size(dsp,1)+1;
                dmap(key) = loc;
                dsp(loc,:) = delta.';
            end
            p = p+1;
            row(p) = loc;
            col(p) = a+(q-1)*mb;       % H(a,q) in MATLAB vec ordering
        end
    end
end
if p ~= maxnnz
    error('cqpsos:InternalMapError','Unexpected coefficient-map size.');
end
lsp = size(dsp,1);
Agram = sparse(row,col,val,lsp,mb^2);

%% Laurent coefficients of the objective on the unit torus
b = complex(zeros(lsp,1));
b(1) = real(trace(Q));          % Q(i,i)*conj(z_i)*z_i = Q(i,i)

% z'*Q*z: coefficient Q(i,j) has exponent delta = e_j-e_i.
for i = 1:n
    for j = 1:n
        if i == j
            continue
        end
        delta = zeros(n,1);
        delta(j) = delta(j)+1;
        delta(i) = delta(i)-1;
        first = find(delta~=0,1,'first');
        if delta(first)>0
            loc = dmap(delta_key(delta));
            b(loc) = b(loc)+Q(i,j);
        end
    end
end

% real(e'*z) = (e'*z + z'*e)/2.  For delta = e_i, its coefficient
% is conjugate(e_i)/2; the coefficient at -e_i is its conjugate.
for i = 1:n
    delta = zeros(n,1);
    delta(i) = 1;
    loc = dmap(delta_key(delta));
    b(loc) = b(loc)+conj(e(i))/2;
end

%% Complete SeDuMi model
% The constant equation is gamma + sum_alpha H(alpha,alpha) = b(1).
gamma_column = sparse(1,1,1,lsp,1);
A = [gamma_column,Agram];
c = sparse(1+mb^2,1);
c(1) = -1;                     % minimize -gamma

K = struct();
K.f = 1;                       % gamma is a real free scalar
K.s = mb;                      % one PSD block of size mb
K.scomplex = 1;                % that PSD block is Hermitian
K.ycomplex = (2:lsp).';        % impose full complex coefficient equations

dAAt = full(sum(abs(A).^2,2));

end

function key = delta_key(delta)
% Exact key for a small integer exponent vector.
key = sprintf('%+d,',delta);
end
