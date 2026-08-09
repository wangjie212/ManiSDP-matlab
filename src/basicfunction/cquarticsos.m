%CQUARTICSOS  Fast order-2 reduced mixed complex SOS relaxation.
%
%   [A,b,c,K] = cquarticsos(Q,n)
%
% outputs SeDuMi data for the second-order complex SOS relaxation
%
%       minimize    [z]_2' * Q * [z]_2
%       subject to  abs(z_i)^2 = 1,  i = 1,...,n,
%
% using the reduced mixed Gram basis
%
%       {z^alpha*conj(z)^beta : |alpha|+|beta| <= 2}
%       \ {z_i*conj(z_i) : i=1,...,n}.
%
% This implementation is algebraically equivalent to the original code,
% but removes containers.Map/string keys and all entrywise mb-by-mb loops.
% The coefficient matrix A depends only on n and is cached between calls.
% Run "clear cquarticsos" to release or rebuild the cache.
%
% NOTE: c(1)=1 is retained from the original implementation.  Under the
% standard SeDuMi convention min c'*x, use c(1)=-1 to maximize gamma.

function [A,b,c,K] = cquarticsos(Q,n)

if nargin < 1
    error('cquarticsos:NotEnoughInputs','Q is required.');
end
if ~isnumeric(Q) || ndims(Q) ~= 2 || size(Q,1) ~= size(Q,2)
    error('cquarticsos:InvalidQ','Q must be a square numeric matrix.');
end

if nargin < 2 || isempty(n)
    % size(Q,1) = (n+1)(n+2)/2 for the holomorphic degree-2 basis.
    ntry = (-3 + sqrt(1 + 8*size(Q,1)))/2;
    if abs(ntry-round(ntry)) > 1e-12 || ntry < 1
        error('cquarticsos:CannotInferDimension', ...
            'Supply n explicitly; size(Q,1) is not nchoosek(n+2,2).');
    end
    n = round(ntry);
end
if ~isscalar(n) || ~isnumeric(n) || ~isreal(n) || n < 1 || n ~= floor(n)
    error('cquarticsos:InvalidDimension','n must be a positive integer.');
end

mh = nchoosek(n+2,2);
if ~isequal(size(Q),[mh,mh])
    error('cquarticsos:InvalidQSize', ...
        'Q must be %d-by-%d because [z]_2 has %d monomials.',mh,mh,mh);
end

hermerr = norm(Q-Q','fro');
qscale = max(1,norm(Q,'fro'));
if hermerr > 1e-12*qscale
    error('cquarticsos:NonHermitianQ', ...
        'Q must be Hermitian so that the objective is real-valued.');
end
Q = (Q+Q')/2;

% A, c, K, and the objective coefficient map depend only on n.
persistent cache
if isempty(cache) || cache.n ~= n
    cache = build_cached_model(n,mh);
end

A = cache.A;
c = cache.c;
K = cache.K;

% Assemble the Laurent coefficients of [z]_2'*Q*[z]_2.  The sparse
% constructor performs the required summation over repeated exponents.
qcoeff = Q(cache.objective_linear_indices);
b = full(sparse(cache.objective_rows, ...
                ones(numel(qcoeff),1), ...
                qcoeff,cache.number_of_rows,1));
b(1) = real(b(1));

end

function model = build_cached_model(n,mh)
%BUILD_CACHED_MODEL  Construct all Q-independent data for a fixed n.

[sp,degs] = build_bases(n,mh);
mb = size(degs,2);

%% Coefficient map for v'*H*v
% Enumerate one unordered pair {a,q}.  Canonicalizing delta and then
% explicitly inserting both symmetric matrix positions gives exactly the
% same matrices as the final symmetrization step in the original code.
[ia,iq] = ordered_upper_pairs(mb);
D = degs(:,iq)-degs(:,ia);
[D,~] = canonicalize_deltas(D);

% Stable ordering keeps delta=0 as the first coefficient equation.
[dsp,~,gram_rows] = unique(D.','rows','stable');
number_of_rows = size(dsp,1);
if any(dsp(1,:) ~= 0)
    error('cquarticsos:InternalMapError', ...
        'The zero Laurent exponent is not the first coefficient row.');
end

is_diagonal = (ia == iq);
is_off_diagonal = ~is_diagonal;
upper_columns = ia + (iq-1)*mb;
lower_columns = iq(is_off_diagonal) + ...
                (ia(is_off_diagonal)-1)*mb;

upper_values = 0.5*ones(numel(ia),1);
upper_values(is_diagonal) = 1;

Agram = sparse([gram_rows; gram_rows(is_off_diagonal)], ...
               [upper_columns; lower_columns], ...
               [upper_values; 0.5*ones(nnz(is_off_diagonal),1)], ...
               number_of_rows,mb^2);

% Every Gram variable occurs in exactly one coefficient equation.
if nnz(Agram) ~= mb^2
    error('cquarticsos:InternalMapError', ...
        'Unexpected number of Gram coefficient entries.');
end

%% Coefficient map for [z]_2'*Q*[z]_2
[oa,oq] = ordered_upper_pairs(mh);
Dobj = sp(:,oq)-sp(:,oa);
[Dobj,flip] = canonicalize_deltas(Dobj);

[found,objective_rows] = ismember(Dobj.',dsp,'rows');
if ~all(found)
    error('cquarticsos:InternalMapError', ...
        'An objective monomial is missing from the SOS coefficient map.');
end

% If the raw exponent has the wrong orientation, use Q(q,a), not Q(a,q),
% so that its coefficient belongs to the canonical exponent delta.
objective_linear_indices = oa + (oq-1)*mh;
objective_linear_indices(flip) = ...
    oq(flip) + (oa(flip)-1)*mh;

%% Complete SeDuMi model
A = [sparse(1,1,1,number_of_rows,1),Agram];

c = sparse(1+mb^2,1);
c(1) = 1;  % Preserves the original code; use -1 to maximize gamma.

K = struct();
K.f = 1;
K.l = 0;
K.q = [];
K.r = [];
K.s = mb;

model = struct();
model.n = n;
model.A = A;
model.c = c;
model.K = K;
model.objective_rows = objective_rows;
model.objective_linear_indices = objective_linear_indices;
model.number_of_rows = number_of_rows;

end

function [sp,degs] = build_bases(n,mh)
%BUILD_BASES  Vectorized exponent construction for both monomial bases.

E = int8(eye(n));

% Holomorphic basis [z]_2 = [1,z_i,z_i*z_j].
sp = zeros(n,mh,'int8');
sp(:,2:n+1) = E;
[hi,hj] = ordered_upper_pairs(n);
sp(:,n+2:end) = E(:,hi)+E(:,hj);

% Laurent exponents of the reduced mixed basis.  On the unit torus,
% z^alpha*conj(z)^beta is represented by alpha-beta.
W = [E,-E];
[mi,mj] = ordered_upper_pairs(2*n);
is_unit_monomial = (mi <= n) & (mj == n+mi);
mi = mi(~is_unit_monomial);
mj = mj(~is_unit_monomial);

degs = [zeros(n,1,'int8'),W,W(:,mi)+W(:,mj)];
expected_mb = nchoosek(2*n+2,2)-n;
if size(degs,2) ~= expected_mb
    error('cquarticsos:InternalBasisError', ...
        'Unexpected reduced mixed basis size.');
end
if size(unique(degs.','rows'),1) ~= expected_mb
    error('cquarticsos:InternalBasisError', ...
        'The reduced mixed basis contains duplicate torus monomials.');
end

end

function [i,j] = ordered_upper_pairs(m)
%ORDERED_UPPER_PAIRS  Return (i,j), 1 <= i <= j <= m.
% Ordering matches the original nested loops: i=1:m, j=i:m.

number_of_pairs = m*(m+1)/2;
i = zeros(number_of_pairs,1);
j = zeros(number_of_pairs,1);
cursor = 0;
for ii = 1:m
    count = m-ii+1;
    ind = cursor + (1:count);
    i(ind) = ii;
    j(ind) = (ii:m).';
    cursor = cursor + count;
end

end

function [D,flip] = canonicalize_deltas(D)
%CANONICALIZE_DELTAS  Select one member of each pair {delta,-delta}.
% The retained member has positive first nonzero coordinate.

number_of_deltas = size(D,2);
[~,first_nonzero] = max(D ~= 0,[],1);
linear_indices = sub2ind(size(D),first_nonzero,1:number_of_deltas);
flip = D(linear_indices) < 0;
D(:,flip) = -D(:,flip);
flip = flip(:);

end
