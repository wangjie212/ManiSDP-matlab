%CQUARTICMOM  Fast reduced-basis order-2 complex moment relaxation.
%
%   [At,b,c,K] = cquarticmom(n,Q)
%
% outputs SeDuMi data for the second-order complex moment relaxation
%
%       minimize    v_2(z)'*Q*v_2(z)
%       subject to  abs(z_i)^2 = 1,  i=1,...,n,
%
% using the reduced mixed basis
%
%       {z^alpha*conj(z)^beta : |alpha|+|beta| <= 2}
%          \ {z_i*conj(z_i) : i=1,...,n}.
%
% The output convention is
%
%       minimize    c'*x
%       subject to  At'*x = b,  x in K,
%
% with x = vec(X) and K.s = mb.
%
% PERFORMANCE NOTES
%   * All Q-independent data are cached for the most recently used n.
%   * Laurent differences are grouped by an exact mixed-radix encoding,
%     avoiding an mb^2-by-n double matrix in the usual range of n.
%   * Moment-consistency equations are assembled in one sparse call.
%   * The objective is accumulated only from the mh^2 nonzero entries of
%     the embedded holomorphic block.
%
% Use
%
%       clear cquarticmom
%
% to release the cached model data.

function [At,b,c,K] = cquarticmom(n,Q)

%% Validate input
if nargin < 2
    error('cquarticmom:NotEnoughInputs','Both n and Q are required.');
end
if ~isscalar(n) || ~isnumeric(n) || ~isreal(n) || ...
        n < 1 || n ~= floor(n)
    error('cquarticmom:InvalidDimension','n must be a positive integer.');
end
if ~isnumeric(Q) || ndims(Q) ~= 2
    error('cquarticmom:InvalidQ','Q must be a numeric matrix.');
end

mh = nchoosek(n+2,2);
if ~isequal(size(Q),[mh,mh])
    error('cquarticmom:InvalidQSize', ...
        'Q must be %d-by-%d because length(v_2(z))=%d.',mh,mh,mh);
end

hermerr = norm(Q-Q','fro');
qscale = max(1,norm(Q,'fro'));
if hermerr > 1e-12*qscale
    error('cquarticmom:NonHermitianQ', ...
        'Q must be Hermitian so that v_2(z)''*Q*v_2(z) is real-valued.');
end
Q = 0.5*(Q+Q');

%% Build or reuse the Q-independent affine model
persistent cache
if isempty(cache) || cache.n ~= n
    cache = build_cached_model(n,mh);
end

At = cache.At;
b  = cache.b;
K  = cache.K;

%% Objective coefficients
% Only the holomorphic principal block is nonzero before moment-group
% averaging. sparse(...) sums all Q entries belonging to the same group.
group_sum = full(sparse(cache.objective_group, ...
                        ones(mh^2,1),Q(:),cache.ngroup,1));

c0 = group_sum(cache.entry_group) ./ cache.entry_group_count;

% Match the canonicalization performed by the original implementation:
% first Hermitianize C, then symmetrize its real SDP coefficient matrix.
cH = 0.5*(c0 + conj(c0(cache.transpose_order)));
c  = real(0.5*(cH + cH(cache.transpose_order)));

end

function model = build_cached_model(n,mh)
%BUILD_CACHED_MODEL Construct all data depending only on n.

%% Reduced full mixed basis, preserving the original ordering
mb = nchoosek(2*n+2,2)-n;
alpha = zeros(n,mb,'int8');
beta  = zeros(n,mb,'int8');
p = 1;                                  % constant monomial

% Degree-one monomials w=[z;conj(z)].
for i = 1:2*n
    p = p+1;
    if i <= n
        alpha(i,p) = 1;
    else
        beta(i-n,p) = 1;
    end
end

% Degree-two monomials, omitting z_i*conj(z_i).
for i = 1:2*n
    for j = i:2*n
        if i <= n && j == n+i
            continue
        end

        p = p+1;
        if i <= n
            alpha(i,p) = alpha(i,p)+1;
        else
            beta(i-n,p) = beta(i-n,p)+1;
        end
        if j <= n
            alpha(j,p) = alpha(j,p)+1;
        else
            beta(j-n,p) = beta(j-n,p)+1;
        end
    end
end

if p ~= mb || any(sum(double(alpha)+double(beta),1) > 2)
    error('cquarticmom:InternalBasisError', ...
        'Unexpected reduced mixed basis.');
end

degs = alpha-beta;
if size(unique(degs.','rows'),1) ~= mb
    error('cquarticmom:InternalBasisError', ...
        'The reduced basis contains duplicate Laurent monomials.');
end

holind = find(all(beta==0,1));
if numel(holind) ~= mh
    error('cquarticmom:InternalBasisError', ...
        'Could not identify the holomorphic degree-two subbasis.');
end

%% Group X(a,q) by degs(:,a)-degs(:,q)
% A base-9 digit represents one difference coordinate in {-4,...,4}.
% Blocks of at most 16 coordinates remain exact IEEE doubles because
% 9^16 < 2^53. Multiple blocks provide an exact collision-free key for
% arbitrary n while using only ceil(n/16) columns instead of n columns.
nentry = mb^2;
block_size = 16;
nblock = ceil(n/block_size);
entry_key = zeros(nentry,nblock);
D = double(degs.');                       % mb-by-n

for blk = 1:nblock
    first = (blk-1)*block_size+1;
    last  = min(blk*block_size,n);
    len   = last-first+1;

    weights = 9.^(0:len-1).';
    code = D(:,first:last)*weights;
    offset = 4*sum(weights);

    key_block = offset + code-code.';     % row a, column q
    entry_key(:,blk) = key_block(:);      % MATLAB vec ordering
end

if nblock == 1
    [~,representative,entry_group] = unique(entry_key,'sorted');
else
    [~,representative,entry_group] = unique(entry_key,'rows','sorted');
end

ngroup = max(entry_group);
group_count = accumarray(entry_group,1,[ngroup,1]);

% Determine which member of each conjugate pair {delta,-delta} is kept.
rep_a = mod(representative-1,mb)+1;
rep_q = floor((representative-1)/mb)+1;
group_delta = degs(:,rep_a)-degs(:,rep_q);
is_zero_group = all(group_delta==0,1);

nz = (group_delta~=0);
[has_nonzero,first_nonzero] = max(nz,[],1);
linear_first = sub2ind(size(group_delta),first_nonzero,1:ngroup);
is_positive_group = has_nonzero & (group_delta(linear_first)>0);
retained_group = is_zero_group | is_positive_group;

zero_group = find(is_zero_group);
if numel(zero_group) ~= 1 || representative(zero_group) ~= 1
    error('cquarticmom:InternalGroupError', ...
        'The zero Laurent group is inconsistent.');
end

%% Vectorized moment-consistency equations
% representative(g) is the first vec(X) entry in group g. Retain all
% entries in the zero group and in one orientation of each nonzero pair,
% except for the representative itself.
all_entry = (1:nentry).';
rep_for_entry = representative(entry_group);
keep = retained_group(entry_group).' & (all_entry ~= rep_for_entry);
current = all_entry(keep);
reference = rep_for_entry(keep);

ncons = 1+numel(current);
cons_col = (2:ncons).';

transpose_order = reshape(1:nentry,mb,mb).';
transpose_order = transpose_order(:);

% Directly form symmetric coefficient matrices. Duplicate triplets on
% diagonal entries are intentionally summed by sparse(...).
rows = [1; current; transpose_order(current); ...
           reference; transpose_order(reference)];
cols = [1; cons_col; cons_col; cons_col; cons_col];
vals = [1; 0.5*ones(numel(current),1); ...
          0.5*ones(numel(current),1); ...
         -0.5*ones(numel(current),1); ...
         -0.5*ones(numel(current),1)];

At = sparse(rows,cols,vals,nentry,ncons);
b = zeros(ncons,1);
b(1) = 1;

%% Objective-group lookup
% Q(:) is embedded at C0(holind,holind) in the same column-major order.
objective_entry = holind(:) + (holind(:).'-1)*mb;
objective_group = entry_group(objective_entry(:));

%% SeDuMi cone
K = struct();
K.f = 0;
K.l = 0;
K.q = [];
K.r = [];
K.s = mb;

model = struct();
model.n = n;
model.At = At;
model.b = b;
model.K = K;
model.ngroup = ngroup;
model.entry_group = entry_group;
model.entry_group_count = group_count(entry_group);
model.objective_group = objective_group;
model.transpose_order = transpose_order;

end
