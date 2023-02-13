function [A,b,c,K] = generate_hamming(k,d)
% PROTOTYPE:
%    [A,b,c,K] = generate_hamming(k,d)
%
% PARAMETERS:
%    k: number of bits per vector.
%    d: vector of allowed hamming distances
%
% PURPOSE:
%    Create an SDP in SeDuMi format that computes the theta function
%    for the Hamming graph H_{k,d}.
%
% RETURN VALUES:
%    A,b,c,K: SDP problem in SeDuMi format

% 
% REVISION INFORMATION:
% $Id: generate_hamming.m,v 1.1 2000/05/01 15:42:46 schmieta Exp $
% $Log: generate_hamming.m,v $
% Revision 1.1  2000/05/01 15:42:46  schmieta
% Initial revision
%
%
n = 2^k;
if n > flintmax
  error('Can''t generate Hamming graph. k to large');
end
% Generate bit patterns for neighbors with distances in d
bitpat = [];
for i = d
  b = bitset(0,nchoosek(1:k,i)',1);
  if size(b,1) ~= 1
    b = sum(b);
  end
  bitpat = [bitpat, b];
end
Ai = [];
Aj = [];
C = sparse(n,n);
start = 1;
for i = 0:(n-1)
  neighbors = bitxor(i, bitpat);
  neighbors = neighbors(find(neighbors > i));
  if length(neighbors) > 0
    C(i+1,neighbors + 1) = 1;
    C(neighbors+1,i+1) = 1;
    % Add a row to A for each neighbor
    rows = start:(start+length(neighbors)-1);
    Ai = [Ai; rows'; rows'];
    Aj = [Aj; neighbors'*n + i + 1];
    Aj = [Aj; neighbors' + i*n + 1];
    start = start + length(neighbors);
  end
end
A = [sparse(vec(speye(n,n))'); sparse(Ai, Aj, ones(length(Aj),1),start-1,n^2)];
c = sparse(-vec(1 - C));
b = zeros(size(A,1),1);
b(1) = 1;
b = sparse(b);
K.s = n;
end