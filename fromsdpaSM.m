function [At,b,c,n,m]=fromsdpaSM(fname)
% Reads SDP problem of Single Matrix from sparse SDPA-formatted input file.
% 

fid = fopen(fname, 'r');
if fid == -1
    error('File not found.')
end
%The number of equality constraints is the first number in the file
m='';
while(isempty(m))
    m = fscanf(fid,'%d',1);
    if fgetl(fid) == -1;
        fclose(fid);
        error('Invalid SDPA file. Number of equality constraints not found.')
    end
end
%The number of semidefinite blocks
nblocks = fscanf(fid, '%d', 1);
if isempty(nblocks)
    fclose(fid);
    error('Invalid SDPA file. Number of semidefinite blocks not found.')
end
fgetl(fid);
%The dimension of the blocks
%Negative means diagonal block, we convert that to nonnegative variables
%.,(){} are omitted
dims = sscanf(regexprep(fgetl(fid), '[\.\,(){}]', ' '), '%d', nblocks)';
%Dimensions cannot be 0, and their number must be nblocks
if any(dims == 0) || length(dims) ~= nblocks
    fclose(fid);
    error('Invalid SDPA file. Invalid semidefinite block dimensions.')
end
%nblocks = length(dims);
cdims = cumsum(dims);
n = cdims(end);
cdims = [0, cdims(1:end-1)];
N = n^2;

%Vector b
%,(){} are omitted
b = sscanf(regexprep(fgetl(fid),'[\,(){}]',' '),'%f',m);
if length(b) ~= m
    fclose(fid);
    error('Invalid SDPA file. The right-hand side vector is not of the right dimension.')
end
%If b is very sparse then we store it as sparse
if nnz(b)/m < 0.1
    b=sparse(b);
end

%Coefficients
%It is much faster to get all the numbers at once than to read the file
%line by line
try
    E = fscanf(fid,'%d %d %d %d %f',[5. inf]);
catch
    fclose(fid);
    error('Invalid SDPA file. Error reading the coefficients.')
end
fclose(fid);
%We are done with the file

%Extract the objective
cE = E(:,E(1,:)==0);

%Repeating indices in the sparse matrix structure create the sum of
%elements, we need to clear one for the diagonals
data2 = cE(5,:);
data2(cE(3,:)==cE(4,:)) = 0;
%we need the minus sign because the SDPA format assumes maximization while
%SeDuMi uses minimization by default
%This magic reshuffles the coefficients to the right place
c = -sparse([(cE(3,:)+cdims(cE(2,:))-1)*n+(cE(4,:)+cdims(cE(2,:))),(cE(4,:)+cdims(cE(2,:))-1)*n+(cE(3,:)+cdims(cE(2,:)))],...
    1,...
    [cE(5,:),data2],...
    N,1);
clear cE data2


%Get rid of the objective coefficients from E
AtE = E(:,E(1,:)~=0);
clear E
%Take all the coefficients
data2 = AtE(5,:);
%The coefficients for the diagonal elements in the blocks
data2(AtE(3,:)==AtE(4,:)) = 0;
At = sparse([((AtE(3,:)-1)+cdims(AtE(2,:)))*n+AtE(4,:)+cdims(AtE(2,:)),((AtE(4,:)-1)+cdims(AtE(2,:)))*n+AtE(3,:)+cdims(AtE(2,:))],...
    [AtE(1,:),AtE(1,:)],...
    [AtE(5,:),data2],...
    N,m);
clear AtE data2

%Finally, let us correct the sparsity
%c and At are always stored as sparse
[i,j,v] = find(c);
c = sparse(i,1,v,N,1);
[i,j,v] = find(At);
At = sparse(i,j,v,N,m);