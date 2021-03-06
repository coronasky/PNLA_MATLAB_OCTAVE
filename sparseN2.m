function [N, nM, tol]=sparseN(polysys,d)
% [N, nM, tol]=sparseN(polysys,d)
% -------------------------------
% Computes a sparse basis for the null space of the Macaulay matrix for a
% given polynomial system polysys and degree d.
%
% N             =   matrix, sparse basis for null space of M(d)
%
% nM            =   vector, 
%
% tol           =   scalar, the tolerance used for the rank-test,
%                   max(size(M(d))) * eps(s), where s is an upper bound on
%                   the largest singular value of M(d)
%
% polysys       =   cell containing coefficients and monomials exponents of the
%                   set of polynomial equations.
%
% d             =	scalar, degree for which M(d) has to be constructed from polysys
%
% CALLS
% -----
%
% getM.m
%
% Kim Batselier, december 2013

s=size(polysys,1);
n=size(polysys{1,2},2);

% normalize each coefficient vector to unit 1-norm such that the largest
% singular value of M(d) is upper bounded by sqrt(s)
for i=1:s
    polysys{i,1}=polysys{i,1}/norm(polysys{i,1},1);
end

% if M is sparse column compressed than accessing columns
% is faster than rows
initialM=getM(polysys,d,1); 
[p, q]=size(initialM);

tol=max(p,q)*eps(sqrt(s));
nM=zeros(1,p);
nM(1)=norm(initialM(1,:));

% indices of nonzero coefficients
nindex=find(initialM(1,:));
% indices of zero coefficients
zindex=setdiff(1:q,nindex);

%% initalize sparse N, first time we need to make it explicit
% column indices
colI=zeros(1,length(zindex)+2*(length(nindex)-1));
colI(1:length(zindex))=1:length(zindex); % need to finish remaining indices

% row indices
rowI=zeros(1,length(zindex)+2*(length(nindex)-1));
rowI(1:length(zindex))=zindex; % need to finish the remaining indices

% values
values=zeros(1,length(zindex)+2*(length(nindex)-1));
values(1:length(zindex))=ones(1,length(zindex));

pivotI=nindex(find(abs(initialM(1,nindex))==max(abs(initialM(1,nindex))),1));
initialM(1,:)=initialM(1,:)/initialM(1,pivotI);
tempI=setdiff(nindex,pivotI);

% make the new M
M=spalloc(p-1,q-1,nnz(initialM));
M(:,1:length(zindex))=initialM(2:end,zindex);

% finish M and remaining indices/values for initialN
for j=1:length(tempI)
    % remaining columns of M
    M(:,length(zindex)+j)=initialM(1,tempI(j))*initialM(2:end,pivotI)-initialM(1,pivotI)*initialM(2:end,tempI(j));
    
    % remaining indices that determine new N
    colI(length(zindex)+2*(j-1)+1:length(zindex)+2*j)=(length(zindex)+j)*ones(1,2);
    rowI(length(zindex)+2*(j-1)+1:length(zindex)+2*j)=[pivotI tempI(j)];
    values(length(zindex)+2*(j-1)+1:length(zindex)+2*j)=[initialM(1,tempI(j)) -initialM(1,pivotI)];    
end

% now we can make the first N
N=sparse(rowI,colI,values,q,q-1,length(zindex)+2*(length(nindex)-1));
[p,q]=size(M);

for i=1:p
    
    nM(i+1)=norm(M(i,:));
    if nM(i+1) > tol
        nindex=find(M(i,:));
        zindex=setdiff(1:q,nindex);
        
        pivotI=nindex(find(abs(M(i,nindex))==max(abs(M(i,nindex))),1));
        tempI=setdiff(nindex,pivotI);
        M(i,:)=M(i,:)/M(i,pivotI);
          
        tempM=M;
        M=spalloc(p,q-1,nnz(tempM));
        M(:,1:length(zindex))=tempM(:,zindex);
        
        tempN=N;
        N=spalloc(size(N,1),q-1,nnz(tempN));
        N(:,1:length(zindex))=tempN(:,zindex);
        
        for j=1:length(tempI)
            % remaining columns of M            
            M(:,length(zindex)+j)=tempM(i,tempI(j))*tempM(1:end,pivotI)-tempM(i,pivotI)*tempM(1:end,tempI(j));
            
            N(:,length(zindex)+j)=tempM(i,tempI(j))*tempN(:,pivotI)-tempM(i,pivotI)*tempN(:,tempI(j));    
        end
        [p, q]=size(M);
%     else
       % first row of M is numerically zero, remove it
%        M=M(2:end,:);
%        [p,q]=size(M); 
    end
end