% test file of rcandecomp
clear all;
addpath('C:\matlab\R2012a\toolbox\suiteSparse\SPQR\maTLAB\')
% http://www.math.uic.edu/~jan/Demo/quadfor2
% Jan Verschelde and Karin Gatermann:
%`Symmetric Newton Polytopes for Solving Sparse Polynomial Systems',
%Adv. Appl. Math., 16(1): 95-127, 1995.
% [x1 x2 w1 w2]
%  w1 + w2 - 1;
polysys{1,1} = [1 1 -1];
polysys{1,2} = [0 0 1 0;0 0 0 1; 0 0 0 0];
%  w1*x1 + w2*x2;
polysys{2,1} = [1 1];
polysys{2,2} = [1 0 1 0;0 1 0 1];
%  w1*x1**2 + w2*x2**2 - 2/3;
polysys{3,1} = [1 1 -2/3];
polysys{3,2} = [2 0 1 0;0 2 0 1;0 0 0 0];
%  w1*x1**3 + w2*x2**3;
polysys{4,1} = [1 1];
polysys{4,2} = [3 0 1 0;0 3 0 1];

%=================== begin test
getD0(polysys)
getDorig(polysys)
tic
M=getM(polysys,getD0(polysys));
N=sparse(null(M));
c2(getD0(polysys))=size(N,2);
for i=getD0(polysys)+1:getD0(polysys)+8,
    N=updateN(N,getMex(polysys,i,i-1,1),1);   % sparse qr
    c2(i)=size(N,2);
end,toc
c2
% see roots
% qdsparf(polysys)
[root d c ns check cr digits] = qdsparf(polysys)
[a, b, R] = candecomp(polysys,4);
[a_, b_, R_] = rcandecomp(polysys,4)

update=vec2polysys(R_,4);
polysys=[polysys;update];

getD0(polysys)
getDorig(polysys)
tic
M=getM(polysys,getD0(polysys));
N=sparse(null(M));
c2(getD0(polysys))=size(N,2);
for i=getD0(polysys)+1:getD0(polysys)+8,
    N=updateN(N,getMex(polysys,i,i-1,1),1);   % sparse qr
    c2(i)=size(N,2);
end,toc
c2
