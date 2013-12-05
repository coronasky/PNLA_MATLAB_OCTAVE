% test file of rcandecomp
clear all;
addpath('C:\matlab\R2012a\toolbox\suiteSparse\SPQR\maTLAB\')
% http://www.math.uic.edu/~jan/Demo/quad or2
% Jan Verschelde and Karin Gatermann:
%`Symmetric Newton Polytopes for Solving Sparse Polynomial Systems',
% %Adv. Appl. Math., 16(1): 95-127, 1995.
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


% % Haotian's equations, N=5, simplified, forced pr1 != 0
% strings{1}='-30*a[1]*a[4]^2*a[6]-15*a[1]*a[2]^2*a[4]-30*a[1]*a[2]*a[3]*a[4]-15*a[1]*a[6]^3-30*a[1]*a[2]^2*a[6]-30*a[1]*a[3]^2*a[6]-15*a[1]*a[5]^2*a[6]+5*a[1]*a[6]-25*a[1]^2*a[5]+a[5]';
% strings{2}='-9*a[1]*a[4]^3-18*a[1]*a[4]*a[6]^2-18*a[1]*a[2]^2*a[4]+18*a[1]*a[2]*a[4]*a[5]-9*a[1]*a[3]^2*a[4]-18*a[1]*a[4]*a[5]^2+3*a[1]*a[4]-9*a[1]*a[2]^2*a[6]-18*a[1]*a[2]*a[3]*a[6]-9*a[1]^2*a[3]+a[3]';
% strings{3}='-3*a[1]*a[4]^2*a[6]-3*a[1]*a[2]^2*a[4]+6*a[1]*a[2]*a[4]*a[5]-6*a[1]*a[3]*a[4]*a[5]-6*a[1]*a[2]*a[3]*a[6]-a[1]^2*a[2]+a[2]+3*a[1]*a[3]^2*a[6]';
% strings{4}='15*a[1]*a[2]*a[4]^2-30*a[1]*a[4]^2*a[5]-15*a[1]*a[5]*a[6]^2+25*a[1]^2*a[6]-a[6]-15*a[1]*a[2]^2*a[3]-30*a[1]*a[2]^2*a[5]-15*a[1]*a[2]*a[3]^2-30*a[1]*a[3]^2*a[5]-15*a[1]*a[5]^3+5*a[1]*a[5]';
% strings{5}='-9*a[1]*a[3]*a[4]^2-18*a[1]*a[2]*a[4]*a[6]+9*a[1]^2*a[4]-a[4]-18*a[1]*a[3]*a[6]^2-3*a[1]*a[2]^3-18*a[1]*a[2]^2*a[3]-9*a[1]*a[2]^2*a[5]-18*a[1]*a[2]*a[3]*a[5]-9*a[1]*a[3]^3-18*a[1]*a[3]*a[5]^2+3*a[1]*a[3]';
% strings{6}='6*a[2]*a[4]^2-3*a[4]^2*a[5]+6*a[2]*a[4]*a[6]+6*a[3]*a[4]*a[6]+6*a[2]*a[6]^2+3*a[2]^3+3*a[2]^2*a[3]+6*a[2]*a[3]^2+6*a[2]*a[3]*a[5]+6*a[2]*a[5]^2-a[2]+3*a[3]^2*a[5]';
% strings{7}='a[2]*a[7]-1';   
% polysys=lti2polysys(strings,[],[],'na',7);
%=================== begin test
% getD0(polysys)
% getDorig(polysys)
% tic
% M=getM(polysys,getD0(polysys)); 
% N=sparse(null(M));
% c2(getD0(polysys))=size(N,2);
for i=getD0(polysys)+1:getD0(polysys)+4,
    i
    N=updateN(N,getMex(polysys,i,i-1,1),1);   % sparse qr
    c2(i)=size(N,2);
end,toc
% c2
% see roots
% qdsparf(polysys)
% [root, d, c, ns, check, cr, digits] = qdsparf(polysys)


% [a, b, R] = candecomp(polysys,4);
% [a_, b_, R_] = rcandecomp(polysys,6);
[a_, b_, R_] = rcandecomp(polysys,5);      % inc order will make more extra eqns, eliminate more inf roots.
% update=vec2polysys(R_,7);     % n variables
update=vec2polysys(R_,4);
% polysys=[polysys;update(1,:)];   % add condecomp together
% polysys=[polysys;update];

getD0(polysys)
getDorig(polysys)

M=getM(polysys,getD0(polysys));
figure;spy(M);
N=sparse(null(M));
c2(getD0(polysys))=size(N,2);
for i=getD0(polysys)+1:getD0(polysys)+5,
    tic
    N=updateN(N,getMex(polysys,i,i-1,1),1);   % sparse qr
    c2(i)=size(N,2);
%     [i,toc]
end
c2
% [root, d, c, ns, check, cr, digits] = qdsparf(polysys)

