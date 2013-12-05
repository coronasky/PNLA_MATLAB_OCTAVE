% HB for polynomial operational matrix
% example: bioch example w/ state immersion
% x'=-x-x/(1+10x)-u u=cos(wt)
% Forced system, NO phase assumption
% v2: add multivariable x to make a larger system

clear all;
syms w; 
w=2*pi;
n=1;   % single side, m=2n+1 harmonics in all
nsp=3; % # state space = 2+2 extra linear variables
m=2*n+1;
D=1j*w*Dm(n);   % build differential operational matrix
P0=sym('cr%d_0',[nsp 1]);
Pr=sym('cr%d_%d',[nsp n]);  % unknown variables real part      e.g., cr3_2: second harmonic for 3rd variable
Pi=sym('ci%d_%d',[nsp n]);  % unknown variables imag part

Pr=sym(Pr,'real');Pi=sym(Pi,'real');P0=sym(P0,'real');  % all variables are real
P=[fliplr(Pr) P0 Pr] +1j*[-fliplr(Pi) zeros(nsp,1) Pi];       % unknown tones in complex representation


% bioch example w/ state immersion

C0=eye(nsp);C0(nsp,nsp)=0;  % etra one
G1=[1 0 1;-1 1 0;-1 0 1];
% G2=[0 0 0 0; 0 10 0 0];
G2=zeros(nsp,nsp*nsp);G2(nsp,nsp)=10;
B=[1;0;0];
U=zeros(1,m);U(n)=1/2;U(n+2)=1/2; % input harmonic


Pm=C0*P*D+G1*P+B*U; % (C0x'+G1x+Bu) term, under Phi(m) basis
x2p=G2*kron(P,P)*S(m,m);  % G2x*x term, under Phi(2m-1) basis
x2p=x2p(:,n+1:3*n+1);  % shrink x2p to Phi(m) fourier basis




Pmall=Pm+x2p;
%%%%% complex representation
sr=real(expand(Pmall));       % get real part of equations
si=imag(expand(Pmall));       % get imag part of equations
sr=sr(:,1:n+1);si=si(:,1:n);    % truncate duplicated equations Phi(m) basis


eqn_all=[reshape(sr,[],1);reshape(si,[],1)];
vpa(eqn_all,11)    % generate eqns for our solver
% eqn_all*1000    % generate eqns for Maple (Grobner basis)

syms p;
var=[Pr(1),w,p];