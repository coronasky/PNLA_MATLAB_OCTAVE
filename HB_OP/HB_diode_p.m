% HB for polynomial operational matrix
% example: Diode example w/ state immersion
% x'=-x-exp(4x)+1+u u=cos(wt)
% Forced system, NO phase assumption
clear all;
syms w; 
w=2*pi;
n=2;   % single side, m=2n+1 harmonics in all
nsp=2; % # state space = 2;
m=2*n+1;
D=1j*w*Dm(n);   % build differential operational matrix
P0=sym('cr%d_0',[nsp 1]);
Pr=sym('cr%d_%d',[nsp n]);  % unknown variables real part      cr3_2: second harmonic for 3rd variable
Pi=sym('ci%d_%d',[nsp n]);  % unknown variables imag part

Pr=sym(Pr,'real');Pi=sym(Pi,'real');P0=sym(P0,'real');  % all variables are real
P=[fliplr(Pr) P0 Pr] +1j*[-fliplr(Pi) zeros(nsp,1) Pi];       % unknown tones in complex representation


% Diode example w/ state immersion
% x'+x+exp(4x)-1-u=0 u=cos(wt)
% G1=-[1 1;4 4];
% G2=-[0 0 0 0; 0 4 0 4];
% D1=[0 0;0 4];
% B=[1;4];
G1=[-1 1;4 -4];
G2=[0 0 0 0; 0 -4 0 4];
D1=[0 0;0 4];
B=[1;-4];
U=zeros(1,m);U(n)=1/2;U(n+2)=1/2; % input harmonic
% one=zeros(m,1);one(n+1)=1;    % constant '1'

Pm=-P*D+G1*P; % (-x'+G1x) term, under Phi(m) basis
x2p=G2*kron(P,P)*S(m,m);  % -G2x*x term, under Phi(2m-1) basis
x2p=x2p(:,n+1:3*n+1);  % shrink x2p to Phi(m) fourier basis
xdp=D1*kron(P,U)*S(m,m);  % -D1x*u term, under Phi(2m-1) basis
xdp=xdp(:,n+1:3*n+1);  % shrink xdp to Phi(m) fourier basis



Pmall=Pm+x2p+xdp+B*U;
%%%%% complex representation
sr=real(expand(Pmall));       % get real part of equations
si=imag(expand(Pmall));       % get imag part of equations
sr=sr(:,1:n+1);si=si(:,1:n);    % truncate duplicated equations Phi(m) basis


eqn_all=[reshape(sr,[],1);reshape(si,[],1)];
vpa(eqn_all,11)    % generate eqns for our solver
% eqn_all*1000    % generate eqns for Maple (Grobner basis)

syms p;
var=[Pr(1),w,p];