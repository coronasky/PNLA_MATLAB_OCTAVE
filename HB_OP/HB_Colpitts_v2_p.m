% HB for polynomial operational matrix
% example: Colpitts ocillator
% assumption: all ci0=0, ci1_1=0
% v2 Taylor expansion instead of state immersion

syms w; 
n=1;   % single side, m=2n+1 harmonics in all
nsp=3; % # state space = 3;
n_lower=3;
m=2*n+1;
D=1j*w*Dm(n);   % build differential operational matrix
P0=sym('cr%d_0',[nsp 1]);
Pr=sym('cr%d_%d',[nsp n]);  % unknown variables real part      cr3_2: second harmonic for 3rd variable
Pi=sym('ci%d_%d',[nsp n]);  % unknown variables imag part

Pr=sym(Pr,'real');Pi=sym(Pi,'real');P0=sym(P0,'real');w=sym(w,'real');  % all variables are real
Pi(1,1)=0;  % manually set ci1_1=0 to fix the phase
P=[fliplr(Pr) P0 Pr] +1j*[-fliplr(Pi) zeros(nsp,1) Pi];       % unknown tones in complex representation


% Colpitts Oscillator
% x'=G1x+G2x*x
G1=[0 1.5 1.5;0 0 1.5;-1/3 -1/3 -0.5];
G2=[0 0 0 0 -3/4 0 0 0 0;...
    0 0 0 0 0 0 0 0 0;...
    0 0 0 0 0 0 0 0 0];
 

Pm=P*D-G1*P; % (x'-G1x) term, under Phi(m) basis
x2p=-G2*kron(P,P)*S(m,m);  % -G2x*x term, under Phi(2m-1) basis
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
% eqn_lower=[eqn_lower;p*Pr(1)-1]
% num=solve(eqn_lower,var)