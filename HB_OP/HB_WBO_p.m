% HB for polynomial operational matrix
% example: WBO
% v2, add several assumptions

syms w pr0; 
n=2;   % single side, m=2n+1 harmonics in all
n_lower=3;
m=2*n+1;
D=1j*w*Dm(n);   % build differential operational matrix
Pr=sym('pr',[1 n]);  % unknown variables real part
Pi=sym('pi',[1 n]);  % unknown variables imag part

Pr=sym(Pr,'real');Pi=sym(Pi,'real');pr0=sym(pr0,'real');w=sym(w,'real');  % all variables are real
Pi(1)=0;pr0=0;%Pr(4)=0;Pi(4)=0;   % manually set these variables to be zero !!!! added in v2
% for i=2:2:n,
%     Pi(i)=0;
%     Pr(i)=0;
% end
Pi(2:2:end)=0;Pr(2:2:end)=0; % even order coef. all 0
P=[fliplr(Pr) pr0 Pr] +1j*[-fliplr(Pi) 0 Pi];       % unknown tones in complex representation


% WBO
% x''-.234x'+x+6.585x'x^2-3.334x'x^4=0
Pm=P*D*D-0.234*P*D+P; % (x''-0.234x'+x) term, under Phi(m) basis

x2p=kron(kron(P,P)*S(m,m),P*D)*S(2*m-1,m);
x2p=x2p(2*n+1:4*n+1);  % shrink P3m to Phi(m) fourier basis
P3m=6.585*x2p;         % (6.585x'x^2) term, under Phi(3m-2) basis
x4p=kron(kron(P,P)*S(m,m),x2p)*S(2*m-1,m);
x4p=x4p(2*n+1:4*n+1);  % shrink P3m to Phi(m) fourier basis
P5m=-3.334*x4p;

Pmall=Pm+P3m+P5m;
%%%%% complex representation
sr=real(expand(Pmall.'));       % get real part of equations
si=imag(expand(Pmall.'));       % get imag part of equations
% sr=sr(1:n+1);si=si(1:n);        % truncate duplicated equations Phi(m) basis

eqn2=[sr;si];
% for Pl and Pu test
P_high=[Pi,Pr];
n_lower=7;
Pi(end-n+n_lower+1:end)=0;Pr(end-n+n_lower+1:end)=0; % higher order coef. all 0
P_low7=[Pi,Pr];   % lower order replacement
n_lower=5;
Pi(end-n+n_lower+1:end)=0;Pr(end-n+n_lower+1:end)=0; % higher order coef. all 0
P_low5=[Pi,Pr];   % lower order replacement
n_lower=3;
Pi(end-n+n_lower+1:end)=0;Pr(end-n+n_lower+1:end)=0; % higher order coef. all 0
P_low3=[Pi,Pr];   % lower order replacement
n_lower=1;
Pi(end-n+n_lower+1:end)=0;Pr(end-n+n_lower+1:end)=0; % higher order coef. all 0
P_low1=[Pi,Pr];   % lower order replacement

sr1=sr(2-mod(n,2):2:n);
si1=si(2-mod(n,2):2:n);

% eqn_all=[sr1;si1];
eqn_all=sym(zeros(n+1,1));
eqn_all(1:2:end)=flipud(sr1); 
eqn_all(2:2:end)=flipud(si1);
eqn_lower7=expand(subs(eqn_all,P_high,P_low7)); % replace by lower order only
eqn_lower5=expand(subs(eqn_all,P_high,P_low5)); % replace by lower order only
eqn_lower3=expand(subs(eqn_all,P_high,P_low3)); % replace by lower order only
eqn_lower1=expand(subs(eqn_all,P_high,P_low1)); % replace by lower order only
eqn_upper7=expand(eqn_all - eqn_lower7);
eqn_upper5=expand(eqn_lower7-eqn_lower5);
eqn_upper3=expand(eqn_lower5-eqn_lower3);
eqn_upper1=expand(eqn_lower3-eqn_lower1);
eqn_upper0=eqn_lower1;

% P_low
syms p;
var=[Pr(1),w,p];
% eqn_lower=[eqn_lower;p*Pr(1)-1]
% num=solve(eqn_lower,var)