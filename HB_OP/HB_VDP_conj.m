% function eqn=HB_VDP_conj(n)

% HB for polynomial operational matrix
% example: VDP
% v2, add several assumptions
% conjugate unknowns version

syms w pr0; 
n=9;   % single side, m=2n+1 harmonics in all
n_lower=5;   % all to n, lower terms to n_lower
m=2*n+1;
D=1j*w*Dm(n);   % build differential operational matrix
Pr=sym('pr',[1 n]);  % unknown variables real part
Pi=sym('pi',[1 n]);  % unknown variables imag part

Pr=sym(Pr,'real');Pi=sym(Pi,'real');pr0=sym(pr0,'real');w=sym(w,'real');  % all variables are real
Pi(1)=0;Pi(2)=0;pr0=0;Pr(2)=0;%Pr(4)=0;Pi(4)=0;   % manually set these variables to be zero !!!! added in v2
% Pi(5)=0;Pr(5)=0; % eqv test
% Pi(6)=0;Pr(6)=0; % eqv test 
% Pi(7)=0;Pr(7)=0; % eqv test
% P=[fliplr(Pr) pr0 Pr] +1j*[-fliplr(Pi) 0 Pi];       % unknown tones in complex representation
P1=sym('cp',[1 n]);       % unknown tones normal representation
P2=sym('cn',[1 n]);       % unknown tones normal representation conjugate
P1(2:2:end)=0;P2(2:2:end)=0; % even order coef. all 0
P=[fliplr(P1),0,P2];

n_lower=7;
P1(end-n+n_lower+1:end)=0;P2(end-n+n_lower+1:end)=0; % higher order coef. all 0
P_lower7=[fliplr(P1),0,P2];   % lower order replacement
n_lower=5;
P1(end-n+n_lower+1:end)=0;P2(end-n+n_lower+1:end)=0; % higher order coef. all 0
P_lower5=[fliplr(P1),0,P2];   % lower order replacement
n_lower=3;
P1(end-n+n_lower+1:end)=0;P2(end-n+n_lower+1:end)=0; % higher order coef. all 0
P_lower3=[fliplr(P1),0,P2];  % lower order replacement
n_lower=1;
P1(end-n+n_lower+1:end)=0;P2(end-n+n_lower+1:end)=0; % higher order coef. all 0
P_lower1=[fliplr(P1),0,P2];   % lower order replacement

% VDP
% x''-x'+x+3x'x^2=0
Pm=P*D*D-P*D+P; % (x''-x'+x) term, under Phi(m) basis
% Pm=[zeros(1,2*n), Pm, zeros(1,2*n)];   % extend Pm to Phi(3m-2) fourier basis
P3m=3*kron(kron(P,P)*S(m,m),P*D)*S(2*m-1,m);  % (3x'x^2) term, under Phi(3m-2) basis
P3m=P3m(2*n+1:4*n+1);  % shrink P3m to Phi(m) fourier basis


Pmall=Pm+P3m;
%%%%% conjugate representation


eqn_all=expand(Pmall(1:2:end).');   % even order eqn all 0, omit
eqn_all=flipud(eqn_all(1:(n+1)/2));

eqn_lower7=expand(subs(eqn_all,P,P_lower7)); % replace by lower order only
eqn_lower5=expand(subs(eqn_all,P,P_lower5)); % replace by lower order only
eqn_lower3=expand(subs(eqn_all,P,P_lower3)); % replace by lower order only
eqn_lower1=expand(subs(eqn_all,P,P_lower1)); % replace by lower order only

eqn_upper7=expand(eqn_all - eqn_lower7);
eqn_upper5=expand(eqn_lower7 - eqn_lower5);
eqn_upper3=expand(eqn_lower5 - eqn_lower3);
eqn_upper1=expand(eqn_lower3 - eqn_lower1);