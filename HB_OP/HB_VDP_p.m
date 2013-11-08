% HB for polynomial operational matrix
% example: VDP

syms w pr0; 
n=2;   % single side, m=2n+1 harmonics in all
m=2*n+1;
D=1j*w*Dm(n);   % build differential operational matrix
Pr=sym('pr',[1 n]);  % unknown variables real part
Pi=sym('pi',[1 n]);  % unknown variables imag part
Pr=sym(Pr,'real');Pi=sym(Pi,'real');pr0=sym(pr0,'real');w=sym(w,'real');  % all variables are real
P=[fliplr(Pr) pr0 Pr] +1j*[-fliplr(Pi) 0 Pi];       % unknown tones in complex representation


% VDP
% x''-x'+x+3x'x^2=0
Pm=P*D*D-P*D+P; % (x''-x'+x) term, under Phi(m) basis
Pm=[zeros(1,2*n), Pm, zeros(1,2*n)];   % extend Pm to Phi(3m-2) fourier basis
P3m=3*kron(kron(P,P)*S(m,m),P*D)*S(2*m-1,m);  % (3x'x^2) term, under Phi(3m-2) basis

Pmall=Pm+P3m;
%%%%% complex representation
sr=real(expand(Pmall.'));       % get real part of equations
si=imag(expand(Pmall.'));       % get imag part of equations
sr=sr(1:3*n+1);si=si(1:3*n);  % truncate duplicated equations

eqn=[sr;si];
