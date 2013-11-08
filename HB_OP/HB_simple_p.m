% HB for polynomial operational matrix
% example: simple x'+x-x*x-coswt=0;
global b;

syms w pr0; 
n=3;   % single side, m=2n+1 harmonics in all
m=2*n+1;
D=1j*w*Dm(n);   % build differential operational matrix
Pr=sym('pr',[1 n]);  % unknown variables real part
Pi=sym('pi',[1 n]);  % unknown variables imag part
Pr=sym(Pr,'real');Pi=sym(Pi,'real');pr0=sym(pr0,'real');%w=sym(w,'real');
% P=[fliplr(Pr) pr0 Pr] +1j*[-fliplr(Pi) 0 Pi];       % unknown tones in complex representation
P=sym('p',[1 m]);       % unknown tones normal representation

% simple
% x'+x-x*x-coswt=0;
U=zeros(m,1);U(n)=1/2;U(n+2)=1/2;
PD=P*D;
Pm=PD+P-b*U.'; % x'+x-coswt term, Phi(m) basis
Pm=[zeros(1,1*n), Pm, zeros(1,1*n)];   % extend Pm to Phi(2m-1) fourier basis
P3m=-1*kron(P,P)*S(m,m);  % x^2 term, Phi(2m-1) basis
Pm=Pm.';P3m=P3m.';PD=PD.';
Pmall=Pm+P3m;


%%%%% compared with HB
re=tones;
repl=[whb fliplr(re(2:n+1)') re(1) re(2:n+1).'];
% repl=[whb Xhb.'];
diff=subs(Pmall,[w P],repl);
% diff1=subs(PD,[w P],repl)
max(abs(diff));
Pmall;diff