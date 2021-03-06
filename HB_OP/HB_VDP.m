% HB for polynomial operational matrix
% example: VDP

syms w pr0; 
n=1;   % single side, m=2n+1 harmonics in all
m=2*n+1;
D=1j*w*Dm(n);   % build differential operational matrix
Pr=sym('pr',[1 n]);  % unknown variables real part
Pi=sym('pi',[1 n]);  % unknown variables imag part
Pr=sym(Pr,'real');Pi=sym(Pi,'real');pr0=sym(pr0,'real');w=sym(w,'real');
P=[fliplr(Pr) pr0 Pr] +1j*[-fliplr(Pi) 0 Pi];       % unknown tones in complex representation
P=sym('p',[1 m]);       % unknown tones normal representation

% VDP
% x''-x'+x+3x'x^2=0
Pm=P*D*D-P*D+P; % x''-x'+x term, , Phi(m) basis
Pm=[zeros(1,2*n), Pm, zeros(1,2*n)];   % extend Pm to Phi(3m-2) fourier basis
P3m=3*kron(kron(P,P)*S(m,m),P*D)*S(2*m-1,m);  % 3x'x^2 term, Phi(3m-2) basis
Pm=Pm.';P3m=P3m.';
Pmall=Pm+P3m;


%%%% test (don't comment!!!!!)
re=tones;

%%%%% complex representation
% sr=real(expand(Pmall));
% si=imag(expand(Pmall));
% sr=sr(1:3*n+1);si=si(1:3*n+1);   % truncate duplicated equations
% %%%%%%%% test for time domain solution
% repl=[2*pi/f real(re(1:n+1))' imag(re(2:n+1))']
% err_re=subs(sr,[w pr0 Pr Pi],repl);
% err_im=subs(si,[w pr0 Pr Pi],repl);
% diff=[err_re;err_im];
% max(abs(diff))

%%%%% normal representation
% repl=[2*pi/f fliplr(re(2:n+1)') re(1) re(2:n+1).'];
%%%%% compared with HB
repl=[whb Xhb'];
diff=subs(Pmall,[w P],repl);
max(abs(diff))
Pmall;diff