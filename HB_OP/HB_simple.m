% original HB
% for x'+x-x*x-coswt
global b;
b=3;
n=3; % 3 tones

w=2*pi;whb=w;
X=zeros(n,1);
U=zeros(n,1);U(floor(n/2))=1/2;U(floor(n/2)+2)=1/2;
Omega=1j*w*diag([-floor(n/2):floor(n/2)]);

idftM=zeros(n,n);
for i=-floor(n/2):floor(n/2),
    for k=-floor(n/2):floor(n/2),
       idftM(k+ceil(n/2),i+ceil(n/2))=exp(1j*2*pi*i/n*k)/n; 
    end
end
dftM=idftM'*n;


Fnew=Omega*X+dftM*(idftM*X-1*(idftM*X).^2)-b*U;
Xnew=X-(Omega+eye(n))\Fnew;
while (norm(Fnew)>1e-7)
    X=Xnew;
    Fnew=Omega*X+dftM*(idftM*X-1*(idftM*X).^2)-b*U;
    Xnew=X-(Omega+eye(n))\Fnew;
end


Xhb=X;
%---------td

f=1;
ntone=100; % 100 sample frequency
fs=f/ntone;
t=[0:fs:100*f-fs];
% w=2*pi/f;
% alpha=0.1;
[T,Y1] = ode23(@(t,y) simple_ode(t,y) ,t,0);   % already PSS
 

% figure;hold on;grid on;
plot(Y1,'b');

xt=Y1(end-ntone+1:end);
plot(Y1(end-ntone:end));
ln=length(xt);
y_spectrum=fft(xt,ln)/ln;
tones=y_spectrum;