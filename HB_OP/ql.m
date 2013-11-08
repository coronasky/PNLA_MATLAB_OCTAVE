% VDP fft test
% clc;
clear;
f=6.6655;
ntone=100; % 100 sample frequency
fs=f/ntone;
t=[0:fs:f-fs];
w=2*pi/f;
% alpha=0.1;
[T,Y1] = ode23(@(t,y) ql_ode(t,y) ,t,[0,1.2546]);   % already PSS
 

figure;hold on;grid on;
plot(Y1(:,2),Y1(:,1),'b');

xt=Y1(:,1);

ln=length(xt);
y_spectrum=fft(xt,ln);
% fprintf('x0 of original system is %f\n',abs(y_spectrum(1)/ln));
% fprintf('x1 of original system is %f\n',abs(y_spectrum(2)/ln));
% fprintf('x2 of original system is %f\n',abs(y_spectrum(3)/ln));
% fprintf('x3 of original system is %f\n',abs(y_spectrum(4)/ln));
tones5=fft(xt(1:20:end),5)/5;

tones=y_spectrum/ln;
y_spectrum(end-2:end)/ln;

% test

x=[-0.0132 + 0.0042i, -0.0521 - 0.0444i, 0.0662 - 0.5776i, 0.0662 + 0.5776i, -0.0521 + 0.0444i, -0.0132 - 0.0042i]*[exp(5*1j*w*t);exp(3*1j*w*t);exp(1j*w*t);exp(-1j*w*t);exp(-3*1j*w*t);exp(-5*1j*w*t)];

% plot(abs(xt-x'));
% diff=xt-x';
% df=fft(diff)/ln;