clc;
t=[0:0.01:100];


[T,Y1] = ode45(@WBO_ode ,t,[0.1;0]);


figure;hold on;grid on;% 
% plot(Y1(:,2),Y1(:,1),'LineWidth',1);
plot(t,Y1(:,1),'LineWidth',1);


% t=[-1.5:0.1:1.5];
% plot(t,-t+t.^3);
% plot(t,-t.*(-t+t.^ 3));