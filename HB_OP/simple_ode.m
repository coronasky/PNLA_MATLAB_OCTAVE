function dy = simple_ode(t,y)
global b;

dy=-y+1*y^2+b*cos(2*pi*t);
