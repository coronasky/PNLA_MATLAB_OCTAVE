% Jacobian of polynomial equations
syms w pr1 pr3 pi3
vars=[w,pr1, pr3, pi3];
eqn=eqn_all;
% eqn(2)=eqn(2)*100;
% eqn(1)=eqn_all(3);eqn(2)=eqn_all(6);eqn(3)=eqn_all(2);eqn(4)=eqn_all(5);eqn(5)=eqn_all(1);eqn(6)=eqn_all(4);
Jac=jacobian(eqn,vars);

% load('save_Jacob')

w = 0.94294236004556073569267041920248;
c1 = 0.58163632516652939737707591650654;
c3r = -0.024523254236563325095870940201166;
c3i = 0.064145783801961427976067869239161;
% c5i = -0.009478181088859126891376836651277;
% c5r = -0.011010754430848876951384997121712;


init=[w,c1,0,0];
x=init.';
beta=inv(subs(Jac,vars,x));
eta=beta*subs(eqn,vars,x);
for i=1:10,
    f=subs(eqn,vars,x);
    j=subs(Jac,vars,x);
    x=x-j\f;
end

%-------Jacobian of Jacobian
n=size(Jac,1);
F=cell(0);Fsym=cell(0);
F{1}=zeros(n);Fsym{1}=sym(F{1});
for i=2:4,   % 0 - 3rd order
    % i-1 order
    combin=nchoosek(n+i-2,n-1);
    F{i}=zeros(n,n*combin);
    Fsym{i}=sym(F{i});
end


for i=0:3   % 3rd for VDP
    s=0;
    for i1=0:i,
        for i2=0:i-i1
            for i3=0:i-i1-i2
%                 for i4=0:i-i1-i2-i3
%                     for i5=0:i-i1-i2-i3-i4
                        i4=i-i1-i2-i3;                        
                        % eval coef of element in Jacobian
                        tempJ=zeros(n);                        
                        for j=1:n,
                            for k=1:n,
                               p=feval(symengine,'poly',Jac(j,k),'[w,pr1, pr3, pi3]'); 
                               tempJ(j,k)=feval(symengine,'coeff',p,'[w,pr1, pr3, pi3]',...
                                           strcat('[',int2str(i1),',',int2str(i2),',',...
                                           int2str(i3),',',int2str(i4),']'));               
                            end
                        end
                        F{i+1}(:,s*n+1:(s+1)*n)=tempJ;
%                         Fsym{i+1}(:,s*n+1:(s+1)*n)=w^i1*pr1^i2*pr3^i3*pi3^i4*pr5^i5*pi5^i6;
                        s=s+1;
%                     end
%                 end               
            end
        end        
    end
    
%     norm(F{i+1})
end

% init=[w,c1,c3r,c3i,0,0];
% init=[1,1/sqrt(3),0,0];
init=[1,1/sqrt(3),0,0];
x=init.';

for i=1:1,
    f=subs(eqn,vars,x);
    j=subs(Jac,vars,x);
    beta=inv(j);
    eta=beta*f;
    K=norm(F{2},2)+norm(F{3},2)*(4*norm(eta,2)+2*norm(x,2))+3*norm(F{4},2)*(2*norm(eta,2)+norm(x,2))^2;
    K=K*2;
    [norm(eta,2), norm(beta,2), K, norm(eta,2)*norm(beta,2)*K]
    x=x-j\f;
end
% err=x-init.';