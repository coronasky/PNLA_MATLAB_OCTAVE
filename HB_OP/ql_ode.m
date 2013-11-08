function dy = ql_ode(t,y)

% u=sin(2*pi*f*t);
dy=zeros(2,1);
% dy(1)=-y(2)-y(1)+u;
% dy(2)=2*(-y(2)-y(1)+u);


dy(1)=y(2);    % x'=-x-y+u; y=-x^2+1
dy(2)=-y(1)+y(2)-3*y(1)*y(1)*y(2);
% dy(3)=alpha*y(3);

