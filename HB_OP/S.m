function x = S(m1,m2)
% construct compact basis matrix for multiplication
% input m1,m2
x=zeros(m1*m2,m1+m2-1);
for i=1:m1
   x((i-1)*m2+1:i*m2,i:i+m2-1)=eye(m2); 
end

end