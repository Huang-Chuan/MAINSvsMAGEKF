function J=numjacobian(f,x)

n=length(x);
fx=feval(f,x); % evaluate the function at point x
%step=0.0031; % difference step
step=12*1024/(4e6); % difference step

for i=1:n
   xstep = x;
   % compute the i-th component partial derivative
   % numerically using first order forward difference approximation
   xstep(i)=x(i)+step;
   J(:,i)=(feval(f,xstep)-fx)/step;
end
end
