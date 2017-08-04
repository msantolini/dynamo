function J = compJacobian(f,x)
% numerical computation of the jacobian of f around x
% using simple Euler method

n = length(x);
fx = feval(f,x); % evaluate the function at point x

for i=1:n
    xstep = x;
    % compute the i-th component partial derivative
    % numerically using first order forward difference approximation
    
    if (x(i) == 0)
        step = 1e-12; % difference step you choose, can be 1e-10 if you like
        xstep(i)=x(i)+step;
    else
        xstep(i)= (1 + 1e-3) * x(i);
        step = xstep(i) - x(i);
    end
    J(i,:)=(feval(f,xstep)-fx)/step;
end;
end
