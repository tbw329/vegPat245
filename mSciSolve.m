function [x,converged,J]=mSciSolve(f,x0,tol,maxit,jach)

%f - function to solve for f(x) = 0
%x0 - initial guess
%df - Jacobian at x0
%jach: finite difference to use for calculating Jacobian
%tol
%maxIt: max number of newton iterations b4 giving up

%initialise
xold = x0;
J = mSciJacobian(f,xold,jach);
converged = 0;

for i = 1:maxit   
    x = xold - (J^-1)*f(xold);%Newton iteration
    J = mSciJacobian(f,x,jach); %Jacobian at new x
    if (all(abs(x-xold)) < tol) && (all(f(x)) < tol) %Break if close to ans.
        converged = 1;
        break
    else    
        xold = x; %repeat if not converged yet
    end
end
    