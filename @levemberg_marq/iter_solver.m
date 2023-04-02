function [sol r_jac r_sol]= iter_solver(self,fun,jac,x,fitData)
         
% this function is called multiple times to solve iterative values of anonymous function and its Jacobian.
         if nargin>3
            r_jac = feval(jac,x,fitData(:,1));
            sol = feval(fun,fitData(:,2),fitData(:,1),x);
                   
         else
            sol = feval(fun,fitData(:,2),fitData(:,1),x);
         end 
         r_sol = (sol'*sol)/2;
end 