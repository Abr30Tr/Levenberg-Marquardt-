    function [sol r_sol A_res G_res meu norm_G] = solve_inc(self,fun,jac,x,fitData) 
      % this function solves the initial values of the anonymous function, Jacobian, and meu. 
    if nargin > 3
            [sol r_jac r_sol]= iter_solver(self,fun,jac,x,fitData);
        end 
%          tr_jac = r_jac';
         A_res = r_jac'*r_jac;   
         G_res = r_jac'*sol;
         norm_G = norm(G_res,inf);  
         meu = self.tau * max(diag(A_res));
         
    end