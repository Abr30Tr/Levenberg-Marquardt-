    function [x lm_result  hist] = save_result (self,fun,jac,sol,r_sol,A_res,G_res,norm_G,meu,x,fitData)
             [X hist meu norm_H A_res norm_G k_r_sol i stop] = solve_levemberg(self,fun,jac,sol,r_sol,A_res,G_res,norm_G,meu,x,fitData);       
           if  nargin > 2
               ii = 1 : i+1;  
               x = X(:,ii);   
               hist = struct('sol',hist(1,ii), 'norm_G',hist(2,ii), 'meu',hist(3,ii));
           else 
               x = X; 
           end
           if  stop 
               k_sol = NaN; 
               k_r_sol = NaN; 
           end
               lm_result = [norm_G  k_r_sol  norm_H];
    end 