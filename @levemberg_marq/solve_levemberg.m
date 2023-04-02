function [X hist meu norm_H A_res norm_G k_r_sol i stop] =solve_levemberg(self,fun,jac,sol,r_sol,A_res,G_res,norm_G,meu,x,fitData)
         % this solves the iterative value of updates (xnew). it also solves values of stopping criteria. 
         h = 0;   k = 0;  i = 0; stop = 0; 

            while ~k 
                if norm_G <= self.eps1
                   k=1;
                
                   else 
                      % solve value for the input of if condition 
                       chalk_1 = self.eps2*(norm(x) + self.eps2);
                       [meu h ] = decent_solve(self,A_res, G_res, meu);
                       norm_H = norm(h);
                    if  norm_H <= chalk_1
                      k =1
                    end
                end
              if ~k 
                      x_new = x + h; 
                      % soving the value of sigma 
                      h = x_new - x;
                      sub_2 = (meu * h - G_res);
                      h_tran = h';
                      sigma_f = 0.5 * h_tran * sub_2;
                      [k_sol k_jac k_r_sol] = iter_solver(self,fun,jac,x_new,fitData);
                      r_sub =  sol - k_sol;
                      r_add =  sol + k_sol;
                      r_sub =  r_sub';
                      sigma_l = 0.5 * r_sub *r_add;
                      sigma =  sigma_f/sigma_l;
                      if sigma > 0
                          x = x_new;
                          i = i + 1;
                          % solve value for the input of if condition 
                          sol = k_sol;
%                           r_sol=k_r_sol;
                          A_res = k_jac'*k_jac;      
                          G_res = k_jac'*k_sol;
                          norm_G = norm(G_res,inf);   
                          norm_X = norm(x);
               
                          % solving for new mue 
                          c_max = 2*sigma - 1;
                          c_max = c_max^3;
                            
                          meu = meu * max(1/3,(1-c_max));
                          v = 2;
                              if  nargout > 1
                                  X(:,i+1) = x_new;   
                                  hist(:,i+1) = [k_r_sol norm_G meu]';
                              end 
                      else 
                          meu = meu * v;
                          v = 2* v ;  

                      end 
                          if  k > self.max_iter    
                              stop = 1; 
                          end
                end 
            
            end 
end 
 
