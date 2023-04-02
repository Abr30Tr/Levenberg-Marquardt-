function [meu h ] = decent_solve(self,A_res,G_res,meu)
          % solve (JTJ + µI) hlm = -JTr depending on meu h can be steepest
          % desent or newton raphson
          % update the value of meu
      
          f_max = max(abs(A_res(:)));
          if meu == 0
             h=((-1)*G_res)/A_res;
          end 
          mat = size(A_res); flag = 1; 
          while flag
                [R,flag] = chol(A_res +meu*eye(mat)) ;
          end
          if flag 
              meu = max(self.tau*meu, self.eps*max(f_max));
          end 
         % b = R*R';
          h = R \ (R' \ (-G_res));
         % h = -1*(G_res)/b;        
end 