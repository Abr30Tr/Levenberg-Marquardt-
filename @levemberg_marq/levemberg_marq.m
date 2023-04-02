%
classdef levemberg_marq < handle
  properties  (SetAccess = private, Hidden = true)
    fun
    jac
    eps1
    eps2
    tau
    max_iter
  end
  methods
    %
    %  Class contructor

    %
    function self = levemberg_marq( varargin )
      % set default values
      self.eps1      = 1e-8;
      self.eps2      = 1e-14;
      self.max_iter  = 100;
      self.tau      = 1e-4;
      if nargin > 0
        %
        % in a true class require to check the input
        % to be correct, 1 real number > 0 etc.
        %
        self.eps1 = varargin{1};
      end
      if nargin > 1
        self.eps2 = varargin{2};
      end
      if nargin > 2
        self.max_iter = varargin{3};
      end
      if nargin > 3
         self.tau = varargin{3};
      end
    end
     %
     function set_tol_esp1( self, new_tol )
      self.eps1  = new_tol;
    end
    %
    %
    function set_tol_esp2( self, new_tol_1 )
      self.esp2 = new_tol_1;
    end
    %
    function set_max_iter( self, max_iter )
      self.max_iter = max_iter;
    end
    %
    % store the function fun and its jacobian jac
    % in the class
    %
    function setup( self, f, jac )
      self.fun   = f;
      self.jac = jac;
    end
    %
    % Declare the methods solver without inplementation.
    % The implementation is in the file levemberg_marq and  in the
    % same directory
    %
    [sol r_sol A_res G_res meu norm_G] = solve_inc(self,fun,jac,x,fitData)
    [x lm_result  hist] = save_result(self,fun,jac,sol,r_sol,A_res,G_res,norm_G,meu,x,fitData)
  end
end
