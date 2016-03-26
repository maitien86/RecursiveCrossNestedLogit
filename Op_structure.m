% Optimization structure
classdef Op_structure
   properties
      Optim_Method = OptimizeConstant.LINE_SEARCH_METHOD; % or 'BTR'
      Switching = OptimizeConstant.NONE;
      Hessian_approx = OptimizeConstant.BFGS;     
      Selected_Hessian_approx = OptimizeConstant.BFGS;
      
      n = 5;          % size of vector
      m = 4;
      alpha = 0.0001; % Paramters for linesearch algorithm
      gamma = 0.9 ;   % Paramters for linesearch algorithm
      maxIter = 200;
      nFev = 0;       % Number of function evaluation
      tol = 0.00001;
      
      % For switching algorithms
      prev_x = [];
      prev_value = 0;
      nSwitch = 0;
      prevHind = 1;
      
      
      x = [];
      grad = [];
      deltaGrad = [];
      step = [];
      stepLength = 0;
      k = 0;
      value = 0;
      deltaValue = 0;
      H = [];   
      Ak = [];% Hessian approximation.
      radius = 1;   
      % for switching
      nH = 2; % number of Hessian approximation
      Ai = [];
      set_hes = [];
      Hi = [];
      
      % trust region radius.   
      delta = 1;
      ETA1 = 0.2;
      ETA2 = 0.75;
   end
end