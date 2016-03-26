%% MAIN function

function [] = RCNLoptimizer(LS)

    %% Credits
    Credits;
    globalVar;
    global resultsTXT; 
    global Zsave;

    %% Data 

    file_linkIncidence = './Input/linkIncidence.txt';
    file_AttEstimatedtime = './Input/TravelTime.txt';
    file_turnAngles = './Input/TurnAngles.txt';
    file_observations = './simulatedData/ObservationsAll.txt';

    %% Set estimation type
    global isFixedMu;
    isFixedMu = 0;

    %% Initialize the optimizer structure
     %Op_old = Op;
     if strcmp(LS,'LS')
         isLinkSizeInclusive = true;
     else
         isLinkSizeInclusive = false;
     end
     isFixedUturn = false;
     loadData;
     Op = Op_structure;
     initialize_optimization_structure();
    
     if isLinkSizeInclusive == false
         Op.x = [-2.494,-0.933,-0.411,-4.459, -0.0, 0.0, -0.0,-0.00, 0.000, -0.00]';
%         Op.x = [-1.378421e+00;-5.172495e-01;-6.508789e-02;-2.906635e+00;6.365889e-01;-1.918225e-01;-2.654382e-02;-2.626613e-04;-2.847599e-05;4.745908e-01];
     else
         Op.x = [-3.060,-1.057,-0.353,-4.431, -0.227, -0.0, 0.0, -0.0,0.001, 0.001, -0.001]';
%         Op.x = [-1.566538e+00;-5.684770e-01;-7.248798e-02;-2.964155e+00;-1.149496e-01;4.434236e-01;-1.567649e-01;-2.059576e-02;1.835406e-08;-1.737376e-07;-4.832981e-01];
     end
     
    %% Relax Att
    getAtt();
    for i = Op.m+1: Op.n
        u = sparse(zeros(size(incidenceFull)));
        Atts(i).value = (u);
    end
    initialize_switching_structure();
    Op.Switching = OptimizeConstant.SW_RETRO;
    Op.Optim_Method = OptimizeConstant.TRUST_REGION_METHOD;
    Op.Hessian_approx = OptimizeConstant.BFGS;
    Gradient = zeros(nbobs,Op.n);
    Zsave = objArray(nbobs);
    for i = 1:nbobs
        Zsave(i).value = [];
    end

    %% Starting optimization
    tic ;
    disp('Start Optimizing ....');
    [Op.value, Op.grad ] = LL(Op.x);
    Op.delta = norm(Op.grad) * 0.1;
    PrintOut(Op);
    header = [sprintf('%s \n',file_observations) Op.Optim_Method];
    header = [header sprintf('\nNumber of observations = %d \n', nbobs)];
    header = [header sprintf('Hessian approx methods = %s \n', OptimizeConstant.getHessianApprox(Op.Hessian_approx))];
    resultsTXT = header;

    %% Loop
    while (true)    
      Op.k = Op.k + 1;
      if strcmp(Op.Optim_Method,OptimizeConstant.LINE_SEARCH_METHOD);
        ok = line_search_iterate();
        if ok == true
            PrintOut(Op);
        else
            disp(' Unsuccessful process ...')
            break;
        end
      else
        ok = btr_interate();
        PrintOut(Op);
      end
      [isStop, Stoppingtype, isSuccess] = CheckStopping(Op);  
      %----------------------------------------
      % Check stopping criteria
      if(isStop == true)
          isSuccess
          fprintf('The algorithm stops, due to %s \n', Stoppingtype);
          resultsTXT = [resultsTXT sprintf('The algorithm stops, due to %s \n', Stoppingtype)];
          break;
      end
    end
    %% Compute variance - Covariance matrix
    disp(' Calculating VAR-COV ...');
    getCov;
    %Finishing ...
    
    ElapsedTtime = toc
    resultsTXT = [resultsTXT sprintf('\n Number of function evaluation %d \n', Op.nFev)];
    resultsTXT = [resultsTXT sprintf('\n Estimated time %d \n', ElapsedTtime)];
    dlmwrite('output.txt',resultsTXT);
    %% Send email nitification    
end
