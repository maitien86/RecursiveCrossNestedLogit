function [] = RCNLoptimizer_retro(LS)

    %% Credits
    Credits;
    globalVar;
    global resultsTXT; 
    global Zsave;
    % Initialize email notification
    % notifyMail('set','amyeuphich@gmail.com','sntal2908');

    %% Data 

    file_linkIncidence = './Input/linkIncidence.txt';
    file_AttEstimatedtime = './Input/ATTRIBUTEestimatedtime.txt';
    file_turnAngles = './Input/ATTRIBUTEturnangles.txt';
    file_observations = './Input/observationsForEstimBAI.txt';
    %file_observations = './simulatedData/ObservationsAll_NoLS.txt';

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
     %nbobs = 5;
     Op = Op_structure;
     initialize_optimization_structure();

     %Op.x = Op_old.x(1:Op.n);
     %Op.delta = 0.05;
     %Op.H = getHessian();
     if isLinkSizeInclusive == false
         Op.x = [-2.494,-0.933,-0.411,-4.459, -0.0, 0.0, -0.0,-0.001, 0.0010, -0.001]';
     else
         Op.x = [-3.060,-1.057,-0.353,-4.431, -0.227, -0.0, 0.0, -0.0,0.001, 0.001, -0.001]';
     end
     %Op.x = [-3.060,-2.057,-3.353,-4.431, -0.0, 0.0, -0.0]';
     %Op.x = [-2,-1,-1,-4, -0.2, -0.1, -1]';
     %Op.x = [-2.494;-0.933;-0.411;-4.459;0;0;0];
     %mu = Scale(:,1)+0.2;
     %Op.x = [-1;-1;-1;-1];
     % Op.x = [-1.340446e+00,-7.776641e-01,-8.111989e-01,-5.534458e+00,-2.35543e+00]';

    %% Relax Att
    getAtt();
    for i = Op.m+1: Op.n
        u = sparse(zeros(size(incidenceFull)));
        Atts(i).value = (u);
    end
    initialize_switching_structure();
    Op.Switching = OptimizeConstant.SW_RETRO;
    Op.Optim_Method = OptimizeConstant.TRUST_REGION_METHOD;
    Op.Hessian_approx = OptimizeConstant.BHHH;
    Gradient = zeros(nbobs,Op.n);
    Zsave = objArray(nbobs);
    for i = 1:nbobs
        Zsave(i).value = [];
    end

    % Generate Observations
    %   createSimulatedObs;
    %   file_observations = './simulatedData/ObservationsAll.txt';
    %  loadData;

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
        ok = btr_swretro_interate();
        %ok = btr_swpred_iterate();
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
    %PrintOut(Op);
    %disp(' Calculating VAR-COV ...');
    getCov;
    %Finishing ...
    ElapsedTtime = toc
    resultsTXT = [resultsTXT sprintf('\n Number of function evaluation %d \n', Op.nFev)];
    resultsTXT = [resultsTXT sprintf('\n Estimated time %d \n', ElapsedTtime)];
    %% Send email nitification    
end
