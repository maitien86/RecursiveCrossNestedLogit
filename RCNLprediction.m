%% RCNL prediction

function [] = RCNLprediction(LS, indexSample)

    %% Credits
    Credits;
    globalVar;
    global resultsTXT; 
    global Zsave;
    global SampleObs;
    %% Data 

    file_linkIncidence = './Input/linkIncidence.txt';
    file_AttEstimatedtime = './Input/ATTRIBUTEestimatedtime.txt';
    file_turnAngles = './Input/ATTRIBUTEturnangles.txt';
    file_observations = './Input/observationsForEstimBAI.txt';
    %file_observations = './simulatedData/ObservationsAll_NoLS.txt';

    %% Set estimation type
    %global isFixedMu;
    %isFixedMu = 0;
    
    % Inilialize prediction sample
    PredSample = spconvert(load('../PredSampleObs.1.40.txt'));
    TXT = 'Mixed prediction:';
    TXT = [TXT sprintf('Sample:%s,%s |', indexSample,LS)];
    index = str2num(indexSample); 
    train = PredSample(index*2-1,:);
    test  = PredSample(index*2,:);
    SampleObs = train; 

    %% Initialize the optimizer structure
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
     if isLinkSizeInclusive == false
         Op.x = [-2.494,-0.933,-0.411,-4.459, -0.0, 0.0, -0.0,-0.001, 0.0010, -0.001]';
         Op.x = [-1.378421e+00;-5.172495e-01;-6.508789e-02;-2.906635e+00;6.365889e-01;-1.918225e-01;-2.654382e-02;-2.626613e-04;-2.847599e-05;4.745908e-01];
     else
         Op.x = [-3.060,-1.057,-0.353,-4.431, -0.227, -0.0, 0.0, -0.0,0.001, 0.001, -0.001]';
         Op.x = [-1.566538e+00;-5.684770e-01;-7.248798e-02;-2.964155e+00;-1.149496e-01;4.434236e-01;-1.567649e-01;-2.059576e-02;1.835406e-08;-1.737376e-07;-4.832981e-01];
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
    resultsTXT = header;
    PrintOut(Op);
    TXT = [TXT sprintf('Gradient: %f |',norm(Op.grad))];
    % Compute predicted LL
    SampleObs = test;
    PredLL = LL(Op.x);
    
    %% Finishing, 
    TXT = [TXT sprintf('PredLL: %f\n',PredLL)];
    fileID = fopen('ResultsRCNLPred.txt','at+');
    fprintf(fileID,TXT);
    reTxt = resultsTXT;
   
    
end
