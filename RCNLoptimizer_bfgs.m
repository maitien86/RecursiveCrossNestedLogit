function [] = RCNLoptimizer(LS)

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
         Op.x = [-1.378421e+00;-5.172495e-01;-6.508789e-02;-2.906635e+00;6.365889e-01;-1.918225e-01;-2.654382e-02;-2.626613e-04;-2.847599e-05;4.745908e-01];
     else
         Op.x = [-3.060,-1.057,-0.353,-4.431, -0.227, -0.0, 0.0, -0.0,0.001, 0.001, -0.001]';
         Op.x = [-1.566538e+00;-5.684770e-01;-7.248798e-02;-2.964155e+00;-1.149496e-01;4.434236e-01;-1.567649e-01;-2.059576e-02;1.835406e-08;-1.737376e-07;-4.832981e-01];
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
    Op.Hessian_approx = OptimizeConstant.BFGS;
    Gradient = zeros(nbobs,Op.n);
    Zsave = objArray(nbobs);
    for i = 1:nbobs
        Zsave(i).value = [];
    end
    %% Starting optimization
    tic ;
    %nbobs = 2;
    disp('Start Optimizing ....');
    
    %options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton','GradObj','on');
    %[x,fval,exitflag,output,grad] = fminunc(@LL,Op.x,options)
    
    options_const =  optimoptions(@fmincon,'Display','iter','Algorithm','interior-point','GradObj','on');
    [x,fval,exitflag,output,lambda,grad] = fmincon(@LL,Op.x,[],[],[],[],[],[],[],options_const)
    toc
    %% Send email nitification    
end
