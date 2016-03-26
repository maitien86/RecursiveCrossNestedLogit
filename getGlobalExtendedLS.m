%   Compute Link Size attribute from data
%   
%%
function ExpV_is_ok = getGlobalExtendedLS()
    global incidenceFull; 
    global Obs;     % Observation
    global LSatt;
    global EstimatedTime;
    global TurnAngles;
    global LeftTurn;
    global Uturn;    
    beta = [-2.5,-1,-0.4,-4]';
    %beta = [-1,0,0,0]';
    % ----------------------------------------------------
    mu = 1; % MU IS NORMALIZED TO ONE
    Obs_temp = Obs;
    Obs_temp(find(Obs)) =  Obs(find(Obs)) + 1;
    % Save data
    IF = (incidenceFull);
    ET = (EstimatedTime);
    TA = (TurnAngles);
    LT = (LeftTurn);
    UT = (Uturn);
    
    incidenceFull = transferMarix(incidenceFull);
    EstimatedTime = transferMarix(EstimatedTime);
    TurnAngles = transferMarix(TurnAngles);
    LeftTurn = transferMarix(LeftTurn);
    Uturn = transferMarix(Uturn);    
    lastIndexNetworkState = size(incidenceFull,1);
    dest = unique(Obs_temp(:,1));
    orig = unique(Obs_temp(:,2));
    incidenceFull(1,orig') = 1;
    EstimatedTime(1,orig') = 0.05;
    incidenceFull(dest,lastIndexNetworkState) = 1;
    EstimatedTime(dest,lastIndexNetworkState) = 0.05;
    
    M = getM(beta,false);
    [expV, expVokBool] = getExpV(M);
    if (expVokBool == 0)
       ExpV_is_ok = false;
       disp('ExpV is not fesible')
       return; 
    end  
    P = getP(expV, M);
    G = sparse(zeros(size(expV)));
    G(1) = 1;
    I = speye(size(P));
    F = (I-P')\G;                        
    if (min(F) < 0)                
        ToZero = find(F <= 0);
        for i=1:size(ToZero,1)
            F(ToZero(i)) = 1e-9;
        end
    end
    ips = F;
    ips = ips(1:lastIndexNetworkState); 
    ips(size(incidenceFull,2),1) = 0;   
    I = find(incidenceFull);
    [nbnonzero, c] = size(I);
    ind1 = zeros(nbnonzero,1);
    ind2 = zeros(nbnonzero,1);
    s = zeros(nbnonzero,1);
    for i = 1:nbnonzero
        [k a] = ind2sub(size(incidenceFull), I(i));
        ind1(i) = k;
        ind2(i) = a;
        s(i) =  ips(k) * P(k,a);
    end    
    LSatt = sparse(ind1, ind2, s, size(incidenceFull,1), size(incidenceFull,2));
    % Recover real size
    incidenceFull = IF;
    EstimatedTime = ET;
    TurnAngles = TA;
    LeftTurn = LT;
    Uturn = UT;
    LSatt = LSatt(2:size(IF,1)+1,2:size(IF,2)+1);
    LSatt = 100 * (LSatt .* ET);
    ExpV_is_ok = true;
end