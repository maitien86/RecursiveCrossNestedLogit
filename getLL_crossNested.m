% Compute the loglikelohood value and its gradient.
% Based on V -  Relaxing scales mu.
%% NOTE: ROW: U * diag,  COL: diag * U.
function [LL, grad] = getLL_crossNested()

    global incidenceFull_DM; 
    global Gradient;
    global Op;
    global Mfull;
    global Ufull;
    global Atts;
    global Obs;     % Observation
    global nbobs;  
    global isLinkSizeInclusive;
    global lastIndexNetworkState;
    global mu;
    global LSatt_DM;
    global LinkSize;
    global prevLink;
    global nDlinks; % number of dummy links
    global nRlinks; % number of real links
    global SampleObs;


    %% Setting sample
    if isempty(SampleObs)
        sample = 1:nbobs;
    else
        sample = SampleObs;
        nbobs = size(find(SampleObs),2);
    end
  
    
    %% Get M, U
    [lastIndexNetworkState, maxDest] = size(incidenceFull_DM);     
    LL = 0;
    grad = zeros(1, Op.n);
    mu = getMu();
    N = lastIndexNetworkState + 1;
     %% Compute gradient of Mu
    gradMu = getGradMu(mu);
    %%
    %% Compute phi(a|k)=mu_a / mu_k;    
    %Mfull = getM(Op.x, false); % matrix with exp utility for given beta
    M = incidenceFull_DM(1:lastIndexNetworkState,1:lastIndexNetworkState);            
    M(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
    M(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));      
    MI = sparse(M); 
    MI(find(M)) = 1;
    dPhi = objArray(Op.n - Op.m);
    a = mu;
    k = 1 ./ mu;
    GrdMuCol = objArray(Op.n);
    kM  = spdiags(k,0,N,N); %bsxfun(@times,k', MI')';
    aM =  spdiags(a,0,N,N); %bsxfun(@times,a', MI')';
    for i = 1: Op.n
        GrdMuCol(i).value = spdiags(gradMu(:,i),0,N,N); %bsxfun(@times, gradMu(:,i)', MI);
    end
%     for i = Op.m + 1: Op.n
%        dPhi(i - Op.m).value = spdiags(-Scale(:,i - Op.m),0,N,N) * MI + MI * spdiags(Scale(:,i - Op.m),0,N,N);  
%     end        
    phi = (kM * MI) .* (MI * aM);  

    e = ones(size(M,1),1); 
    %% Compute gradient of phi(.)
    gradPhi = objArray(Op.n);
    for i = 1: Op.n
        if i <= Op.m
           gradPhi(i).value = sparse(MI * 0); 
        else
           gradPhi(i).value =  (MI * spdiags(gradMu(:,i),0,N,N)) .* (kM * MI) - (MI * aM) .*  (spdiags(k .* k .* gradMu(:,i),0,N,N) * MI);
           gradPhi(i).value = sparse(gradPhi(i).value);
        end
    end    
    %% 
    if isLinkSizeInclusive == true
        sizeOfParams = Op.m - 1;
    else
        sizeOfParams = Op.m;
    end 
    AttLc = objArray(Op.n);
    for i = 1 : sizeOfParams
        AttLc(i).value =  sparse(Atts(i).value(1:lastIndexNetworkState,1:lastIndexNetworkState));
        AttLc(i).value(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
        AttLc(i).value(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    end
    for i = Op.m+1: Op.n
        AttLc(i).value =  sparse(AttLc(1).value * 0);
    end
    
    gradVd = sparse(zeros(N ,Op.n));    
    %%Create temp matrices
    
    %% Loop over observations
    %fprintf(' n = ');
    if isLinkSizeInclusive == false
        prevOD = [0];
        ii = 1;
    else
        prevOD = [0,0];
        ii = 2;
    end
    for n = 1:nbobs
        %tic;
%         fprintf('%d - ',n);
%         if mod(n,30) == 0
%             fprintf('\n');
%         end
        dest = Obs(sample(n), 1);
        if ~isequal(prevOD,Obs(sample(n),1:ii)) ;
            prevOD = Obs(sample(n),1:ii);
            %% set Link Size attributes        
            if isLinkSizeInclusive == true
                LinkSize = LSatt_DM(sample(n)).value;
                Atts(Op.m).value = LinkSize;
                AttLc(Op.m).value =  (LinkSize(1:lastIndexNetworkState,1:lastIndexNetworkState));
                AttLc(Op.m).value(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
                AttLc(Op.m).value(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
            end
            %% Compute M matrix
            Mfull = getM(Op.x, kM, isLinkSizeInclusive); % matrix with exp utility for given beta
            M = Mfull(1:lastIndexNetworkState,1:lastIndexNetworkState);            
            addColumn = Mfull(:,dest);
            M(:,lastIndexNetworkState+1) = addColumn;
            M(lastIndexNetworkState+1,:) = zeros(1,lastIndexNetworkState+1);
            M = sparse(M);
            Ufull = getU(Op.x, isLinkSizeInclusive); % matrix with exp utility for given beta
            U = Ufull(1:lastIndexNetworkState,1:lastIndexNetworkState);            
            addColumn = Ufull(:,dest);
            U(:,lastIndexNetworkState+1) = addColumn;
            U(lastIndexNetworkState+1,:) = zeros(1,lastIndexNetworkState+1);
            U = sparse(U);

            %% get Z 
            MI = sparse(M); 
            MI(find(M)) = 1;
            phi =  (kM * MI) .* (MI * aM);
            phi = sparse(phi);
            [Z, expVokBool] = getZ(M, MI, phi);
            %[Z, expVokBool] = getZ_NK(M, MI, phi);
            %[Z, expVokBool] = getV_NK(U, M, MI, phi);
            if (expVokBool == 0)
                LL = OptimizeConstant.LL_ERROR_VALUE;
                grad = ones(Op.n,1);
                disp('The parameters not fesible')
                return; 
            end
            %% Compute V
           V = log(Z);
           V = mu .* V; %(bsxfun(@times,mu,V)); 
            %% Compute gradient of V - respect to attributes parameters   
            gradV = objArray(Op.n);    
            for i = 1:Op.n
                % Compute gradient of M
                h = gradMu(:,i) .* log(Z);
                %U1 = sparse(bsxfun(@times,k, AttLc(i).value));
                %U1 = kM .* AttLc(i).value;
                U1 = kM * AttLc(i).value;
                %U2 = sparse(bsxfun(@times,k .* k .* gradMu(:,i) , U)); 
                %U2 = kM .* kM .*  GrdMuCol(i).value .* U;
                U2 = (kM .* kM .*  GrdMuCol(i).value) * U;

                gradM = sparse(M .* (U1 - U2));
                % Compute gradient of V
                %Zd = sparse(bsxfun(@times,Z',MI));
                Zd = MI * spdiags(Z,0,N,N);
                X = MI;
                X(find(MI)) =  Zd(find(MI)) .^ (phi(find(MI)));
                %X = bsxfun(@ldivide,Z',X')';
                X = spdiags(1 ./ Z,0,N,N) * X;

                H = M .* X;
                %U1 = bsxfun(@times, log(Z'), H);
                U1 = H * spdiags(log(Z),0,N,N);

                %U2 = aM .* U1 .* gradPhi(i).value;
                U2 = (aM * U1) .* gradPhi(i).value;

                U3 = U1 * GrdMuCol(i).value ;
                %S =  aM .* gradM .* X + U2 - U3;
                S =   aM * (gradM .* X) + U2 - U3;
                gradV(i).value = sparse((speye(size(M)) - H)\( S * e + h));   

            end
        end
        %% Compute gradient of Z    
        lnPn = 0;
        for i = 1: Op.n
             Gradient(n,i) = 0;
             gradVd(:,i) =  gradV(i).value;
        end
        path = Obs(sample(n),:);
        lpath = size(find(path),2);
        % Compute log-likelihood and gradient
        lnPpath = 0;
        gradPpath = zeros(1, Op.n);
        
        for i = 2:lpath - 1
            fLink = path(i);
            tLink = path(i + 1);
            eLink = min(path(i+1),lastIndexNetworkState + 1);
            nextL = find(incidenceFull_DM(fLink,:));
            nd =  size(nextL,2); % number of dummies to add
            P = 0;
            gradP = zeros(1, Op.n);
            for j = 1:nd
                dLink = nextL(j);
                temp = (Ufull(fLink,dLink) + V(dLink) - V(fLink))/mu(fLink) + (Ufull(dLink,tLink) + V(eLink) - V(dLink))/mu(dLink) ; 
                %temp =  Ufull(dLink,tLink);% + Ufull(dLink,tLink);%V(dLink) - V(fLink) + V(eLink) - V(dLink);
                P = P + exp(temp);             
                for t = 1:Op.n
                    p1 =   (AttLc(t).value(fLink,dLink) + gradVd(dLink,t) - gradVd(fLink,t))/mu(fLink);
                    p2 =   (Ufull(fLink,dLink) + V(dLink) - V(fLink))/(mu(fLink)^2) * gradMu(fLink,t) ;
                    p3 =   (AttLc(t).value(dLink,eLink) + gradVd(eLink,t) - gradVd(dLink,t))/mu(dLink);
                    p4 =   (Ufull(dLink,tLink) + V(eLink) - V(dLink))/(mu(dLink)^2) * gradMu(dLink,t) ;
                    gradP(t) = gradP(t) + (p1 - p2 + p3 - p4) * exp(temp);
%                    gradP(t) =  gradP(t) + (AttLc(t).value(dLink,eLink)) * exp(temp);% % gradP(t) + ( gradVd(dLink,t) - gradVd(fLink,t) + gradVd(eLink,t) - gradVd(dLink,t)) * exp(temp);
                end
                %gradP = gradP / P;
            end
            P;
            lnPpath = lnPpath + log(P);
            gradPpath =  gradPpath + gradP / P;
        end
        
        Gradient(n,:) = Gradient(n,:) + gradPpath;
        lnPn = lnPn + lnPpath;  
        LL =  LL + (lnPn - LL)/n;
        grad = grad + (Gradient(n,:) - grad)/n;
        Gradient(n,:) = - Gradient(n,:);
        %toc
    end
 %   fprintf('\n');
    LL = -1 * LL; % IN ORDER TO HAVE A MIN PROBLEM
    grad =  - grad';
end
