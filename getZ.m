function [Z, expVokBool] = getZ(M, MI, phi)
    global Op;
    N = size(MI,1);
    b = zeros(N,1);
    b(N) = 1;
    A = speye(N) - M;
    Z = A\b;
    Z = full(Z);
    %Z = ones(N,1);
    j = 0;
     while(1)
        j = j+1 ;     
        Zprev = Z;
        Z = getNextZ(Z, M, MI, b,phi);
        %norm(log(Z) - log(Zprev))
        if mod(j,10) == 0
           %j         
           if norm(log(Z) - log(Zprev)) < 0.0001;% norm(Op.grad)*10 %0.0001
             break;
           end
           %norm(log(Z) - log(Zprev))
           %Z(100)
           if (j > 1500)             
            break ;
           end
        end
     end
     residual = norm(log(Z) - log(Zprev));
     % Check feasible
     minele = min(Z(:));
     expVokBool = 1;
     if minele == 0 || minele < OptimizeConstant.NUM_ERROR
       expVokBool = 0;
       fprintf('min zero');
     end 
     
     Zabs = abs(Z);
     D = Z - Zprev;
     %resNorm = norm(D(:));
     if residual > 10 || (~isreal(Z))
       expVokBool = 0;
     end    
 end

function Znext = getNextZ(Z, M, MI, b, phi)
    n = size(M,1);
    e = ones(n,1);
    I = speye(n);
    %U = sparse(bsxfun(@times,Z',MI));
    
    Zsp =  spdiags(Z,0,n,n);
%      U = MI * (Zsp);
%      X = sparse(MI * 0);
%      T = find(MI);
% % %     tic;
%       X(T) =  U(T) .^ phi(T);
%      V1 = X;
      
    U = MI * Zsp + MI * realmin;
    [~,~,sMI] = find(MI * realmin);
    [~,~,sU] = find(U);
    [i,j,sPhi] = find(phi);
    s = (sU - sMI) .^ sPhi;
    X = sparse(i,j,s,n,n);
    
    Z =  (M .* X) * e;
    Znext = Z + b;   
end
