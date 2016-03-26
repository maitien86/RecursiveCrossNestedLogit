%   Get gradient of mu
%%
function gradMu = getGradMu(Mu)   
    global Scale;
    global prevLink;
    global nDlinks; % number of dummy links
    global nRlinks; % number of real links
    global Op;
    t = (Op.n - Op.m)/2;
    theta = Op.x(Op.m+1: Op.m + t);
    lambda = Op.x(Op.m+t+1: Op.n);
    N = size(Scale,1);
    n1 = nRlinks+1; n2 = nRlinks + nDlinks; 
   
    
     %% Compute gradient of Mu
    gradMu = (zeros(N,Op.n));
    for i = Op.m+1 : Op.m+3
         gradMu(1:nRlinks,i) = Mu(1:nRlinks) .* Scale(1:nRlinks,i - Op.m);
         gradMu(n1:n2,i) = Mu(n1:n2) .* Scale(prevLink(n1:n2), i-Op.m);
    end
    
    for i = Op.m+t+1 : Op.n
         gradMu(1:nRlinks,i) = Mu(1:nRlinks) * 0;
         gradMu(n1:n2,i) = - 2 * Op.x(i) * (Mu(n1:n2) .* Scale(n1:n2, i - Op.m - t));
    end 
    
%    Mu = exp(Scale(:,1:t) * theta);
%    Mu(n1:n2) = exp(Scale(prevLink(n1:n2),1:t) * theta - Scale(n1:n2,1:t) * (lambda.^2));
end