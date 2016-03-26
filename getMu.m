%   Get mu
%%
function Mu = getMu()   
    global Scale;
    global prevLink;
    global nDlinks; % number of dummy links
    global nRlinks; % number of real links
    global Op;
    t = (Op.n - Op.m)/2;
    theta = Op.x(Op.m+1: Op.m + t);
    lambda = Op.x(Op.m+t+1: Op.n);
    Mu = exp(Scale(:,1:t) * theta);
    n1 = nRlinks+1; n2 = nRlinks + nDlinks; 
    Mu(n1:n2) = exp(Scale(prevLink(n1:n2),1:t) * theta - Scale(n1:n2,1:t) * (lambda.^2));
end