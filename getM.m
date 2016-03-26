%   Get MUtility
%%
function Mfull = getM(x,kM,isLS)   

    global incidenceFull_DM;
    global EstimatedTime_DM;
    global TurnAngles_DM;
    global LeftTurn_DM;
    global Uturn_DM;    
    global LinkSize;
    global isFixedUturn;  
    global ConstantAtts;
    global mu;
    u1 = x(1) * EstimatedTime_DM;
    u2 = x(2) * TurnAngles_DM;
    u3 = x(3) * LeftTurn_DM;
    if isFixedUturn == false
        u4 = x(4) * Uturn_DM;
        if isLS == true
            u5 = x(5) * LinkSize;
        else
            u5 = 0 * LeftTurn_DM;
        end
    else 
        u4 = -20 * Uturn_DM;
        if isLS == true
            u5 = x(4) * LinkSize;
        else
            u5 = 0 * LeftTurn_DM;
        end
    end
    u = sparse(u1 + u2 + u3 + u4 + u5 + ConstantAtts);
%     a = 1 ./ mu;
%     if size(a,1) > size(u,1)
%        a = a(1:size(u,1));
%     else
%        a(size(a,1)+1:size(u,1),1) = ones(size(u,1) - size(a,1),1); 
%     end
%     u = (bsxfun(@times,a,u));
    sizeU = size(u,1);
    kM = kM(1:sizeU,1:sizeU); 
    u = kM * u;
    expM = u ;
    expM(find(incidenceFull_DM)) = exp(u(find(incidenceFull_DM)));
    %expM = exp(u);
    Mfull = incidenceFull_DM .* expM;

end