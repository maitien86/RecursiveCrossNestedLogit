%   Get Utility
%%
function Ufull = getU(x, isLS)

    global incidenceFull_DM;
    global EstimatedTime_DM;
    global TurnAngles_DM;
    global LeftTurn_DM;
    global Uturn_DM;
    global ConstantAtts;
    global LinkSize;
    global isFixedUturn;
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
    u = sparse(u1 + u2 + u3 + u4 + u5 + ConstantAtts) ;
    Ufull = incidenceFull_DM .* u;
end


