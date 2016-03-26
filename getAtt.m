%   Get attribute 
function [] = getAtt()

    global incidenceFull_DM;
    global EstimatedTime_DM;
    global TurnAngles_DM;
    global LeftTurn_DM;
    global Uturn_DM;
    global Atts;
    global Op;  
    Atts  = objArray(Op.n);
    
    Incidence = incidenceFull_DM;
    Atts(1).value = (Incidence .* EstimatedTime_DM);
    Atts(2).value = (Incidence .* TurnAngles_DM);
    Atts(3).value = (Incidence .* LeftTurn_DM);
    Atts(4).value = (Incidence .* Uturn_DM);
    %Atts(5).value = Matrix2D([]);
end
