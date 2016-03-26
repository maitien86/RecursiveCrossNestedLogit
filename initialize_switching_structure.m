% Initialie optimization structure for predictive model
%%
%
function [] = initialize_switching_structure()
    global Op;
    Op.Switching = OptimizeConstant.SWITCHING;
    if strcmp(Op.Optim_Method,OptimizeConstant.TRUST_REGION_METHOD)
        Op.nH = 2;
    else
        Op.nH = 2;
    end
    %Op.set_hes = [OptimizeConstant.BHHH;OptimizeConstant.CB_BFGS; OptimizeConstant.CB_SR1;  OptimizeConstant.SSA_BFGS; OptimizeConstant.SSA_SR1; ];
    Op.set_hes = [OptimizeConstant.BHHH;OptimizeConstant.SR1;OptimizeConstant.SSA_SR1];
    Op.Hi = objArray(Op.nH);
    Op.Ai = objArray(Op.nH);
    for i = 1:Op.nH
        Op.Hi(i).value = eye(Op.n);
        if Op.set_hes(i) == OptimizeConstant.CB_SR1 || Op.set_hes(i) == OptimizeConstant.CB_BFGS
            Op.Ai(i).value = zeros(Op.n);
        end
    end
end