%
%   Predictive trust region iterate implementation
%   MAI ANH TIEN - 9.August.2013
%
%%
function isSuccess = btr_swpred_iterate()    
    global Op;
    Op.step = btr_step_steihaug_toint(Op);
    if isnan(Op.step)
       Op.H = Op.Hi(1).value; % back to BHHH     
    end
    Op.step = btr_step_steihaug_toint(Op);
    xold = Op.x;
    [newValue gradplus] = LL(Op.x + Op.step);
    Op.deltaGrad = gradplus - Op.grad;    
    % compute quadratic model
    md = Op.value + Op.grad'*Op.step + 0.5*Op.step' * Op.H * Op.step;
    rho = (Op.value - newValue) / (Op.value - md);
    if btr_accept_candidate(rho)
        Op.deltaValue = newValue - Op.value;
        oldValue = Op.value;
        oldGrad = Op.grad;
        Op.value = newValue;
        Op.grad = gradplus;
        % Select next Hessian ::
        ind = swpred_select_hessian(oldValue, Op.value, oldGrad, Op.step);
        update_multi_hessians_approx();
        Op.H = Op.Hi(ind).value;
        Op.Selected_Hessian_approx = Op.set_hes(ind);
        isSuccess = true;
    else
        isSuccess = false;
        Op.x = xold;
    end    
    % Update radius
    btr_update_radius(rho);  
end
