% select hessian for predictive model
%%
%
function ind = swpred_select_hessian(oldValue, newValue, grad, step)
    global Op;
    ind = 0;
    mError = bitmax;
    for i =1:Op.m
       error = abs(newValue - get_quadratic_model(oldValue, grad, step, Op.Hi(i).value))
       if (mError >= error)
           ind = i;
           mError = error;
       end
    end
    if ind ~= Op.prevHind 
        Op.nSwitch = Op.nSwitch + 1;
    end    
    Op.prevHind = ind;
    
end
%% 
%% Compute: f + g'p + 0.5p'Hp
function model = get_quadratic_model(f,g,p,H)
     model = f + g'*p + 0.5*p' * H * p;
end