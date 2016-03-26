function [f g] = LL(x)
    global Op;
    Op.x = x;
    [f g] = getLL_crossNested();
    Op.nFev  = Op.nFev + 1;
end