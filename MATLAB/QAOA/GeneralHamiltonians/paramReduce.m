function param = paramReduce(param, sym)
%paramReduce remove degeneracies of parameters due to symmetries
%
%   paramOut = paramReduce(param, sym)
%
%   sym is a string containing descriptions of symmetries
%     accepted ones are:
%       TR = Time-Reversal symmetry
%       Z2 = Z2 symmetry from X^{\otimes N} (Beta is periodic up to Pi/2)
%       Pi_gamma = Gamma is periodic up to Pi


p = length(param)/2;
bTR = contains(sym, 'TR');
bZ2 = contains(sym, 'Z2');

if bTR && param(1) < 0 % time reversal symmetry
    param = -param;
end


if contains(sym, 'Pi_gamma')
    param(1:p) = mod(param(1:p),pi);
    if bTR && param(1) > pi/2
        param(1:p) = pi-param(1:p);
        param(p+1:end) = -param(p+1:end);
    end
    
    for ind = 2:p % try to impose smoothness
        delta =  param(ind)-param(ind-1);
        if delta > pi/2
            param(ind) = param(ind) - pi;
        elseif delta < -pi/2
            param(ind) = param(ind) + pi;
        end
    end
end

if bZ2 % Z2 symmetry for all X rotation
    betas = mod(param(p+1:end), pi/2);

    for ind = 2:p % try to impose smoothness
        delta =  betas(ind)-betas(ind-1);
        if delta > pi/4
            betas(ind) = betas(ind) - pi/2;
        elseif delta < -pi/4
            betas(ind) = betas(ind) + pi/2;
        end
    end
    
    param(p+1:end) = betas;
end
