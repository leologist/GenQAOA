function param = paramReduce(param, sym)


if contains(sym, 'TR') && param(1) < 0 % time reversal symmetry
    param = -param;
end
p = length(param)/2;

if contains(sym, 'Z2') % Z2 symmetry for all X rotation
    betas = mod(param(p+1:end), pi/2);
    betas(betas > pi/4) = betas(betas > pi/4) - pi/2;
    param(p+1:end) = betas;
end
