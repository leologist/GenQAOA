function paramOut = interpParam(paramIn, p)

    pIn = length(paramIn)/2;
    xpIn = linspace(0,1,pIn);
    xpOut = linspace(0,1,p);
    paramOut = [interp1(xpIn, paramIn(1:pIn), xpOut, 'linear'), ...
                 interp1(xpIn, paramIn(pIn+1:end), xpOut, 'linear')];

end