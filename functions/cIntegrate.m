function t = cIntegrate(time, curve, fraction, normMode, centering)
    % Nalezne pozici idx, ve ktere je integral I(0,idx) = fraction*I(0,end)
    % Hleda se integral +curve, nebo -curve v zavislosti na tom, jestli
    % ma krivka vzdalenejsi maximum nebo minimum od hodnoty centering.
    % Hodnota normMode urcuje sekundarni transformaci krivky pred integraci:
    %   0: krivka se nijak netransformuje
    %   1: odstrani se zaporne hodnoty (tj. ty na druhou stranu od centering)
    %   2: minimum krivky se posune do nuly

    testMax = max(curve);
    testMin = min(curve);
    
    if abs(testMax-centering) > abs(centering-testMin)
        curve = curve - centering;
    else
        curve = centering - curve;
    end

    if normMode == 1  % odstrani zaporne hodnoty
        curve(curve < 0) = 0;
    elseif normMode == 2  % posune minimum krivky do nuly
        curve = curve - min(curve);
    end
    
    cint = cumtrapz(time, curve);
    
    idx = find(cint < fraction*cint(end), 1, 'last');
    t = time(idx);

end
