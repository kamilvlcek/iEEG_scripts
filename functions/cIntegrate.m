function t = cIntegrate(time, curve, fraction, normMode)
    % Nalezne pozici idx, ve ktere je integral I(0,idx) = fraction*I(0,end)

    testMax = max(curve);
    testMin = min(curve);
    
    if abs(testMax-0.5) > abs(0.5-testMin)
        curve = curve - 0.5;
    else
        curve = 0.5 - curve;
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
