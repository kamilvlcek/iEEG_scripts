function [t,aint] = cIntegrate(time, curve, fraction, normMode, centering,iTime)
    % Nalezne pozici idx, ve ktere je integral I(0,idx) = fraction*I(0,end)
    % Hleda se integral +curve, nebo -curve v zavislosti na tom, jestli
    % ma krivka vzdalenejsi maximum nebo minimum od hodnoty centering.
    % Hodnota normMode urcuje sekundarni transformaci krivky pred integraci:
    %   0: krivka se nijak netransformuje
    %   1: odstrani se zaporne hodnoty (tj. ty na druhou stranu od centering)
    %   2: minimum krivky se posune do nuly
    % iTime - kamil 13.9.2019 - indexy krivky, ze kterych se ma fraction pocitat

    if ~exist('iTime','var') || ~any(iTime) %pokud neni casovy index pouzit nebo pokud je prazny (same false)
        iTime = true(size(curve));
    end
    fiTime = find(iTime);     
    
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
    if sum(iTime)>1 %funkce cumtrapz musi mit alespon 2 prvni v array
        cint = cumtrapz(time(iTime), curve(iTime));   %area under significant parts of the curve 
        idx = find(cint < fraction*cint(end), 1, 'last');
        idx = fiTime(idx);  %z relativnich indexu v ramci iTime udelam absolutni v ramci curve
        if isempty(idx)
            t = 0; % t cant be empte
        else
            t = time(idx);            
        end
        aint = cint(end); %cint is cummulative vector
    else
        aint = 0; %area is zero
        t = 0; %nic se nenaslo
    end

end
