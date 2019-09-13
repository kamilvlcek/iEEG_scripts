function [amax, idx, idxFrac] = cMax(curve, fraction, centering,iTime)
    % CMAX Nalezne maximum / minimum krivky curve, jeho hodnotu a
    % pozici, a pozici poloviny maxima. 
    % iTime - kamil 13.9.2019 - indexy krivky, ze kterych se ma max/min pocitat

    if ~exist('fraction','var'), fraction = 0.5; end
    if ~exist('centering','var'), centering = 0; end
    if ~exist('iTime','var') || ~any(iTime) %pokud neni casovy index pouzit nebo pokud je prazny (same false)
        iTime = true(size(curve));
    end    
    
    [testMax, idxMax] = max(curve(iTime));
    [testMin, idxMin] = min(curve(iTime));
    fiTime = find(iTime); 
    idxMax = fiTime(idxMax); %z relativnich indexu v ramci iTime udelam absolutni v ramci curve
    idxMin = fiTime(idxMin);
    if abs(testMax-centering) > abs(centering-testMin)
        amax = testMax;
        idx = idxMax;
        %TODO: Zobecnit na 2D vstup - nekolik krivek curve ve sloupcich (find v tom pripade nefunguje)
        idxFrac = find(curve(iTime) > fraction*(amax-centering)+centering, 1, 'first'); % prvni vyskyt poloviny maxima
    else
        amax = testMin;
        idx = idxMin;
        idxFrac = find(curve(iTime) < fraction*(amax-centering)+centering, 1, 'first'); % prvni vyskyt poloviny maxima
    end
    idxFrac = fiTime(idxFrac); %z relativnich indexu v ramci iTime udelam absolutni v ramci curve    
end
