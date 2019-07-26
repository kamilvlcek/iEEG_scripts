function [amax, idx, idxFrac] = cMax(curve, fraction, centering)
    % CMAX Nalezne maximum / minimum krivky curve, jeho hodnotu a
    % pozici, a pozici poloviny maxima. 

    switch nargin
        case 1
            centering = 0;
            fraction = 0.5;
        case 2
            centering = 0;
        otherwise
    end
    
    [testMax, idxMax] = max(curve);
    [testMin, idxMin] = min(curve);
    if abs(testMax-centering) > abs(centering-testMin)
        amax = testMax;
        idx = idxMax;
        %TODO: Zobecnit na 2D vstup - nekolik krivek curve ve sloupcich (find v tom pripade nefunguje)
        idxFrac = find(curve > fraction*(amax-centering)+centering, 1, 'first'); % prvni vyskyt poloviny maxima
    else
        amax = testMin;
        idx = idxMin;
        idxFrac = find(curve < fraction*(amax-centering)+centering, 1, 'first'); % prvni vyskyt poloviny maxima
    end
end
