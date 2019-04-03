function [amax, idx, idxHalf] = cMax(curve)
    % Nalezne maximum / minimum krivky curve, jeho hodnotu a
    % pozici, a pozici poloviny maxima
    [testMax, idxMax] = max(curve);
    [testMin, idxMin] = min(curve);
    if abs(testMax-0.5) > abs(0.5-testMin)
        amax = testMax;
        idx = idxMax;
        idxHalf = find(curve > (amax+0.5)/2, 1, 'first'); % prvni vyskyt poloviny maxima
    else
        amax = testMin;
        idx = idxMin;
        idxHalf = find(curve < (amax+0.5)/2, 1, 'first'); % prvni vyskyt poloviny maxima
    end
end
