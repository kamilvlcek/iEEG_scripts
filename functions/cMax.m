function [amax, idx, idxFrac] = cMax(curve, fraction)
    % Nalezne maximum / minimum krivky curve, jeho hodnotu a
    % pozici, a pozici poloviny maxima
    [testMax, idxMax] = max(curve);
    [testMin, idxMin] = min(curve);
    if abs(testMax-0.5) > abs(0.5-testMin)
        amax = testMax;
        idx = idxMax;
        idxFrac = find(curve > fraction*(amax-0.5)+0.5, 1, 'first'); % prvni vyskyt poloviny maxima
    else
        amax = testMin;
        idx = idxMin;
        idxFrac = find(curve < fraction*(amax-0.5)+0.5, 1, 'first'); % prvni vyskyt poloviny maxima
    end
end
