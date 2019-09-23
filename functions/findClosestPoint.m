function [minIndex, minDistance] = findClosestPoint(p1, p2, points, tolerance)
%findClosestPoint(p1, p2, points, tolerance) 
% Vraci dvojici [nejblizsi kanal, vzdalenost od primky]. Pokud je vzdalenost
% vetsi, nez nastavena tolerance, vraci [0 inf]

  minIndex = 0;     % index nejblizsiho kanalu
  ray = p1 - p2;    % raycast - primka mezi p1 a p2
  numPoints = size(points, 2);
  distancesRay = zeros(1, numPoints);   % vzdalenosti kanalu od raycast primky
  for i = 1:numPoints
    testPoint = points(:, i)';  % pozice testovaneho kanalu
    pdir = testPoint - p2;      % smerovy vektor k bodu na primce
    distancesRay(i) = norm(cross(ray, pdir)) / norm(ray);   % vzdalenost od primky
  end
  pointsNearLine = find(distancesRay < tolerance);  % vyberu pouze body blizko primky (ve vzdalenosti < tolerance)
  minDistance = inf;
  for i = pointsNearLine    % z bodu lezicich blizko primce najdu ten "nejbliz k obrazovce"
    testPoint = points(:, i)';  % pozice testovaneho kanalu
    distanceP1 = norm(p1 - testPoint);  % vzdalenost od mista kliknuti v "rovine obrazovky"
    if distanceP1 < minDistance    % pokud je bod nejlepsi nalezeny, ulozim jeho index
        minDistance = distanceP1;
        minIndex = i;
    end
  end
end