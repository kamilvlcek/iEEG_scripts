function indexes = findZeros(d, m, pl)
% Funkcia na vyhladanie postupky nul v EEG datach
% d = povodne EEG data
% m = kolko minimalne nul za sebou chceme najst
% pl = plot true/false
%
% indexes = vrati indexy pozicii, ktorym predchadza (m-1) nul za sebou 
% pre kazdy channel zvlast

    % d = time x electrodes
    [row, col] = find(d==0);

    t = zeros(size(d));
    t(sub2ind(size(t),row, col)) = 1;

    channels = size(t,2);
    orig = t;
    result = t;
    for i = 1:(m-1)
        result = [result; zeros(1,channels)]; % pridam na koniec jeden riadok nul, aby som mohla pouzit maticovy and
        orig = [zeros(1,channels); orig]; % povodnu maticu posuvam o jedno a porovnavam s result
        result = result & orig;
    end
    
    [r, c] = find(result == 1);
    indexes = [r, c];
    
    if pl
        figure;
        plot(1:channels, sum(result)); 
        xlabel('Channels'); ylabel('Occurences');

        figure;
        scatter(c, r, 'x');
        xlabel('Channels'); ylabel('Time points');
    end
end
