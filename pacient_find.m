function [ p ] = pacient_find( pacienti,nick )
%PACIENT_FIND hleda pacienta ve strukture pacienti podle nick
%   vrati jeho index

nalezen = false;
for p = 1:numel(pacienti)
    if strfind(pacienti(p).folder,nick)
        nalezen = true;
        break; %nasel jsem pacienta
    end        
end
if ~nalezen  
    p = -1;
end

end

