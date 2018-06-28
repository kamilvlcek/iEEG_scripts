function [str]= basename(filename)
    %BASENAME vraci filename s koncem path pro identifikaci pacienta
    %
    if isempty(filename)
        str = filename;
        return;
    end
    [~,basen,ext] = fileparts(filename); %takhle to nechci
    %fslash = strfind(filename,'\');
    %str = filename(fslash(end-2)+1:end); %dve casti path pred basename      
    str = [basen ext];
    
end

