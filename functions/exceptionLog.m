function [ errorMessage ] = exceptionLog( exception )
%EXCEPTIONLOG vrati detailni popis chyby z exception


    errorMessage = sprintf('** Error Message:%s, identifier: %s \n', exception.message, exception.identifier);
    for s = 1:numel(exception.stack)
        errorMessage = [ errorMessage sprintf('** Error in function %s() at line %d of file %s\n', exception.stack(s).name, exception.stack(s).line,strrep(exception.stack(s).file,'\','\\')) ]; %#ok<AGROW>
    end

end

