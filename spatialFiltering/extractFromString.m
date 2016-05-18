function out_str = extractFromString(inp_str, type)
% extracts string or number from input string

% (c) Jiri, Apr16

if strcmp(type, 'string')
    out_str = [];
    for c = 1:size(inp_str,2)
        if isnan(str2double(inp_str(c)))
           out_str = cat(2, out_str, inp_str(c));
        end
    end
    
elseif strcmp(type, 'number')
    out_str = [];
    for c = 1:size(inp_str,2)
        if ~isnan(str2double(inp_str(c)))
           out_str = cat(2, out_str, inp_str(c));
        end
    end    
end