classdef CPacients < handle
    %CPACIENTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        Property1
    end
    
    methods (Access = public)
        function obj = CPacients()
            %CPACIENTS Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function List = ListAll(obj,tests)
            %LISTALL returs struct with all pacients across all mentioned tests
            %   Detailed explanation goes here
            for t = 1:numel(tests)
                
            end
        end
    end
end

