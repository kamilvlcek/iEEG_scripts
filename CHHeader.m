classdef CHHeader
    %CHHEADER Trida na praci s headerem od Jirky Hammera
    
    properties (Access = public)
        H; %data od Jirky Hammera
        E; % unikatni jmena elektrod, napr RG, I, GC ...
    end
    
    methods (Access = public)
        function obj = CHHeader(H)
            %konstruktor
            obj.H = H;     
              %nefunguje protoze name je napriklad RG12 a ne jen RG, takze neni unikatni pro kazdou elektrodu
%             E = cell(size(H.electrodes,2),1); %#ok<*PROP>
%             for iE= 1:size(H.electrodes,2)
%                 E{iE}=H.electrodes(iE).name;
%             end
%             [U,idx]=unique(cell2mat(E)); 
%             obj.E = cell(numel(U),1);
%             for iE = 1:numel(idx)
%                 obj.E{iE}= E{idx};
%             end
                
        end       
       
    end
    
end

