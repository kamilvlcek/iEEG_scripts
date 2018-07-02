classdef CSelCh < matlab.mixin.Copyable 
    %CSelCh Trida uchovavajici a poskutujici seznamy vybranych kanalu pro soubory CHIlbertMulti aj
    %   pred frekvence, kontrasty aj
    
    properties (Access = public)
        selCh;
        n; %pocet ulozenych dat
    end
    
    methods (Access = public)
        function obj = CSelCh()
            obj.selCh = cell(1,4); %prvni sloupec je filename, druhy selCh, treti katstr, ctvrty freq 
            obj.n = 0;
        end
        function obj = SetSelCh(obj,selCh,filename,katstr,freq)            
            %ulozi vyber kanalu            
            if ~exist('katstr','var'), katstr = ''; end
            if ~exist('freq','var'), freq = ''; end
            obj.selCh(obj.n+1,:) = {filename,selCh,katstr,freq};
            obj.n = obj.n +1;
        end
        function selCh = GetSelCh(obj,filename)
            s = find(~cellfun(@isempty,strfind(obj.selCh(:,1),filename)),1); %najdu pouze prvni vyhovujici soubor
            selCh = obj.selCh{s,2};
        end
    end
    
end

