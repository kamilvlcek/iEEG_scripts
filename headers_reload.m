function [ nahrazeno ] = headers_reload( testname,filename )
%HEADERS_RELOAD funkce vymeni headery ve vsech souborech za aktualni
%   zkontroluje, jestli je header stejny
[ pacienti, setup  ] = pacienti_setup_load( testname );
nahrazeno = {};
for p = 1:numel(pacienti) % cyklus pacienti
    if pacienti(p).todo       
        headerfile = [setup.basedir pacienti(p).folder '\' pacienti(p).header];
        if(exist(headerfile,'file')~=2)
            msg = ['Header neexistuje: ' pacienti(p).folder '\\' pacienti(p).header];
            warning(msg);
        else
            E = pacient_load(pacienti(p).folder,testname,filename);
            load(headerfile); %nacte H
            if ~isfield(E.CH.H,'filename') || isempty(E.CH.H.filename) || ~strcmp(E.CH.H.filename,pacienti(p).header)
                %jen pokud je jmeno header jine
                if isprop(E,'reference') && ~isempty(E.reference) 
                    if strcmp(E.reference,'Bipolar') %u bipolarni reference je jiny pocet kanalu 
                        CH = CHHeader(H);
                        CH.RejectChannels( pacienti(p).rjch); %musim vyradit vyrazene kanaly, protoze ty se vyrazuji v bipolarni referenci
                        CH.ChangeReference('b'); %ostatni reference zatim neresim, nepouzivam
                        H = CH.H;
                        RjCh = CH.RjCh;
                    elseif strcmp(E.reference,'original') 
                        CH = CHHeader(H);
                        CH.RejectChannels( pacienti(p).rjch);
                        H = CH.H;
                        RjCh = CH.RjCh;
                    else
                        warning(['neznama reference ' E.reference]);
                        continue; %nechci nacitat novy header
                    end
                end                
                E.GetHHeader(H,pacienti(p).header);
                E.RejectChannels(RjCh);
                E.Save();
                disp(['**** ' pacienti(p).folder ' - header nahrazen: ' ,pacienti(p).header]);
                nahrazeno = [ nahrazeno; [pacienti(p).folder ':' pacienti(p).header] ]; %#ok<AGROW>
            end
        end
    end
end
end

