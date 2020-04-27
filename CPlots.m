classdef CPlots < matlab.mixin.Copyable
    %CPLOTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Eh; %handle to CiEEGData object       
        PlotRChMean; %data for  PlotResponseChMean        
    end
    methods (Access = public)
        function obj = CPlots(E) %constructor
            obj.Eh = E; %save handle to main object
        end 
        function PlotResponseChMean(obj,kategories,channels)
            if ~exist('kategories','var') || isempty(kategories)                 
                WpA = obj.Eh.WpActive; %jen zkratka                
                if ~isempty(obj.Eh.Wp) && isfield(obj.Eh.Wp(WpA), 'kats')
                    kategories = obj.Eh.Wp(WpA).kats; %pokud jsou kategorie v parametru, prvni volba je pouzit je ze statistiky
                else
                    kategories = obj.Eh.PsyData.Categories();
                end
            end
            if ~exist('channels','var') || isempty(channels) 
                channels = obj.Eh.CH.plotCh2D.chshow;                
                figuretitle = [obj.Eh.CH.plotCh2D.chshowstr ' chns: ' num2str(numel(channels))];                
            else
                figuretitle = ['channels: ' num2str(numel(channels))];                
            end
            
            obj.PlotRChMean.channels = channels;
            obj.PlotRChMean.kategories = kategories;
            
            chnshow = mat2str(channels(1:(min(20,numel(channels)))));
            if numel(channels) > 20, chnshow = [chnshow ' ...']; end
            popis = ['ChShow:  ' obj.Eh.CH.plotCh2D.chshowstr '=' chnshow];
            
            hue = 0.8;
            colorsErrorBars = cellfun(@(a) min(a+hue, 1), obj.Eh.colorskat, 'UniformOutput', false);
            katlinewidth = 2;
            if isfield(obj.PlotRChMean,'fh') && (verLessThan('matlab','9.0') || isvalid(obj.PlotRChMean.fh)) %isvalid je od verze 2016
                figure(obj.PlotRChMean.fh); %pouziju uz vytvoreny graf
                clf(obj.PlotRChMean.fh); %graf vycistim
            else
                obj.PlotRChMean.fh = figure('Name','PlotResponseChMean','CloseRequestFcn', @obj.tearDownFigCallbackPlotResponseChMean);
            end
                        
            T = linspace(obj.Eh.epochtime(1),obj.Eh.epochtime(2),size(obj.Eh.d,1)); %od podnetu do maxima epochy. Pred podnetem signifikanci nepocitam
            ymax = 0;
            for k = 1 : numel(kategories) %index 1-3 (nebo 4)
                katnum = kategories(k);
                colorkatk = [obj.Eh.colorskat{katnum+1} ; colorsErrorBars{katnum+1}]; %dve barvy, na caru a stderr plochu kolem                
                [katdata,~,RjEpCh] = obj.Eh.CategoryData(katnum,[],[],channels);%katdata now time x x channels epochs                
                CHM = zeros(numel(channels),size(katdata,1));
                for ich = 1:numel(channels)
                    CHM(ich,:) = mean(katdata(:,channels(ich),~RjEpCh(1,:)),3);
                end
                M = mean(CHM,1);               
                E = std(CHM,[],1)/sqrt(size(CHM,1)); %std err of mean                
                ymax = max(ymax,max(M+E));
                ciplot(M+E, M-E, T, colorkatk(2,:)); %funguje dobre pri kopii do corelu, ulozim handle na barevny pas 
                hold on;
                plot(T,M,'LineWidth',katlinewidth,'Color',colorkatk(1,:));  %prumerna odpoved,  ulozim si handle na krivku                      
            end            
            xlim(obj.Eh.epochtime(1:2)); 
            title(figuretitle);
            text(0,0.99*ymax,popis, 'FontSize', 10);
            obj.PlotRChMean.filterListener = addlistener(obj.Eh.CH, 'FilterChanged', @obj.filterChangedCallbackPlotResponseChMean);
        end
        function filterChangedCallbackPlotResponseChMean(obj,~,~)   %update chart if the filter is changed         
            if ~isequal(obj.Eh.CH.plotCh2D.chshow, obj.PlotRChMean.channels) %jinak se to z nejakeho duvodu po zmene filtru vola porad dokolecka
                obj.PlotResponseChMean(obj.PlotRChMean.kategories);
            end
        end
        function tearDownFigCallbackPlotResponseChMean(obj,src,~)
            delete(obj.PlotRChMean.filterListener);
            delete(src);
        end  
    end
    methods (Static,Access = public)       
        function PlotElectrodeEpiEvents(els,RjCh,DE,objtabs,tabs_orig,epochs,samples,epochtime,timeax,timeaxy,sec,time_n,elmaxmax,shift,iD)
            %27.4. - nahrazuju cast funkce CiEEGData.PlotElectrode toutu funkci, kvuli zkraceni, ale predavam strasne moc argumentu
            %treba bych mohl uz driv predat hodne z nich? A mit je jako properties tehle tridy?
            epieventsum = 0;  
            weights = [];
            DE.Clear_iDEtabs();
            for ch = els
                if ~ismember(ch,RjCh)
                    if epochs <= 1 %neepochovana data
                        [epitime, weight] = DE.GetEvents( [objtabs(iD(1)) objtabs(iD(2))],ch,tabs_orig(1)); 
                    else
                        epochy = sec : min(sec+ceil(time_n/samples)-1 , epochs ); %cisla zobrazenych epoch
                        tabs = [ objtabs(1,epochy)' objtabs(end,epochy)' ]; %zacatky a konce zobrazenych epoch
                        [epitime, weight] = DE.GetEvents( tabs,ch,tabs_orig(1)) ;                            
                        if numel(epitime) > 0
                            epitime(:,1) = epitime(:,1) + epochtime(1);
                        end
                    end
                    if size(epitime,1) > 0
                        plot(epitime(:,1),shift(elmaxmax-ch+1,1),'o', 'MarkerEdgeColor','r','MarkerFaceColor','r', 'MarkerSize',6);
                        epieventsum = epieventsum + size(epitime,1);
                        weights = [weights; round(weight*100)/100]; %#ok<AGROW>
                    end
                end
            end             
            text(timeax(1)+((timeax(end)-timeax(1))/8),timeaxy,['epileptic events: ' num2str(epieventsum) ' (' num2str(min(weights)) '-' num2str(max(weights)) ')']);         
        end
        function spatne = EpochsEpi(RjEpochCh,objels,CHH)
            figure('Name','Rejected epochs in individual channels');
            pocty = sum(RjEpochCh,2)/size(RjEpochCh,2); %pocty epoch v jednotlivych kanalech / celkovym poctem epoch
            plot(pocty,'.-');
            for el = 1:numel(objels)-1
                line([objels(el) objels(el)]+1,[0 1],'Color',[0.5 0.5 0.5]);
            end
            line([1 CHH.H.selCh_H(end)],[0.30 0.30],'Color','red');
            spatne = find(pocty >= 0.30); 
            disp(['kanaly s podilem vyrazenych epoch >= 0.3:' mat2str(spatne)]);                
        end
        function PlotEpochsKatPPA(E,ch,s)
            %udaje z PPAconditionsEEG_info.xlsx
            %unikatni cisla 325 obrazku, tak jak sly za sebou
            %ch - kanal k vykresleni
            %s - sample k vykresleni (okamzik v case)
            jpgnum = [168 46 121 231 303 4 102 178 203 283 70 181 149 199 317 99 2 76 44 226 210 175 81 64 80 165 143 47 247 126 240 190 185 180 318 306 251 26 272 96 69 273 59 278 195 188 289 218 184 41 74 48 161 3 63 205 182 235 8 100 141 248 298 255 157 284 109 106 40 145 65 164 67 130 51 228 103 319 111 155 139 287 73 146 34 62 5 253 68 75 214 71 325 232 250 20 204 276 224 305 91 201 233 55 163 200 60 260 77 294 206 215 245 236 266 194 187 11 85 110 220 54 167 186 316 312 118 137 95 135 323 144 254 264 244 92 241 58 173 208 45 324 134 290 90 9 209 160 97 222 18 14 219 300 93 280 243 38 36 176 52 17 147 56 309 170 191 35 119 33 249 216 279 322 292 72 258 89 275 133 50 117 108 107 13 189 152 39 193 177 263 285 22 88 87 24 66 313 131 304 230 301 153 297 259 125 286 113 239 140 198 267 227 308 257 84 223 183 269 171 94 83 132 86 158 142 211 150 154 129 21 270 237 42 242 25 234 28 217 15 197 23 192 256 202 29 315 225 123 1 101 311 265 321 16 114 7 179 57 151 53 27 196 271 320 282 159 31 274 37 212 172 299 262 19 207 112 120 238 295 82 116 12 156 246 302 79 61 261 104 281 30 138 314 32 148 310 122 124 136 115 128 268 213 98 162 221 291 293 127 288 169 49 6 252 229 10 277 166 78 174 307 105 43 296 ];
            
            %unikatni cisla umelych scen, v nespravnem poradi
            artificial = [231 226 247 240 251 272 273 235 248 255 228 253 232 250 233 260 245 236 266 254 264 244 241 243 249 258 275 263 230 259 239 267 227 257 269 270 237 242 234 256 265 271 274 262 238 246 261 268 252 229 ];
            iArt = ismember(jpgnum,artificial); %indexy v jpgnum, kde jsou cisla z artificial
            %unikatni cisla prirozenych scen
            natural = [303 283 317 318 306 278 289 298 284 319 287 325 276 305 294 316 312 323 324 290 300 280 309 279 322 292 285 313 304 301 297 286 308 315 311 321 320 282 299 295 302 281 314 310 291 293 288 277 307 296  ];
            iNat = ismember(jpgnum,natural);
            
            figure;
            d = squeeze(E.d(s,ch,1:325)); %jen prvnich 325 epoch
            x = 1:325; %fixne potech epoch
            plot(x,d);
            hold on;
            plot(x(iNat),d(iNat),'og');
            nmean = mean(d(iNat));
            nerr = std(d(iNat))/sqrt(sum(iNat));
            plot(x(iArt),d(iArt),'*r');
            amean = mean(d(iArt));
            aerr = std(d(iArt))/sqrt(sum(iArt));
            line([1 325],[nmean nmean],'Color',[0 1 0]);
            line([1 325],[nmean+nerr nmean+nerr],'Color',[0 0.5 0]);
            line([1 325],[amean amean],'Color',[1 0 0]);
            line([1 325],[amean+aerr amean+aerr],'Color',[0.5 0 0]);
            ylim([-2 2]);
            fprintf('Art: %f +- %f, Nat %f +- %f\n',amean,aerr,nmean,nerr);
        end

    end
    
end

