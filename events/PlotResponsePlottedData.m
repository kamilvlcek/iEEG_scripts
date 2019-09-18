classdef (ConstructOnLoad) PlotResponsePlottedData < event.EventData
   properties
      plottedChannel
   end
   
   methods
      function data = PlotResponsePlottedData(ch)
         data.plottedChannel = ch;
      end
   end
end
