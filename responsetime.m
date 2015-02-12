function [ latencies ] = responsetime( ALLEEG,dataset )
%RESPONSETIME Summary of this function goes here
%   Detailed explanation goes here


d = dataset;
latencies = zeros(size(ALLEEG(d).epoch,2),1);
    for e = 1:size(ALLEEG(d).epoch,2)
        for ev = 1: size(ALLEEG(d).epoch(e).eventtype,2);
            if strcmp(ALLEEG(d).epoch(e).eventtype{ev},'response')==1
                latencies(e)=ALLEEG(d).epoch(e).eventlatency{ev};
            end
        end
    
    end
    %plot(latencies);
end

