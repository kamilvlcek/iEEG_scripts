function [ output_args ] = ppalog( casy )
%PPALOG zobrazi vystup z ppalog.php - ITI z PPA testu
%   Detailed explanation goes here
d = casy(2:end)-casy(1:end-1);
d(d>2) = [];
%figure();
%hist(d,[0.5:0.01:1.5]);
figure();
plot(d,'.');
ylim([1 1.6]);
end

