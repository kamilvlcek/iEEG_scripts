%Allo
Interval1=[575 625]; %ms
Interval2=[650 700];
SetAllo = 2;

M1 = means(ALLEEG(1,SetAllo).data,floor((Interval1+200)/1000*512));
M2 = means(ALLEEG(1,SetAllo).data,floor((Interval2+200)/1000*512));
Allo = [M1 M2];

SetEgo = 3;
M1 = means(ALLEEG(1,SetEgo).data,floor((Interval1+200)/1000*512));
M2 = means(ALLEEG(1,SetEgo).data,floor((Interval2+200)/1000*512));
Ego = [M1 M2];

figure('Name','575 625ms')
plot( ([Allo(:,1) Ego(:,1)]));

figure('Name','650 700ms')
plot( ([Allo(:,2) Ego(:,2)]));

