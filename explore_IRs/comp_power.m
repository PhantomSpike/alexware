function [Lp,power_ratio] = comp_power(data1,data2)
%Compare the total power of two signals
data1 = data1(:);
data2 = data2(:);
p1 = sum(data1.^2);
p2 = sum(data2.^2);
Lp = 10*log10(p1/p2);
power_ratio = 100*(p1/p2);
