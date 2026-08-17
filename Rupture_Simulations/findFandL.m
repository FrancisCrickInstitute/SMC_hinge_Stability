function F = findFandL(t)
% this function finds F1, F2, L1 and L2 that satisfy all conditions

%               o               -
%               |               ^
%               |               x0
%             / o \
%       /----      ----\
%   o--                 --o
%   |<--------- D ------->|

global L10 L20 D x0

L1 = t(1);
L2 = t(2);
F1 = t(3);
F2 = t(4);

F = zeros(4,1);

F(1) = F2 - 0.0828*(0.25/(1-L2/L20)/(1-L2/L20)-0.25+L2/L20);
F(2) = F1 - 0.0828*(0.25/(1-L1/L10)/(1-L1/L10)-0.25+L1/L10);
% 0.0828 is kT/Lp = 4.14/50 pN for Lp=50 nm (dsDNA)
alpha = atan(2*(x0-L1)/D);
F(3) = F1 - 2*F2*sin(alpha);
F(4) = L2*sin(alpha)/2 + L1 - x0;

