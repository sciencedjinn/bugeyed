function d = andist(a, b)
% calculates the difference of the angles (a-b) in degrees

dplus = mod((a-b), 360);
D = [dplus dplus-360];
[~, i] = min(abs(D));
d = D(i);