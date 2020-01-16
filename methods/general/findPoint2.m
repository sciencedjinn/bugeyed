function outNumber=findPoint2(azi, ele, G1, G2)
%FINDPOINT(azi, ele, G1, G2) finds the point (azi, ele) in the array given by azi-values in
%G1 and ele-values in G2.
%if there is no Point in a 0.001-neighborhood, findPoint return NaN,
%otherwise the position number(s) of the point(s)

%uses 'spdist'

%used by 'findNext'

[c, outNumber]=min(spdist(azi, ele, G1, G2));
%if c>0.1
%    outNumber=NaN;
%    disp(['** Warning: Point ', num2str(azi), ' ', num2str(ele), ' not found in array!']);
%end