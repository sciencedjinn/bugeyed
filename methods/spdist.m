function d = spdist(az1, el1, az2, el2)
% SPDIST(azi1, ele1, azi2, ele2) computes the spherical distance of two points on a unit sphere
% Input and output are in degrees, and can be multidimensional

% used by 'findPoint'
% This is equivalent to calling distance(el1, az1, el2, az2) from the Mapping toolbox

eta1 = pi/2-pi/180.*el1;
eta2 = pi/2-pi/180.*el2;
phi1 = pi/180.*az1;
phi2 = pi/180.*az2;

dist = acos(sin(eta1).*sin(eta2).*cos(phi1).*cos(phi2) + sin(eta1).*sin(eta2).*sin(phi1).*sin(phi2)+cos(eta1).*cos(eta2));
d = 180/pi.*dist;
