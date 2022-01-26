function outPoints = deprojectV(inPoints, inView)
% DEPROJECTV(inPoints, inView) deproject points from the plane back onto the sphere.
% Viewing Angle inView in spherical coordinates (degrees)
% inPoints in spherical (degrees), outPoints in cartesian coordinates

% used by 'findNext'

[X, Y, Z] = sph2cart(pi/180*inPoints(:, 1), pi/180*inPoints(:, 2), 1);
% rotate to make View = [0 0]

% 1. rotate around Z by azi to make azi = 0
t = pi/180*inView(1);
RotMatrixZ = [cos(t) sin(t) 0; -sin(t) cos(t) 0; 0 0 1];

% 2. then rotate around y by -ele to make ele = 0
t2 = -pi/180*inView(2);
RotMatrixY = [cos(t2) 0 -sin(t2); 0 1 0; sin(t2) 0 cos(t2)];

% 3. do it
for i = 1:length(X)
    RotatedPoints(i, :) = (RotMatrixZ*RotMatrixY*[X(i); Y(i); Z(i)])';
end

% View = [0 0]

[outPoints(:, 1), outPoints(:, 2), R] = cart2sph(RotatedPoints(:, 1), RotatedPoints(:, 2), RotatedPoints(:, 3));
outPoints = 180/pi*outPoints;
