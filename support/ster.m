function st=ster(ang)
%STER(ang) computes the steradiant of a cone with cone angle ang (in degrees)

st=pi*(2-2*cos(pi/360*ang)+sin(pi/360*ang).^2);