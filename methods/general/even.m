function r=even(s)
%EVEN(s) returns 1 if the input is an even integer, 0 if any other number (odd
%integer or floating point). Returns 1, as well, for even integers in
%floating point representation

%test for correct input format
if nargin==0
    error('*** Error: Not enough input arguments');
end
if ~isa(i, 'numeric')
    error('*** Error: Function "even" is only defined for numerics');
end

%Main program
if 2*round(s/2)==s
    r=1;
else
    r=0;
end