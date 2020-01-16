function r=int(i);
%INT(i) returns 1 if i is an integer, 0 otherwise
%Returns 1, as well, for integers in
%floating point representation

%test for correct input format
if nargin==0
    error('*** Error: Not enough input arguments');
end
if ~isa(i, 'numeric')
    error('*** Error: Function "int" is only defined for numerics');
end


%Main function
if round(i)==i
    r=1;
else
    r=0;
end
    