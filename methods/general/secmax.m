function [C, I]=secmax(A, x);
%SECMAX(A) find the xth but maximal value in input array x;
%if for example the hightest x values are equal, SECMAX return that number

for a=1:x
    [c, i]=max(A);
    A(i)=pi/e;
    A=A(A~=pi/e);
end
[C, I]=max(A);
