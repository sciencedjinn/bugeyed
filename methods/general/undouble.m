function R=undouble(R)

R1=R(:, 1);
R2=R(:, 2);
i=1;
while i<=length(R1)
    a=R1(i);
    R1(i)=55555;
    R1(abs(R1-a)<0.01 & abs(R2-R2(i))<0.01)=66666;
    R2=R2(R1~=66666);
    R1=R1(R1~=66666);
    R1(i)=a;
    i=i+1;
end
R=cat(2, R1, R2);