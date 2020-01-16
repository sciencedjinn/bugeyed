function r=crosscor(inX, inY);

mx=mean(inX);
my=mean(inY);
for d=0:length(inX)-1
    x=circshift(inX, [0 d]);
    y=inY;
    r(d+1)=(sum((x-mx).*(y-my))) / (sqrt(sum((x-mx).^2)).*sqrt(sum((y-my).^2)));
end
