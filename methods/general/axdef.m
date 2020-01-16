%AXDEFPAT: properties of the current axes are set to default values (only upper right
%quadrant is shown)


axis equal;
axis([0 180 0 90]);
set(gca,'xtick',[0 30 60 90 120 150 180]);
set(gca,'ytick',[0 30 60 90]);
xlabel('azimuth [\circ]'); 
ylabel('elevation [\circ]');