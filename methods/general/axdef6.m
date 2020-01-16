%AXDEF6: properties of the current axes are set to default values
%FLOW FIELDS

axis equal;
axis([-195 195 -90 90]);
set(gca,'xtick',[-180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180]);
set(gca,'ytick',[-90 -60 -30 0 30 60 90]);
xlabel('azimuth [\circ]'); 
ylabel('elevation [\circ]');