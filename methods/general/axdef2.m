%AXDEF2: properties of the current axes are set to default values

axis equal;
set(gca, 'xtick', [-180 -150 -120 -90 -60 -30 0 30 60 90 120 150 180]);
set(gca, 'ytick', [-90 -60 -30 0 30 60 90]);
axis([-180 180 -90 90]);
xlabel('azimuth [\circ]'); 
ylabel('elevation [\circ]');