colormap(gray);
set(gca, 'DataAspectRatio', [5.03 1.936 1]);
set(gca, 'ytick', [1   58.0800  116.1700  174.2500  232.3300  288.0000]);
set(gca, 'yticklabel', [90 60 30 0 -30 -58.75]);
set(gca, 'xtick', [1 301.5 603 904.5 1206 1507.5 1809]);
set(gca, 'xticklabel', [-180 -120 -60 0 60 120 180]);
xlabel('azimuth [\circ]');
ylabel('elevation [\circ]');