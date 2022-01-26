function makePatches(inGeomName)

vFile = bugeyed_fileName(inGeomName, 'voronoi');
pFile = bugeyed_fileName(inGeomName, 'patches');

%% load Voronoi and convert to a quick-to-plot format
load(vFile, 'outVOR');

polSize = cellfun(@(x) size(x, 1), outVOR); % number of vertices for each polygon
for np = min(polSize):max(polSize)
    % find all polygons with np vertices
    theseVerts = find(polSize==np);

    Fx{np} = cell2mat(cellfun(@(x) x(:, 1), outVOR(theseVerts), 'UniformOutput', false)); %#ok<*AGROW>
    Fy{np} = cell2mat(cellfun(@(x) x(:, 2), outVOR(theseVerts), 'UniformOutput', false));
    Id{np} = theseVerts;
end

save(pFile, 'Fx', 'Fy', 'Id')