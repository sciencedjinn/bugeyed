function outVOR = makeVoronoi(inGeomName, overwrite)
% makeVORONOI(rootFolder, inGeomName, overwrite) calculates and saves a Voronoi Diagram of the Pattern given by the name string 'inGeomName', 
% if variable overwrite==1, all previous files will automatically be
% overwritten; if 2, a dialog will ask the user (default)

% used by 'makeAllVoronoi', 'showVoronoi'

%% paths
vFolder = fullfile(rootFolder, 'be_input', 'voronoi', inGeomName);
vFile = fullfile(vFolder, [inGeomName, 'Voronoi.mat']);
gFolder = fullfile(rootFolder, 'be_input', 'geom', inGeomName);
gFile = fullfile(gFolder, [inGeomName, 'Points.txt']);

if nargin<2
    overwrite=2;
end

if exist(vFolder, 'file')
    if overwrite==2
        button = questdlg(['Voronoi diagram of ', inGeomName, ' already exists. Overwrite?'], 'File already exists','Yes','No','No');
    elseif overwrite==1
        button='Yes';
    elseif overwrite==0
        button='No';
    else
        error('Invalid value for variable overwrite');
    end
    if strcmp(button, 'No')
        disp('** Warning: Voronoi diagram was NOT created. Process aborted by user');
        return;
    end
else
    mkdir(vFolder);
end

disp(['* Calculating Voronoi diagram ', inGeomName, ' ...']);
Geom = dlmread(gFile);
outVOR = voronoiPattern(Geom);

%and save the result
save(vFile, 'outVOR');
fprintf('\bdone\n');

function VOR = voronoiPattern(inGeom)
    % VORONOIPATTERN(inGeom) calculates the voronoi cells for all points in the input array inGeom
    h = waitbar(0, 'Calculating Voronoi diagram...(0%)','Name', 'Progress');
    
    VOR = cell(1, size(inGeom, 1));
    for i = 1:size(inGeom, 1)
        VOR{i} = findVoronoiCell(i, inGeom(:, 1), inGeom(:, 2));
        waitbar(i/size(inGeom, 1), h, ['Calculating Voronoi diagram....(', sprintf('%d', round(100*i/size(inGeom, 1))), '%)']);
    end

    close(h);


function outPoints = findVoronoiCell(inPos, inGeom1, inGeom2)

    H1 = inGeom1;
    H2 = inGeom2;
    
    sel = spdist(H1(inPos), H2(inPos), H1, H2)<=10;
    H1 = H1(sel);
    H2 = H2(sel);
    
    % find position of p in H1 and H2
    pos2 = findPoint(inGeom1(inPos), inGeom2(inPos), H1, H2);

    % project Points into Plane
    voroPoints = projectV(cat(2, H1, H2), [inGeom1(inPos) inGeom2(inPos)]);

    [v, c] = voronoin([voroPoints(:, 1) voroPoints(:, 2)]);


    v1 = v(:, 1);
    v2 = v(:, 2);
    C = c{pos2};
    outPoints(:, 1) = v1(C);
    outPoints(:, 2) = v2(C);

    outPoints = deprojectV(outPoints, -[inGeom1(inPos) inGeom2(inPos)]);

    %% plot
    %figure(1); clf;
    %scatter(inGeom1, inGeom2, 1, 'b');
    %hold on; scatter(inGeom1(inPos), inGeom2(inPos), 'r');

    %figure(2); clf;
    %scatter(H1, H2, 1, 'g');
    %hold on; scatter(H1(pos2), H2(pos2), 'r');

    %figure(3); clf;
    %scatter(voroPoints(:, 1), voroPoints(:, 2), 1, 'c');
    %hold on; scatter(voroPoints(pos2, 1), voroPoints(pos2, 2), 'r');

    %figure(4); clf;
    %v1=v(:, 1);v2=v(:, 2);
    %for i=1:length(c)
    %plot(v1(c{i}), v2(c{i}), 'k');hold on;
    %end
    %scatter(voroPoints(:, 1), voroPoints(:, 2), 1, 'c');
    %hold on; scatter(voroPoints(pos2, 1), voroPoints(pos2, 2), 'r');

    %figure(5); clf
    %plot(outPoints(:, 1), outPoints(:, 2), 'y');
    %hold on; scatter(voroPoints(pos2, 1), voroPoints(pos2, 2), 'r');

    %figure(6); clf
    %plot(outPoints(:, 1), outPoints(:, 2), 'y');
    %hold on; scatter(inGeom1(inPos), inGeom2(inPos), 'r');

    %plot
    %plot(outPoints(:, 1), outPoints(:, 2));
    %hold on;
    %scatter(H1(pos2), H2(pos2));
    

