function outVOR = makeVoronoi(inGeomName, overwrite, maxDist)
% MAKEVORONOI(inGeomName, overwrite) calculates and saves a Voronoi diagram of the pattern given by the name string or path 'inGeomName', 
% if variable overwrite==1, all previous files will automatically be overwritten; if 2, a dialog will ask the user (default)
% maxDist is the maximum (great-circle) distance, in degree, that will be searched for neighbours; lower maxDist means faster computation

% used by 'makeAllVoronoi', 'showVoronoi'

%% inputs
if nargin<3, maxDist = 180; end
if nargin<2, overwrite = 2; end
assert(islogical(overwrite) || isnumeric(overwrite), 'Input argument overwrite must be logical or numeric');
assert(ismember(overwrite, [0,1,2]), 'Invalid value for variable overwrite: %g', overwrite);

%% paths
[vFile, vFolder] = bugeyed_fileName(inGeomName, 'voronoi');
gFile = bugeyed_fileName(inGeomName, 'points');

%% check whether folders and files exist, and whether to overwrite them

if ~isfolder(vFolder)
    mkdir(vFolder);
else
    if isfile(vFile)
        switch overwrite
            case 2
                button = questdlg(['Voronoi diagram of ', inGeomName, ' already exists. Overwrite?'], 'File already exists', 'Yes', 'No', 'No');
            case 1
                button = 'Yes';
            case 0
                button = 'No';
        end
        if strcmp(button, 'No')
            disp('** Warning: Voronoi diagram was NOT created. Process aborted by user');
            return;
        end
    end
end

disp(['* Calculating Voronoi diagram ', inGeomName, ' ...']);
geom   = dlmread(gFile);
outVOR = voronoiPattern(geom, maxDist);

%and save the result
save(vFile, 'outVOR');
fprintf('\bdone\n');
end

%% sub functions
function vor = voronoiPattern(inGeom, maxDist)
    % VORONOIPATTERN(inGeom) calculates the voronoi cells for all points in the input array inGeom
    assert(size(inGeom, 2)==2, 'inGeom needs to be a 2-column matrix');
    h = waitbar(0, 'Calculating Voronoi diagram...(0%)','Name', 'Progress');
    
    vor = cell(1, size(inGeom, 1));
    for i = 1:size(inGeom, 1)
        vor{i} = findVoronoiCell(i, inGeom, maxDist);
        waitbar(i/size(inGeom, 1), h, ['Calculating Voronoi diagram....(', sprintf('%d', round(100*i/size(inGeom, 1))), '%)']);
    end

    close(h);
end

function outPoints = findVoronoiCell(centreIndex, allPoints, maxDist)
    % FINDVORONOICELL(inPos, inGeom, maxDist) calculates the voronoi cell for a single point in the input array inGeom

    centrePoint = allPoints(centreIndex, :); % the point around which the voronoi cell should be found
    sel = spdist(centrePoint(1), centrePoint(2), allPoints(:, 1), allPoints(:, 2))<=maxDist; % reduce the point set to reduce computation time
    
    if nnz(sel)<3
        warning('Not enough points for point no. %d. (Increase beInfo.maxVoronoiDist!)', centreIndex);
        outPoints = [NaN NaN];
    else
        % find position of p in reduced set
        allInds     = 1:length(sel);
        ind2        = find(allInds(sel)==centreIndex); % position of our target point in the new, reduced point set
        % project points into plane
        voroPoints  = sub_rotateSphereD(allPoints(sel, :), centrePoint);
        [v, c]      = voronoin(voroPoints);
        outPoints   = v(c{ind2}, :);    %#ok<FNDSB> 
        outPoints   = sub_unRotateSphereD(outPoints, -centrePoint);
    end

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
end


function outPoints = sub_rotateSphereD(inPoints, targetView)
    % SUB_ROTATESPHD(inPoints, inView) roates a sphere so that it's coordinate targetView is at [0 0] (useful to project points on a sphere onto a plane)
    % targetView in spherical coordinates (degrees)
    % inPoints and outPoints in spherical coordinates (degrees, [az el])
    
    [X, Y, Z] = sph2cart(deg2rad(inPoints(:, 1)), deg2rad(inPoints(:, 2)), 1);
    RotatedPoints = sub_rotate(sub_rotate([X Y Z], [0 0 -targetView(1)]), [0 targetView(2) 0]);
    [outPoints(:, 1), outPoints(:, 2)] = cart2sph(RotatedPoints(:, 1), RotatedPoints(:, 2), RotatedPoints(:, 3));
    outPoints = rad2deg(outPoints);
end

function outPoints = sub_unRotateSphereD(inPoints, targetView)
    % SUB_UNROTATESPHD(inPoints, inView) undoes the rotation of sub_rotateSphereD
    % targetView in spherical coordinates (degrees)
    % inPoints and outPoints in spherical coordinates (degrees, [az el])
    
    [X, Y, Z] = sph2cart(deg2rad(inPoints(:, 1)), deg2rad(inPoints(:, 2)), 1);
    RotatedPoints = sub_rotate([X Y Z], [0 targetView(2) -targetView(1)]);    
    [outPoints(:, 1), outPoints(:, 2)] = cart2sph(RotatedPoints(:, 1), RotatedPoints(:, 2), RotatedPoints(:, 3));
    outPoints = rad2deg(outPoints);
end

function P = sub_rotate(P, R)
    % SUB_ROTATE rotates all points in Nx3 matrix P by rotation vector R in degrees, first around x, then y, then z
    % P = sub_rotate(P, R)
    %
    % P - Nx3 matrix of points [x y z]
    % R - 1x3 vector of rotation angles, in degrees [x y z]

    assert(size(P, 2)==3, 'P needs to be a Nx3 array');
    assert(length(R)==3, 'R needs to be a 1x3 or 3x1 vector');
    P = (sub_rotMat(R)*P')';
end

function R = sub_rotMat(r)
    % ROTMAT creates a rotation matrix R from rotation vector r (in degrees, [x y z])
    r  = deg2rad(r);
    Rx = [1 0 0; 0 cos(r(1)) -sin(r(1)); 0 sin(r(1)) cos(r(1))];
    Ry = [cos(r(2)) 0 sin(r(2)); 0 1 0; -sin(r(2)) 0 cos(r(2))];
    Rz = [cos(r(3)) -sin(r(3)) 0; sin(r(3)) cos(r(3)) 0; 0 0 1];
    R  = Rz*Ry*Rx;
end