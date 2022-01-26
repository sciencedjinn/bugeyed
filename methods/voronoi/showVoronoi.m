function outPh = showVoronoi(inGeomName, plotType, rgb, h, inPh, forceRecalc, maxDist)
% SHOWVORONOI (inGeomName) creates planar/spherical plot of the Voronoi diagram of 
% the input pattern in folder 'inGeomName'

% uses 'makeTriang'

if nargin<7, maxDist = 180; end
if nargin<6, forceRecalc = false; end
if nargin<2, plotType = 'plane'; end

[vFile, vPath, shortName] = bugeyed_fileName(inGeomName, 'voronoi');
pFile = bugeyed_fileName(inGeomName, 'patches');

if nargin<6, forceRecalc = false; end

if forceRecalc
    makeVoronoi(inGeomName, 1, maxDist);
    makePatches(inGeomName);
else
    if ~exist(vPath, 'file') || ~exist(vFile, 'file')
        button = questdlg(['Voronoi diagram of ', shortName, ' does not exist. Do you want to create it now?'], 'No Voronoi diagram found', 'Yes', 'No', 'Yes');
        if strcmp(button, 'No')
            disp('** Warning: Voronoi diagram does not exist. Plotting cancelled');
            return;
        else
            makeVoronoi(inGeomName, 0, maxDist);
        end
    end
    
    
    if ~exist(pFile, 'file')
        makePatches(inGeomName);
    end
end

load(pFile, 'Fx', 'Fy', 'Id')

%% plot
if nargin < 4
    figure('Name', 'Voronoi diagram', 'Position', [0 0 1024 700]);
    h = axes;
end
       
if nargin>4 && ~isempty(inPh)
    outPh = inPh;
else
    outPh = [];
end
hold(h, 'off');

switch plotType
    case {'p', 'pl', 'plane', 'planar'}
        for i = 1:length(Fx)
            fx = Fx{i}; fy = Fy{i}; fc = rgb(Id{i}, :);
            % delete segments that cross from one side to the other
            sel = max(fx, [], 1) - min(fx, [], 1) > 30 | max(fy, [], 1) - min(fy, [], 1) > 30;
            
            if ~all(sel)
                fx = fx(:, ~sel); fy = fy(:, ~sel); 
                fc = fc(~sel, :);
                fc = reshape(fc, [size(fc, 1) 1 size(fc, 2)]);
                if nargin>4 && ~isempty(inPh)
                    set(inPh(i), 'CData', fc);
                elseif nargin>2
                    outPh(i) = patch(h, fx, fy, fc, 'edgecolor', 'none');
                else
                    outPh(i) = patch(h, fx, fy, 'b');
                end
                hold(h, 'on');
            end
        end
    case {'sp', 's', 'sphere', 'spherical'}
        hold(h, 'on')
        sphere;
        colormap([1 1 1]);
        for i = 1:length(Fx)
            fx = Fx{i}; fy = Fy{i}; fc = rgb(Id{i}, :);
            fc = reshape(fc, [size(fc, 1) 1 size(fc, 2)]);
            [x, y, z] = sph2cart(pi/180*fx, pi/180*fy, 1.01);
            patch(x, y, z, fc)
        end
        [Angle(1), Angle(2), Angle(3)] = sph2cart(deg2rad(-15), deg2rad(30), 3);
        [x, y, z] = sph2cart(0, 0, 1);
        arrow = quiver3(x, y, z, x/2, y/2, z/2, 'k');
        for i = 1:size(arrow, 1)
            set(arrow(i), 'LineWidth', 4);
        end

        axis equal;
        %scatter3(1.01, 0, 0, 'r');
%         set(h, 'CameraPosition', Angle);
%         set(h, 'CameraViewAngle', 50);
%         set(h, 'Visible', 'off');
        title(['Voronoi diagram of ', shortName]);
end

