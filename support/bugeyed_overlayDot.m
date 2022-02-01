function im = bugeyed_overlayDot(im, visField, pos, diam, opacity)

if nargin<1
    im = rand([3600 1800 3]);
    visField = [-180 180 -90 90];
    diam = 120;
    pos = [0 0];
    opacity = 0.5;
    im_ori = im;
end

pos         = mod(pos+180, 360)-180;

azvec         = linspace(visField(1), visField(2), size(im, 1)+1); % vector of azimuthal positions in input image pixel edges
azvec         = azvec(1:end-1)+diff(azvec)/2;                                % vector of azimuthal positions in input image pixel centres
elvec         = linspace(visField(3), visField(4), size(im, 2)+1); % vector of elevation positions in input image pixel edges
elvec         = elvec(1:end-1)+diff(elvec)/2;                                % vector of azimuthal positions in input image pixel centres
[azi, ele]    = ndgrid(azvec, elvec);


% dpp_a       = median(abs(diff(azi)));   % degrees per pixel along azimuth of original images
% dpp_e       = median(abs(diff(ele)));   % degrees per pixel along elevation of original images
% Area        = (dpp_e*ones(length(azi), 1)) * (dpp_a*cosd(ele)); % The area (in degrees^2) covered by each pixel in the image


r = diam/2;
sel = spdist(pos(1), pos(2), azi(:), ele(:))<=r;
for c = 1:3
    thisC = im(:, :, c);
    thisC(sel) = thisC(sel)/10^opacity;
    im(:, :, c) = thisC;
end

% for iEle = 1:length(ele)
%     e = ele(iEle);
%     eDist = min([r abs(pos(2)-e)]);
%     aRad = sqrt(r^2-eDist^2); % 
% %     aWidth = aWidth*cosd(e); % correct for flat monitor distortion %% TODO: There is also azimuthal distortion!
% %     aWidth = aWidth/cosd(e); % correct for equirectangular projection 
% 
%     iAzi = azi>=pos(1)-aRad & azi<=pos(1)+aRad;
%     im(iAzi, iEle, :) = im(iAzi, iEle, :)/10^opacity;
% end

% figure(8); imagesc(im_ori); axis image
% figure(9); imagesc(im); axis image
end
