function maxfig(handle)

if nargin<1
    handle = gcf;
end
handle.WindowState = 'maximized';

%% Legacy Code:
% warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
% pause(0.01); % makes this work even if the figure has only just been opened 
% jf = get(handle,'JavaFrame');
% set(jf,'Maximized',1);
% % jf.getFigurePanelContainer.getComponent(0).getTopLevelAncestor.setMaximized(1)