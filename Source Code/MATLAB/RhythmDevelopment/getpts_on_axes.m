%% getpts_on_axes
% For a target axes, allows the user to specify points on a single asxes 
% with mouse-clicks until a shift-, right- double- or alt-click adds a 
% final point and resumes program execution, or until the user presses 
% "return". "Backspace" or "Delete" removes the most-recently selected 
% point. 
%
% Usage: 
%   [x,y] = getpts_on_axes()
%   [x,y] = getpts_on_axes(hax)
% Where: 
%   hax = a scalar handle to an axes. If omitted, gca is assumed.
%   [x,y] = Nx1 vectors of the x and y coordinates, in points relative to
%       the axes origin, of the selected points. See notes in comments to
%       convert to other units. 
%
% Usage notes: 
%   Although hax should be the handle to an axes, it may be a handle to any
% graphics object which can be unambiguously related to an axes. If hax is 
% a figure, that figure's "CurrentAxes" is targeted; if the figure's 
% "CurrentAxes" is not set, the function throws an error. If h is a 
% descendent of an axes (such as a line object, an image, an annotation, 
% etc.), the ancestor axes is targeted.    

function [x,y] = getpts_on_axes(hax)
% Specify points with mouse. 

% Check the inputs:
if ~nargin
    hax = gca;
else
    % Don't bother trying to process non-handle/non-scalar input:
    if ~ishandle(hax) || ~isscalar(hax)
        error('Input must be a scalar handle to an axes or figure.');
    end
    % Check the type of handle that was passed. Allow figures with a
    % current axes, axes, or children of axes:
    switch get(hax,'type')
        case 'figure',
            hax = get(hax, 'CurrentAxes');
            % We could make this like ginupt, which adds an axes if the
            % figure doesn't have one, but this is cleaner:
            if isempty(hax), error('Specified figure does not contain an axes.'); end;
        case 'axes' % Do nothing; this is what we want
        otherwise
            % User may have passed the handle returned by plot, line,
            % etc--check to see if it has an axes ancestor:
            hax = ancestor(hax,'axes');
            if isempty(hax)
                error('Input must be a scalar handle to an axes or figure.');
            end
    end
end
hfig = handle(ancestor(hax,'figure'));
hax = handle(hax);

% Store the old functions, if any, to restore at end of point selection:
oldButtonMotion = hfig.WindowButtonMotionFcn;
oldButtonUp = hfig.WindowButtonUpFcn;
oldKeyPress = hfig.WindowKeyPressFcn;

% Set the callbacks:
hfig.WindowButtonMotionFcn = {@buttonMotionFcn, hfig, hax};
hfig.WindowButtonUpFcn = {@buttonUpFcn, hfig, hax};
hfig.WindowKeyPressFcn = {@keyPressFcn, hfig};

% Wait until one of the callbacks clears the WindowButtonMotionFunction:
% Note: The documentation for `waitfor` indicates that it throws an error
% if the figure is deleted while it is executing; empirically, that's not
% the case, so no need to wrap in try-catch.
waitfor(hfig,'WindowButtonMotionFcn', '');

if ishandle(hfig)
    pts = getappdata(hfig, 'SelectedPoints');
    % Make the points relative to the axes origin:
    % Note: to leave realitve to the figure, delete the subtraction.
    scrn_size = get(0,'ScreenSize');
    % Figure values are ratios, multiply by screen pixel size to get pixels
    fig_size = hfig.Position.*[scrn_size(3:4) scrn_size(3:4)];
    % Axes values are also ratios, multiply by figure pixels to convert 
    ax_size = hax.Position.*[fig_size(3:4) fig_size(3:4)];
    % Grab axes limits
    ax_XLim = hax.XLim;
    ax_YLim = hax.YLim;
    % Calculate converstion from screen location to axis value
% %     ax_Xspace = ax_XLim(2)/ax_size(3);    
% %     ax_Yspace = ax_YLim(2)/ax_size(4);
    ax_Xspace = (ax_XLim(2)-ax_XLim(1))/ax_size(3);    
    ax_Yspace = (ax_YLim(2)-ax_YLim(1))/ax_size(4);
    % Subtract axis start position from axis figure position to get
    % location in axis coordinates (Y axis start top instead of bottom so
    % subtract value from maximum)
    x = (pts(:,1) - ax_size(1))*ax_Xspace+ax_XLim(1);
    y = (ax_size(4)-(pts(:,2) - ax_size(2)))*ax_Yspace+ax_YLim(1);
    
    % %     x = pts(:,1);
    % %     y = pts(:,2);
    % %     x = (x/hax.Position(3))*diff(hax.XLim);
    % %     y = (y/hax.Position(4))*diff(hax.YLim);

    % To convert to normalized, divide x and y by hax.Position(3) and
    % hax.Position(4), respectively; to convert to data units, multiply the
    % normalized result by diff(hax.XLim) and diff(hax.YLim), respectively.
    % For this version, simply return in points relative to the axes
    % origin.

    % Cleanup: Return the figure to exactly the same state as when we
    % started (or as close as possible):
    rmappdata(hfig,'SelectedPoints');
    oldPointer = getappdata(hfig,'oldPointer');
    if isempty(oldPointer), oldPointer = 'arrow'; end;
    hfig.Pointer = oldPointer;
    rmappdata(hfig,'oldPointer');
    hfig.WindowButtonMotionFcn = oldButtonMotion;
    hfig.WindowButtonUpFcn = oldButtonUp;
    hfig.WindowKeyPressFcn = oldKeyPress;
else
    x = [];
    y = [];
end

function buttonMotionFcn(~, ~, hfig, hax)
% Change the cursor to "crosshair" when it is over the target axes, and
% back to the original cursor when it is not. 

% Store the current units properties to restore at the end of this callback:
axOldUnits = hax.Units;
figOldUnits = hfig.Units;

% Change Units to "points" for easy calculation:
hax.Units = 'points';
hfig.Units = 'points';

% Position convert position to [left, bottom, right, top]:
xLimits = cumsum([hax.Position(1), hax.Position(3)]);
yLimits = cumsum([hax.Position(2), hax.Position(4)]);

pt = hfig.CurrentPoint;

% If the point is on the current axes, update the pointer as required:
if xLimits(1) <= pt(1) && pt(1) <= xLimits(2) && ...
        yLimits(1) <= pt(2) && pt(2) <= yLimits(2)
    if ~strcmp(hfig.Pointer, 'crosshair')
        setappdata(hfig,'oldPointer',hfig.Pointer);
        hfig.Pointer = 'crosshair';
    end
else
    % If the point is not on the current axes, set the pointer to the
    % original value (or "arrow" if the original value is unavailable):
    oldPointer = getappdata(hfig,'oldPointer');
    if isempty(oldPointer), oldPointer = 'arrow'; end;
    hfig.Pointer = oldPointer;
end

% Return the units to the prior values:
hax.Units = axOldUnits;
hfig.Units = figOldUnits;

function buttonUpFcn(~, ~, hfig, hax)
% Add selected point to the list of points in the figure appdata. If
% SelectionType is not "normal", resume program execution. 

% Store the current units properties to restore at the end of this callback:
axOldUnits = hax.Units;
figOldUnits = hfig.Units;

% Change Units to "points" for easy calculation:
hax.Units = 'points';
hfig.Units = 'points';

% Position convert position to [left, bottom, right, top]:
xLimits = cumsum([hax.Position(1), hax.Position(3)]);
yLimits = cumsum([hax.Position(2), hax.Position(4)]);

pt = hfig.CurrentPoint;

% Only process the selection if it was on the target axes:
if xLimits(1) <= pt(1) && pt(1) <= xLimits(2) && ...
        yLimits(1) <= pt(2) && pt(2) <= yLimits(2)
    % Add the new point to the figure appdata:
    pts = getappdata(hfig, 'SelectedPoints');
    if isempty(pts), pts = pt;
    else pts = [pts; pt];
    end
    setappdata(hfig, 'SelectedPoints', pts);
    % If the user double-clicked, alt-clicked, right-clicked, etc., end
    % collecting points by clearing the WindowButtonMotionFunction
    % (triggering the end of `waitfor`):
    if ~strcmp(hfig.SelectionType,'normal')
        hfig.WindowButtonMotionFcn = '';
    end
end

% Restore the original values of units:
hax.Units = axOldUnits;
hfig.Units = figOldUnits;


function keyPressFcn(~, event, hfig)
% Resume program execution if the user presses "return".
% Remove the last point (if any) if the user presses "backspace" or
% "delete". Note: this does not limit the user to removing a single point,
% but it would be trivial to do so.

switch event.Key
    case 'return'
        hfig.WindowButtonMotionFcn = '';
    case {'backspace','delete'}
        pts = getappdata(hfig, 'SelectedPoints');
        if ~isempty(pts)
            pts(end,:) = [];
            setappdata(hfig,'SelectedPoints',pts);
        end
end