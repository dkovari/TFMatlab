function h = textborder(x, y, string, text_color, border_color, varargin)
%TEXTBORDER Display text with border.
%   TEXTBORDER(X, Y, STRING)
%   Creates text on the current figure with a one-pixel border around it.
%   The default colors are white text on a black border, which provides
%   high contrast in most situations.
%   
%   TEXTBORDER(X, Y, STRING, TEXT_COLOR, BORDER_COLOR)
%   Optional TEXT_COLOR and BORDER_COLOR specify the colors to be used.
%   
%   Optional properties for the native TEXT function (such as 'FontSize')
%   can be supplied after all the other parameters.
%   Since usually the units of the parent axes are not pixels, resizing it
%   may subtly change the border of the text out of position. Either set
%   the right size for the figure before calling TEXTBORDER, or always
%   redraw the figure after resizing it.
%
% Output:
%   h = handle to text object. The text objects for the outlines are stored
%   in h.UserData.hOutlines
%
%   A listener is added to the String property of h so that anytime you
%   change h.String =... it will automatically change the outlines
%
%   
%   
%   Author: João F. Henriques, April 2010
%   Modified by Daniel T. Kovari, 1/9/2017

	if isempty(string), return; end
	
	if nargin < 5, border_color = 'k'; end  %default: black border
	if nargin < 4, text_color = 'w'; end  %default: white text
	
	%border around the text, composed of 4 text objects
	offsets = [0 -1; -1 0; 0 1; 1 0];
    hOL = gobjects(4,1);
	for k = 1:4
		hOL(k) = text(x, y, string, 'Color',border_color, varargin{:});
		
		%add offset in pixels
		set(hOL(k), 'Units','pixels')
		pos = get(hOL(k), 'Position');
		set(hOL(k), 'Position', [pos(1:2) + offsets(k,:), 0])
		set(hOL(k), 'Units','data')
	end
	
	%the actual text inside the border
	h = text(x, y, string, 'Color',text_color, varargin{:});
	
	%same process as above but with 0 offset; corrects small roundoff
	%errors
	set(h, 'Units','pixels')
	pos = get(h, 'Position');
	set(h, 'Position', [pos(1:2), 0])
	set(h, 'Units','data')
    
    %Add outlines to userdata
    ud.hOutlines = hOL;
    set(h,'UserData',ud);
    
    %% String Change Fcn
    function StringChange(~,~)
        for j=1:4
            hOL(j).String = h.String;
        end
    end
    addlistener(h,'String','PostSet',@StringChange);
    
    %% Resize Fcn
    function SizeChange(~,~)
        orig_units = h.Units;
        h.Units = 'pixels';
        h_pos = h.Position;
        for j=1:4
            set(hOL(k), 'Units','pixels')
            set(hOL(k), 'Position', [h_pos(1:2) + offsets(k,:), 0])
            set(hOL(k), 'Units',orig_units)
        end
        h.Units = orig_units;
    end
    %addlistener(h.Parent,'SizeChanged',@SizeChange);
    addlistener(h,'Position','PostSet',@SizeChange);
    %% Delete fcn
    function DelFcn(~,~)
        try
            delete(hOL);
        catch
        end
    end
    addlistener(h,'ObjectBeingDestroyed',@DelFcn);
end

