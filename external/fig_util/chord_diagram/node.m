classdef node < handle
% NODE Helper class for circularGraph. Not intended for direct user manipulation.
%%
% Copyright 2016 The MathWorks, Inc.
  properties (Access = public)
    Label = '';             % String
    Connection = line(0,0); % Array of lines
    Position;               % [x,y] coordinates
    start_pos;
    end_pos;
    arc_xpos;
    arc_ypos;
    Color = [0 0 0];        % [r g b]
    Visible = true;         % Logical true or false
  end
  
  properties (Access = public, Dependent = true)
    Extent; % Width of text label
  end
  
  properties (Access = private)
    TextLabel;    % Text graphics object
    NodeMarker;   % Line that makes the node visible
    Marker = '-'; % Marker symbol when the node is 'on'
  end
  
  properties (Access = private, Constant)
    labelOffsetFactor = 1.1;
  end
  
  methods
    function this = node(x,y,t,dt)
      % Constructor
      this.Position = [x,y];
      this.start_pos = [cos(t-dt/2) sin(t-dt/2)];
      this.end_pos = [cos(t+dt/2) sin(t+dt/2)];
      ddt = linspace(-dt/2, dt/2, 50);
      this.arc_xpos = cos(t+ddt);
      this.arc_ypos = sin(t+ddt);
      this.Connection = line(0,0);
      makeLine(this);
    end
    
    function makeLine(this)
      % Make the node's line graphic object
      this.NodeMarker = line(...
        this.arc_xpos,...
        this.arc_ypos,...
        repmat(2, size(this.arc_xpos)),...
        'LineWidth', 3, ...
        'Color',this.Color,...
        'LineStyle','-',...
        'PickableParts','all',...
        'ButtonDownFcn',@node.ButtonDownFcn,...
        'UserData',this);
    end
    
    function set.Visible(this,value)
      this.Visible = value;
      updateVisible(this);
    end
    
    function set.Color(this,value)
      this.Color = value;
      updateColor(this);
    end
    
    function set.Label(this,value)
      this.Label = value;
      updateTextLabel(this);
    end
    
    function value = get.Extent(this)
      value = this.TextLabel.Extent(3);
    end
    
    function updateVisible(this)
      if this.Visible
        %this.NodeMarker.Marker = 'o';
        this.NodeMarker.Marker = '-';
        set(this.Connection,'Color',this.Color);
        for i = 1:length(this.Connection)
          this.Connection(i).ZData = ones(size(this.Connection(i).XData));
        end
      else
        this.NodeMarker.Marker = 'x';
        set(this.Connection,'Color',0.9*[1 1 1]);
        for i = 1:length(this.Connection)
          this.Connection(i).ZData = zeros(size(this.Connection(i).XData));
        end
      end
    end
    
    function updateColor(this)
      this.NodeMarker.Color = this.Color;
      set(this.Connection,'Color',this.Color);
    end
    
    function updateTextLabel(this)
      delete(this.TextLabel);
      
      x = this.Position(1);
      y = this.Position(2);
      t = atan2(y,x);
      
      this.TextLabel = text(0,0,this.Label);
      
      this.TextLabel.Position = node.labelOffsetFactor*this.Position;
      if abs(t) > pi/2
        this.TextLabel.Rotation = 180*(t/pi + 1);
        this.TextLabel.HorizontalAlignment = 'right';
      else
        this.TextLabel.Rotation = t*180/pi;
      end
    end
    
    function delete(this)
      % Destructor
      delete(this.Connection(:))
      delete(this.TextLabel);
      delete(this.NodeMarker);
      delete(this);
    end

  end
  
  methods (Static = true)
    function ButtonDownFcn(this,~)
      n = this.UserData;
      if n.Visible
        n.Visible = false;
      else
        n.Visible = true;
      end
    end
  end
end