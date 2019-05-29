function [hl,hb]=jdPlotBounded(varargin)
    
    % [hl,hb]=jdPlotBounded(varargin)
    % Note that you can change the properties of the line and the bounds
    % after plotting using the output handles 'hl' and 'hb' respectively
    
    p = inputParser;   % Create an instance of the inputParser class.
    p.addParamValue('x',[],@(x)isvector(x));
    p.addParamValue('y',[],@(x)isvector(x));
    p.addParamValue('eu',[],@(x)isvector(x));
    p.addParamValue('ed',[],@(x)isvector(x));
    p.addParamValue('Color','',@(x)isempty(x)||ischar(x)||(isnumeric(x)&&numel(x)==3));
    p.addParamValue('FaceColor','b',@(x)isempty(x)||ischar(x)||(isnumeric(x)&&numel(x)==3));
    p.addParamValue('LineColor','b',@(x)isempty(x)||ischar(x)||(isnumeric(x)&&numel(x)==3));
    p.addParamValue('FaceAlpha',0.18,@(x)isnumeric(x));
    p.addParamValue('LineWidth',1.5,@(x)isnumeric(x));
    p.addParamValue('LineStyle','-',@(x)ischar(x));
    p.addParamValue('axes',gca,@(x)isa(x,'matlab.graphics.axis.Axes'));
    p.parse(varargin{:});
    %
    x=p.Results.x(:)';
    y=p.Results.y(:)';
    eu=p.Results.eu(:)';
    ed=p.Results.ed(:)';
    assert(~any(diff([numel(x) numel(y) numel(eu) numel(ed)])),'X Y EU ED must have the same dimensionality');
    if ~isempty(p.Results.Color)
        lineCol=p.Results.Color;
        faceCol=lineCol;
    else
        lineCol=p.Results.LineColor;
        faceCol=p.Results.FaceColor;
    end
    %
  %  if sum(abs([ed ud]))>0
    PVX=[x x(end:-1:1)]; % Patch vertices X
    PVY=[y-ed y(end:-1:1)+eu(end:-1:1)]; % Patch vertices Y
    hb=patch(p.Results.axes,PVX(:),PVY(:),faceCol,'FaceAlpha',p.Results.FaceAlpha,'LineStyle','none');
    wasHolding=ishold;
    if ~wasHolding
        hold(p.Results.axes,'on');
    end
    hl=plot(p.Results.axes,x,y,'-','Color',lineCol,'LineWidth',p.Results.LineWidth,'LineStyle',p.Results.LineStyle);
    if ~wasHolding
        hold(p.Results.axes,'off');
    end
end
