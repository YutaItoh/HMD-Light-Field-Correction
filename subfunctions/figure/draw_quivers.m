%% DRAW_QUIVERS draws 2D error vectors
%
% NOTE: This function turns on "hold" for the current figure pane
%
% DRAW_QUIVERS(src,dist)
% src:  Source 2D points (2xN)
% dist: Distination 2D points (2xN)
% quiverscale: length of quivers, default is 1.0

% Copyright (c) Yuta Itoh 2014

function draw_quivers(src,dist,quiverscale)
hold on;

if nargin<=2
    quiverscale=1.0;
end

meter2pixel=1;
x = src(1,:)*meter2pixel;
y = src(2,:)*meter2pixel;
x2 = dist(1,:);
y2 = dist(2,:);
%x2 = x0 + w/2-0.5;
%y2 = h - (y0 + h/2-0.5);
% optical flow between the 1st and 2nd ST planes
px=x2-x;
py=y2-y;
quiverc(x,y,px,py,quiverscale);
end