function [cols rows] = eyeDva2Pixels(x,y,framePix,frameDva);
%
% [cols rows] = eyeDva2Pixels(x,y,framePix,frameDva);
%
% Transforms X and Y eye position coordinates from dva to pixel row and
% column indices. For plotting eye traces overlaid onto movie frames.
%
% Input arguments
% X and Y must be vectors of the same length specifying eye position in
% degrees of visual angle, relative to center of the screen. FRAMESIZE
% spefices the size of the movie display in pixels, [WIDTH HEIGHT].
%
% Output arguments
% [COLS ROWS] are the points in the image frame (in pixels) corresponding
% to locations in X and Y.
%
% last modified 2013-jul-16
% dbtm

% convert dva to pixel values
ppdX = framePix(1)/frameDva(1);
ppdY = framePix(2)/frameDva(2);

ctrX = framePix(1)/2;
ctrY = framePix(2)/2;

cols = round(x*ppdX + ctrX);
rows = round(-y*ppdY + ctrY);

% censor off-screen eye positions.
outX = find(cols<1 | cols>framePix(1));
cols(outX) = NaN;
outY = find(rows<1 | rows>framePix(2));
rows(outY) = NaN;