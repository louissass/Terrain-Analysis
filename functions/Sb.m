function [ slopebreak ] = Sb(Z, cellsize, near_distance, far_distance, direction )
%UNTITLED Adam Winstral's slopebreak term, Sb, where breaks from the
%prevailing wind are defined as:

% the near field Sx minus the far field Sx

% inputs: 
% Z matrix of DEM Z values. This should be the first output from geotiffread or similar. Note the DEM must have square cells
% cellsize from the second output of geotiffread, xextent=yextent=cellsize
% near_distance pair of numbers, the first should =cellsize, and the second = the outside radius for the near field Sx calculation
% far_distance pair of numbers, the first should =the minimum distance for the far field shelter calculation. So it should be > to near_distance(2). the second number is the maximum range of the far distance calculation
% direction of the prevailing wind in degrees, 0-360% 
 
Sx_near = Sx(Z, cellsize, near_distance(1), near_distance(2), direction);
Sx_far = Sx(Z, cellsize, far_distance(1), far_distance(2), direction);
slopebreak = Sx_near - Sx_far;

end

