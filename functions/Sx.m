function [ shelter ] = Sx(Z, cellsize, min_distance, max_distance, direction )
%Adam Winstral's shelter term, Sx, where shelter from the
%prevailing wind is defined as:

% max of the tan dz/run for some predefined distance along a predefined
% direction
%   
% -I am stacking shifted dems for efficiency rather than looping
%    through dems.
% -I round dx and dy to the nearest whole cellsize, so it isn't
% interpolating at fractions of a cellsize. 
%     
% inputs: 
% Z matrix of DEM Z values. This should be the first output from geotiffread or similar. Note the DEM must have square cells
% cellsize from the second output of geotiffread, xextent=yextent=cellsize
% min_distance should be set to the cellsize for the standard use of Sx, can be set higher for using Sx as a function to calculate Slopebreak
% max_distance = the total range for the Sx calculation
% direction of the prevailing wind in degrees, 0-360

shelter = zeros(size(Z));
for i = 1:3
 
    adjust = (i-2)*0;
dx = sind(direction + adjust);
dy = cosd(direction + adjust);

[r,c] = size(Z);

if min_distance<cellsize 
    error('the min_distance needs to be greater than the cellsize')
end
if max_distance<min_distance 
    error('the max_distance needs to be greater than the min_distance')
end

loop_end = max_distance/cellsize;
loop_start = min_distance/cellsize;

for n = loop_start:loop_end;
    
    addX = round(abs(dx)*n);
    addY = round(abs(dy)*n);
    zerosX = zeros(r,addX);
    zerosY = zeros(addY,c+addX);
    
    if dx<=0
        Zi = [Z,zerosX];
        Zd = [zerosX,Z];
        
    else
        Zi = [zerosX,Z];
        Zd = [Z,zerosX];
    end
    
    if dy<=0
        Zi = [zerosY;Zi];
        Zd = [Zd;zerosY];
    else
        Zi = [Zi;zerosY];
        Zd = [zerosY;Zd];
    end
    
 temp = atan((Zd - Zi)./(n*cellsize));
 
    if dx<=0
        temp = temp(:,1:c);
    else
        temp = temp(:,addX+1:end);
    end
    
    if dy<=0
        temp = temp(addY+1:end,:);
    else
        temp = temp(1:r,:);
    end
 
 shelt(:,:,n) = temp;
    
    
    
end
    
pshelter = max(shelt,[],3);


    
shelter = pshelter + shelter;
end
shelter = shelter./3;
end

