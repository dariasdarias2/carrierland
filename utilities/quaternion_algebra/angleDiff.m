function output = angleDiff( eulerA, eulerB )
%**************************************************************************
%function output = angleDiff( eulerA, eulerB )
%
% Date: 13.11.2020
% DLR - Institute of Communications and Navigation
% Daniel Arias Medina
%
% Estimate the smallest angle difference between two angles (or vector of angles)
% output = euler1 - euler2 
%
% Input Parameters:
%   eulerA - input angle(s) [rad] (if vector, it shall be a column vector)
%   eulerB - input angle(s) [rad] (if vector, it shall be a column vector)
%   dim - direction in which the angle difference is estimated
%
% Output Parameters:
%   output - smallest difference [rad]
%
% Reference:
%       - The smallest difference between 2 angles, StackOverflow: 
%          https://stackoverflow.com/questions/1878907/the-smallest-difference-between-2-angles
%
%************************************************************************** 
 
if size(eulerA) ~= size(eulerB)
    error('Error: wrong inputs!')
end

output = nan( size(eulerA) );
for i=1:size(eulerA,1)
    for j=1:size(eulerA,2)
        output(i,j) = atan2( sin( eulerA(i,j)-eulerB(i,j) ), cos(eulerA(i,j)-eulerB(i,j)) );
    end
end


end