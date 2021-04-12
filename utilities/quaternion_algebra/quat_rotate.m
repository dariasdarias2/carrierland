function aRotated = quat_rotate ( q, a )
%**************************************************************************
%function aRotated = quat_rotate ( p, a )
%
% Date: 19.09.2017
% DLR - Institute of Communications and Navigation
% Daniel Arias Medina
%
% Rotation of the vector "a" using the quaternion rotation:
%   aRot = q \otimes [0; a] \otimes \quat_inv(q)
%
% 
% The quaternion is a column vector which follows the Hamilton convention,
% where the real part is the first component of the vector:
% q = [ qw, qx, qy, qz  ]' 
% 
% Input Parameters:
%   q - quaternion -> should be [4,1]
%   a - vector to rotate, should be [3,1]
%
% Output Parameters:
%   aRotated - vector [3,1] rotated using the quaternion
%
% Reference:
%       - Joan Sola 2017: Quaternion Dynamics for Error State Kalman Filter
%
%**************************************************************************

if max(size(q)) ~= 4 || min(size(q)) ~= 1
    error('Error: wrong size of the input quaternion')
end

if max(size(a)) ~= 3 || min(size(a)) ~= 1
    error('Error: wrong size of the input vector')
end


% In case the quaternion is given as row vector -> turn it into column
if size(q,1) == 1 && size(q,2) == 4
    q = reshape(q, [4,1]);
end
% In case the quaternion is given as row vector -> turn it into column
if size(a,1) == 1 && size(a,2) == 3
    a = reshape(a, [3,1]);
end

if ~(quat_norm(q) > 1-10^-4 && quat_norm(q) < 1+10^-4)
    warning('Warning: not a proper quaternion rotation')
    q = quat_normalize(q);
end

tmp = quat_mult( q, [0; a], quat_conj(q) );

aRotated = tmp(2:4);

end