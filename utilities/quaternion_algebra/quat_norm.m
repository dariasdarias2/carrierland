function qNorm = quat_norm (q )
%**************************************************************************
%function qOut = quat_norm( q )
%
% Date: 19.09.2017
% DLR - Institute of Communications and Navigation
% Daniel Arias Medina
%
% Calculate the norm of the quaternion:
% qOut = sqrt( qw^2 + qx^2 + qy^2 + qz^2 )
%
% 
% Input Parameters:
%   q - input quaternion 
%
% Output Parameters:
%   qNorm - output norm of the quaternion 
%
% Reference:
%       - Joan Sola 2017: Quaternion Dynamics for Error State Kalman Filter
%
%**************************************************************************
 
if max(size(q)) ~= 4 || min(size(q)) ~= 1
    error('Error: wrong size of the input quaternion')
end


qNorm = sqrt(  q(1)^2 + q(2)^2 + q(3)^2 + q(4)^2  );

end