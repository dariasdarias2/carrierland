function qOut = quat_normalize (q )
%**************************************************************************
%function qOut = quat_normalize( q )
%
% Date: 19.09.2017
% DLR - Institute of Communications and Navigation
% Daniel Arias Medina
%
% Normalize the input quaternion
% qOut = q/||q||
%
%
% Input Parameters:
%   q - input quaternion 
%
% Output Parameters:
%   qOut - normalized quaternion 
%
% Reference:
%       - Joan Sola 2017: Quaternion Dynamics for Error State Kalman Filter
%
%**************************************************************************
 
if max(size(q)) ~= 4 || min(size(q)) ~= 1
    error('Error: wrong size of the input quaternion')
end

% In case the quaternion is given as row vector -> turn it into column vectors
if size(q,1) == 1 && size(q,2) == 4
    q = reshape(q, [4,1]);
end

qOut = q/quat_norm(q);

end