function qOut = quat_conj (q )
%**************************************************************************
%function qOut = quat_conj( q )
%
% Date: 19.09.2017
% DLR - Institute of Communications and Navigation
% Daniel Arias Medina
%
% Calculate the conjugate of the quaternion:
% qOut = q* = [qw; -qu]
%
% The quaternion is a column vector which follows the Hamilton convention,
% where the real part is the first component of the vector:
% p = [ pw, px, py, pz  ]' 
% 
% The size of p should be [4,1]. This function do not support
% multiple quaternion conjugation, just to keep the function as simple
% as possible.
%
% Input Parameters:
%   q - input quaternion 
%
% Output Parameters:
%   qOut - output quaternion 
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

qw = q(1);
qu = q(2:4);

qOut = [ qw; -qu ];

end