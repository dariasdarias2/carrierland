function R = quat_2_DCM( q )
%**************************************************************************
% function R = quat_2_DCM( q )
%
% Date: 19.09.2017
% DLR - Institute of Communications and Navigation
% Daniel Arias Medina
%
% Calculate the rotation matrix from the input quaternion
% R{q} = (qw^2 - qu' qu) I + 2 qu qu' + 2 qw [qu]_x
%
% The quaternion is a column vector which follows the Hamilton convention,
% where the real part is the first component of the vector:
% q = [ qw, qx, qy, qz  ]' 
% 
% The size of p should be [4,1]. This function do not support
% multiple quaternion input, just to keep the function as simple
% as possible.
%
% Input Parameters:
%   q - input quaternion 
%
% Output Parameters:
%   R - rotation matrix, size [3x3]
%
% Reference:
%       - Joan Sola 2017: Quaternion Dynamics for Error State Kalman Filter
%
% NOTE: One should also check that the inverse and the transpose matrices of the
% estimated rotation matrix R are roughly identical. If not, the matrix is
% not a proper rotation matrix
%
%**************************************************************************
 
if max(size(q)) ~= 4 || min(size(q)) ~= 1
    error('Error: wrong size of the input quaternion')
end

% In case the quaternion is given as row vector -> turn it into column vectors
if size(q,1) == 1 && size(q,2) == 4
    q = reshape(q, [4,1]);
end

R = zeros(3);
qw = q(1);
qu = q(2:4);

R = (qw^2 - qu' * qu) * eye(3) + 2*qu*qu' + 2*qw*skewMatrix(qu);

if det(R) == -1 
    warning('Warning: improper rotation matrix')
elseif ~( det(R) < 1 + 10^-3 &&  det(R) > 1 - 10^-3 ) % the determinant is roughly 1
    warning('Warning: not a rotation matrix')
end

end