function prodMat = quat_productMatrix ( q, side )
%**************************************************************************
%function prodMat = quat_productMatrix ( q, side )
%
% Date: 19.09.2017
% DLR - Institute of Communications and Navigation
% Daniel Arias Medina
%
% Calculate the product matrix for a quaternion multiplication:
%    p \otimes q = [p]_R  q
%    p \otimes q = [q]_L   p
%
% where both [p]_R and [q]_L are [4x4 matrices]
% 
% The quaternion is a column vector which follows the Hamilton convention,
% where the real part is the first component of the vector:
% p = [ pw, px, py, pz  ]' 
% 
% The size of p should be [4,1]. This function do not support
% multiple quaternion input, just to keep the function as simple
% as possible.
%
% Input Parameters:
%   q - input quaternion 
%  side - string, either "right" or "left". 
%
% Output Parameters:
%   prodMat - output product matrix
%
% Reference:
%       - Joan Sola 2017: Quaternion Dynamics for Error State Kalman Filter
%
%**************************************************************************
 
if max(size(q)) ~= 4 || min(size(q)) ~= 1
    error('wrong size of the input quaternion')
end

% In case the quaternion is given as row vector -> turn it into column vectors
if size(q,1) == 1 && size(q,2) == 4
    q = reshape(q, [4,1]);
end

if strcmp( side, 'left' ) == 0   &&   strcmp( side, 'right' ) == 0
    error('wrong size of the input quaternion')
end

qw = q(1);
qu = q(2:4);

if strcmp( side, 'left' )
    prodMat = qw*eye(4) + [ 0,  -qu'; qu, skewMatrix(qu) ];
    
elseif strcmp( side, 'right' )
    prodMat = qw*eye(4) + [ 0,  -qu'; qu, -skewMatrix(qu) ];
end

end