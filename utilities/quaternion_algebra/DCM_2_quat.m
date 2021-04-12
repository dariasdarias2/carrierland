function qOut = DCM_2_quat( R )
%**************************************************************************
% function qOut = DCM_2_quat( q )
%
% Date: 19.09.2017
% DLR - Institute of Communications and Navigation
% Daniel Arias Medina
%
% Calculate the quaternion rotation from the input rotation matrix,
% given that the definition of the quaternion is as follows:
%   q = [  cos(phi/2);  u sin(phi/2) ]
%   
%   phi = arccos(  (trace(R)-1)/2  )
%      u = ( (R - R') )^{v} / (2*sin(phi))
%
% where the operation ( )^{v} is the inverse of [ ]_x   

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
 
if size(R,1) ~= 3 || size(R,2) ~= 3
    error('Error: wrong size of the input rotation matrix')
end

phi = acos( 0.5*(trace(R) -1)   );
u = [ -R(2,3), R(1,3), -R(1,2) ]';

qOut = zeros(4,1);
qOut(1) = cos(phi/2);
qOut(2:4) = u*sin(phi/2);

qOut = quat_normalize(qOut);



%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
[v,d] = eig(R-eye(3));

% The following code assumes the eigenvalues returned are not necessarily
% sorted by size. This may be overcautious on my part.
d = diag(abs(d));   % Extract eigenvalues
[s, ind] = sort(d); % Find index of smallest one
if d(ind(1)) > 0.001   % Hopefully it is close to 0
    warning('Rotation matrix is dubious');
end

axis = v(:,ind(1)); % Extract appropriate eigenvector

if abs(norm(axis) - 1) > .0001     % Debug
    warning('non unit rotation axis');
end

% Now determine the rotation angle
twocostheta = trace(R)-1;
twosinthetav = [R(3,2)-R(2,3), R(1,3)-R(3,1), R(2,1)-R(1,2)]';
twosintheta = axis'*twosinthetav;

theta = atan2(twosintheta, twocostheta);

q = [cos(theta/2); axis*sin(theta/2)];
qOut = q;

end