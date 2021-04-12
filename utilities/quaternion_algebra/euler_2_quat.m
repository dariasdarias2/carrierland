function qOut = euler_2_quat( euler )
%**************************************************************************
% function qOut = euler_2_quat( euler )
%
% Date: 19.09.2017
% DLR - Institute of Communications and Navigation
% Daniel Arias Medina
%
% Calculate the quaternion from the input Euler angles
% 
% The quaternion is a column vector which follows the Hamilton convention,
% where the real part is the first component of the vector:
% q = [ qw, qx, qy, qz  ]' 
% 
% Input Parameters:
%   euler - euler angles vector
%
% Output Parameters:
%   qOut - quaternion 
%
% Reference:
%       - Wikipedia...  https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
%
%**************************************************************************
 
phi = euler(1);
theta = euler(2);
psi = euler(3);

cosPhi = cos(phi/2);
sinPhi = sin(phi/2);
cosTheta = cos(theta/2);
sinTheta = sin(theta/2);
cosPsi = cos(psi/2);
sinPsi = sin(psi/2);

q0 = cosPhi * cosTheta * cosPsi + sinPhi * sinTheta * sinPsi;
q1 = sinPhi * cosTheta * cosPsi - cosPhi * sinTheta * sinPsi;
q2 = cosPhi * sinTheta * cosPsi + sinPhi * cosTheta * sinPsi;
q3 = cosPhi * cosTheta * sinPsi - sinPhi * sinTheta * cosPsi;

qOut = [ q0, q1, q2, q3]';

end