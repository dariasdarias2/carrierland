function eulerOut = quat_2_euler( q )
%**************************************************************************
% function eulerOut = quat_2_euler( q )
%
% Date: 19.09.2017
% DLR - Institute of Communications and Navigation
% Daniel Arias Medina
%
% Calculate the Euler angles from the input quaternion
% [ phi; theta; psi ] 
% 
% phi = atan2(  2*(q0*q1 + q2*q3), 1-2*(q1^2+q2^2)  ) 
% theta = asin( 2*(q0*q2 - q3*q1)  )
% psi = atan2(  2*(q0*q3 + q1*q2) , 1-2*(q2^2+q3^2)   ) 
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
%   eulerOut - euler angles
%
% Reference:
%       - Wikipedia...  https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
%
%**************************************************************************

% rigthInputSize = (max(size(q)) == 4) && ( min(size(q)) == 1 );
% if rigthInputSize == false
%     error('Error: wrong size of the input quaternion')
% end


q0 = q(1);
q1 = q(2);
q2 = q(3);
q3 = q(4);

phi = atan2(  2*(q0*q1 + q2*q3), 1-2*(q1^2+q2^2)  );
theta = asin( 2*(q0*q2 - q3*q1)  );
psi = atan2(  2*(q0*q3 + q1*q2) , 1-2*(q2^2+q3^2)   ) ;

eulerOut = [phi, theta, psi]';

end