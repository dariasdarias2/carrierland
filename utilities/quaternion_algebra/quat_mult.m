function qOut = quat_mult (p, q, r )
%**************************************************************************
%function qOut = quat_mult(p,q)
%
% Date: 19.09.2017
% DLR - Institute of Communications and Navigation
% Daniel Arias Medina
%
% Calculate the following quaternion multiplication:
% qOut = p \otimes q
%
% If you input 3 quaternions, the solution would be:
% qOut = p \otimes q \otimes r
%
% The quaternion is a column vector which follows the Hamilton convention,
% where the real part is the first component of the vector:
% p = [ pw, px, py, pz  ]' 
% 
% The size of p and q should be [4,1]. This function do not support
% multiple quaternion multiplication, just to keep the function as simple
% as possible.
%
% Input Parameters:
%   p - input 1st quaternion / quaternion matrix (column-wise)
%   q - input 2nd quaternion / quaternion matrix (column-wise)
%   ( r - not mandatory )
%
% Output Parameters:
%   qOut - output quaternion / quaternion matrix (column-wise)
%
% Reference:
%       - Joan Sola 2017: Quaternion Dynamics for Error State Kalman Filter
%
%**************************************************************************
if nargin == 2
    if size(p,1) ~= size(q,1)  || size(p,2) ~= size(q,2)
        error('Error: the input quaternions have a different size')
    end
    
%     if max(size(p)) ~= 4 || min(size(p)) ~= 1
%         error('Error: wrong size of the input quaternion')
%     end
    
    % In case the quaternions are given as row vectors -> turn them into column
    % vectors
    if size(p,1) == 1 && size(p,2) == 4
        p = reshape(p, [4,1]);
        q = reshape(q, [4,1]);
    end
    
    pw = p(1);
    pu = p(2:4);
    qw = q(1);
    qu = q(2:4);
    
    qOut = [ pw*qw - pu'*qu; ...
                    pw*qu + qw*pu + cross( pu, qu )];
                
    
end

if nargin == 3
        if size(p,1) ~= size(q,1)  || size(p,2) ~= size(q,2) ||   size(p,1) ~= size(r,1)  || size(p,2) ~= size(r,2)
        error('Error: the input quaternions have a different size')
    end
    
    if max(size(p)) ~= 4 || min(size(p)) ~= 1
        error('Error: wrong size of the input quaternion')
    end
    
    % In case the quaternions are given as row vectors -> turn them into column
    % vectors
    if size(p,1) == 1 && size(p,2) == 4 
        p = reshape(p, [4,1]);
        q = reshape(q, [4,1]);
        r = reshape(r, [4,1]);
    end
    
    pw = p(1);
    pu = p(2:4);
    qw = q(1);
    qu = q(2:4);    
    qOut1 = [ pw*qw - pu'*qu; ...
                     pw*qu + qw*pu + cross( pu, qu )];
        
    pw = qOut1(1);
    pu = qOut1(2:4);
    qw = r(1);
    qu = r(2:4);
    qOut = [ pw*qw - pu'*qu; ...
                    pw*qu + qw*pu + cross( pu, qu )];    
end

end