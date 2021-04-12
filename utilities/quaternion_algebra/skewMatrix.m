function [ out_matrix ] = skewMatrix( in_vec )
%**************************************************************************
%
% Date: 21.06.2017
% DLR Neustrelitz
% Author: Daniel Arias Medina
%
% Function to build the skew-symmetric matrix from the input vector. Skew
% symmetric matrices come in handy for some rotation operations or general
% algebra. A skew-symmetric operation can also substitute a cross-product:
% a x b = [a]_x b --> where a,b are vectors and []_x represents the
% skew-matrix of the vector.
% 
% Properties of skew-symmetric matrices: 
%  -A = A'
%  if A is skew symmetric -> Ab = -bA'
%
%  a = [x,y,z] -> [a]_x = [ 0    -z  y ;
%                                       z    0   -x ;
%                                       -y   x   0  ]
%
%**************************************************************************
out_matrix = zeros(3);
out_matrix = [0,                -in_vec(3),         in_vec(2); ...
                        in_vec(3),    0,                    -in_vec(1); ...
                        -in_vec(2),   in_vec(1),        0];


end

