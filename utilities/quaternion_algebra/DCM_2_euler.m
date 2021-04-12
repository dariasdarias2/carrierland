function euler_out = DCM_2_euler(dcm_in)
%**************************************************************************
%function euler_out = DCM2Euler(dcm_in)
%
% Date: 06.06.2017
% DLR Neustrelitz
% Daniel Arias Medina / Michailas Romanovas
%
% Calculates the Euler angles from given array (3D) of DCM matrices ->
% Adapted from DCM_2_Euler from Michailas (possible error in MyTheta??)
%
% Input Parameters:
%   dcm_in - input Direction Cosine Matrix / Array of Direction cosine
%   matrices (3 x 3 x M)
%
% Output Parameters:
%   euler_out - output Euler angles / matrix of Euler angles (3 x M)
%
% Remarks:
%   The Euler angles are formed columnwise as [\phi \theta \psi] (similar
%   to quaternions) and are in radians. The DCM matrix for M rotations has
%   the form [3 x 3 x M]
%
%   The mathematical expression is based on Groves book, p. 28 (Eq. 2.17)
%
%**************************************************************************
MyPhi = zeros(1,size(dcm_in,3));
MyTheta = zeros(1,size(dcm_in,3));
MyPsi = zeros(1,size(dcm_in,3));

for iCounter = 1:size(dcm_in,3)                                                          % Counting DCM matrices
    MyPhi(1, iCounter) = atan2(dcm_in(3,2,iCounter),dcm_in(3,3,iCounter));               % Phi Euler Angle
    MyTheta(1, iCounter) = -asin(dcm_in(3,1,iCounter));                                  % Theta Euler Angle
    MyPsi(1, iCounter) = atan2(dcm_in(2,1,iCounter),dcm_in(1,1,iCounter));               % Psi Euler Angle
end
euler_out = [MyPhi; MyTheta; MyPsi];


