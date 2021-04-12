function dcm_out = euler_2_DCM(euler_in)
%**************************************************************************
%function euler_out = Euler2DCM(dcm_in)
%
% Date: 06.06.2017
% DLR Neustrelitz
% Daniel Arias Medina
%
% Calculates the Direct Cosine Direction (DCM) matrix from the input Euler
% angles
%
% Input Parameters:
%   euler_in - input Euler Angles (3 x 1 x M)
%
% Output Parameters:
%   dcm_out  - output DCM matrix (3 3 x M)
%
% Remarks:
%   The Euler angles are formed columnwise as [\phi \theta \psi], [pitch,
%   roll, yaw].
%
%   The mathematical expression is based on Groves Book 2008, pag. 28
%
%**************************************************************************
if size(euler_in,2)==1 % One single Euler angle as input
    dcm_out = zeros(3,3);
else
    dcm_out = zeros(3,3,size(euler_in,2));
end

% Groves book
for iCounter = 1:size(euler_in,2)                                                          % Counting DCM matrices
    
    R_yaw = [ cos(euler_in(3, iCounter)),   sin(euler_in(3, iCounter)),     0;  ...
                     -sin(euler_in(3, iCounter)),    cos(euler_in(3, iCounter)),     0;...
                     0,                                             0,                                           1];
    
    R_pitch = [ cos(euler_in(2, iCounter)),         0,      -sin(euler_in(2, iCounter)); ...
                        0,                                                1,      0;...
                       sin(euler_in(2, iCounter)),         0,      cos(euler_in(2, iCounter))];
                   
    R_roll = [ 1,                                                   0,                      0;...
                    0,          cos(euler_in(1, iCounter)),     sin(euler_in(1, iCounter));...
                    0,          -sin(euler_in(1, iCounter)),    cos(euler_in(1, iCounter))];
                
%     R_groves = R_roll * R_pitch * R_yaw;
    R_groves = R_yaw * R_pitch * R_roll;
end

X = euler_in(1);
Y = euler_in(2);
Z = euler_in(3);

cosX = cos(X);
sinX = sin(X);
cosY = cos(Y);
sinY = sin(Y);
cosZ = cos(Z);
sinZ = sin(Z);
        
% Wikipedia...         Link: https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions
% This agrees with Groves (page 28, eq. 2.15 )
R_wiki = [ cosY * cosZ,          -cosX * sinZ + sinX * sinY * cosZ,          sinX * sinZ + cosX * sinY * cosZ; ...
                cosY * sinZ,           cosX * cosZ + sinX * sinY * sinZ,           -sinX * cosZ + cosX * sinY * sinZ; ...
                -sinY,                     sinX * cosY,                                              cosX * cosY];        

% dcm_out = R_groves;                
% dcm_out = R_tu;                
dcm_out = R_wiki;                
        
end        



