function positionEnu = ecef2enu( positionEcef, referencePositionEcef )
%**************************************************************************
%
% Date: 10.08.2015
% DLR 
% Author: Daniel Arias Medina
%
% This functions returns the position/s in ENU frame, given a vector of
% position(s) in the ECEF frame and the reference ECEF position.
% 
% As latitude and longitude values belong to the ellipsoid coordinate
% system, we need to get the values of latitude and longitude using the
% function "ECEF_2_Ellipsoid". 
%  
% The procedure is explained on the page 188 of the ESA Book: "GNSS Data
% Processing, Vol I: Fundamentals and Algorithms" 
%
%
%**************************************************************************

[latitude, longitude, ~] = ECEF_2_Ellipsoid(referencePositionEcef);


Rotation_Matrix_ECEF2ENU = [           -sin(longitude),                    cos(longitude),                   0         ;...
                             -cos(longitude) * sin(latitude),     -sin(longitude) * sin(latitude),  	cos(latitude) ;...
                              cos(longitude) * cos(latitude),      sin(longitude) * cos(latitude),      sin(latitude) ];
                          
      
nPosition = size( positionEcef, 2 );
positionEnu = NaN( 3, nPosition );
for iPosition = 1:nPosition
    positionEnu(:, iPosition) = Rotation_Matrix_ECEF2ENU * ( positionEcef(:,iPosition) - referencePositionEcef );
end