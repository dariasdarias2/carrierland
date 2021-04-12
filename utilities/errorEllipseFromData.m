function ellipse = errorEllipseFromData( data, confidenceInterval )
    % When wanting to get the right values for the chi-squared
    % distribution, you can do the following:
%     x = 0:0.001:20;
%     yy= chi2cdf(x,degreeOfFreedom);
%     [~,ind]= min( abs( (confidenceLevel) - yy ) );
%     chisquare_val = sqrt( x(ind) );
    
    covariance = cov( data );
    [eigenvec, eigenval ] = eig(covariance);
    % Get the index of the largest eigenvector
    [largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
    largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);
    % Get the largest eigenvalue
    largest_eigenval = max(max(eigenval));
    % Get the smallest eigenvector and eigenvalue
    if(largest_eigenvec_ind_c == 1)
        smallest_eigenval = max(eigenval(:,2));
        smallest_eigenvec = eigenvec(:,2);
    else
        smallest_eigenval = max(eigenval(:,1));
        smallest_eigenvec = eigenvec(1,:);
    end
    % Calculate the angle between the x-axis and the largest eigenvector
    angle = atan2(largest_eigenvec(2), largest_eigenvec(1));
    if(angle < 0)
        angle = angle + 2*pi;
    end
    %Define a rotation matrix
    R = [ cos(angle) sin(angle); -sin(angle) cos(angle) ];
    theta_grid = linspace(0,2*pi);
    % Confidence interval 
    switch confidenceInterval
        case 68
            chisquare_val = sqrt(2.2960);
        case 90
            chisquare_val = sqrt(4.605);
        case 95
            chisquare_val = sqrt(5.991);
        case 99
            chisquare_val = sqrt(9.210);
        case 9973
            chisquare_val = sqrt(11.8290);
        case 999
            chisquare_val = sqrt(13.816);
        otherwise
            chisquare_val = sqrt(4.605);
    end
    % Build the semiaxis
    a = chisquare_val*sqrt(largest_eigenval);
    b = chisquare_val*sqrt(smallest_eigenval);
    % the ellipse in x and y coordinates
    ellipse_x  = a*cos( theta_grid );
    ellipse_y  = b*sin( theta_grid );
    % rotated ellipse
    ellipse = [ellipse_x;ellipse_y]' * R;
end