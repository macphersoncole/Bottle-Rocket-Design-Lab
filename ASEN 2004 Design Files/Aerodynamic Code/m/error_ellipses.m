function [] = error_ellipses(data,ta)
 
    %% Replace this section of code with your real data
    % Simulate and plot 100 data points, (YOU SHOULD USE REAL DATA HERE!)
    x = data(:,1);
    y = data(:,2);
    
    %%

    % Calculate covariance matrix
    P = cov(x,y);
    mean_x = mean(x);
    mean_y = mean(y);
    figure; 
    plot(mean_x,mean_y,'y.','markersize',24)
    axis equal; 
    grid on; 
    xlabel('x, Downrange Distance [m]'); 
    ylabel('y, Crossrange Distance [m]'); 
    title('Landing Points with Error Ellipses');
    hold on;
    plot(x,y,'k.','markersize',6)
    if ta == 1
        plot(40.5,-8.61,'m.','markersize',24)
    end

    % Calculate the define the error ellipses
    n=100; % Number of points around ellipse
    p=0:pi/n:2*pi; % angles around a circle

    [eigvec,eigval] = eig(P); % Compute eigen-stuff
    xy_vect = [cos(p'),sin(p')] * sqrt(eigval) * eigvec'; % Transformation
    x_vect = xy_vect(:,1);
    y_vect = xy_vect(:,2);

    % Plot the error ellipses overlaid on the same figure
    plot(1*x_vect+mean_x, 1*y_vect+mean_y, 'b', 'Linewidth', 1.5)
    plot(2*x_vect+mean_x, 2*y_vect+mean_y, 'g', 'Linewidth', 1.5)
    plot(3*x_vect+mean_x, 3*y_vect+mean_y, 'r', 'Linewidth', 1.5)
    if ta == 1
        legend('Center','Landing Spots','Actual Landing Spot','1\sigma','2\sigma','3\sigma','location','southeast');
    else
        legend('Center','Landing Spots','1\sigma','2\sigma','3\sigma');
        semimajor = max(3*x_vect);
        semiminor = max(3*y_vect);
        area = pi*semimajor*semiminor;
        fprintf(['The semimajor axis is ' num2str(semimajor) ' [m]. The semiminor axis is ' num2str(semiminor) ' [m]. The area is ' num2str(area) ' [m^2]. \n']);
    end

end
