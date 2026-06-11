function radius = radiusToAchieveOverlap(overlap, spacing)

    overlap = min([1 max([0.001 overlap]) ]);

    a = overlap * pi /2;
    theta = (1.5 * a)^(1/3);
    maxSteps = 9;
    step = 0;
    keepGoing = true;

    while (step < maxSteps) && (keepGoing)
        step = step + 1;
        d = (a + 0.5*sin(2*theta) -theta)/(1-cos(2*theta));
        theta = theta + d;
        if (theta < 1e-5)
            keepGoing = false;
        end
    end

    x = cos(theta);
    y = sin(theta);
    achievedOverlap = 2*(theta-x * y)/pi;

    distance = 2*x;
    radius = 1;
    f = spacing / (2*x);
    distance = distance * f;
    radius = radius * f;

    debug = false;
    if (debug)
        figure(1); clf;
        xx1 = -distance/2 + radius*cosd(0:360);
        yy1 = radius*sind(0:360);
        xx2 = distance/2 + radius*cosd(0:360);
        yy2 = yy1;
    
        plot(xx1, yy1, 'r-');
        hold on;
        plot(xx2, yy2, 'b-');
        plot(-distance/2, 0, 'ro');
        plot(distance/2, 0, 'bo');
        axis 'equal'
    end

end