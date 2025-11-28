function [XB, YB] = Locations(m,p,t,c,currentnaca,r,lin,ifplot)
% Cosine spacing along the chord
theta = linspace(0, pi, lin/2);
x = c/2 * (1 - cos(theta));  

% Calculate thickness distribution
yt = (t / 0.2) * c *(0.2969 *sqrt(x/c) - .1260*(x/c) - .3515*(x/c).^2 + .2843*(x/c).^3 - .1036*(x/c).^4); 

% Preallocate camber line and angles
yc = zeros(1,length(x));
xzi = zeros(1,length(x));

% Calculate camber line and surface angles
for i = 1:length(x) 
    if x(i) < (p*c)
        yc(i) = m*x(i)*(2*p-(x(i)/c))./p^2;
        xzi(i) = atan((m*(2*p - x(i)/c))./p.^2 - (m*x(i))./(c*p.^2));
    else
        yc(i) = m*(c-x(i))*(1+(x(i)/c)-(2*p))./(1-p)^2;
        xzi(i) = atan((m*(c - x(i)))./(c*(p - 1).^2) - (m*(x(i)./c - 2*p + 1))./(p - 1)^2);
    end
end

% Calculate upper and lower surface positions
xu = x - yt.*sin(xzi);
xl = x + yt.*sin(xzi);
yu = yc + yt.*cos(xzi);
yl = yc - yt.*cos(xzi);

% Remove first point from upper surface to avoid duplication at trailing edge
xu2 = xu(2:end);
yu2 = yu(2:end);

% Flip lower surface to go from trailing edge to leading edge
xl2 = fliplr(xl);
yl2 = fliplr(yl);

XB = [xl2, xu2]; 
YB = [yl2, yu2];

% Plotting
str = currentnaca;
if ifplot == 1
    figure(r)
    hold on;
    plot(XB,YB,'LineWidth',1.5,'Color','black');
    plot(x,yc,'LineWidth',1.5,'Color','red')
    axis equal;
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    title(currentnaca)
    grid on;
    legend("Airfoil Surface","Camber Line")
    print(str,'-dpng','-r300');
    hold off;
end
end