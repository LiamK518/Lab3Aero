
clc 
clear all
close

alpha = linspace(-4,4,2000);
r = 1;
ifplota = 1;
ifplotb = 0;

%% Task 1 

NacaStringT1 = ["NACA 0018","NACA 2418"];
Alpha = 5;
for j = 1:length(NacaStringT1)
[m,p,t,currentnaca] = AirfoilComps(NacaStringT1(j));
c = 100;
pannum = c;
% m,p,t,c all naca comp
%current naca is naca name 
%r is a variable tracked for figures
%pan num- number of pannels 
%ifplot 1 = on 0 = off
[XB,YB]= Locations(m,p,t,c,currentnaca,r,pannum,ifplota);
Cl2(j) = Vortex_Panel(XB,YB,Alpha);
r = r+1;
end

%% Task 2 Part 1
Naca_Strings_T2P1 = ["NACA 0012"];
for e = 1:length(Naca_Strings_T2P1)
[m,p,t,currentnaca2] = AirfoilComps(Naca_Strings_T2P1(e));
c = 100;
[XB2,YB2]= Locations(m,p,t,c,currentnaca2,r,c,ifplotb);
r = r+1;
CLE = Vortex_Panel(XB2,YB2,Alpha); % Exact CL 

%Varried number of pannels
for f = 3:99
[XB2,YB2]= Locations(m,p,t,c,currentnaca2,r,f,ifplotb);
Cl_T2P1(f-2) = Vortex_Panel(XB2,YB2,Alpha);
r = r+1;
end
for d = 1:length(Cl_T2P1)
    if Cl_T2P1(d) >= (.99*CLE)
        numpan = d;
        break
    end
end
hold off;
figure(r)
hold on
b = linspace(3,c,length(Cl_T2P1));
plot(b,Cl_T2P1,'LineWidth',1.5,'Color','blue');
yline(CLE(e),LineWidth=1.3)
xline(numpan)
ylim([CLE(e)-.1 CLE(e)+.1]);
xlim([4,c])
xlabel('Number of Panels');
ylabel('Lift Coefficient');
title(['Lift Coefficient vs Number of Panels for ', currentnaca]);
%Legend
legend('Lift Coefficient for each number of pannels', 'Exact CL', 'Pannel Number of 1% of Exact CL');
grid on;
axis tight
print(currentnaca,'-dpng','-r300'); 
end







%% Task 2 part 2
 a = 1;
colors = lines(length(NacaStrings)); 
for j = 1:length(NacaStrings)
[m,p,t] = AirfoilComps(NacaStrings(j));
c = 100;
[XB2,YB2]= Locations(m,p,t,c,currentnaca,r,c,ifplotb);
for h = 1:length(alpha)
Cl2(a,h) = Vortex_Panel(XB2,YB2,alpha(h));
end
a= a+1;
end
for q =1:length(NacaStrings)
    hold on
figure(30)
plot(alpha, Cl2(q,:), 'color', colors(q,:))
end
xlabel('Angle of Attack (degrees)');
ylabel('Lift Coefficient (Cl)');
title('Lift Coefficient vs Angle of Attack for Various Airfoils');
yline(0);
xline(0);
grid on
legend1 = legend("NACA 0066","NACA 0012","NACA 0018","","");
set(legend1,...
'Position',[0.563654402336621 0.359636302437642 0.117667121418827 0.0799319727891155]);
fontsize(legend1,20,"points")

Cl4 = Cl2(1,:);
Cl5 = Cl2(2,:);
Cl6 = Cl2(3,:);
% Calculate the lift coefficient for the last airfoil and find the index for the desired lift coefficient
for s = 1:length(Cl4)
    if Cl4(s) >= 0
        disp(Cl4(s))
        index2 = s;
        break;
    end
end

zlift = alpha(index2);
xline(zlift,'LineStyle','--')

print("Cl vs Alpha for Various Airfoils",'-dpng','-r300');








%% Functions

function [m,p,t,currentnaca] = AirfoilComps(NacaStrings)
currentnaca = char(NacaStrings);
m = str2double(currentnaca(6)) / 100;
p = str2double(currentnaca(7))/ 10;
t = str2double(currentnaca(8:end)) / 100;
end


function [XB, YB] = Locations(m,p,t,c,currentnaca,r,lin,ifplot)
x = linspace(0,c,lin);
yt = (t / 0.2) * c *(0.2969 *sqrt(x/c)- .1260*(x/c)-.3515*(x/c).^2+.2843*(x/c).^3 -.1036*(x/c).^4); 
for i = 1:length(x) 
    if x(i)<(p*c)
    yc(i)=m*x(i)*(2*p-(x(i)/c))./p^2;
    xzi(i) = atan((m*(2*p - x(i)/c))./p.^2 - (m*x(i))./(c*p.^2));
    else
    yc(i)=m*(c-x(i))*(1+(x(i)/c)-(2*p))./(1-p)^2;
    xzi(i) =atan((m*(c - x(i)))./(c*(p - 1).^2) - (m*(x(i)./c - 2*p + 1))./(p - 1)^2);
    end

end

%postions
xu = x - yt .*sin(xzi);
xl = x + yt .*sin(xzi);
yu = yc + yt.*cos(xzi);
yl = yc - yt.*cos(xzi);
xu2 = xu(2:end);
yu2 = yu(2:end);

xl2 = fliplr(xl);
yl2 = fliplr(yl);
%combine 
% Combine upper and lower surfaces into a single array for vortex panel method
XB = [xl2,xu2]; 
YB = [yl2,yu2];
str = currentnaca;
if ifplot == 1
figure(r)
hold on;
plot(XB,YB,'LineWidth',1.5,'Color','black');
plot(x,yc)
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




function [CL] = Vortex_Panel(XB,YB,ALPHA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:                           %
%                                  %
% XB  = Boundary Points x-location %
% YB  = Boundary Points y-location %
% ALPHA = AOA in degrees           %
%                                  %
% Output:                          %
%                                  %
% CL = Sectional Lift Coefficient  %
% improves efficiency by preallocating matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
% Convert to Radians %
%%%%%%%%%%%%%%%%%%%%%%

ALPHA = ALPHA*pi/180;

%%%%%%%%%%%%%%%%%%%%%
% Compute the Chord %
%%%%%%%%%%%%%%%%%%%%%

CHORD = max(XB)-min(XB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the Number of Panels %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = max(size(XB,1),size(XB,2))-1;
MP1 = M+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preallocate Matrices for Efficiency %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = zeros(1,M);
Y = zeros(1,M);
S = zeros(1,M);
THETA = zeros(1,M);
SINE = zeros(1,M);
COSINE = zeros(1,M);
RHS = zeros(1,M);
CN1 = zeros(M);
CN2 = zeros(M);
CT1 = zeros(M);
CT2 = zeros(M);
AN = zeros(M);
AT = zeros(M);
V = zeros(1,M);
CP = zeros(1,M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intra-Panel Relationships:                                  %
%                                                             %
% Determine the Control Points, Panel Sizes, and Panel Angles %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:M
    IP1 = I+1;
    X(I) = 0.5*(XB(I)+XB(IP1));
    Y(I) = 0.5*(YB(I)+YB(IP1));
    S(I) = sqrt( (XB(IP1)-XB(I))^2 +( YB(IP1)-YB(I))^2 );
    THETA(I) = atan2( YB(IP1)-YB(I), XB(IP1)-XB(I) );
    SINE(I) = sin( THETA(I) );
    COSINE(I) = cos( THETA(I) );
    RHS(I) = sin( THETA(I)-ALPHA );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inter-Panel Relationships:             %
%                                        %
% Determine the Integrals between Panels %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:M
    for J = 1:M
        if I == J
            CN1(I,J) = -1.0;
            CN2(I,J) = 1.0;
            CT1(I,J) = 0.5*pi;
            CT2(I,J) = 0.5*pi;
        else
            A = -(X(I)-XB(J))*COSINE(J) - (Y(I)-YB(J))*SINE(J);
            B = (X(I)-XB(J))^2 + (Y(I)-YB(J))^2;
            C = sin( THETA(I)-THETA(J) );
            D = cos( THETA(I)-THETA(J) );
            E = (X(I)-XB(J))*SINE(J) - (Y(I)-YB(J))*COSINE(J);
            F = log( 1.0 + S(J)*(S(J)+2*A)/B );
            G = atan2( E*S(J), B+A*S(J) );
            P = (X(I)-XB(J)) * sin( THETA(I) - 2*THETA(J) ) ...
              + (Y(I)-YB(J)) * cos( THETA(I) - 2*THETA(J) );
            Q = (X(I)-XB(J)) * cos( THETA(I) - 2*THETA(J) ) ...
              - (Y(I)-YB(J)) * sin( THETA(I) - 2*THETA(J) );
            CN2(I,J) = D + 0.5*Q*F/S(J) - (A*C+D*E)*G/S(J);
            CN1(I,J) = 0.5*D*F + C*G - CN2(I,J);
            CT2(I,J) = C + 0.5*P*F/S(J) + (A*D-C*E)*G/S(J);
            CT1(I,J) = 0.5*C*F - D*G - CT2(I,J);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inter-Panel Relationships:           %
%                                      %
% Determine the Influence Coefficients %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:M
    AN(I,1) = CN1(I,1);
    AN(I,MP1) = CN2(I,M);
    AT(I,1) = CT1(I,1);
    AT(I,MP1) = CT2(I,M);
    for J = 2:M
        AN(I,J) = CN1(I,J) + CN2(I,J-1);
        AT(I,J) = CT1(I,J) + CT2(I,J-1);
    end
end
AN(MP1,1) = 1.0;
AN(MP1,MP1) = 1.0;
for J = 2:M
    AN(MP1,J) = 0.0;
end
RHS(MP1) = 0.0;

%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for the gammas %
%%%%%%%%%%%%%%%%%%%%%%%%

GAMA = AN\RHS';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for Tangential Veloity and Coefficient of Pressure %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for I = 1:M
    V(I) = cos( THETA(I)-ALPHA );
    for J = 1:MP1
        V(I) = V(I) + AT(I,J)*GAMA(J);
    end
    CP(I) = 1.0 - V(I)^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for Sectional Coefficient of Lift %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CIRCULATION = sum(S.*V);
CL = 2*CIRCULATION/CHORD;

end
