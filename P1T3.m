clc
clear
close all

%% Inputs
NacaStrings = ["NACA 0006","NACA 0012","NACA 0018"];
%NacaStrings = ["NACA 0012","NACA 2412","NACA 4412"];
%NacaStrings = ["NACA 0018","NACA 2418"];
numberPannels=100;
c=100;
alphas=linspace(-6,12,100); %degrees

%% Default Figures
%Default figure settings
set(0, 'defaultFigureUnits', 'inches', 'defaultFigurePosition', [1 1 8 5]);
% figures are 8" wide and 5" tall, with the bottom left corner of the figure beginning 1" up, and 1" to the right from the bottom left corner of your screen -- adjust the size of the figure to your liking
set(0,'defaultLineLineWidth',2.5) % sets all line widths of plotted lines
set(0,'DefaultaxesLineWidth', 1.5) % sets line widths of axes
set(0,'DefaultaxesFontSize', 14)
set(0,'DefaultTextFontSize', 14)
set(0,'DefaultaxesFontName', 'Times new Roman')
set(0,'DefaultlegendFontName', 'Times new Roman')
set(0,'defaultAxesXGrid','on')
set(0,'defaultAxesYGrid','on')

%% Get CL, alphaL=0, a0 and Plot

for i=1:length(NacaStrings)

    %Get NACA details,locations
    [m,p,t,Naca]=AirfoilComps(NacaStrings(i),c);
    [XB,YB]=Locations(m,p,t,c,Naca,i,numberPannels,1);
    difference=100;

    %Get CL for Alphas
    for j=1:length(alphas)
        Cl(i,j)= Vortex_Panel(XB,YB,alphas(j));

        %Find alphaL=0
        if abs(Cl(i,j))<=difference
            alphaL0(i)=alphas(j);
            difference=abs(Cl(i,j));
        end
    end

    %Find a0
    changeCl=diff(Cl(i,:));
    liftSlope=changeCl./(alphas(2)-alphas(1));
    a0(i)=mean(liftSlope);

    %alphaL=0 thin airfoil
    fun1= @(theta) ((2*m/p)-(m*(1-cos(theta))/p/p)).*(cos(theta)-1);
    fun2= @(theta) ((m/((1-p)^2))*((2*p)-1+cos(theta))).*(cos(theta)-1);
    intBreak=acos(1-(2*p));
    if intBreak==0
        alphaL0thin(i)=(-1/pi)*(integral(fun2,intBreak,pi))*(180/pi);
    else
        alphaL0thin(i)=(-1/pi)*(integral(fun1,0,intBreak)+integral(fun2,intBreak,pi))*(180/pi);
    end

    %a0 thin airfoil
    a0thin(i)=2*pi*pi/180;
end

%Plot Results
figure
plot(alphas,Cl)
title('Cl vs Angle Of Attack Vortex Panel Method')
xlabel('\alpha (^{o})')
ylabel('Cl')
legend(NacaStrings,'Location','north')
print('CLAOAVortexPannel','-dpng','-r300');
a0
a0thin
alphaL0
alphaL0thin

%% Functions
function [m,p,t,currentnaca] = AirfoilComps(NacaStrings,c)
%Use 
currentnaca = char(NacaStrings);
m = str2double(currentnaca(6)) / c;
p = str2double(currentnaca(7))/ (c/10);
t = str2double(currentnaca(8:end)) / c;
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
