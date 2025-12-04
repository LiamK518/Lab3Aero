clc
clear
close all

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

%% P3T1 
addpath(genpath("./Darpan PLLT/"));
%Cessna 180
s=36;
cr=5+(4/12);
ct=3+(7/12);
%via P1T3
a0t=0.118720677642972*180/pi;
aerot=0*pi/180;
a0r=0.118489454792571*180/pi;
aeror=-2.181818181818182*pi/180;
alphast=linspace(-10,15,50)*pi/180;
alphasr=alphast+(2*pi/180);
N=50;

for i=1:length(alphasr)
    [e(i), c_L(i), c_Di(i)] = PLLT(s,a0t,a0r,ct,cr,aerot,aeror,alphast(i),alphasr(i),N);
end
figure
plot(alphast*180/pi,c_L)
title('CL vs AoA For Cessna 180')
xlabel('Angle of Attack (degrees)')
ylabel('CL')
print('CLvAoA','-dpng','-r300');


%% P3T2
%clvsalpha
data = load('Data.mat');
CLVSCD = load('Clvscd2.mat');

%alpha
alpha1 = data.data(1:24,1)';

%cl
cl_1 = data.data(1:24,2)';

%alpha linspaced
alpha_lin1 = linspace(min(alpha1),max(alpha1),10000);

%interp cl values for alpha
clinterp = interp1(alpha1,cl_1,alpha_lin1,"linear");

% CL vs CD Data 

cl = CLVSCD.data(:,1)';
cd = CLVSCD.data(:,2)';
%Linspace the cl values 
cl_lin = linspace(min(cl),max(cl),10000);

%Interp CD values using cl_lin
cdinterp= interp1(cl,cd,cl_lin);

%Plot of CL vs CD




%Find cd at values of alpha Using cl cd and cl that was used for cl vs alpha 
cd_at_alpha = interp1(cl, cd, cl_1, 'spline');

%New alpha to match size of cd
alpha2 = linspace(min(alpha1),max(alpha1),length(cd));

%Find coefficents of variables 
p = polyfit(alpha1,cd_at_alpha,3);


%New alpha for smooth curve
alpha3 = linspace(min(alpha1),max(alpha1), 100);
cd_at_alpha_interp = polyval(p, alpha3);

figure
plot(alpha3,cd_at_alpha_interp)
hold on
xlim([-10 10])
% Now plot alpha vs CD
plot(alpha1, cd_at_alpha, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6)
legend( 'Mathmatical Model','Experimental Data', Location='best')
title('Sectional Drag Coefficent vs Alpha Models')
xlabel('Angle of Attack (degrees)')
ylabel('Sectional Drag Coefficent')
print('SecDvAoAModel','-dpng','-r300');

grid on

%% P3T3
figure
plot(alphast*180/pi,c_Di)
hold on
cdo= interp1(alpha3,cd_at_alpha_interp,alphast*180/pi,'linear','extrap');
plot(alphast*180/pi,cdo)
plot(alphast*180/pi,c_Di+cdo)
legend('Induced Drag','Profile Drag','Total Drag','Location','best')
title('Drag Components vs AoA For Cessna 180')
xlabel('Angle of Attack (degrees)')
ylabel('Drag Coefficient Magnitude')
print('Drag Components','-dpng','-r300');

totalDrag=c_Di+(cdo);

%% P3T4
% w=2500; %lb
% rho=0.001756;%slug/ft3
% area=160.5;%ft^2
% velocities=linspace(0,1000,1000); %ft/s
% ReqCL=w.*2./rho./velocities./velocities./area;
% ReqAoA=interp1(c_L,alphast*180/pi,ReqCL);
% ReqCD=interp1(alphast*180/pi,totalDrag,ReqAoA);
% ReqT=ReqCD.*area.*velocities.*velocities.*rho./2;

w = 2500;              % weight [lb]
rho = 0.001756;        % density [slug/ft^3]
area = 160.5;          % wing area [ft^2]

% Start velocities slightly above 0 to avoid division by zero
velocities = linspace(10, 1000, 1000);  % ft/s (start at 10, not 0)

% Required lift coefficient for level flight
ReqCL = (w * 2) ./ (rho * velocities.^2 * area);

% Required angle of attack
% Make sure c_L and alphast are your actual variable names
ReqAoA = interp1(c_L, alphast*180/pi, ReqCL, 'linear', 'extrap');

% Required drag coefficient  
ReqCD = interp1(alphast*180/pi, totalDrag, ReqAoA, 'linear', 'extrap');

% Required thrust
ReqT = 0.5 * ReqCD * area * rho .* velocities.^2;

% Plot results
figure;
subplot(2,2,1);
plot(velocities, ReqCL);
xlabel('Velocity [ft/s]'); ylabel('Required C_L');
title('Required Lift Coefficient'); grid on;

subplot(2,2,2);
plot(velocities, ReqAoA);
xlabel('Velocity [ft/s]'); ylabel('Required AoA [deg]');
title('Required Angle of Attack'); grid on;

subplot(2,2,3);
plot(velocities, ReqCD);
xlabel('Velocity [ft/s]'); ylabel('Required C_D');
title('Required Drag Coefficient'); grid on;

subplot(2,2,4);
plot(velocities, ReqT);
xlabel('Velocity [ft/s]'); ylabel('Required Thrust [lb]');
title('Thrust Required'); grid on;
figure
plot(velocities/1.68781,ReqT)
title('Thrust Required vs Velocity for Cessna 180')
xlabel('Velocity (knots)')
ylabel('Thrust Required (lbs)')

