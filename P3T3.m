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

%% P3T2
figure
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
hold on;
plot(alpha_lin1,clinterp);
hold off
% CL vs CD Data 

cl = CLVSCD.data(:,1)';
cd = CLVSCD.data(:,2)';
%Linspace the cl values 
cl_lin = linspace(min(cl),max(cl),10000);

%Interp CD values using cl_lin
cdinterp= interp1(cl,cd,cl_lin);
figure
hold on;

%Plot of CL vs CD
plot(cl_lin,cdinterp);



%Find cd at values of alpha Using cl cd and cl that was used for cl vs alpha 
cd_at_alpha = interp1(cl, cd, cl_1, 'spline');

%New alpha to match size of cd
alpha2 = linspace(min(alpha1),max(alpha1),length(cd));

%Find coefficents of variables 
p = polyfit(alpha2,cd,3);


%New alpha for smooth curve
alpha3 = linspace(min(alpha1),max(alpha1), 100);
cd_at_alpha_interp = polyval(p, alpha3);

figure;
hold on;
plot(alpha3,cd_at_alpha_interp)
% Now plot alpha vs CD
plot(alpha1, cd_at_alpha, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Alpha (degrees)');
ylabel('CD');
title('Alpha vs CD');
xlim([-10 10])
grid on;

%% P3T3
figure
plot(alphast*180/pi,c_Di)
hold on
cdo= interp1(alpha3,cd_at_alpha_interp,alphast*180/pi,'linear','extrap');
plot(alphast*180/pi,cdo)
plot(alphast*180/pi,c_Di+cdo)
legend('Induced Drag','Profile Drag','Total Drag','Location','best')
title('CD vs AoA For Cessna 180')
xlabel('Angle of Attack (degrees)')
ylabel('Drag Coefficient Magnitude')