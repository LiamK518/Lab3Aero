clc;
close all;
clear;

%clvsalpha
data = load('Data.mat');
CLVSCD = load('CLVSCD.mat');
x = data.data(1:25,1)';
y = data.data(1:25,2)';
xi = linspace(min(x),max(x),10000);
q = interp1(x,y,xi,"linear");
hold on;
plot(xi,q);
hold off
cl = CLVSCD.data(:,1);
cd = CLVSCD.data(:,2);
cl_lin = linspace(min(cl),max(cl),10000);
cdinterp= interp1(cl,cd,cl_lin);
figure(2)
hold on;
plot(cl_lin,cdinterp);

alpha = linspace(-15,12,length(y2));
cd_at_alpha = interp1(x2, y2, y, 'spline');
p = polyfit(alpha,y2,3);
alpha2 = linspace(min(x),max(x), 100);
cd_at_alpha_interp = polyval(p, alpha2);
figure;
hold on;
plot(alpha2,cd_at_alpha_interp)

% Now plot alpha vs CD
plot(x, cd_at_alpha, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Alpha (degrees)');
ylabel('CD');
title('Alpha vs CD');
xlim([-10 10])
grid on;