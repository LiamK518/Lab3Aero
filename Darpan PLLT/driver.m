clear; clc; close all;

N_terms = 50;
a0      = 2*pi;
aero    = 0;
alpha   = 5*pi/180;      % 5 degrees and also no twist

AR_list    = [4 6 8 10];
lambda = linspace(0,1,100);
b = 1;

difference = zeros(length(AR_list),length(lambda)); % preallocate vector for difference between the ar and the eliptical wing

for ia = 1:length(AR_list)
    AR = AR_list(ia);
    for k = 1:length(lambda)
        lam = lambda(k);

        cr = 2*b/(AR*(1+lam));
        ct = lam*cr;
        

        [e,~,~] = PLLT(b,a0,a0,ct,cr,aero,aero,alpha,alpha,N_terms);
        difference(ia,k) = 1/e - 1;
    end
end

figure; hold on;
for ia = 1:length(AR_list)
    plot(lambda, difference(ia,:), LineWidth=1.5);   %need to add labels and colors
    % doesn't look right yet but i can't figure out why??
end
grid on;
xlabel("Taper Ratio");
ylabel("Induced Drag Factor")
legend("AR = 4", "AR = 6", "AR = 8", "AR = 10")
hold off;
