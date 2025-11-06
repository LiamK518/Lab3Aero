clc 
clear 
close all
NACA=1327;
m=1;
p=3;
thickness=27;
c=100;

t = thickness / 100;
yt = (t / 0.2) * c *(.2969 *sqrt(x/c)- .1260*(x/c)-.3515*(x/c)^2+.2843*(x/c)^3 -.1036*(x/c)^4);

if i<p*c
    yc(i)=m*i*(2*p-(i/c))/p/p;
else
    yc(i)=m*(c-x)*(1+(x/c)-(2*p))/(1-p)/(1-p);
end

%postions
