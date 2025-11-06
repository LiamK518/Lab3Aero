clc
clear
close all

NACA=1327;
m=1
p=3
t=27
c=100


if i<p*c
    yc(i)=m*i*(2*p-(i/c))/p/p;
else
    yc(i)=m*(c-x)*(1+(x/c)-(2*p))/(1-p)/(1-p);
end