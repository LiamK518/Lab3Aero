function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
%{
=============== outputs ===============
e - spanwise efficiency
c_L - coefficient of lift
c_D - coefficient of drag

=============== inputs ===============
b - span (ft)

a0_t - cross sectional lift slope at tips (rads)
a0_r - cross sectional lift slope at root (rads)

c_t - chord  at the tips (ft)
c_r - chord at the roots (ft)

aero_t - zero lift angle of attack at tips (rads)
aero_r - zero lift angle of attack at roots (rads)

geo_t - geometric angle of attack at tips (rads)
geo_r - geometric angle of attack at roots (rads)

N - number of odd terms to include in the series expansion for circulation
%}

%% Calculate theta in terms of y
theta = acos( -(2y)/b );

alpha = (2Gamma)/(a0_t*y0*Vinf*)
end
