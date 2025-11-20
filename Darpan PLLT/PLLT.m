function [e, c_L, c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)


i  = (1:N)';
th = i*pi/(2*N);          % find the theta spots along the half span
mu = cos(th);

%calculations
a0   = a0_r   + (a0_t   - a0_r)   .* mu;
c    = c_r    + (c_t    - c_r)    .* mu;
aL0  = aero_r + (aero_t - aero_r) .* mu;
alpha= geo_r  + (geo_t  - geo_r)  .* mu;

d = alpha - aL0;

M = zeros(N,N);
for p = 1:N
    thp = th(p);
    a0p = a0(p);
    cp  = c(p);
    for q = 1:N
        n = 2*q - 1;              % only want the odds
        M(p,q) = (4*b/(a0p*cp))*sin(n*thp) + n*sin(n*thp)/sin(thp);
    end
end
% from anderson 5.2 and the slides
A = M\d;

S  = 0.5*b*(c_r + c_t);
AR = b^2/S;

A1  = A(1);
c_L = A1*pi*AR;

nvec = (2*(1:N)-1)';

if N > 1
    r    = A(2:end)/A1;
    nsub = nvec(2:end);
    delta = sum(nsub .* (r.^2));
else
    delta = 0;
end

e    = 1/(1 + delta);
c_Di = c_L^2/(pi*AR)*(1 + delta);
end
