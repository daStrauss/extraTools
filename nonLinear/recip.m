% hrm. I've never investigated methods for solving the ode:
% y' = -1/y
% a simple explicit method is easy to come by
% but it appears that an implicit method requires solving a
% quadratic equation.

dt = 0.1;

maxIts = 1000;


y0 = 10;
ye(1) = y0;
yip(1) = y0;
yin(1) = y0;


for ix = 2:maxIts
    ye(ix) = ye(ix-1) - dt/ye(ix-1);
    
    yip(ix) = (yip(ix-1) + sqrt(yip(ix-1)^2 - 4*dt))/2;
    yin(ix) = (yin(ix-1) - sqrt(yin(ix-1)^2 - 4*dt))/2;
    
end


figure(1);
plot((1:1000)*dt,[ye' yip'])
xlabel('time')
ylabel('y(t)')
legend('explicit scheme','implicit scheme')


