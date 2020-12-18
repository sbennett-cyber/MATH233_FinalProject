function [out] = Plot(Start, End, StpInX, u, BC, IC)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

Min=min(u);
Max=max(u);
x=linspace(Start, End, StpInX)';
u1=[BC(1);u(1:StpInX-2);IC];
u12=[BC(2);u(StpInX-1:2*(StpInX-2));IC];
u2=[IC;u(2*(StpInX-2)+1:3*(StpInX-2));BC(3)];
subplot(1,3,1)
plot(x,u1,'r', 'LineWidth',3)
axis([0 End -3 7])
title('Upper channel <0>')
xlabel('X')
ylabel('Population Density')

subplot(1,3,2)
plot(x,u12,'r', 'LineWidth',3)
axis([0 End -3 7])
title('Channel <1>')
xlabel('X')
ylabel('Population Density')

subplot(1,3,3)
plot(x,u2,'b', 'LineWidth',3)
axis([0 End -3 7])
title('Lower channel r')
xlabel('X')
ylabel('Population Density')
drawnow
out=0;
end

