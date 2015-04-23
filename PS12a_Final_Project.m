% Deven Hurt, Raul Jordan, Ha Le, Maille Radford
% Final Project
close all
clear all

%% Part 1a

C = 1.1;

% Height of person divided by two for crouching
L = 1.8./2; 

% Average width of shoulders is 18 inches. We halve it.
R = 0.23;
A = 2*R*L;
p = 1.2;

% Coefficient of static friction
mu = 0.05;

% m/s^2
g = 9.8;
% kg
m = 70; 

% Loops over Theta to figure out best value using evolver
theta = 60; 

%m (olympic standards)
x = 101.6; 
h = x*sind(theta); 

v1 = sqrt((2.*m.*g.*h - 2.*mu.*m.*g.*cosd(theta).*(h./sind(theta))) ...
    ./(m+p.*A.*C.*(h./sind(theta))))

%% Part 1b

% L is the integral from 0 to any constant we want that we'll solve for of
% 20/sqrt(1+4x^2)... The constant is 20/tan(theta) (whatever we make the
% second height, h2... we just arbitrarily assigned it as 20)
h2 = 20;
fun = @(x) h2./sqrt(1+4.*x.^2);
L = integral(fun, 0, h2./tan(theta));
% L = integral(h2./sqrt(1+4.*x.^2)), 0, (h2./tan(theta)))

v2 = sqrt((2.*m.*g.*h2 - 2.*mu.*m.*g.*cosd(theta).*L)./(m+p.*A.*C.*L))+v1
% Reconsider this equation... considering the way we do conservation of
% potential energy... where would we consider v1 cuz there is kinetic
% energy @ this switch

%% Part 1 c

theta2 = 45; 

num = m.*v2.^2 - 2.*m.*g.*h2 - C.*p.*A.*v2.^2.*(h2./sind(theta2)) ...
    - 2.*mu.*m.*g.*cosd(theta).*(h2./sind(theta))

v3 = sqrt(num./m)



%% Part 2

% kg
m = 70; 
% m
h = 1.8; 
% Using the BSA equation, except divided by 2
A = (sqrt(m.*h/3600))./2; 
% consider what is the best area to use
C = 0.6;
p = 1.2;
D = (p.*C.*A)./2;
% m/s^2, gravity
g = 9.8;   
% time step
dt = 0.01;   

% define initial conditions
x0 = 0;
y0 = 5; 
v3 = v3; 
theta2 = 30;
vx0 = v3.*cosd(theta2);
vy0 = v3.*sind(theta2);

% iterations
N = 1000; 

% define data arrays
x = zeros(1,N+1); 
x(1) = x0;
y = zeros(1,N+1);
y(1) = y0;
vx = zeros(1,N+1); 
vx(1) = vx0;
vy = zeros(1,N+1); 
vy(1) = vy0;

i = 1;
j = 1;

while i < N
    ax = -(D./m).*v3*vx(i);
    vx(i+1) = vx(i) + ax.*dt;
    x(i+1) = x(i) + vx(i).*dt + 0.5.*ax.*dt.^2;
    ay = -g - ((D./m).*v3*vy(i));
    vy(i+1) = vy(i) + ay.*dt;
    y(i+1) = y(i) + vy(i).*dt + 0.5.*ay.*dt.^2;
    
    if y(i+1) < 0 
         i = N;
         j = j + 1;
    else 
         i = i + 1;
         j = j + 1;
    end
end

plot(x(1:N),y(1:N),'r-')

max_x = max(x(1:N))
    
    