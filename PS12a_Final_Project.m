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

% Height of person
h = 1.8; 
dt = 0.01;  
% iterations
N = 100; 

x_distance = 101.6; 

% Loops over Theta to figure out best value using evolver
theta = linspace(10,80, 100);
theta2 = linspace(10,80,100);

% Initialize the best x 
bestx = 0;

% Best indeces for theta and theta2 respectively
bestn = 1;
bestm = 1;

for n = 1:length(theta)
    for m = 1:length(theta2)
   
    h = x_distance*sind(theta(n)); 

    v1 = sqrt((2*m*g*h - 2*mu*m*g*cosd(theta(n))*(h/sind(theta(n)))) ...
        /(m+p*A*C*(h/sind(theta(n)))));

    %% Part 1b

    % L is the integral from 0 to any constant we want that we'll solve for of
    % 20/sqrt(1+4x^2)... The constant is 20/tan(theta) (whatever we make the
    % second height, h2... we just arbitrarily assigned it as 20)
    h2 = 20;
    fun = @(q) h2./sqrt(1+4.*q.^2);
    L = integral(fun, 0, h2/tan(theta(n)));

    v2 = sqrt((2*m*g*h2 - 2*mu*m*g*cosd(theta(n))*abs(L))./(m+p*A*C*L)) + v1;
    % Reconsider this equation... considering the way we do conservation of
    % potential energy... where would we consider v1 cuz there is kinetic
    % energy @ this switch

    %% Part 1 c


    num = m*v2^2 - 2*m*g*h2 - C*p*A*v2^2*(h2/sind(theta2(m))) ...
        - 2*mu*m*g*cosd(theta2(m))*(h2./sind(theta2(m)));

    v3 = sqrt(num/m);

    
    % Using the BSA equation, except divided by 2
    A = (sqrt(m*h/3600))/2; 
   
    D = (p*C*A)/2;
   
    % time step

    % define initial conditions
    x0 = 0;
    y0 = 5; 
    vx0 = v3*cosd(theta2(m));
    vy0 = v3*sind(theta2(m));


    % define data arrays
    
    x = zeros(N,N); 
    x(n, m) = x0;
    y = zeros(N,N);
    y(n,m) = y0;
    vx = zeros(N); 
    vx(1) = vx0;
    vy = zeros(N); 
    vy(1) = vy0;

    i = 1;
    j = 1;

    while i < N
        ax = -(D/m)*v3*vx(i);
        vx(i+1) = vx(i) + ax.*dt;
        x(i+1) = x(i,i) + vx(i).*dt + 0.5.*ax.*dt.^2;
        ay = -g - ((D/m)*v3*vy(i));
        vy(i+1) = vy(i) + ay.*dt;
        y(i+1,i+1) = y(i) + vy(i).*dt + 0.5.*ay.*dt.^2;

        if y(i+1, i+1) < 0 
             i = N;
             j = j + 1;
        else 
             i = i + 1;
             j = j + 1;
        end
        i = i + 1;
    end
    
    if (x(n,m) > bestx)
        bestx = x(n, m);
        bestn = n;
        bestm = m;
    end

    end
end


best_theta = theta(bestn); 
best_theta2 = theta2(bestm);