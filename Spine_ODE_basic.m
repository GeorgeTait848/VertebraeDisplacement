% basic test for Spine ODE

%ODE = @(t,y) [m*y(1) + C*(2*y(1)-y(2)) - K*(y())]

clear;
clf;

% symbolic stuff for testing

M = 10;
%spring constants
C1 = 100;
C2 = 100;
%damping constants
K1 = 100;
K2 = 100;

hold on;

syms t y1(t) y2(t) y3(t);

h = @(t) (5*sin((pi/2)*(t-3))) .* ((t>3) & (t<5));

ODE1 = M*diff(y1,t,2) == - C1*(2*y1-y2) - K1*(2*diff(y1,t)-diff(y2));

ODE2 = M*diff(y2,t,2) == - C2*(2*y2-y1-y3) - K2*(2*diff(y2,t)-diff(y3,t)-diff(y1,t));

ODE3 = M*diff(y3,t,2) == - C2*(2*y3-y2) - K2*(2*diff(y3,t)-diff(y2,t));

[V,S] = odeToVectorField([ODE1, ODE2, ODE3]);
V(2,1) = V(2,1) + h;

Meqn = matlabFunction(V,'vars', {'t','Y'});

cond = [0 0 0 0 0 0]; 
sol = ode45(Meqn,[0 10],cond);

x = linspace(0, 10, 10000);

g2 = @(x) deval(sol,x,2) + 1;
g4 = @(x) deval(sol,x,4) + 2;
g6 = @(x) deval(sol,x,6) + 3;

g2x = g2(x);
g4x = g4(x);
g6x = g6(x);

fplot(g2,'LineWidth', 2, 'color', "red", 'LineStyle',"-");
fplot(g4,'LineWidth', 2, 'color', "blue", 'LineStyle',"-");
fplot(g6,'LineWidth', 2, 'color', "black", 'LineStyle',"-");

plot(x,h(x),'LineWidth', 2, 'color', "green", 'LineStyle', "-");


axis([0 10 -10 10]);