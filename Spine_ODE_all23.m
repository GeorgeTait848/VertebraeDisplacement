
clear;
clf;

% changeable parameters:

tStart = 0; % Start time
tFinal = 50; % End time
nPoints = 10000; % number of evaluation points
yLimPos = 25; % positive limit of y axis on plot
yLimNeg = -10; % negative limit of y axis on plot

V = 23; % number of Vertibrae
InitConds = zeros(2*V,1); % Initial conditions

M = ones(V,1); % Mass vector
C = ones(V,1); % spring constant vector
K = ones(V,1); % damping constant vector

crestStart = 25; % guess of when teh crest wave starts
crestEnd = 100; % guess of when the crest wave ends

M(1,1) = 9;
M(23,1) = 25;

% Initialising:

% - Symbolic ODE:

syms t;

T = sprintfc('Y%d(t)', 1:V);
syms(T{:});

h = @(t) (10*sin((pi/20)*(t-30))).* ((t>30) & (t<50)); % single bump sine curve



ODE(1) = M(1)*diff(str2sym(T{1}),t,2) == - K(1)*(2*diff(str2sym(T{1}),t) - diff(str2sym(T{2}),t)) + h;
%- C(1)*(2*str2sym(T{1}) - str2sym(T{2}))  +

for i = 2:V-1
    
    ODE(i) = M(i)*diff(str2sym(T{i}),t,2) == - C(i)*(2*str2sym(T{i}) - str2sym(T{i-1}) - str2sym(T{i+1})) - K(i)*(2*diff(str2sym(T{i}),t) - diff(str2sym(T{i-1}),t) - diff(str2sym(T{i+1}),t));

end

ODE(V) = M(V)*diff(str2sym(T{V}),t,2) == - C(V)*(2*str2sym(T{V}) - str2sym(T{V-1})) - K(V)*(2*diff(str2sym(T{V}),t) - diff(str2sym(T{V-1}),t));

% - ODE to system of linear ODEs

[ODEVF,S] = odeToVectorField(ODE);

% - System to Matlab functions

Meqn = matlabFunction(ODEVF,'vars', {'t','Y'});

% solve ODE

sol = ode45(Meqn,[0 tFinal],InitConds);

% evaluating and plotting ODE

figure(1);

hold on;

x = linspace(0,tFinal,nPoints);

% - plotting each vertibrae y-displacement

for j = 1:V

    g = @(x) deval(sol,x,2*j) + j;

    gx = g(x);

    gfig = fplot(g,'LineWidth', 2, 'color', 'red' , 'LineStyle',"-");

    gfig.ShowPoles = 'off';

end

% - Plotting bump

plot(x,h(x),'LineWidth', 2, 'color', "green", 'LineStyle', "-");

% - Setting axes

axis([0 tFinal yLimNeg yLimPos]);

hold off;

figure(2);

hold on;

x = linspace(0,tFinal,nPoints);

% - plotting each vertibrae y-displacement

gV = @(x) deval(sol,x,2*V);

gVx = g(x);

gVfig = fplot(gV,'LineWidth', 2, 'color', 'red' , 'LineStyle',"-");

gVfig.ShowPoles = 'off';


% - Setting axes

axis([0 tFinal yLimNeg yLimPos]);

hold off;

% finding the gradient of wave crests

gradWave = zeros(V,1);
gradIndex = zeros(V,1);

for k = 1:V
    
    g = @(x) deval(sol,x,2*k);

    gx = g(x);

    [gradWavetmp, gradIndextmp] = max(gx);
        %(nPoints/tFinal)*crestStart):((nPoints/tFinal)*(crestEnd))));

    gradWave(k) = gradWavetmp;
    gradIndex(k) = gradIndextmp;

end

waveTime1 = gradIndex(1)/(nPoints/tFinal);
waveTime2 = gradIndex(V)/(nPoints/tFinal);

gradOfWave = (gradWave(1) - gradWave(V))/(waveTime2 - waveTime1);



