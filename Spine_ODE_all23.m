
clear;
clf;

timeStart = tic;

% changeable parameters:

tStart = 0; % Start time
tFinal = 500; % End time
nPoints = 10000; % number of evaluation points
yLimPos = 25; % positive limit of y axis on plot
yLimNeg = -10; % negative limit of y axis on plot

V = 24; % number of Vertibrae should be 24 (V is the two top vertibrae and the head, 2:V-1 are the pre-sacral (normal) vertibrae, 1 is the fused sacrum and coccyx)
InitConds = zeros(2*V,1); % Initial conditions

% Measurements
% Length = mm (Milimeters)
% mass = g (Grams)
% force = N (Newtons)

% SPINE MASSES

M = zeros(24,1); % Mass vector
M(1,1) = 10.6*9; % Sacrum and Coccyx mass (g) (fused) (Average) (Lawrence and Latimer)
M(2:6,1) = 17.9; % Lumbar Masses (g) (Average) (Lawrence and Latimer)
M(7:18,1) = 8.7; % Thoracic Masses (g) (Average) (Lawrence and Latimer)
M(19:23,1) = 6.3; % Cervical Masses (g) (without top two which dont have a disc) (Average) (Lawrence and Latimer)
M(24,1) = 470 + (6.3*2); % Mass of Head region (g) (top two Vertibrae and skull, brain ect) (Yoganandan et. al (for the 4700))

% All spring constants are for the spring underneath it
% SPINE STIFFNESS (SPRING CONSTANT)

K = zeros(24,1); % spring constant vector
K(1,1) = 1; % Sacrum and Coccyx Spring Constants (N/mm) (This doesnt matter as it is not used) (Qiao and Rahmatallla)
K(2:6,1) = 36.526; % Lumbar Spring Constants (N/mm) (Qiao and Rahmatallla) (ks3/1000 in there model as it is the closest to the Lumbar region)
K(7:18,1) = 36.341; % Thoracic Spring Constants (N/mm) (Qiao and Rahmatallla) (ks2/1000 in there model as it is the closest to the Thoracic region)
K(19:24,1) = 20.662; % Cervivcal Spring Constants (N/mm) (Qiao and Rahmatallla) (ks1/1000 in there model as it is the closest to the Cervical region)

% SPINE DAMPING

C = zeros(V,1); % spring damping vector
C(1,1) = sqrt(95.4); %crtitical damping between bottom of spine and vibration
C(2:6,1) = 262.720; % Lumbar + Sacrum and Coccyx Damping Constants (g/s) (Qiao and Rahmatallla) (cs3*1000 in there model as it is the closest to the Lumbar region)
C(7:18,1) = 44.590; % Thoracic Damping Constants (g/s) (Qiao and Rahmatallla) (cs2*1000 in there model as it is the closest to the Thoracic region)
C(19:24,1) = 39.660; % Cervivcal Damping Constants (g/s) (Qiao and Rahmatallla) (cs1*1000 in there model as it is the closest to the Cervical region)

% wave speed paramters

crestStart = 25; % guess of when teh crest wave starts
crestEnd = 100; % guess of when the crest wave ends

%M(1,1) = 9;
%M(23,1) = 25;

% Initialising:

% - Symbolic ODE:

syms t;

T = sprintfc('Y%d(t)', 1:V);
syms(T{:});

h = @(t) (10*sin((pi/10)*(t-10))).* ((t>10) & (t<20)); % single bump sine curve



ODE(1) = M(1)*diff(str2sym(T{1}),t,2) == - C(1)*(2*diff(str2sym(T{1}),t) - diff(str2sym(T{2}),t)) - K(1)*(2*str2sym(T{1}) - str2sym(T{2}))  + h;
%
%

for i = 2:V-1
    
    ODE(i) = M(i)*diff(str2sym(T{i}),t,2) == - K(i)*(2*str2sym(T{i}) - str2sym(T{i-1}) - str2sym(T{i+1})) - C(i)*(2*diff(str2sym(T{i}),t) - diff(str2sym(T{i-1}),t) - diff(str2sym(T{i+1}),t));

end

ODE(V) = M(V)*diff(str2sym(T{V}),t,2) == - K(V)*(2*str2sym(T{V}) - str2sym(T{V-1})) - C(V)*(2*diff(str2sym(T{V}),t) - diff(str2sym(T{V-1}),t));

% - ODE to system of linear ODEs

[ODEVF,S] = odeToVectorField(ODE);

% - System to Matlab functions

Meqn = matlabFunction(ODEVF,'vars', {'t','Y'});

% timing stuff

timeODEBuild = toc(timeStart);
disp("ODE Initialised in " + num2str(timeODEBuild));

% solve ODE

sol = ode45(Meqn,[0 tFinal],InitConds);

% evaluating and plotting ODE

figure(1);

hold on;

x = linspace(0,tFinal,nPoints);

% - plotting each vertibrae y-displacement

for j = 1:V

    g = @(x) deval(sol,x,2*j-1) + j;

    gx = g(x);

    gfig = fplot(g,'LineWidth', 2, 'color', 'red' , 'LineStyle',"-");

    gfig.ShowPoles = 'off';

end

% - Plotting bump

plot(x,h(x),'LineWidth', 2, 'color', "green", 'LineStyle', "-");

% - Setting axes

axis([0 tFinal yLimNeg yLimPos]);

% - axis titles and shit

ax.FontSize = 20;
title("Spinal Displacement of each vertibrae",'interpreter','latex','fontsize',40);
ylabel("Displacement",'interpreter','latex','fontsize',30);
xlabel("Time $(t)$",'interpreter','latex','fontsize',30);

% timing stuff

timeODEPlot = toc(timeStart);
disp("Solution plotted in " + num2str(timeODEPlot));

return
% head displacement

figure(2);

x = linspace(0,tFinal,nPoints);

% - plotting each vertibrae y-displacement

gV = @(x) deval(sol,x,2*V);

gVx = g(x);

gVfig = fplot(gV,'LineWidth', 2, 'color', 'red' , 'LineStyle',"-");

gVfig.ShowPoles = 'off';

% - Setting axes

axis([0 tFinal yLimNeg yLimPos]);

% - axis titles and shit

ax.FontSize = 20;
title("23rd vertibrae (head)",'interpreter','latex','fontsize',40);
ylabel("Displacement",'interpreter','latex','fontsize',30);
xlabel("Time $(t)$",'interpreter','latex','fontsize',30);

hold off;

% distance bewteen each vertibrae

figure(3)

hold on

gjdif = zeros(V-1,nPoints);

for j = 1:V-1

    gj = @(x) deval(sol,x,2*j);

    gj1 = @(x) deval(sol,x,(2*j)+2);

    gjdif(j,:) = gj1(x) - gj(x);

    gjdiffig = plot(x,gjdif(j,:),'LineWidth', 2, 'color', 'red' , 'LineStyle',"-");

end

% - Plotting bump

plot(x,h(x),'LineWidth', 2, 'color', "green", 'LineStyle', "-");

% - Setting axes

axis([0 tFinal yLimNeg yLimPos]);

% - axis titles and shit

ax.FontSize = 20;
title("All 23 vertibrae",'interpreter','latex','fontsize',40);
ylabel("Displacement",'interpreter','latex','fontsize',30);
xlabel("Time $(t)$",'interpreter','latex','fontsize',30);

hold off

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


% notes

%driving forces vs head forces
% plot for different vertibrae distances



