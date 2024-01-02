
clear;
clf;

timeStart = tic;

% changeable parameters:

tStart = 0; % Start time
tFinal = 2500; % End time
nPoints = 10000; % number of evaluation points
yLimPos = 700; % positive limit of y axis on plot
yLimNeg = -100; % negative limit of y axis on plot

V = 24; % number of Vertibrae should be 24 (V is the two top vertibrae and the head, 2:V-1 are the pre-sacral (normal) vertibrae, 1 is the fused sacrum and coccyx)
InitConds = zeros(2*V,1); % Initial conditions

% Measurements
% Length = mm (Milimeters)
% mass = g (Grams)
% force = N (Newtons)

% SPINE MASSES

M = zeros(24,1); % Mass vector
M(1,1) = 54; % Sacrum and Coccyx mass (g) (fused) (Average) (Lowrence and Latimer)

M(2,1) = 19.6; % Lumbar Masses (g) (Average) (Lowrence and Latimer)
M(3,1) = 19.4;
M(4,1) = 19.0;
M(5,1) = 17.0;
M(6,1) = 14.5;

M(7,1) = 12.2; % Thoracic Masses (g) (Average) (Lowrence and Latimer)
M(8,1) = 11.1;
M(9,1) = 10.4;
M(10,1) = 9.5;
M(11,1) = 8.8;
M(12,1) = 8.2;
M(13,1) = 7.8;
M(14,1) = 7.1;
M(15,1) = 6.8;
M(16,1) = 6.9;
M(17,1) = 7.6;
M(18,1) = 8.3;

M(19,1) = 6.6; % Cervical Masses (g) (without top two which dont have a disc) (Average) (Lowrence and Latimer)
M(20,1) = 5.8;
M(21,1) = 5.4;
M(22,1) = 5.2;
M(23,1) = 4.9;

M(24,1) = 4700 + 8.3 + 7.6; % Mass of Head region (g) (top two Vertibrae and skull, brain ect) (Yoganandan et. al (for the 4700))

% All spring constants are for the spring underneath it
% SPINE STIFFNESS (SPRING CONSTANT)

K = zeros(24,1); % spring constant vector
K(1,1) = 1; % Sacrum and Coccyx Spring Constants (N/mm) (these arent discs and as such need to be calculated spereately i ahve made them have critical damping) (Qiao and Rahmatallla)
K(2:6,1) = 363.4073; %365.2599; % Lumbar Spring Constants (N/mm) (Qiao and Rahmatallla) (ks3/1000 in there model as it is the closest to the Lumbar region)
K(7:18,1) = 0.5*(363.4073 + 206.6199); % Thoracic Spring Constants (N/mm) (Qiao and Rahmatallla) (ks2/1000 in there model as it is the closest to the Thoracic region)
K(19:24,1) = 206.6199; % Cervivcal Spring Constants (N/mm) (Qiao and Rahmatallla) (ks1/1000 in there model as it is the closest to the Cervical region)
K = K.*1;
% SPINE DAMPING

C = zeros(V,1); % spring damping vector
C(1,1) = 2*sqrt(M(1,1)*K(1,1)); %crtitical damping between bottom of spine and vibration (conceptually it like the buttcheeks or seat)
C(2:6,1) = 0.4459; %2.6276; % Lumbar + Sacrum and Coccyx Damping Constants (Ns/mm) (Qiao and Rahmatallla) (cs3/1000 in there model as it is the closest to the Lumbar region)
C(7:18,1) = 0.5*(0.4459 + 0.3966); % Thoracic Damping Constants (Ns/mm) (Qiao and Rahmatallla) (cs2/1000 in there model as it is the closest to the Thoracic region)
C(19:24,1) = 0.3966; % Cervivcal Damping Constants (Ns/mm) (Qiao and Rahmatallla) (cs1/1000 in there model as it is the closest to the Cervical region)
C = C.*1;
% wave speed paramters

crestStart = 25; % guess of when teh crest wave starts
crestEnd = 100; % guess of when the crest wave ends

% bump function

syms t;

x = linspace(0,tFinal,nPoints);

% - "random" noise function ( sine interpolation )
% function must be continuous for the code to work so interpolating random points is necessary

nEqns = 100;

k = randn(nEqns,1);

p = zeros(nEqns,1);

for i = 1:nEqns

    P(i,1) = k(i,1)*sin(0.01*i*t);

end

%h = 0.5*sum(P);

%h = 300*exp(-0.01.*((t-50).^2));

%h = @(t) (100*sin((pi/50)*(t-10))) .* ((t>10) & (t<60)); % single bump sine curve

%h = 1*-cos((68/1000)*pi*t); % nat resonance 68/1000?

h = @(t) (100) .* ((t>10) & (t<60));

% Initialising:

% - Symbolic ODE:

T = sprintfc('Y%d(t)', 1:V);
syms(T{:});

ODE(1) = M(1)*diff(str2sym(T{1}),t,2) == - C(1)*(2*diff(str2sym(T{1}),t) - diff(str2sym(T{2}),t)) - K(1)*(2*str2sym(T{1}) - str2sym(T{2})) + h;

for i = 2:V-1
    
    ODE(i) = M(i)*diff(str2sym(T{i}),t,2) == - K(i)*(2*str2sym(T{i}) - str2sym(T{i-1}) - str2sym(T{i+1})) - C(i)*(2*diff(str2sym(T{i}),t) - diff(str2sym(T{i-1}),t) - diff(str2sym(T{i+1}),t));

end

ODE(V) = M(V)*diff(str2sym(T{V}),t,2) == - K(V)*(2*str2sym(T{V}) - str2sym(T{V-1})) - C(V)*(2*diff(str2sym(T{V}),t) - diff(str2sym(T{V-1}),t));

% - ODE to system of linear ODEs

[ODEVF,S] = odeToVectorField(ODE);
%ODEVF(1,1) = ODEVF(1,1) +h;
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

% - plotting each vertibrae y-displacement

VertDist = [12.6; 12.4; 12.1; 12.2; 14.3; 16.1; 16.7; 17.5; 17.9; 17.6; 18.1; 19.3; 19.8; 21.3; 22.2; 24.2; 25.5; 24.5; 24.7; 25.5; 24.1; 25.3; 106; 1];
VertDist = flip(VertDist,1);

for j = 1:V

    g = @(x) deval(sol,x,2*j-1) + sum(VertDist(1:j));

    gx = g(x);

    gfig = fplot(g,'LineWidth', 2, 'color', 'red' , 'LineStyle',"-");

    gfig.ShowPoles = 'off';

end

% - Plotting bump

%plot(x,subs(h,t,x),'LineWidth', 2, 'color', "green", 'LineStyle', "-");
plot(x,h(x),'LineWidth', 2, 'color', "green", 'LineStyle', "-");

% - Setting axes

axis([0 tFinal yLimNeg yLimPos]);

% - axis titles and shit

ax.FontSize = 20;
title("Spinal Displacement of each vertibrae",'interpreter','latex','fontsize',40);
ylabel("Displacement (mm)",'interpreter','latex','fontsize',30);
xlabel("Time (ms)",'interpreter','latex','fontsize',30);

% timing stuff

timeODEPlot = toc(timeStart);
disp("Solution plotted in " + num2str(timeODEPlot));

% head displacement
figure(2);

% - plotting head displacement

gV = @(x) deval(sol,x,2*V-1);

gVx = g(x);

gVfig = fplot(gV,'LineWidth', 2, 'color', 'red' , 'LineStyle',"-");

gVfig.ShowPoles = 'off';

% - Setting axes

axis([0 tFinal yLimNeg yLimPos]);

% - axis titles and shit

ax.FontSize = 20;
title("23rd vertibrae (head)",'interpreter','latex','fontsize',40);
ylabel("Displacement (mm) ",'interpreter','latex','fontsize',30);
xlabel("Time $(s)$",'interpreter','latex','fontsize',30);

hold off;

% distance bewteen each vertibrae

figure(3)

hold on

gjdif = zeros(V-1,nPoints);

VertDist = [12.6; 12.4; 12.1; 12.2; 14.3; 16.1; 16.7; 17.5; 17.9; 17.6; 18.1; 19.3; 19.8; 21.3; 22.2; 24.2; 25.5; 24.5; 24.7; 25.5; 24.1; 25.3; 106; 1];
VertDist = flip(VertDist,1);

for j = 1:V-1

    gj = @(x) deval(sol,x,(2*j)-1);

    gj1 = @(x) deval(sol,x,(2*j)+1);

    gjdif(j,:) = - (gj1(x) - gj(x) ) + sum(VertDist(1:j));

    gjdiffig = plot(x,gjdif(j,:),'LineWidth', 2, 'color', 'red' , 'LineStyle',"-");

end

% - Plotting bump

%plot(x,subs(h,t,x),'LineWidth', 2, 'color', "green", 'LineStyle', "-");
plot(x,h(x),'LineWidth', 2, 'color', "green", 'LineStyle', "-");

% - Setting axes

axis([0 tFinal yLimNeg yLimPos]);

% - axis titles and shit

ax.FontSize = 20;
title("All 23 vertibrae distance between",'interpreter','latex','fontsize',40);
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



