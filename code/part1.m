clearvars
clearvars -GLOBAL
close all
set(0,'DefaultFigureWindowStyle', 'docked')

%CONSTANTS AND VARIABLES

global C
global Em T
global BoundX BoundY
global Pxp Px Pyp Py Vx Vy
global Vtherm
global nElectrons
global t_mn



C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665;                      % metres (32.1740 ft) per s²


Em = 0.26 * C.m_0;                    % Mass of the Electron
BoundX = 200e-9;                    % X boundary
BoundY = 100e-9;                    % Y boundary
T = 300;                            % Semiconductor temperature
            
t_mn = 0.2e-12;                     % Mean time between collisions
                    
TimeSteps = 100;                   % Number of time steps

nElectrons = 1000;                   % Number of electrons

dt = 1e-14;                         % Time Step


Efield = 0.5/BoundX;
Force = Efield * C.q_0;
Accel = Force/Em;

Pxp(1: nElectrons) = rand(nElectrons, 1) * BoundX;
Pyp(1: nElectrons) = rand(nElectrons, 1) * BoundY;

Vtherm = sqrt(2 * C.kb * T/Em);

Vx(1: nElectrons) = randn(nElectrons, 1) * (Vtherm/sqrt(2));
Vy(1: nElectrons) = randn(nElectrons, 1) * (Vtherm/sqrt(2));

myColors = ['r' 'b' 'g' 'y' 'm' ];
myColorTyp = 1;

Pscat = 1 - exp(-(dt/t_mn));

a = randi(nElectrons,5,1);

% Current Density
J = zeros(1,TimeSteps - 1);
% I(1) = nElectrons * mean(abs(Vx)) * C.q_0 * BoundX * BoundY;
TAvgp = 300;
for i=2:TimeSteps

   Px(1: nElectrons) = Pxp(1: nElectrons) + (Vx .* dt);
   Py(1: nElectrons) = Pyp(1: nElectrons) + (Vy .* dt);

      %Apply 0.1V at one end
   Vx = Vx + ((1/2) * Accel * dt);
   
   if(Pscat > rand())
     Vx(1: nElectrons) = randn(nElectrons, 1) * (Vtherm/sqrt(2));
     Vy(1: nElectrons) = randn(nElectrons, 1) * (Vtherm/sqrt(2));
   end 

   Vy((Py>BoundY) | (Py<0)) = -Vy((Py>BoundY) | (Py<0));
   
   for j=1:5
       subplot(2,1,1);
       plot([Pxp(a(j)) Px(a(j))], [Pyp(a(j)) Py(a(j))],myColors(j));
       xlim([0 BoundX]);
       ylim([0 BoundY]);
   end
   pause(0.1)
   hold on
   title('2-D plot of particle trajectories');
   
   VxAbs = abs(Vx);
   VyAbs = abs(Vy);
   
   TAvg = (mean((VxAbs.^2)+ (VyAbs.^2)) * Em)/(2 * C.kb);
 
   subplot(2,1,2);
   plot([i-1 i],[TAvgp TAvg],'r');
   xlim([0 TimeSteps]);
   ylim([0 800]);
   pause(0.1)
   hold on
   title('Average Temperature');
   
   Px(Px>BoundX) = Px(Px>BoundX)-BoundX;
   Px(Px<0) = BoundX;
   
   Pxp = Px;
   Pyp = Py;
   TAvgp = TAvg;
   J(i) = nElectrons * mean(abs(Vx)) * C.q_0;
end
figure
plot(linspace(2,TimeSteps,TimeSteps),J);
title('Current Plot');

%Temperature Map

Vmag = sqrt(Vx.^2 + Vy.^2);
Tmap = (Em * (Vmag.^2))./(2 * C.kb);
mapX = linspace(min(Px), max(Px), 100);
mapY = linspace(min(Py), max(Py), 50);
[X,Y] = meshgrid(mapX, mapY);
Tsurf = griddata(Px,Py,Tmap,X,Y);
figure
surf(Tsurf);

% Density Map
Jx = linspace(0, BoundX, 100);
Jy = linspace(0, BoundY, 50);
CurrJ = histcounts2(Py, Px, Jy, Jx);
figure
surf(CurrJ);