%{
% Author: Maharshi Gurjar
% ELEC 4700 - Modeling of Integrated Devices
% Assignment 3
% Part 1
%}
%%
clc; close all; clear;
set(0, 'DefaultFigureWindowStyle', 'docked')
%Define simulation envrionment and constants
M0 = 9.10938356e-31; %Rest mass of electron
Mass_n = 0.26*M0; %Effective mass of electron
Qe = -1.60217662e-19; %Electron Charge
Vx = 0.1; % Voltage along the x-axis
Vy = 0; %Voltage along the y-axis
EConcentration = 1e15 * 100^2; %Electron concentration (1/m^2)
T = 300; % Simulation envrionment temperature (K)
k = 1.38064852e-23; % Boltzmans constant
V_thermal = sqrt(2*k*T/Mass_n); %Thermal Veleocity
Height = 100e-9; % The height of the simulation environment
Length = 200e-9; % The lengthof the simulation environment
nElectrons = 30e3; % Total number of electrons to simulate
nPlotted_Electrons = 10; %Total number of electrons displayed 
Time_Step = Height/V_thermal/100; % Time step of simulation
Iterations = 1000; % Number of iternations to simulate
Show_Movie = 1; %Display steps control
% The mean free path is determined by multipling the thermal velocity 
... by the mean time between collisions: 
MFP = V_thermal * 0.2e-12 %Mean free path 
%The state of the electron  (postion and velocity) is stored in a array
... where each index refers to [x-position y-position v-in-x v-in-y]
Electron_State = zeros(nElectrons,4);
%Setting the top/bottom of the boxes specularity
Top_Specular = 1;
Bottom_Specular = 1;
%Temperature will be recorded in the array below
Temperature = zeros(Iterations,1);
%Create a scattering probability
P_Scatterieng =1 - exp(-Time_Step/0.2e-12);
%Create a distribution using the matlab makedist function
Velocity_PDF = makedist('Normal', 'mu', 0, 'sigma', sqrt(k*T/Mass_n));
%Generate a random inital population postion and velocity
%%
% The electric field compoennts (assuming uniform fields) given by:
Ex = Vx/Length;
Ey = Vy/Height;
%%
% Force on individual force given by:
Fx = Qe*Ex
Fy = Qe*Ey
%%
Temperature = zeros(Iterations,1);
J = zeros(Iterations,2);
%%
for i = 1:nElectrons
   Electron_State(i,:) = [Length*rand() Height*rand() random(Velocity_PDF) random(Velocity_PDF)];
end
%%
figure(1);
subplot(3,1,1)
plot([],[]);
axis ([0 Length/1e-9 0 Height/1e-9]);
title("Trajectories of electron velocities",'interpreter','latex');
xlabel('x (nm)');
ylabel('y (nm)');
grid on;
subplot(3,1,2)
Temp_Plot = animatedline;
title("Semiconducter region temperature",'interpreter','latex');
xlabel('Time (s)');
ylabel('Temperature (K)');
grid on;
axis tight;
subplot(3,1,3)
Current_Plot = animatedline;
title("Jx - Drift Current Density along x-axis",'interpreter','latex');
xlabel('Time (s)');
ylabel('Current Density (A/m)');
grid on;
axis tight;
%We will now move (iterate) over time, updating the positions and direction
...while plotting the state
%%
for i = 1:Iterations
    %%
    %The line below updates Vx and Vy by calculating the change in
    %veloctity under the influence of the electric field
    Electron_State(:,3) = Electron_State(:,3) + Fx*Time_Step/Mass_n;
    Electron_State(:,4) = Electron_State(:,4) + Fy*Time_Step/Mass_n;
    
    %The line below updates the x,y position by moving it to a new position
    ... using its current position + the velocity*(time step)
    Electron_State(:,1:2) = Electron_State(:,1:2) + Time_Step.*Electron_State(:,3:4);
    
    %Checking boundary conditions using Matlab matrix equations 
    
    %Check if and move all electrons at X=200nm Bound:
    Electron_State((Electron_State(:,1)>Length),1) = Electron_State((Electron_State(:,1)>Length),1) - Length;
 
    %Check if and move all electrons at X=0nm Bound:
    Electron_State((Electron_State(:,1)<0),1) =Electron_State((Electron_State(:,1)<0),1) + Length;
    
    %Check (if) and move all electrons at Y Bounds and if specular or diffusive:
    if (Top_Specular == 1)
       Electron_State((Electron_State(:,2)>Height),4) = -1*Electron_State((Electron_State(:,2)>Height),4) ;
       Electron_State((Electron_State(:,2)>Height),2) = 2*Height - Electron_State((Electron_State(:,2)>Height),2);
    else
       %Electron_State((Electron_State(:,2)>Height),2) = Height;
       Electron_State((Electron_State(:,2)>Height),4) = -random(Velocity_PDF);
       Electron_State((Electron_State(:,2)>Height),3) = random(Velocity_PDF);
    end
    if (Bottom_Specular == 1)
        Electron_State((Electron_State(:,2)<0),4) = -1*Electron_State((Electron_State(:,2)<0),4) ;
        Electron_State((Electron_State(:,2)<0),2) = -Electron_State((Electron_State(:,2)<0),2);
    else
       %Electron_State((Electron_State(:,2)<0),2) = 0;
       Electron_State((Electron_State(:,2)<0),4) = random(Velocity_PDF);
       Electron_State((Electron_State(:,2)<0),3) = random(Velocity_PDF);
    end
    %%
    %Add scattering
     j = rand(nElectrons,1) < P_Scatterieng;
    Electron_State(j,3:4) = random(Velocity_PDF,[sum(j),2]);
    
    % Stores the Electron [x y] posistions in the Trajectories vector
    ... for each different electron in a new coloum
    for j = 1: nPlotted_Electrons
       Trajectories_x(i,j) = Electron_State(j,1); 
       Trajectories_y(i,j) = Electron_State(j,2);
    end
    %To calcuatle the themal energy, Maxwell's principle of equipartion 
    ... is used,  where the final equation then becomes;
    Temperature(i) = ( sum (Electron_State(:,3).^2) + sum(Electron_State(:,4).^2)) * Mass_n / k / 2 / nElectrons;
    
    %To calculate the current density
    J(i,1) = Qe.*EConcentration.*mean(Electron_State(:,3));
    J(i,2) = Qe.*EConcentration.*mean(Electron_State(:,4));
    % Add the temperature and current data to the respective plots using
    % addpoints command
    addpoints(Temp_Plot,Time_Step.*i, Temperature(i));
    addpoints(Current_Plot,Time_Step.*i,J(i,1));
    
    %Shows the pathing of the electron, as well as the updating trajectory
    if Show_Movie && mod(i,50)
       figure(1)
       subplot(3,1,1);
       hold off;
       plot(Electron_State(1:nPlotted_Electrons,1)./1e-9,Electron_State(1:nPlotted_Electrons,2)./1e-9,'o');
       grid on;
       axis([0 Length/1e-9 0 Height/1e-9]);
       xlabel('x (nm)');
       ylabel('y (nm)');
       title("Trajectories of electron velocities",'interpreter','latex');
       hold on;
    end
end
%%
figure(1)
subplot(3,1,1)
hold on;
plot(Trajectories_x(:,1:nPlotted_Electrons)./1e-9, Trajectories_y(:,1:nPlotted_Electrons)./1e-9,'.');
hold off;
axis([0 Length/1e-9 0 Height/1e-9]);
xlabel('x (nm)');
ylabel('y (nm)');
grid on;
title("Trajectories of electron velocities",'interpreter','latex');
%saveas(gcf,fullfile('D:\School Work\ELEC 4700\My 4700 Code\Assignment 3\Simulation Results','[Part1]TrajResults.png'),'png')

%%
Density = hist3(Electron_State(:,1:2),[200 100])';
N = 20;
sigma = 3;
%Creating a Gaussian filtering matrix
[x,y]=meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
G=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
G=G./sum(G(:));

figure("name","Electron Density (with Gaussian filtering)")
Density = conv2(Density,G,'same') / (Height./size(Density,1)*Length./size(Density,2));
surf(conv2(Density,G,'same'));
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('Bin count')
axis tight
title("Electron density region",'interpreter','latex');
%saveas(gcf,fullfile('D:\School Work\ELEC 4700\My 4700 Code\Assignment 3\Simulation Results','[Part1]ElectronDensity.png'),'png')
%For the temperature density plot the electrons positions were binned
%Used the same process I did in the first assignment
Temperature_Sum_X = zeros(200,100);
Temperature_Sum_Y = zeros(200,100);
Temperature_Bin_Count = zeros(200,100);
for i=1:nElectrons
    x = floor(Electron_State(i,1)/1e-9);
    y = floor(Electron_State(i,2)/1e-9);
    if(x==0)
        x = 1;
    end
    if(y==0)
        y = 1;
    end
    Temperature_Sum_X(x,y) = Temperature_Sum_X(x,y) + Electron_State(i,3)^2;
    Temperature_Sum_Y(x,y) = Temperature_Sum_Y(x,y) + Electron_State(i,4)^2;
    Temperature_Bin_Count(x,y) = Temperature_Bin_Count(x,y) + 1;
end
% Now, with the velocities added up, calculate the temperatures:
TemperatureDensity = (Temperature_Sum_X + Temperature_Sum_Y).*Mass_n./k./2./Temperature_Bin_Count;
%If somewhere in the calculation the density beame nan then make it 0
TemperatureDensity(isnan(TemperatureDensity)) = 0;
%Transpose the matrix
TemperatureDensity = TemperatureDensity';
[x,y] = meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
G=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
G=G./sum(G(:));
figure("Name","Electron heat denstiy (with Guassian filtering)")
surf(conv2(TemperatureDensity,G,'same'))
set(gca,'YDir','normal');
title('Temperature Map','interpreter','latex');
xlabel('x (nm)');
ylabel('y (nm)');
zlabel('Temperature');
axis tight
%saveas(gcf,fullfile('D:\School Work\ELEC 4700\My 4700 Code\Assignment 3\Simulation Results','[Part1]ElectronHeatDensity.png'),'png')
