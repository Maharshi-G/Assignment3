%{
Author : Maharshi Gurjar
Elec 4700 Assignment 3 -  Finite Difference Method
Part 3
%}
%%
clc; close all; clear;
set(0, 'DefaultFigureWindowStyle', 'docked')
%%
%Define simulation envrionment and constants
M0 = 9.10938356e-31; %Rest mass of electron
Mass_n = 0.26*M0; %Effective mass of electron
Qe = -1.60217662e-19; %Electron Charge
Vx = -0.4; % Voltage along the x-axis
Vy = 0; %Voltage along the y-axis
EConcentration = 1e15 * 100^2; %Electron concentration (1/m^2)
T = 300; % Simulation envrionment temperature (K)
k = 1.38064852e-23; % Boltzmans constant
V_thermal = sqrt(2*k*T/Mass_n); %Thermal Veleocity
Height = 100e-9; % The height of the simulation environment
Length = 200e-9; % The lengthof the simulation environment
nElectrons = 30e3; % Total number of electrons to simulate
nPlotted_Electrons = 20; %Total number of electrons displayed
Time_Step = Height/V_thermal/100; % Time step of simulation
Iterations = 500; % Number of iternations to simulate
Show_Movie = 0; %Display steps control
%%
% The mean free path is determined by multipling the thermal velocity
... by the mean time between collisions:
    MFP = V_thermal * 0.2e-12; %Mean free path
%%
%The state of the electron  (postion and velocity) is stored in a array
... where each index refers to [x-position y-position v-in-x v-in-y]
    Electron_State = zeros(nElectrons,4);
%%
%Temperature will be recorded in the array below
Temperature = zeros(Iterations,1);
%%
%Create a scattering probability
P_Scattering =1 - exp(-Time_Step/0.2e-12);
%%
%Create a distribution using the matlab makedist function
Velocity_PDF = makedist('Normal', 'mu', 0, 'sigma', sqrt(k*T/Mass_n));
%%
%Setting the top/bottom of the boudnary specularity
Top_Specular = 1;
Bottom_Specular = 1;
%Set the box characteristics (0 = diffusive)
Top_Box_Specular = 1;
Bottom_Box_Specular = 1;
%Create Box-positions [x1, x2, y1,y2]
%Box_pos = 1e-9.*[80 120 0 32.5; 80 120 67.5 100]; % [Gap size 35nm]
%Box_pos = 1e-9.*[80 120 0 35; 80 120 65 100]; % [Gap size 30nm]
%Box_pos = 1e-9.*[80 120 0 37.5; 80 120 62.5 100]; % [Gap size 25nm]
Box_pos = 1e-9.*[80 120 0 40; 80 120 60 100]; % [Gap size 20nm] [Control ]
%Box_pos = 1e-9.*[80 120 0 42.5; 80 120 57.5 100]; % [Gap size 15nm]
%Box_pos = 1e-9.*[80 120 0 45; 80 120 55 100]; % [Gap size 10nm]
%Box_pos = 1e-9.*[80 120 0 47.5; 80 120 52.5 100]; % [Gap size 5nm]
%%
% The electric field compoennts (assuming uniform fields) given by:
Ex = Vx/Length;
Ey = Vy/Height;
%%
% Force on individual force given by:
Fx = Qe*Ex;
Fy = Qe*Ey;
%%
% Current Density vector
J = zeros(Iterations,2);
%% Current plot and Temperature plot stuff
Temp_Plot = zeros(Iterations,1);
Current_Plot = zeros(Iterations,1);
%%
%Create the state of the box (specular or diffusive)
%Generate a random inital population postion and velocity
for i = 1:nElectrons
    Electron_State(i,:) = [Length*rand() Height*rand() random(Velocity_PDF) random(Velocity_PDF)];
    Electron_State( (Electron_State(:,1)>Box_pos(1,1) & Electron_State(:,1)<Box_pos(1,2) & ...
        (Electron_State(:,2)<Box_pos(1,4))) ,1) = Length*rand();
    Electron_State( (Electron_State(:,1)>Box_pos(2,1) & Electron_State(:,1)<Box_pos(2,2) & ...
        (Electron_State(:,2)>Box_pos(2,3))) ,1) = Length*rand();
end
%%
%Figure below used to test inital postion of electrons, and if movie is off
%shows the the same.
figure("Name","Electron positions")
plot(Electron_State(1:nPlotted_Electrons,1)./1e-9,Electron_State(1:nPlotted_Electrons,2)./1e-9,'o');
grid on;
axis([0 Length/1e-9 0 Height/1e-9]);
xlabel('x (nm)');
ylabel('y (nm)');
title(sprintf("Plotting (%d/%d) electron ",nPlotted_Electrons,nElectrons));
%saveas(gcf,'Part_Three_Boxes.png')
hold on;
for j=1:size(Box_pos,1)
    plot([Box_pos(j, 1) Box_pos(j, 1) Box_pos(j, 2) Box_pos(j, 2) Box_pos(j, 1)]./1e-9,...
        [Box_pos(j, 3) Box_pos(j, 4) Box_pos(j, 4) Box_pos(j, 3) Box_pos(j, 3)]./1e-9, 'k-');
end
hold off;
%%
%We will now move (iterate) over time, updating the positions and direction
...while plotting the state
    for i = 1:Iterations
    %The line below updates the x,y position by moving it to a new position
    ... using its current position + the velocity*(time step)
    Electron_Prev_State = Electron_State(:,1:2);
    Electron_State(:,3) = Electron_State(:,3) + Fx*Time_Step/Mass_n;
    Electron_State(:,4) = Electron_State(:,4) + Fy*Time_Step/Mass_n;
    Electron_State(:,1:2) = Electron_State(:,1:2) + Time_Step*Electron_State(:,3:4);
    %Checking boundary conditions using Matlab matrix equations/indexing
    %Check (if) and move all electrons at X=200nm Bound:
    Electron_State((Electron_State(:,1)>Length),1) = Electron_State((Electron_State(:,1)>Length),1) - Length;
    
    %Check (if) and move all electrons at X=0nm Bound:
    Electron_State((Electron_State(:,1)<0),1) =Electron_State((Electron_State(:,1)<0),1) + Length;
    
    %Check (if) and move all electrons at Y Bounds and if specular or diffusive:
    Electron_State((Electron_State(:,2)>Height & Top_Specular==1),4) = -1*Electron_State((Electron_State(:,2)>Height & Top_Specular),4);
    Electron_State((Electron_State(:,2)>Height & Top_Specular==1),2) = 2*Height - Electron_State((Electron_State(:,2)>Height & Top_Specular),2);
    Electron_State((Electron_State(:,2)>Height & Top_Specular==0),4) = -random(Velocity_PDF);
    Electron_State((Electron_State(:,2)>Height & Top_Specular==0),3) = random(Velocity_PDF);
    Electron_State((Electron_State(:,2)<0 & Bottom_Specular==1),4) = -Electron_State((Electron_State(:,2)<0 & Bottom_Specular==1),4);
    Electron_State((Electron_State(:,2)<0 & Bottom_Specular==1),2) = -Electron_State((Electron_State(:,2)<0 & Bottom_Specular==1),2);
    Electron_State((Electron_State(:,2)<0 & Bottom_Specular==0),4) = random(Velocity_PDF);
    Electron_State((Electron_State(:,2)<0 & Bottom_Specular==0),3) = random(Velocity_PDF);
    %%
    %BOX SPECULAR CASE :
    %Check if bottom Box is specular and if so then bounce
    Electron_State( (Bottom_Box_Specular==1 & Electron_State(:,1)>Box_pos(1,1) & Electron_State(:,1)<Box_pos(1,2) & ((Electron_State(:,2)<Box_pos(1,4) & Electron_Prev_State(:,2)<Box_pos(1,4)))),3)...
        =-Electron_State( (Bottom_Box_Specular==1 & Electron_State(:,1)>Box_pos(1,1) & Electron_State(:,1)<Box_pos(1,2) & ((Electron_State(:,2)<Box_pos(1,4) & Electron_Prev_State(:,2)<Box_pos(1,4)))),3);
    
    %Check if top Box is specular and if so then bounce
    Electron_State( (Top_Box_Specular==1 & Electron_State(:,1)>Box_pos(2,1) & Electron_State(:,1)<Box_pos(2,2) & ((Electron_State(:,2)>Box_pos(2,3) & Electron_Prev_State(:,2)>Box_pos(2,3)))),3)...
        =-Electron_State( (Top_Box_Specular==1 & Electron_State(:,1)>Box_pos(2,1) & Electron_State(:,1)<Box_pos(2,2) & ((Electron_State(:,2)>Box_pos(2,3) & Electron_Prev_State(:,2)>Box_pos(2,3)))),3);
    
    %Top Horizontal
    Electron_State( (Top_Box_Specular==1 & Electron_State(:,1)<Box_pos(2,2) & Electron_State(:,1)>Box_pos(2,1)) & (Electron_State(:,2)>Box_pos(2,3)) &...
        (Electron_Prev_State(:,1)>Box_pos(2,1) & Electron_Prev_State(:,1)<Box_pos(2,2)),4) = -Electron_State( (Top_Box_Specular==1 & Electron_State(:,1)<Box_pos(2,2) & Electron_State(:,1)>Box_pos(2,1)) & (Electron_State(:,2)>Box_pos(2,3)) &...
        (Electron_Prev_State(:,1)>Box_pos(2,1) & Electron_Prev_State(:,1)<Box_pos(2,2)),4);
    %Bottom Hoizontal
    Electron_State( (Bottom_Box_Specular==1 & Electron_State(:,1)<Box_pos(1,2) & Electron_State(:,1)>Box_pos(1,1)) & (Electron_State(:,2)<Box_pos(1,4)) &...
        (Electron_Prev_State(:,1)>Box_pos(1,1) & Electron_Prev_State(:,1)<Box_pos(1,2)),4) = -Electron_State( (Bottom_Box_Specular==1 & Electron_State(:,1)<Box_pos(1,2) & Electron_State(:,1)>Box_pos(1,1)) & (Electron_State(:,2)<Box_pos(1,4)) &...
        (Electron_Prev_State(:,1)>Box_pos(1,1) & Electron_Prev_State(:,1)<Box_pos(1,2)),4);
    %BOX SPECULAR CASE END :
    %%
    %BOX DIFFUSIVE CASE:
    %Check if bottom Box is specular and if so then bounce
    Electron_State( (Bottom_Box_Specular==0 & Electron_State(:,1)>Box_pos(1,1) & Electron_State(:,1)<Box_pos(1,2) & ((Electron_State(:,2)<Box_pos(1,4) & Electron_Prev_State(:,2)<Box_pos(1,4)))),3)...
        =random(Velocity_PDF);
    Electron_State( (Bottom_Box_Specular==0 & Electron_State(:,1)>Box_pos(1,1) & Electron_State(:,1)<Box_pos(1,4) & ((Electron_State(:,2)<Box_pos(1,4) & Electron_Prev_State(:,2)<Box_pos(1,4)))),4)...
        =random(Velocity_PDF);
    
    %Check if Top Box is specular and if so then bounce
    Electron_State( (Top_Box_Specular==0 & Electron_State(:,1)>Box_pos(2,1) & Electron_State(:,1)<Box_pos(2,2) & ((Electron_State(:,2)>Box_pos(2,3) & Electron_Prev_State(:,2)>Box_pos(2,3)))),3)...
        =random(Velocity_PDF);
    Electron_State( (Top_Box_Specular==0 & Electron_State(:,1)>Box_pos(2,1) & Electron_State(:,1)<Box_pos(2,2) & ((Electron_State(:,2)>Box_pos(2,3)& Electron_Prev_State(:,2)>Box_pos(2,3)))),4)...
        =random(Velocity_PDF);
    
    %Bottom Horizontal
    Electron_State( (Bottom_Box_Specular==0 & Electron_State(:,1)<Box_pos(1,2)& Electron_State(:,1)>Box_pos(1,1)) & (Electron_State(:,2)<Box_pos(1,4)) &...
        (Electron_Prev_State(:,1)>Box_pos(1,1) & Electron_Prev_State(:,1)<Box_pos(1,2)),3) = random(Velocity_PDF);
    Electron_State( (Bottom_Box_Specular==0 & Electron_State(:,1)<Box_pos(1,2) & Electron_State(:,4)<0 & Electron_State(:,1)>Box_pos(1,1)) & (Electron_State(:,2)<Box_pos(1,4)) &...
        (Electron_Prev_State(:,1)>Box_pos(1,1) & Electron_Prev_State(:,1)<Box_pos(1,2)),4) = random(Velocity_PDF);
    %Top Hoizontal
    Electron_State( (Top_Box_Specular==0  & Electron_State(:,1)<Box_pos(2,2) & Electron_State(:,1)>Box_pos(2,1)) & (Electron_State(:,2)>Box_pos(2,3)) &...
        (Electron_Prev_State(:,1)>Box_pos(2,1) & Electron_Prev_State(:,1)<Box_pos(2,2)),3) = random(Velocity_PDF);
    Electron_State( (Top_Box_Specular==0 & Electron_State(:,4)>0 & Electron_State(:,1)<Box_pos(2,2) & Electron_State(:,1)>Box_pos(2,1)) & (Electron_State(:,2)>Box_pos(2,3)) &...
        (Electron_Prev_State(:,1)>Box_pos(2,1) & Electron_Prev_State(:,1)<Box_pos(2,2)),4) = -random(Velocity_PDF);
    %%
    %Add scattering
    j = rand(nElectrons,1) < P_Scattering;
    Electron_State(j,3:4) = random(Velocity_PDF,[sum(j),2]);
    %%
    %Electron shifting
    %Region 1 (Bottom  left)
    Electron_State((Electron_State(:,1)>Box_pos(1,1) & Electron_State(:,1)<Box_pos(1,2) & ((Electron_State(:,2)<Box_pos(1,4) & Electron_Prev_State(:,2)<Box_pos(1,4))) &...
        Electron_Prev_State(:,1)<Box_pos(1,1)),1) = 2*Box_pos(1,1) - Electron_State( (Electron_State(:,1)>Box_pos(1,1) & Electron_State(:,1)<Box_pos(1,2) & ((Electron_State(:,2)<Box_pos(1,4) & Electron_Prev_State(:,2)<Box_pos(1,4))) &...
        Electron_Prev_State(:,1)<Box_pos(1,1)),1);
    %Region 3 (Top left)
    Electron_State((Electron_State(:,1)>Box_pos(2,1) & Electron_State(:,1)<Box_pos(2,2) & ((Electron_State(:,2)>Box_pos(2,3) & Electron_Prev_State(:,2)>Box_pos(2,3))) &...
        Electron_Prev_State(:,1)<Box_pos(2,1)),1) = 2*Box_pos(2,1) - Electron_State( (Electron_State(:,1)>Box_pos(2,1) & Electron_State(:,1)<Box_pos(2,2) & ((Electron_State(:,2)>Box_pos(2,3) & Electron_Prev_State(:,2)>Box_pos(2,3))) &...
        Electron_Prev_State(:,1)<Box_pos(2,1)),1);
    %Region 2 (Bottom right)
    Electron_State((Electron_State(:,1)>Box_pos(1,1) & Electron_State(:,1)<Box_pos(1,2) & ((Electron_State(:,2)<Box_pos(1,4) & Electron_Prev_State(:,2)<Box_pos(1,4))) &...
        Electron_Prev_State(:,1)>Box_pos(1,2)),1) = 2*Box_pos(1,2) - Electron_State( (Electron_State(:,1)>Box_pos(1,1) & Electron_State(:,1)<Box_pos(1,2) & ((Electron_State(:,2)<Box_pos(1,4) & Electron_Prev_State(:,2)<Box_pos(1,4))) &...
        Electron_Prev_State(:,1)>Box_pos(1,2)),1);
    %Region 4 (Top Right)
    Electron_State((Electron_State(:,1)>Box_pos(2,1) & Electron_State(:,1)<Box_pos(2,2) & ((Electron_State(:,2)>Box_pos(2,3) & Electron_Prev_State(:,2)>Box_pos(2,3))) &...
        Electron_Prev_State(:,1)>Box_pos(2,2)),1) = 2*Box_pos(2,2) - Electron_State( (Electron_State(:,1)>Box_pos(2,1) & Electron_State(:,1)<Box_pos(2,2) & ((Electron_State(:,2)>Box_pos(2,3) & Electron_Prev_State(:,2)>Box_pos(2,3))) &...
        Electron_Prev_State(:,1)>Box_pos(2,2)),1);
    
    %Bottom Horizontal (going top down)
    Electron_State( (Electron_State(:,1)<Box_pos(1,2) & Electron_State(:,1)>Box_pos(1,1) & (Electron_State(:,2)<Box_pos(1,4) & Electron_Prev_State(:,2)>Box_pos(1,4)) &...
        (Electron_State(:,4)<0)),2) =  2*Box_pos(1,4) + Electron_State( (Electron_State(:,1)<Box_pos(1,2) & Electron_State(:,1)>Box_pos(1,1) & (Electron_State(:,2)<Box_pos(1,4) & Electron_Prev_State(:,2)>Box_pos(1,4)) &...
        (Electron_State(:,4)<0)),2);
    %Top Horizontal (going bottom up)
    Electron_State( (Electron_State(:,1)<Box_pos(2,2) & Electron_State(:,1)>Box_pos(2,1) & (Electron_State(:,2)>Box_pos(2,3) & Electron_Prev_State(:,2)<Box_pos(2,3)) &...
        (Electron_State(:,4)>0)),2) =  2*Box_pos(2,3) - Electron_State( (Electron_State(:,1)<Box_pos(2,2) & Electron_State(:,1)>Box_pos(2,1) & (Electron_State(:,2)>Box_pos(2,3) & Electron_Prev_State(:,2)<Box_pos(2,3)) &...
        (Electron_State(:,4)>0)),2);
    %Electron shifting end
    %%
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
        % Add the temperature and current data to the respective plots
        Temp_Plot(i,1) = Temperature(i);
        Current_Plot(i,1) = J(i,1)^2 + J(i,2)^2;
        %Shows the pathing of the electron, as well as the updating trajectory
        if Show_Movie && mod(i,2)
            figure(1)
            plot(Electron_State(1:nPlotted_Electrons,1)./1e-9,Electron_State(1:nPlotted_Electrons,2)./1e-9,'o');
            grid on;
            axis([0 Length/1e-9 0 Height/1e-9]);
            xlabel('x (nm)');
            ylabel('y (nm)');
            title(sprintf("Plotting (%d/%d) electron at constant velocity",nPlotted_Electrons,nElectrons));
            hold on;
            for j=1:size(Box_pos,1)
                plot([Box_pos(j, 1) Box_pos(j, 1) Box_pos(j, 2) Box_pos(j, 2) Box_pos(j, 1)]./1e-9,...
                    [Box_pos(j, 3) Box_pos(j, 4) Box_pos(j, 4) Box_pos(j, 3) Box_pos(j, 3)]./1e-9, 'k-');
            end
            hold off;
            pause(0.01)
        end
    end
%%
figure("name","Trajectory, temperature and speed results results")
subplot(2,2,1)
plot(Trajectories_x(:,1:nPlotted_Electrons)./1e-9, Trajectories_y(:,1:nPlotted_Electrons)./1e-9,'.');
hold on;
for j=1:size(Box_pos,1)
    plot([Box_pos(j, 1) Box_pos(j, 1) Box_pos(j, 2) Box_pos(j, 2) Box_pos(j, 1)]./1e-9,...
        [Box_pos(j, 3) Box_pos(j, 4) Box_pos(j, 4) Box_pos(j, 3) Box_pos(j, 3)]./1e-9, 'k-');
end
hold off
axis([0 Length/1e-9 0 Height/1e-9]);
xlabel('x (nm)');
ylabel('y (nm)');
grid on;
title(sprintf("Trajectories of (%.2d/%.2d) electron(s) at constant velocity",nPlotted_Electrons,nElectrons),'interpreter','latex');
subplot(2,2,2)
plot(Time_Step*(0:Iterations-1), Temperature);
grid on;
xlim([0 Time_Step*Iterations])
title(sprintf("Temperature of the region, Average Temperature: %.2f (K)",mean(Temperature)),'interpreter','latex')
xlabel('Time (s)');
ylabel('Temperature (K)');
subplot(2,2,3)
Velocity = sqrt(Electron_State(:,3).^2 + Electron_State(:,4).^2);
histogram(Velocity);
title(sprintf("Electron Velocity, Average Velocity %.2d (m/s)",mean(Velocity)),'interpreter','latex');
xlabel("Speed (m/s)");
ylabel("Number of particles");
grid on;
subplot(2,2,4)
plot(Time_Step*(0:Iterations-1),Current_Plot)
title(sprintf("Jx - Drift Current Density along x-axis, Average : %.2d (A/m)",mean(Current_Plot)),'interpreter','latex');
xlabel('Time (s)');
ylabel('Current Density (A/m)');
grid on;
axis tight;
%saveas(gcf,fullfile('D:\School Work\ELEC 4700\My 4700 Code\Assignment 3\Simulation Results','[Part3]Trajectory_Temperature_DrifCurr.png'),'png')
%%
Density = hist3(Electron_State(:,1:2),[200 100])';
N = 20;
sigma = 3;
%Creating a Gaussian filtering matrix
[x,y]=meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
G=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
G=G./sum(G(:));
Density = conv2(Density,G,'same') / (Height./size(Density,1)*Length./size(Density,2));
%%
%Plot drift current
figure('Name','Drift current (x-axis) and Electron Density')
subplot(2,1,1)
surf(conv2(Density,G,'same'));
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('Bin count')
axis tight
title("Electron density region (3D view)",'interpreter','latex');
%%
% Plot the electron density
subplot(2,1,2)
surf(conv2(Density,G,'same'));
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('Bin count')
axis tight
title("Electron density region (Top down view)",'interpreter','latex');
view(2)
%saveas(gcf,fullfile('D:\School Work\ELEC 4700\My 4700 Code\Assignment 3\Simulation Results','[Part3]ElectronDensity.png'),'png')