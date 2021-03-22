%{
Author : Maharshi Gurjar
Elec 4700 Assignment 3 -  Finite Difference Method 
Part 2
%}
%%
clc; close all; clear;
set(0, 'DefaultFigureWindowStyle', 'docked')
%Define simulation parameters
W = 100e-9; %Width of region
L = 200e-9; %Length of region
L_Box = 40e-9; %Length of Box
W_Box = 40e-9; %Width of Box
MeshSize = 1e-9; % Preset mesh size
nx = ceil(L/MeshSize); %X bins
ny = ceil(W/MeshSize); %Y bins

Conductivity_1 = 1; %Conductivity of region
Conductitivty_2 = 1e-2; %Conductivity in box

Conductivity_Map = zeros(nx,ny);

for i = 1:nx
    for j = 1:ny
        if ((i-1)>0.5*(L-L_Box)/MeshSize && (i-1)<0.5*(L+L_Box)/MeshSize && ((j-1)<W_Box/MeshSize || (j-1)>(W-W_Box)/MeshSize))
            Conductivity_Map(i,j) = Conductitivty_2;
        else
            Conductivity_Map(i,j) = Conductivity_1;
        end 
    end
end
% Numerial issues can happen if the derivatives are too large, thus using
% the imgaussfilt function the matrix is filtered
Conductivity_Map = imgaussfilt(Conductivity_Map,1);


figure(1);
imagesc(linspace(0,L,nx),linspace(0,W,ny),Conductivity_Map);
title('Conductivity of region','interpreter','latex')
xlabel('x (m)')
ylabel('y (m)')
view(2)
axis tight
grid on;
%saveas(gcf,fullfile('D:\School Work\ELEC 4700\My 4700 Code\Assignment 3\Simulation Results','[Part2]ConductivityMap.png'),'png')
G = sparse(ny*nx);
F = zeros(1,ny*nx);


for i = 1:(nx)
    for j = 1:(ny)
        %Create variable for node mapping
        n = j + (i-1)*ny;
        %Calculate changes in x and y 
        if i == 1 % V = 1, x = 0
            G(n,n) = 1;
            F(n) = 0.1;
        elseif i == nx %V=0, x = L
            G(n,n) = 1;
        elseif j == 1 
            xM = j + (i-2)*ny;
            xP = j + i*ny;
            yP = j+1 + (i-1)*ny;
            %Resistance Values 
            rxM = (Conductivity_Map(i,j) + Conductivity_Map(i-1,j))/2;
            rxP = (Conductivity_Map(i,j) + Conductivity_Map(i+1,j))/2;
            ryP= (Conductivity_Map(i,j) + Conductivity_Map(i,j+1))/2;
            %node equations 
            G(n,n) = -(rxM + rxP + ryP);
            G(n,xM) = rxM;
            G(n,xP) = rxP;
            G(n,yP) = ryP; 
        elseif j == ny%BC @ y=W
            xM = j + (i-2)*ny;
            xP = j + i*ny;
            nym = j-1 + (i-1)*ny;
            %Resistance Values 
            rxM = (Conductivity_Map(i,j) + Conductivity_Map(i-1,j))/2;
            rxP = (Conductivity_Map(i,j) + Conductivity_Map(i+1,j))/2;
            rym = (Conductivity_Map(i,j) + Conductivity_Map(i,j-1))/2;
            %node equations 
            G(n,n) = -(rxM + rxP + rym);
            G(n,xM) = rxM;
            G(n,xP) = rxP;
            G(n,nym) = rym;
        %internal nodes
        else
            xM = j + (i-2)*ny;
            xP = j + i*ny;
            nym = j-1 + (i-1)*ny;
            yP = j+1 + (i-1)*ny;
            %Resistor Values
            rxM = (Conductivity_Map(i,j) + Conductivity_Map(i-1,j))/2;
            rxP = (Conductivity_Map(i,j) + Conductivity_Map(i+1,j))/2;
            ryp = (Conductivity_Map(i,j) + Conductivity_Map(i,j+1))/2;
            rym = (Conductivity_Map(i,j) + Conductivity_Map(i,j-1))/2;
            %node equations
            G(n,n) = -(rxM + rxP + rym + ryp);
            G(n,xM) = rxM;
            G(n,xP) = rxP;
            G(n,nym) = rym;
            G(n,yP) = ryp; 
        end
    end
end

V = G\F';
Voltage_Map = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        Voltage_Map(i,j) = V(j + (i-1)*ny);
    end
end

%[X , Y] = meshgrid(0:MeshSize:L,0:MeshSize:W);
figure(2)
subplot(2,1,1)
surf(Voltage_Map) %X',Y',
axis tight
xlabel('x posistion (nm)')
ylabel('y posistion (nm)')
zlabel('Voltage (V)')
title('Voltage across region (3D view)','interpreter','Latex')
view(90,25)
subplot(2,1,2)
surf(Voltage_Map) %X',Y',
axis tight
xlabel('x posistion (nm)')
ylabel('y posistion (nm)')
zlabel('Voltage (V)')
title('Voltage across region (Top down view)','interpreter','Latex')
view(2)
%saveas(gcf,fullfile('D:\School Work\ELEC 4700\My 4700 Code\Assignment 3\Simulation Results','[Part2]VoltageMap.png'),'png')

[Ex, Ey] = gradient(Voltage_Map,MeshSize);
Ex=-Ex;
Ey=-Ey;
figure(3)
quiver(Ex,Ey);
axis tight
xlabel('x posistion (nm)')
ylabel('y posistion (nm)')
title('Electric Field across region','interpreter','Latex')
%saveas(gcf,fullfile('D:\School Work\ELEC 4700\My 4700 Code\Assignment 3\Simulation Results','[Part2]QuiverPlot.png'),'png')