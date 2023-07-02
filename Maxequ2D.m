
%% A Matlab code developed by HUSSEIN A. H. Muhammed March 2023: B.Sc.H and M.Sc. (Honuors).
%% A Finite-Difference Time-Domain program for Electromagnetic wavefields 2D Maxwell's Equations Simulation.
%% Please cite this code as: Hussein Abduelhaleim Hussein Muhammed (2022), Least-Squares ReverseTime Migration
%% in Pseudodepth Domain and Its Application. China University of Petroleum (East China),
%% Master Thesis, School of Geosciences, Dept. of Geophysics, SWPI lab. Library press©.

clc;
clear;
close all;

%% Prepare the movie file
    % 
    % vidObj = VideoWriter('coupled-BE-2d.avi');
    % open(vidObj);

% field constants
c = 299792458;      % Speed of light in vacuum (m/s)
epsilon0 = 8.854e-12;  % Permittivity of free space (F/m)
mu0 = 4*pi*1e-7;    % Permeability of free space (H/m)

% space-time parameters
dx = 0.01;          % Spatial step size (m)
dt = dx/(2*c);      % Time step size (s)
num_steps = 280;    % Number of time steps (notice the numerical disperssion after step 280).

% FDTD grid size
Nx = 100;           % Number of grid points along x-axis
Ny = 100;           % Number of grid points along y-axis

% Electromagnetic wavefield solution arrays
Ex = zeros(Nx, Ny); % x-component of electric field
Ey = zeros(Nx, Ny); % y-component of electric field
Hz = zeros(Nx, Ny); % z-component of magnetic field

% Main simulation loop
for t = 1:num_steps
    % Update electric field
    for i = 2:Nx-1
        for j = 2:Ny-1
            Ex(i, j) = Ex(i, j) + (dt/(epsilon0*dx)) * (Hz(i, j) - Hz(i, j-1));
            Ey(i, j) = Ey(i, j) - (dt/(epsilon0*dx)) * (Hz(i, j) - Hz(i-1, j));
        end
    end
    
    % Update magnetic field
    for i = 1:Nx-1
        for j = 1:Ny-1
            Hz(i, j) = Hz(i, j) + (dt/(mu0*dx)) * (Ey(i+1, j) - Ey(i, j)) - (dt/(mu0*dx)) * (Ex(i, j+1) - Ex(i, j));
        end
    end
    
    % Source function (Gaussian pulse)
    t0 = 40; % Pulse center time
    spread = 12; % Pulse spread
    pulse = exp(-0.5*((t-t0)/spread)^2);
    
    % Add source to Ex field (at the center point)
    Ex(Nx/2, Ny/2) = pulse;
    
    % Visualization (i.e., plotting the electric field magnitude)
    if mod(t, 10) == 0
        figure(1);
        imagesc(sqrt(Ex.^2 + Ey.^2)');
        colormap(jet);
        c = colorbar;
        c.Label.String = 'Electromagnetic wavefield';
        title(['Time Step: ', num2str(t)]);
        xlabel('x');
        ylabel('y');
        drawnow;
    end
     % Write each frame to the file
       % currFrame = getframe(gcf);
       % writeVideo(vidObj,currFrame);

end

% %% Close the file
% close(vidObj);