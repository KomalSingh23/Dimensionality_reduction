% Parameters
L = 1;          % Length of the square domain
Nx = 50;        % Number of spatial grid points in x-direction
Ny = 50;        % Number of spatial grid points in y-direction
T = 1;          % Total simulation time
Nt = 100;       % Number of time steps
alpha = 0.01;   % Thermal diffusivity

% Discretization
dx = L / Nx;
dy = L / Ny;
dt = T / Nt;

% Initialize temperature field with initial conditions
u = 40 * ones(Nx, Ny);  % Set the entire domain to 40°C initially

% Specify different regions with initial temperatures
u(10:20, 10:20) = 90;   % Region with 90°C
u(30:40, 30:40) = 150;  % Region with 150°C
u(20:30, 40:50) = 500;  % Region with 500°C

% Initialize storage for snapshots
snapshots = zeros(Nx * Ny, Nt);

% Measure the simulation running time
tic;

% Main time-stepping loop
for n = 1:Nt
    u_new = u;
    for i = 2:Nx-1
        for j = 2:Ny-1
            % 2D diffusion equation
            u_new(i, j) = u(i, j) + alpha * dt * (u(i+1, j) + u(i-1, j) + u(i, j+1) + u(i, j-1) - 4 * u(i, j)) / (dx^2 + dy^2);
        end
    end
    u = u_new;
    
    % Store snapshots for analysis
    snapshots(:, n) = u(:);
end

% Measure the total simulation time
simulation_time = toc;

% Visualization of the temperature field
[X, Y] = meshgrid(linspace(0, L, Nx), linspace(0, L, Ny));
u_final = reshape(u, Nx, Ny);

% Plot the temperature field
figure;
surf(X, Y, u_final);
xlabel('X');
ylabel('Y');
zlabel('Temperature (°C)');
title('2D Heat Conduction Simulation');

fprintf('Simulation completed in %.2f seconds.\n', simulation_time);
