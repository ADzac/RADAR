addpath('C:\Users\ASUS\Documents\MATLAB');
% Define parameters
a = 1; % Radius of the sphere (in meters)
epsilon_r = 2; % Relative permittivity of the sphere material
mu_r = 1; % Relative permeability of the sphere material
lambda = 1; % Wavelength of the incident wave (in meters)
k = 2*pi/lambda; % Wave number

% Incident and observation angles (in radians)
theta_i = deg2rad(30); % Incident angle
theta_o = deg2rad(45); % Observation angle

% Calculate Mie coefficients (a_n, b_n) using mie scattering functions
[n, a_n, b_n] = plotBistaticRCS(epsilon_r, mu_r, k*a);

% Compute bistatic RCS for a dielectric sphere using Mie coefficients
RCS_bistatic = zeros(size(theta_i)); % Initialize RCS array

for i = 1:length(theta_i)
    % Calculate scattering amplitude matrix element
    S = a_n + b_n; % Scattering amplitude matrix element (sum of Mie coefficients)
    
    % Compute bistatic RCS (absolute value of the scattering amplitude)
    RCS_bistatic(i) = abs(S)^2; % Bistatic RCS for the given angles
end

% Plot RCS vs. incident angle
figure;
plot(rad2deg(theta_i), RCS_bistatic, 'LineWidth', 2);
xlabel('Incident Angle (degrees)');
ylabel('Bistatic RCS');
title('Bistatic RCS of Dielectric Sphere');
grid on;
