clc
close all
clear

fs = linspace(0.5e9,1.5e9,3);
% Define parameters
miu_r = 1; % Permeability of non-magnetic material
eps_0 = 8.854e-12;
epsr = 1;% Permitivity of the material (4 for dielectric 1 for metallic)
sigma = 1e6; % Conductivity of the material (10e-6 for dielectric)
radius = 0.05;
freq = 0.5e3;
k = 2*pi*freq/3e8;  % Wave number
r = 1.2081;  % Distance from the source
phi_SER = 0;  % Azimuthal angle
theta_SER = linspace(0,2*pi,361);  % Polar angle
omega = 2*pi*0.5e3;

er1 = epsr*miu_r/2;
er2_1 = (sigma/(omega*epsr*eps_0)).^2;
er2_2 = sqrt(1 +er2_1);

m_r = sqrt(er1*(er2_2 + 1));
m_i = sqrt(er1*(er2_2 - 1));

index_refraction = m_r + m_i*1i;
% Complex refractive index
x = k*radius;  % Size parameter

% Preallocate arrays
pp = zeros(2, length(theta_SER));
S1 = zeros(size(theta_SER));
S2 = zeros(size(theta_SER));
E_S_theta = zeros(size(theta_SER));

% Compute terms for electric field components
exp_term = exp(1i*k*r);
m1 = real(index_refraction);
m2 = imag(index_refraction);
First_term = abs((exp_term/(-1i*k*r))* cos(phi_SER));
Sec_term = abs((exp_term/(1i*k*r)) * sin(phi_SER));

% Compute scattering functions and electric field component E_theta
for j = 1:length(theta_SER)
    u = cos(theta_SER(j));
    pp(:, j) = Mie_S12(index_refraction, x, u);
    S1(j) = real(pp(1, j)'*pp(1, j));
    S2(j) = real(pp(2, j)'*pp(2, j));
    E_S_theta(j) = First_term * S2(j);
end

E_S_theta = E_S_theta/max(E_S_theta);

% Polar plot E_theta
polarplot(theta_SER, abs(E_S_theta));
title('Electric Field Component E_{\theta}');
