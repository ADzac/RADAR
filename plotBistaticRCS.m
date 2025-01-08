% function plotBistaticRCS(a, epsilon_r, mu_r, lambda)
    clc;
    clear;
    close all;

% Example usage:
a = 0.1;            % Size parameter
epsilon_r = 4;      % Relative permittivity
mu_r = 1;           % Relative permeability
lambda = 500;    % Incident wavelength in nanometers

    % ***INPUT*** Refractive Index (dimensionless in general but could vary with wavelength)
    m = sqrt(epsilon_r * mu_r); 

    % ***INPUT*** Incident Wavelength (nm)
    lambda = lambda * 1e-9; % Convert to meters

    const = 1; % ***INPUT*** ISOINTENSITY DEGREE (Is/Ii)
    n = 360; % ***INPUT*** Smoothness of Plot

    % Calculate angles in degrees and radians
    thetaprime = linspace(0, 180, n+1);
    thetadeg = linspace(0, 180, nmax);

    thetaprimerad = deg2rad(thetaprime);
    anglerad = deg2rad(linspace(0, 360, 2*n+1));

    % Wavenumber
    k = 2 * pi / lambda;

    % Size parameter (dimensionless)
    x = k * a;

    % Number of terms to calculate S1 and S2
    nmax = round(x + 4 * x^(1/3) + 2);

    % Initialize r array
    r = zeros(size(thetaprime));

    % Compute Mie coefficients using the given parameters
    n = (1:nmax);
    nu = (n + 0.5);
    z = m * x;
    m2 = m * m;
    sqx = sqrt(0.5 * pi / x);
    sqz = sqrt(0.5 * pi / z);

    bx = besselj(nu, x) * sqx;
    bz = besselj(nu, z) * sqz;
    yx = bessely(nu, x) * sqx;
    hx = bx + 1i * yx;
    b1x = [sin(x) / x, bx(1:nmax-1)];
    b1z = [sin(z) / z, bz(1:nmax-1)];
    y1x = [-cos(x) / x, yx(1:nmax-1)];
    h1x = b1x + 1i * y1x;
    ax = x * b1x - n .* bx;
    az = z * b1z - n .* bz;
    ahx = x .* h1x - n .* hx;
    an = (m2 .* bz .* ax - bx .* az) ./ (m2 .* bz .* ahx - hx .* az);
    bn = (bz .* ax - bx .* az) ./ (bz .* ahx - hx .* az);

    % Excel-like Arithmetic
    G = bsxfun(@times, B, E .* C + F .* D);
    H = bsxfun(@times, B, E .* D + F .* C);

    % Calculate scattering amplitudes S1 and S2
    S1 = sum(G);
    S2 = sum(H);

    % Calculate the bistatic RCS
    U = abs(S1).^2 + abs(S2).^2;
    v = (1 / k) * sqrt(U / const);
    r = v;

    % Create a symmetric RCS array for polar plot
    R = [r, fliplr(r(2:end-1))];

    % Plot the bistatic RCS in polar coordinates
    polar(anglerad, R);
% end

