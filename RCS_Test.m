clc
clear
close all

% Define array of frequencies to analyze
frequencies = [0.5e9, 1e9, 1.5e9]; % Frequencies in Hz

for freq = frequencies
    disp(['Analyzing frequency: ' num2str(freq/1e9) ' GHz']);

    nsteps = 1001;
    theta_SER = linspace(0, 2*pi, nsteps);
    phi_SER = 0;

    c = 3e8;

    radius = 0.5; % Sphere radius in meters
    k = 2*pi*freq/c;
    x = k*radius;

    miu_r = 1; % Permeability of non-magnetic material
    eps_0 = 8.854e-12;
    epsr = 4; % Permitivity of the material (4 for dielectric)
    sigma = 10e-6; % Conductivity of the material (10e-6 for dielectric)
    omega = 2*pi*freq;

    m_r = sqrt((epsr*miu_r/2)*(sqrt(1 + (sigma/(omega*epsr*eps_0))^2) + 1));
    m_i = sqrt((epsr*miu_r/2)*(sqrt(1 + (sigma/(omega*epsr*eps_0))^2) - 1));

    m = m_r + m_i*1i;

    r = 0.6;
    exp_term = exp(1i*k*r);

    m1 = real(m);
    m2 = imag(m);

    First_term = abs((exp_term/(-1i*k*r))) * cos(phi);
    Sec_term = abs((exp_term/(1i*k*r))) * sin(phi);

    for j = 1:length(theta)
        u = cos(theta(j));
        a(:, j) = Mie_S12(m, x, u);
        S1(j) = real(a(1, j)'*a(1, j));
        S2(j) = real(a(2, j)'*a(2, j));
        E_S_theta(j) = First_term * S2(j);
        E_S_phi(j) = Sec_term * S1(j);
    end

    max_theta = max(E_S_theta);
    max_phi = max(E_S_phi);

    E_S_phi = E_S_phi / max_phi;
    E_S_theta = E_S_theta / max_theta;

    if phi == 0
        E = abs(E_S_theta).^2;
        E_used = "E_S_theta";
    elseif phi == pi/2
        E = abs(E_S_phi).^2;
        E = reshape(E, [nsteps, 1]);
        E_used = "E_S_phi";
    else
        E = abs(E_S_phi) .* abs(E_S_theta);
        E = reshape(E, [nsteps, 1]);
        E_used = "E_S_theta and E_S_phi";
    end

    % Plot
    % Create polar plot of E_S_theta
    figure;
    subplot(1, 2, 1);
    polarplot(theta, E);
    title(['Polar Plot of ' E_used]);

    subplot(1, 2, 2);
    plot(rad2deg(theta), 10*log10(E_S_theta));
    title(['Decibel of ' E_used]);
    xlabel('\phi (degrees)');
    xlim([0, 360]);
end
