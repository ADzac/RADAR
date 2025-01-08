clear;
close all;
clc; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Defining all the parameters
z_coord = 0;

%Numbers of antennas
N_Rx = 5; %Receiver on x-axis
N_Ry = 5; %Receiver on y-axis
N_Tx = 4;
N_Ty = 1;

%Frequency & Spacing
Nf = 1;  %always odd(impair)
c = 3e8;
f_min = 1e9;
f_max = 1e10;
%bande_f = linspace(f_min,f_max,Nf); %If more than 1 frequency
fc = 1e9; %(f_max+f_min)/2;
lc = c/fc;
shannon = lc/2;
k  = 2*pi/lc;

%Length of the line
% Lx = shannon*N_Rx
Lx = 1.5;
Ly = 1;
space_Rx = Lx/N_Rx;
space_Ry = Ly/N_Ry;

%Finding the place of the receiving antennas centered on 0
minRx = (-1)*((N_Rx-1)*shannon)/2;
maxRx = ((N_Rx-1)*shannon)/2;
minRy = (-1)*((N_Ry-1)*shannon)/2;
maxRy = ((N_Ry-1)*shannon)/2;

%Finding the place of the transmitting antennas centered on 0
minTx = minRx;
maxTx = maxRx;
minTy = (-1)*((N_Ty-1)*shannon)/2;
maxTy = ((N_Ty-1)*shannon)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Creating a coordinate for the receiving antennas
X_R = zeros(N_Rx,N_Ry,3);
if N_Ry == 1
    yR_coord = 0;
    for i = 1:N_Ry
        X_R(:,i,1) = linspace(minRx, maxRx,N_Rx); 
        X_R(:,i,2) = yR_coord; 
        X_R(:,i,3) = z_coord;
    end
else
    yR_coord = linspace(minRy,maxRy,N_Ry);
    for i = 1:N_Ry
        X_R(:,i,1) = linspace(minRx, maxRx,N_Rx); 
        X_R(:,i,2) = yR_coord(i); 
        X_R(:,i,3) = z_coord;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating a coordinate for the transmitting antennas
X_T = zeros(N_Tx,N_Ty,3);


% If want to choose the placement of the transmittter
% X_T(1,:) = [0,0,0];
% X_T(2,:) = [0,0,0];

% If want to create a well coordinated transmittter 
% The transmitters are place at the end of the line for even number
% They are always centered on zero for odd number of transmitter
if N_Ty == 1
    yR_coord = 0;
    for i = 1:N_Ty
        X_T(:,i,1) = linspace(minTx-shannon, maxTx+shannon,N_Tx); 
        X_T(:,i,2) = yR_coord; 
        X_T(:,i,3) = z_coord;
    end
else
    yT_coord = linspace(minTy,maxTy,N_Ty);
    for i = 1:N_Ty
        X_T(:,i,1) = linspace(minTx-shannon, maxTx+shannon,N_Tx); 
        X_T(:,i,2) = yT_coord(i); 
        X_T(:,i,3) = z_coord;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Placing a target into the surface z = zc
zc = 3;
N_cible = 1; %Number of target
%List of target(s) with its coordinates
X_C(1,:) = [0,0,3]; 
% X_C(2,:) = [1,0,3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% distances for each Receiver with target
nR = size(X_R);
distancesList_R = zeros(nR(1),nR(2),N_cible);
for h = 1:N_cible
    for i = 1:nR(2)
        for j = 1:nR(1)
            distancesList_R(j,i,h) = sqrt((X_R(j,i,1) - X_C(h,1)).^2 + (X_R(j,i,2) - X_C(h,2)).^2 + (X_R(j,i,3) - X_C(h,3)).^2);
        end
    end
end

% distances for each Receiver with target
nT = size(X_T);
distancesList_T = zeros(nT(1),nT(2),N_cible); % based on the number of transmitter on X and Y and number of targets
for h = 1:N_cible
    for i = 1:nT(2)
        for j = 1:nT(1)
            distancesList_T(j,i,h) = sqrt((X_T(j,i,1) - X_C(h,1)).^2 + (X_T(j,i,2) - X_C(h,2)).^2 + (X_T(j,i,3) - X_C(h,3)).^2);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%X_voxels et Y_voxels
zVoxels_coords = (N_cible);
for i=1:N_cible
    zVoxels_coords(i) = X_C(i,3);
end

L_XVoxels = 6;
XVoxels_min = -1*(L_XVoxels/2);
XVoxels_max = (L_XVoxels/2);
N_XVoxels = 200; % Space between each Voxels is L_Xvoxels/N_XVoxels = 0.4
N_YVoxels = 1;
deltaY = 0.2;
Voxels = zeros(N_XVoxels,N_YVoxels,3);

if N_cible ==1
    if N_YVoxels == 1
        yVoxels_coord = 0;
        for i = 1:N_YVoxels
            Voxels(:,i,1) = linspace(XVoxels_min, XVoxels_max,N_XVoxels); 
            Voxels(:,i,2) = yVoxels_coord;
            Voxels(:,i,3) = zVoxels_coords;
        end
    else
        yVoxels_coord = linspace(-deltaY,deltaY,N_YVoxels);
        for i = 1:N_YVoxels
            Voxels(:,i,1) = linspace(XVoxels_min, XVoxels_max,N_XVoxels); 
            Voxels(:,i,2) = yVoxels_coord(i);
            Voxels(:,i,3) = zVoxels_coords;
        end
    end
else
    % number of targets more than 1 ???
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Distance between Transmitter and Voxels

nTV = size(X_T);
distancesList_TV = zeros(N_XVoxels*N_YVoxels,nTV(1),nTV(2));
for h = 1:N_XVoxels
    for i = 1:nTV(2) %Number of Y lines
        for j = 1:nTV(1) %Number of T on X axis
            distancesList_TV(h,j,i) = sqrt((X_T(j,i,1) - Voxels(h,:,1)).^2 + (X_T(j,i,2) - Voxels(h,:,2)).^2 + (X_T(j,i,3) - Voxels(h,:,3)).^2);
        end
    end
end

%Distance between Receiver and Voxels
nRV = size(X_R);
distancesList_RV = zeros(N_XVoxels*N_YVoxels,nRV(1),nRV(2));
for h = 1:N_XVoxels
    for i = 1:nRV(2) %Number of R on Y
        for j = 1:nRV(1) %Number of R on X
            distancesList_RV(h,j,i) = sqrt((X_R(j,i,1) - Voxels(h,1)).^2 + (X_R(j,i,2) - Voxels(h,2)).^2 + (X_R(j,i,3) - Voxels(h,3)).^2);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%S
S = zeros(Nf);
for l = 1:Nf
    for m = 1:N_Tx*N_Ty %Number of Transmitter
        for n = 1:N_Rx*N_Ry %Number of Receiver
            %S (frequence,which Transmitter,which receiver)
            S(l,m,n) = exp(((-2i*pi*fc)/c)*(distancesList_T(m)+distancesList_R(n))); %bande_f(l)
        end 
    end
end

%Compare addition with sum to see if its the same or not

Somme = sum(S,4);
Compare = false(size(Somme));
% for j = 1:N_cible
    for i = 1:N_Tx*N_Ty*1:N_Rx*N_Ry
        if(S(:,:,i)==Somme(i))
            Compare(i) = true;
        end
    end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Image
% The formula of Image is taken from Microwave Imaging in Security-Two
% Decades in Innovation by Sherif S.Ahmed on page 193
Image = zeros(Nf,N_XVoxels);

for l = 1:Nf
    for k = 1: length(Voxels)
        for m = 1:N_Tx*N_Ty %Number of Transmitter
            for n = 1:N_Rx*N_Ry %Number of Receiver
                Image(l,k) = Image(l,k) + S(l,m,n)*exp(((2i*pi*fc)/c)*(distancesList_TV(k,m)+distancesList_RV(k,n))); %bande_f(l)          
            end
        end
    end
end

% % Initialize variables to store the maximum amplitude and its placement
% max_amplitude = 0;
% max_amplitude_index = 0;
% 
% % Iterate through the list and find the maximum amplitude
% for i = 1:numel(Image)
%     current_complex = Image(i);
%     current_amplitude = abs(current_complex);
% 
%     fprintf('Current Amplitude: %.14f\n', current_amplitude);
% 
%     if current_amplitude > max_amplitude
%         max_amplitude = current_amplitude;
%         max_amplitude_index = Voxels(i,1,1);
%     end
% end

% fprintf('Amplitude: %.14f\n', max_amplitude);
% fprintf('Placement in the list: X-Voxels = %d\n', max_amplitude_index);


% Create a scatter plot
figure;
plot(Voxels(:,:,1),abs(Image) ,'o-', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
xlabel('X-Coordinates of Voxels');
ylabel('Amplitude');
title('Images');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define a grid of points
figure;
[x, y,z] = meshgrid(-Lx:0.1:Lx, -Ly:0.1:Ly,-Lx:0.1:Lx);
plot3(X_T(:,:,1), X_T(:,:,2), X_T(:,:,3), 'x', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
hold on;
plot3(X_R(:,:,1), X_R(:,:,2), X_R(:,:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
plot3(X_C(:,1), X_C(:,2), X_C(:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
% plot3(Voxels(:,:,1), Voxels(:,:,2), Voxels(:,:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
% Add labels and title
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
legend('Transmitter', 'Receiver', 'Target', 'Location', 'Northeast');
% legend('Transmitter', 'Receiver', 'Target', 'Voxels', 'Location', 'Northeast');
view(0,90);