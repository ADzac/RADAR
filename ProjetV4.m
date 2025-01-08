clear;
close all;
clc; 

% Start the timer
tic;

%% Defining all the parameters
z_coord = 0;

%Numbers of antennas
N_Rx = 12; %Receiver on x-axis
N_Ry = 1; %Receiver on y-axis
N_Tx = 1;
N_Ty = 1;

%Frequency & Spacing
Nf = 1;  %always odd(impair)
c = 3e8;
f_min = 0.5e9;
f_max = 1.5e9;
bande_f = linspace(f_min,f_max,Nf); %If more than 1 frequency
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
minTx = -1*(shannon*(N_Tx-1)*N_Rx)/2;
maxTx = (shannon*(N_Tx-1)*N_Rx)/2;
minTy = (-1)*((N_Ty-1)*shannon)/2;
maxTy = ((N_Ty-1)*shannon)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Creating a coordinate for the receiving antennas
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
%% Creating a coordinate for the transmitting antennas
X_T = zeros(N_Tx,N_Ty,3);

% If want to choose the placement of the transmittter
% X_T(1,:) = [0,0,0];
% X_T(2,:) = [0,0,0];

% If want to create a well coordinated transmittter 
% The transmitters are place at the end of the line for even number
% They are always centered on zero for odd number of transmitter
if N_Ty == 1
    for i = 1:N_Ty
        X_T(:,i,1) = linspace(minTx, maxTx,N_Tx); 
        X_T(:,i,2) = yR_coord; 
        X_T(:,i,3) = z_coord;
    end
else
    yT_coord = linspace(minTy,maxTy,N_Ty);
    for i = 1:N_Ty
        X_T(:,i,1) = linspace(minTx, maxTx,N_Tx); 
        X_T(:,i,2) = yT_coord(i); 
        X_T(:,i,3) = z_coord;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Placing a Target into the surface z = zc

N_cible = 1; %Number of target
%List of target(s) with its coordinates
X_C(1,1,:) = [0,0,3]; 
%X_C(2,1,:) = [1,0,5];
y_cible = X_C(1,1,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Distances for each Receiver with target
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

%% Voxels

% L_XVoxels = 6;
XVoxels_min = -3;
XVoxels_max = 3;
N_XVoxels = 51;

% L_ZVoxels = 6;
ZVoxels_min =  2;
ZVoxels_max = 6;
N_ZVoxels = 101; 

YVoxels_min =  -3;
YVoxels_max = 3;
N_YVoxels = 1; 

Voxels = zeros(N_XVoxels,N_ZVoxels,3);
XV = linspace(XVoxels_min, XVoxels_max,N_XVoxels);
ZV = linspace(ZVoxels_min, ZVoxels_max,N_ZVoxels);
YV = linspace(YVoxels_min, YVoxels_max,N_YVoxels);

if N_YVoxels == 1
    yVoxels_coord = y_cible;
    for i=1:N_XVoxels
        for j=1:N_ZVoxels
            Voxels(i,j,1) = XV(i);
            Voxels(i,j,2) = yR_coord;
            Voxels(i,j,3) = ZV(j);
        end
    end
elseif N_ZVoxels==1
elseif N_XVoxels==1
% else
%     yVoxels_coord = X_C(:,:,2);
%     for i=1:N_XVoxels
%         for j=1:N_ZVoxels
%             Voxels(i,j,1) = XV(i);
%             Voxels(i,j,2) = yVoxels_coord(i);
%             Voxels(i,j,3) = ZV(j);
%         end
%     end
end
% Creation of a plane surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Distance between Transmitter and Voxels

nTV = size(X_T);
distancesList_TV = zeros(N_XVoxels,N_ZVoxels,nTV(1),nTV(2));
for h = 1:N_XVoxels
    for i = 1:N_ZVoxels 
        for j = 1:nTV(1) 
            for k=1:nTV(2)
            distancesList_TV(h,i,j,k) = sqrt((X_T(j,k,1) - Voxels(h,i,1)).^2 + (X_T(j,k,2) - Voxels(h,i,2)).^2 + (X_T(j,k,3) - Voxels(h,i,3)).^2);
            end
        end
    end
end

%% Distance between Receiver and Voxels
nRV = size(X_R);
distancesList_RV = zeros(N_XVoxels,N_ZVoxels,nRV(1),nRV(2));
for h = 1:N_XVoxels
    for i = 1:N_ZVoxels 
        for j = 1:nRV(1) 
            for k=1:nRV(2)
            distancesList_RV(h,i,j,k) = sqrt((X_R(j,k,1) - Voxels(h,i,1)).^2 + (X_R(j,k,2) - Voxels(h,i,2)).^2 + (X_R(j,k,3) - Voxels(h,i,3)).^2);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Champ Test
S = zeros(Nf,N_Tx*N_Ty,N_Rx*N_Ry);

for l = 1:Nf
    for m = 1:N_Tx*N_Ty %Number of Transmitter
        for n = 1:N_Rx*N_Ry %Number of Receiver
            for k = 1 :N_cible
                %S (frequence,which Transmitter,which receiver)
                %S(l,m,n) = exp(((-2i*pi*fc)/c)*(distancesList_T(m)+distancesList_R(n)));
                S(l,m,n,k) = exp(-2i*pi*(bande_f(l)/c)*(distancesList_T(m,:,k)+distancesList_R(n,:,k))); %bande_f(l)
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Image with Y constante 
% The formula of Image is taken from Microwave Imaging in Security-Two
% Decades in Innovation by Sherif S.Ahmed on page 193
Image_YCTE = zeros(Nf,N_XVoxels,N_ZVoxels);

for k = 1:Nf
    for l = 1:N_XVoxels
        for m=1:N_ZVoxels
            for n = 1:N_Tx*N_Ty %Number of Transmitter
                for o = 1:N_Rx*N_Ry %Number of Receiver
                    for p=1:N_cible
                        Image_YCTE(k,l,m) = Image_YCTE(k,l,m) + S(k,n,o,p)*exp(((2i*pi*bande_f(k))/c)*(distancesList_TV(l,m,n)+distancesList_RV(l,m,o))); %bande_f(l)
                    end
                end
            end
        end
    end
end

TotalImage_YCTE = sum(Image_YCTE,1)/Nf;

%% Define a grid of points for the Antennas, Voxels and Target(s)
figure;
plot3(X_T(:,:,1), X_T(:,:,2), X_T(:,:,3), 'x', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
hold on;
plot3(X_R(:,:,1), X_R(:,:,2), X_R(:,:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
plot3(X_C(:,1), X_C(:,2), X_C(:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
plot3(Voxels(:,:,1), Voxels(:,:,2), Voxels(:,:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
% Add labels and title
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
% legend('Transmitter', 'Receiver', 'Target', 'Location', 'Northeast');
legend('Transmitter', 'Receiver', 'Target', 'Voxels', 'Location', 'Northeast');
view(0,90);

%% Create a scatter plot
% Reshaping TotalImage to match the dimensions expected by surf
TotalImageReshaped_YCTE = reshape(abs(TotalImage_YCTE), [N_XVoxels,N_ZVoxels]);

% Create grids X and Z
[X,Z] = meshgrid(ZV, XV);

% Plotting the 3D image with a heatmap colorbar
figure;
contourf(X, Z, TotalImageReshaped_YCTE, 'EdgeColor', 'none');
xlabel('Z-axis');
ylabel('X-axis');
zlabel('Y-Axis');
title('3D Image Plot of target on  [' + string(X_C(:,1)) + ',' + string(X_C(:,2)) + ',' + string(X_C(:,3)) + ']' );

% Adjusting view angle
view(45, 45); % Change the view angle as needed

% Adding a colorbar
colorbar;

%% Images on Z voxels lineare and dB
% figure;

% subplot(1,2,1);
% plot(ZVoxels(:,:,3) ,20*log10(abs(TotalImageZ)),'o-', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
% plot(Voxels(:,:,1) ,20*log10(abs(TotalImageX)),'o-', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
% xlabel('Z-Coordinates of Voxels');
% ylabel('dB');
% title('Images of target on  [' + string(X_C(:,1)) + ',' + string(X_C(:,2)) + ',' + string(X_C(:,3)) + ']' );
% grid on;
% view(90,-90);
% 
% subplot(1,2,2);
% plot3(ZVoxels(:,:,3),Voxels(:,:,1),abs(TotalImageZ),abs(TotalImageX),'o-', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
% xlabel('Z-Coordinates of Voxels');
% ylabel('Amplitude');
% title('Images of target on  [' + string(X_C(:,1)) + ',' + string(X_C(:,2)) + ',' + string(X_C(:,3)) + ']' );
% grid on;
% view(90,-90);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Stop the timer
elapsed_time = toc;

% Display the elapsed time
disp(['Elapsed Time: ' num2str(elapsed_time) ' seconds']);