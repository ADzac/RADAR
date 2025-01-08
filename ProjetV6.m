clear;
close all;
clc; 

f = waitbar(0,'Please wait...');

% Start the timer
tic;

%% Defining all the parameters
z_coord = 0;

%Numbers of antennas
N_Rx = 11; %Receiver on x-axis
N_Ry = 2; %Receiver on y-axis
N_Tx = 1;
N_Ty = 1;

Total_T = N_Ty*N_Tx;
Total_R = N_Ry*N_Rx;

%Frequency & Spacing
Nf = 101;  %always odd(impair)
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
    yT_coord = 0;
    for i = 1:N_Ty
        X_T(:,i,1) = linspace(minTx, maxTx,N_Tx); 
        X_T(:,i,2) = yT_coord; 
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
x_cible = X_C(1,1,1);
y_cible = X_C(1,1,2);
z_cible = X_C(1,1,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
waitbar(.2,f,'Loading your data');
%% Distances for each Receiver with target
nR = size(X_R);
distancesList_R = zeros(nR(1)*nR(2),N_cible);
for h = 1:N_cible
    for i = 1:nR(2)
        for j = 1:nR(1)
            distancesList_R(j,i,h) = sqrt((X_R(j,i,1) - X_C(h,1)).^2 + (X_R(j,i,2) - X_C(h,2)).^2 + (X_R(j,i,3) - X_C(h,3)).^2);
        end
    end
end
distancesList_R = reshape(distancesList_R,[Total_R,N_cible]);
% distances for each Receiver with target
nT = size(X_T);
distancesList_T = zeros(nT(1)*nT(2),N_cible); % based on the number of transmitter on X and Y and number of targets
for h = 1:N_cible
    for i = 1:nT(2)
        for j = 1:nT(1)
            distancesList_T(j,i,h) = sqrt((X_T(j,i,1) - X_C(h,1)).^2 + (X_T(j,i,2) - X_C(h,2)).^2 + (X_T(j,i,3) - X_C(h,3)).^2);
        end
    end
end
distancesList_T = reshape(distancesList_T,[Total_T,N_cible]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Paramaters for the voxels

% L_XVoxels = 6;
XVoxels_min = -3;
XVoxels_max = 3;
N_XVoxels = 51;

% L_ZVoxels = 6;
ZVoxels_min = 2;
ZVoxels_max = 4;
N_ZVoxels = 51; 

YVoxels_min =  -3;
YVoxels_max = 3;
N_YVoxels = 51; 

XV = linspace(XVoxels_min, XVoxels_max,N_XVoxels);
ZV = linspace(ZVoxels_min, ZVoxels_max,N_ZVoxels);
YV = linspace(YVoxels_min, YVoxels_max,N_YVoxels);

%% Voxels for X-Z plane with Y cte
Voxels1 = zeros(N_XVoxels,N_ZVoxels,3);

yVoxels_coord = 0;
for i=1:N_XVoxels
    for j=1:N_ZVoxels
        Voxels1(i,j,1) = XV(i);
        Voxels1(i,j,2) = yVoxels_coord;
        Voxels1(i,j,3) = ZV(j);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Voxels for X-Y plane with Z cte

Voxels2 = zeros(N_XVoxels,N_YVoxels,3);

zVoxels_coord = 3;
for i=1:N_XVoxels
    for j=1:N_YVoxels
        Voxels2(i,j,1) = XV(i);
        Voxels2(i,j,2) = YV(j);
        Voxels2(i,j,3) = zVoxels_coord;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Voxels for Y-Z plane with X cte

Voxels3 = zeros(N_YVoxels,N_ZVoxels,3);

xVoxels_coord = 0;
for i=1:N_YVoxels
    for j=1:N_ZVoxels
        Voxels3(i,j,1) = xVoxels_coord;
        Voxels3(i,j,2) = YV(i);
        Voxels3(i,j,3) = ZV(j);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Distance between Transmitter and Voxels Ycte

nTV = size(X_T);
distancesList_TV1 = zeros(N_XVoxels,N_ZVoxels,nTV(1),nTV(2));
for h = 1:N_XVoxels
    for i = 1:N_ZVoxels 
        for j = 1:nTV(1) 
            for k=1:nTV(2)
            distancesList_TV1(h,i,j,k) = sqrt((X_T(j,k,1) - Voxels1(h,i,1)).^2 + (X_T(j,k,2) - Voxels1(h,i,2)).^2 + (X_T(j,k,3) - Voxels1(h,i,3)).^2);
            end
        end
    end
end

waitbar(.4,f,'Loading your data');
%% Distance between Receiver and Voxels Zcte
nRV = size(X_R);
distancesList_RV1 = zeros(N_XVoxels,N_ZVoxels,nRV(1),nRV(2));
for h = 1:N_XVoxels
    for i = 1:N_ZVoxels 
        for j = 1:nRV(1) 
            for k=1:nRV(2)
            distancesList_RV1(h,i,j,k) = sqrt((X_R(j,k,1) - Voxels1(h,i,1)).^2 + (X_R(j,k,2) - Voxels1(h,i,2)).^2 + (X_R(j,k,3) - Voxels1(h,i,3)).^2);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Distance between Transmitter and Voxels Zcte

nTV = size(X_T);
distancesList_TV2 = zeros(N_XVoxels,N_YVoxels,nTV(1),nTV(2));
for h = 1:N_XVoxels
    for i = 1:N_YVoxels 
        for j = 1:nTV(1) 
            for k=1:nTV(2)
            distancesList_TV2(h,i,j,k) = sqrt((X_T(j,k,1) - Voxels2(h,i,1)).^2 + (X_T(j,k,2) - Voxels2(h,i,2)).^2 + (X_T(j,k,3) - Voxels2(h,i,3)).^2);
            end
        end
    end
end

%% Distance between Receiver and Voxels
nRV = size(X_R);
distancesList_RV2 = zeros(N_XVoxels,N_YVoxels,nRV(1),nRV(2));
for h = 1:N_XVoxels
    for i = 1:N_YVoxels 
        for j = 1:nRV(1) 
            for k=1:nRV(2)
            distancesList_RV2(h,i,j,k) = sqrt((X_R(j,k,1) - Voxels2(h,i,1)).^2 + (X_R(j,k,2) - Voxels2(h,i,2)).^2 + (X_R(j,k,3) - Voxels2(h,i,3)).^2);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Distance between Transmitter and Voxels X cte

nTV = size(X_T);
distancesList_TV3 = zeros(N_YVoxels,N_ZVoxels,nTV(1),nTV(2));
for h = 1:N_YVoxels
    for i = 1:N_ZVoxels 
        for j = 1:nTV(1) 
            for k=1:nTV(2)
            distancesList_TV3(h,i,j,k) = sqrt((X_T(j,k,1) - Voxels3(h,i,1)).^2 + (X_T(j,k,2) - Voxels3(h,i,2)).^2 + (X_T(j,k,3) - Voxels3(h,i,3)).^2);
            end
        end
    end
end

%% Distance between Receiver and Voxels
nRV = size(X_R);
distancesList_RV3 = zeros(N_YVoxels,N_ZVoxels,nRV(1),nRV(2));
for h = 1:N_YVoxels
    for i = 1:N_ZVoxels 
        for j = 1:nRV(1) 
            for k=1:nRV(2)
            distancesList_RV3(h,i,j,k) = sqrt((X_R(j,k,1) - Voxels3(h,i,1)).^2 + (X_R(j,k,2) - Voxels3(h,i,2)).^2 + (X_R(j,k,3) - Voxels3(h,i,3)).^2);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
waitbar(.55,f,'Calculating Test field');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Champ Test
S = zeros(Nf,N_Tx*N_Ty,N_Rx*N_Ry,N_cible);

for l = 1:Nf
    for m = 1:N_Tx*N_Ty %Number of Transmitter
        for n = 1:N_Rx*N_Ry %Number of Receiver
                %S (frequence,which Transmitter,which receiver)
                %S(l,m,n) = exp(((-2i*pi*fc)/c)*(distancesList_T(m)+distancesList_R(n)));
                S(l,m,n,:) = exp(-2i*pi*(bande_f(l)/c)*(distancesList_T(m,:)+distancesList_R(n,:))); %bande_f(l)
        end
    end
end

S = sum(S,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
waitbar(.6,f,'Calculating Image');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Image with Y constante 
% The formula of Image is taken from Microwave Imaging in Security-Two
% Decades in Innovation by Sherif S.Ahmed on page 193
Image_YCTE = zeros(Nf,N_XVoxels,N_ZVoxels);

for k = 1:Nf
    for n = 1:N_Tx * N_Ty % Number of Transmitter
        for o = 1:N_Rx * N_Ry % Number of Receiver
            Image_YCTE(k, :, :) = squeeze(Image_YCTE(k, :, :)) + S(k, n, o, :) *exp(((2i * pi * bande_f(k)) / c) * (distancesList_TV1(:, :, n) + distancesList_RV1(:, :, o))); % bande_f(l)
        end
    end
end

TotalImage_YCTE = sum(Image_YCTE,1)/Nf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Image with Z constante 
% The formula of Image is taken from Microwave Imaging in Security-Two
% Decades in Innovation by Sherif S.Ahmed on page 193
Image_ZCTE = zeros(Nf,N_XVoxels,N_YVoxels);

for k = 1:Nf
    for n = 1:N_Tx*N_Ty %Number of Transmitter
        for o = 1:N_Rx*N_Ry %Number of Receiver)
            Image_ZCTE(k, :, :) = squeeze(Image_ZCTE(k, :, :))+ S(k, n, o, :) * exp(((2i * pi * bande_f(k)) / c) * (distancesList_TV2(:, :, n) + distancesList_RV2(:, :, o))); % bande_f(l)
        end
    end
end

TotalImage_ZCTE = sum(Image_ZCTE,1)/Nf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Image with X constante 
% The formula of Image is taken from Microwave Imaging in Security-Two
% Decades in Innovation by Sherif S.Ahmed on page 193
Image_XCTE = zeros(Nf,N_YVoxels,N_ZVoxels);

for k = 1:Nf
    for n = 1:N_Tx*N_Ty %Number of Transmitter
        for o = 1:N_Rx*N_Ry %Number of Receiver
            Image_XCTE(k, :, :) = squeeze(Image_XCTE(k, :, :))+ S(k, n, o, :) * exp(((2i * pi * bande_f(k)) / c) * (distancesList_TV3(:, :, n) + distancesList_RV3(:, :, o)));
        end
    end
end
TotalImage_XCTE = sum(Image_XCTE,1)/Nf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
waitbar(.65,f,'Image construction');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define a grid of points for the Antennas, Voxels and Target(s
%For Y cte
figure;
%subplot(1,3,1);
plot3(X_T(:,:,1), X_T(:,:,2), X_T(:,:,3), 'x', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
hold on;
plot3(X_R(:,:,1), X_R(:,:,2), X_R(:,:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
plot3(X_C(:,1), X_C(:,2), X_C(:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
%plot3(Voxels1(:,:,1), Voxels1(:,:,2), Voxels1(:,:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
% Add labels and title
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
legend('Transmitter', 'Receiver', 'Target', 'Location', 'Northeast');
%legend('Transmitter', 'Receiver', 'Target', 'Voxels', 'Location', 'Northeast');
view(0,90);
% 
% %For Z cte
% subplot(1,3,2);
% plot3(X_T(:,:,1), X_T(:,:,2), X_T(:,:,3), 'x', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
% hold on;
% plot3(X_R(:,:,1), X_R(:,:,2), X_R(:,:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
% plot3(X_C(:,1), X_C(:,2), X_C(:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
% plot3(Voxels2(:,:,1), Voxels2(:,:,2), Voxels2(:,:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
% % Add labels and title
% xlabel('X-axis');
% ylabel('Y-axis');
% zlabel('Z-axis');
% % legend('Transmitter', 'Receiver', 'Target', 'Location', 'Northeast');
% legend('Transmitter', 'Receiver', 'Target', 'Voxels', 'Location', 'Northeast');
% view(0,90);
% 
% %For X cte
% subplot(1,3,3);
% plot3(X_T(:,:,1), X_T(:,:,2), X_T(:,:,3), 'x', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
% hold on;
% plot3(X_R(:,:,1), X_R(:,:,2), X_R(:,:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
% plot3(X_C(:,1), X_C(:,2), X_C(:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
% plot3(Voxels3(:,:,1), Voxels3(:,:,2), Voxels3(:,:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
% % Add labels and title
% xlabel('X-axis');
% ylabel('Y-axis');
% zlabel('Z-axis');
% % legend('Transmitter', 'Receiver', 'Target', 'Location', 'Northeast');
% legend('Transmitter', 'Receiver', 'Target', 'Voxels', 'Location', 'Northeast');
% view(0,90);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
waitbar(.8,f,'Plotting ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reshaping TotalImage to match the dimensions expected by surf
TotalImageReshaped_YCTE = reshape(abs(TotalImage_YCTE), [N_XVoxels,N_ZVoxels]);
max1 = max(max(TotalImageReshaped_YCTE));

TotalImageReshaped_ZCTE = reshape(abs(TotalImage_ZCTE), [N_XVoxels,N_YVoxels]);
max2 = max(max(TotalImageReshaped_ZCTE));

TotalImageReshaped_XCTE = reshape(abs(TotalImage_XCTE), [N_YVoxels,N_ZVoxels]);
max3 =max(max(TotalImageReshaped_XCTE));

Totalmax = [max1,max2,max3];
BigMax =max(Totalmax);
%% Create a scatter plot Y ct
% Create grids X and Z
[X1,Z1] = meshgrid(ZV, XV);

% Plotting the 3D image with a heatmap colorbar
figure;
subplot(1,2,1);
contourf(X1, Z1, TotalImageReshaped_YCTE/BigMax, 'EdgeColor', 'none');
xlabel('Z-axis');
ylabel('X-axis');
title('3D Image Plot with Y cte=' + string(yVoxels_coord) +' of target on  [' + string(X_C(:,1)) + ',' + string(X_C(:,2)) + ',' + string(X_C(:,3)) + ']' );

% Adjusting view angle
view(90, -90); % Change the view angle as needed

% Adding a colorbar
colorbar;
caxis([0 1]);

subplot(1,2,2);
contourf(X1, Z1, 20*log10(TotalImageReshaped_YCTE/BigMax), 'EdgeColor', 'none');
xlabel('Z-axis');
ylabel('X-axis');
colorbar;
caxis([-30 0]);
view(90, -90);
%% Create a scatter plot Z cte
% Create grids X and Z
[X2,Y2] = meshgrid(YV, XV);

% Plotting the 3D image with a heatmap colorbar
figure;
subplot(1,2,1);
contourf(X2, Y2, TotalImageReshaped_ZCTE/BigMax, 'EdgeColor', 'none');
xlabel('Y-axis');
ylabel('X-axis');
title('3D Image Plot with Z cte = '+ string(zVoxels_coord) + ' of target on  [' + string(X_C(:,1)) + ',' + string(X_C(:,2)) + ',' + string(X_C(:,3)) + ']' );

% Adjusting view angle
view(90, -90); % Change the view angle as needed

% Adding a colorbar
colorbar;
caxis([0 1]);

subplot(1,2,2);
contourf(X2, Y2, 20*log10(TotalImageReshaped_ZCTE/BigMax), 'EdgeColor', 'none');
xlabel('Y-axis');
ylabel('X-axis');
colorbar;
caxis([-30 0]);
view(90, -90);
%% Create a scatter plot X cte
% Create grids X and Z
[Y3,Z3] = meshgrid(ZV, YV);

% Plotting the 3D image with a heatmap colorbar
figure;
subplot(1,2,1);
contourf(Y3, Z3, TotalImageReshaped_XCTE/BigMax, 'EdgeColor', 'none');
xlabel('Z-axis');
ylabel('Y-axis');
title('3D Image Plot with X cte = ' + string(xVoxels_coord) + ' of target on  [' + string(X_C(:,1)) + ',' + string(X_C(:,2)) + ',' + string(X_C(:,3)) + ']' );

% Adjusting view angle
view(90, -90); % Change the view angle as needed

% Adding a colorbar
colorbar;
caxis([0 1]);

subplot(1,2,2);
contourf(Y3, Z3, 20*log10(TotalImageReshaped_XCTE/BigMax), 'EdgeColor', 'none');
xlabel('Z-axis');
ylabel('Y-axis');
colorbar;
caxis([-30 0]);
view(90, -90);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
waitbar(1,f,'Finishing');
%% Stop the timer
elapsed_time = toc;
close(f);
% Display the elapsed time
disp(['Elapsed Time: ' num2str(elapsed_time) ' seconds']);