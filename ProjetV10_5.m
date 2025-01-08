clear;
close all;
clc; 

f = waitbar(0,'Please wait...');

% Start the timer
tic;

%% Defining all the parameters
z_coord = 0;
%Structure of T or a square or X
T = 0; % 1 if T 2 if X
%Numbers of antennas
N_Rx = 5; %Receiver on x-axis
N_Ry = 1; %Receiver on y-axis
% for T/X structure that encapsulate the receivers,
% value on Y should be the number of transmitter , 
% if want it to be in middle choose 1
%N_Tx can be left as one or zero for middle case
%ex , 6 => 24  or 4*N_Ty transmitter in total
N_Tx = 2;
N_Ty = 1; 

Total_T = N_Ty*N_Tx;
Total_R = N_Ry*N_Rx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Antennas Radiation
% Define angles (in radians)
theta = linspace(-pi/2, pi/2, 10000); % Angles from 0 to 180 degrees

% Define radiation pattern (cosine function for example)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
disp(k);
%Length of the line
% Lx = shannon*N_Rx
Lx = 0;
Ly = 0;
spaceRx = Lx/N_Rx;
spaceRy = Ly/N_Ry;

if Lx == 0 && Ly == 0
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
else
    %Finding the place of the receiving antennas centered on 0
    minRx = (-1)*Lx/2;
    maxRx = Lx/2;
    minRy = (-1)*Ly/2;
    maxRy = Ly/2;
    %Finding the place of the transmitting antennas centered on 0
    if N_Tx == 1
        minTx = 0;
        maxTx = 0;
    else
        minTx = minRx - spaceRx;
        maxTx = maxRx + spaceRx;
        minTy = minRy - spaceRy;
        maxTy = maxRy + spaceRy;
    end 
end
Y_voxels = 0;
l_phy = c/fc;
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
        %X_R(:,i,1) = 0.3; 
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
    Y_voxels = yT_coord;
    for i = 1:N_Ty
        X_T(:,i,1) = linspace(minTx, maxTx,N_Tx); 
        %X_T(:,i,1) = 0.3;
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
%% T shaped structure of receiving antennas
if T==1
    n = round(N_Rx/2);
    temp = X_R(:,:,:);
    if mod(N_Rx,2) == 0
        disp('Even');
        X_R = temp(n,:,:);
        X_R(:,:,1) = X_R(:,:,1)+(shannon/2);
        temp(:,:,2) = temp(:,:,2)+(shannon/2);
        for g=1:N_Rx
            X_R(2,g,:) = temp(g,n,:);
        end
        X_R(:,:,2) = X_R(:,:,2);
    else
        X_R = temp(n,:,:);
        for g=1:N_Rx
            X_R(2,g,:) = temp(g,n,:);
        end
    end
    r = size(X_R);
    Total_R = r(1)*r(2);
%Additional antennas for T shaped structure
    maxRX = max(max(X_R(:,:,1)));
    minRX = min(min(X_R(:,:,1)));
    maxRY = max(max(X_R(:,:,2)));
    minRY = min(min(X_R(:,:,2)));
    newYT = X_R(1,1,1);
    newXT = X_R(2,2,2);
    tempt = zeros(N_Tx,4,3);
    for i =1:N_Tx
        %X-axis
        maxRX = maxRX+shannon/2;
        tempt(i,1,1) = maxRX;
        tempt(i,1,2) = newYT;
        tempt(i,1,3) = z_coord;
        minRX= minRX-shannon/2;
        tempt(i,2,1) = minRX;
        tempt(i,2,2) = newYT;
        tempt(i,2,3) = z_coord;
        %Y-axis
        maxRY = maxRY+shannon/2;
        tempt(i,3,1) = newXT;
        tempt(i,3,2) = maxRY;
        tempt(i,3,3) = z_coord;
        minRY= minRY-shannon/2;
        tempt(i,4,1) = newXT;
        tempt(i,4,2) = minRY;
        tempt(i,4,3) = z_coord; 
    end
    X_T =tempt;
    newT = size(X_T);
    Total_T = newT(2)*newT(1);
end
%% X-shaped structure
if T==2
    n = round(N_Rx/2);
    if mod(N_Rx,2) == 0
        disp('Even');
        for i = 1 :N_Ry
            temp(i,1,:) = X_R(i,i,:);
            temp(i,2,:) = X_R(i,N_Rx-(i-1),:);
        end
    else
        disp('Odd');
        for i = 1 :N_Ry
            if i == n
                temp(i,1,:) = X_R(i,i,:);
            else
                temp(i,1,:)= X_R(i,i,:) ;
                temp(i,2,:) = X_R(N_Rx-(i-1),i,:);
            end
        end
    end
    X_R = temp;
    r = size(X_R);
    Total_R = r(1)*r(2);
    maxRX = max(max(X_R(:,:,1)));
    minRX = min(min(X_R(:,:,1)));
    maxRY = max(max(X_R(:,:,2)));
    minRY = min(min(X_R(:,:,2)));
    if N_Tx ~= 0
        tempt = zeros(N_Ty,4,3);
        for i =1:N_Ty
            maxRX = maxRX+shannon/4;
            minRX= minRX-shannon/4;
            maxRY = maxRY+shannon/4;
            minRY = minRY - shannon/4;
            %Highest points
            tempt(i,1,1) = maxRX;
            tempt(i,1,2) = maxRY;
            tempt(i,1,3) = z_coord;
            tempt(i,2,1) = minRX;
            tempt(i,2,2) = maxRY;
            tempt(i,2,3) = z_coord;
            %Lowest points
            tempt(i,3,1) = maxRX;
            tempt(i,3,2) = minRY;
            tempt(i,3,3) = z_coord;
            tempt(i,4,1) = minRX;
            tempt(i,4,2) = minRY;
            tempt(i,4,3) = z_coord;
        end
    else
        tempt= zeros(1,1,3);
    end
    X_T =tempt;
    newT = size(X_T);
    Total_T = newT(2)*newT(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Square receiver
if T==3
    maxRX = shannon/4+shannon*(N_Tx/2-1);
    minRX= -shannon/4+shannon*(N_Tx/2-1);
    maxRY =shannon/4+shannon*(N_Tx/2-1);
    minRY = -shannon/4+shannon*(N_Tx/2-1);
    ecart = round(N_Tx-2);
    if N_Rx ~= 0
        tempt = zeros(N_Ry,16,3);
        for i =1:N_Ry
            %Highest points
            tempt(i,1,1) = maxRX+shannon*i;
            tempt(i,1,2) = maxRY+shannon*i;
            tempt(i,1,3) = z_coord;
            tempt(i,2,1) = maxRX+shannon*i;
            tempt(i,2,2) = minRY+shannon*i;
            tempt(i,2,3) = z_coord;
            tempt(i,3,1) = minRX+shannon*i;
            tempt(i,3,2) = maxRY+shannon*i;
            tempt(i,3,3) = z_coord;
            tempt(i,4,1) = minRX+shannon*i;
            tempt(i,4,2) = minRY+shannon*i;
            tempt(i,4,3) = z_coord;

            tempt(i,5,1) = minRX - shannon*(i+ecart);
            tempt(i,5,2) = maxRY + shannon*(i);
            tempt(i,5,3) = z_coord;
            tempt(i,6,1) = minRX - shannon*(i+ecart);
            tempt(i,6,2) = minRY+shannon*(i);
            tempt(i,6,3) = z_coord;
            tempt(i,7,1) = maxRX-shannon*(i+ecart);
            tempt(i,7,2) = maxRY+shannon*(i);
            tempt(i,7,3) = z_coord;
            tempt(i,8,1) = maxRX-shannon*(i+ecart);
            tempt(i,8,2) = minRY+shannon*(i);
            tempt(i,8,3) = z_coord;

            %Lowest points
            tempt(i,9,1) = minRX - shannon*(i+ecart);
            tempt(i,9,2) = maxRY - shannon*(i+ecart);
            tempt(i,9,3) = z_coord;
            tempt(i,10,1) = minRX - shannon*(i+ecart);
            tempt(i,10,2) = minRY-shannon*(i+ecart);
            tempt(i,10,3) = z_coord;
            tempt(i,11,1) = maxRX - shannon*(i+ecart);
            tempt(i,11,2) = minRY-shannon*(i+ecart);
            tempt(i,11,3) = z_coord;
            tempt(i,12,1) = maxRX-shannon*(i+ecart);
            tempt(i,12,2) = maxRY-shannon*(i+ecart);
            tempt(i,12,3) = z_coord;
            
            tempt(i,13,1) = maxRX+shannon*(i);
            tempt(i,13,2) = minRY-shannon*(i+ecart);
            tempt(i,13,3) = z_coord;
            tempt(i,14,1) = maxRX+shannon*(i);
            tempt(i,14,2) = maxRY-shannon*(i+ecart);
            tempt(i,14,3) = z_coord;
            tempt(i,15,1) = minRX+shannon*i;
            tempt(i,15,2) = minRY-shannon*(i+ecart);
            tempt(i,15,3) = z_coord;
            tempt(i,16,1) = minRX+shannon*i;
            tempt(i,16,2) = maxRY-shannon*(i+ecart);
            tempt(i,16,3) = z_coord;
        end
    else
        tempt= zeros(1,1,3);
    end
    X_R =tempt;
    newR = size(X_R);
    Total_R = newR(2)*newR(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Placing a Target into the surface z = zc

N_cible = 1; %Number of target
%List of target(s) with its coordinates
X_C(1,1,:) = [0,0,0.5]; 
%X_C(2,1,:) = [0.6,0.6,4.2];
x_cible = X_C(1,1,1);
y_cible = X_C(1,1,2);
z_cible = X_C(1,1,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
waitbar(.2,f,'Loading your data');
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
distancesList_R = reshape(distancesList_R,[Total_R,N_cible]);
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
distancesList_T = reshape(distancesList_T,[Total_T,N_cible]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Paramaters for the voxels

% L_XVoxels = 6;
XVoxels_min = -3;
XVoxels_max = 3;
N_XVoxels = 51;

% L_ZVoxels = 6;
ZVoxels_min = 0;
ZVoxels_max = 5;
N_ZVoxels = 51; 

YVoxels_min =  -3+Y_voxels;
YVoxels_max = 3+Y_voxels;
N_YVoxels = 51; 

XV = linspace(XVoxels_min, XVoxels_max,N_XVoxels);
ZV = linspace(ZVoxels_min, ZVoxels_max,N_ZVoxels);
YV = linspace(YVoxels_min, YVoxels_max,N_YVoxels);

%% Voxels for X-Z plane with Y cte
Voxels1 = zeros(N_XVoxels,N_ZVoxels,3);

yVoxels_coord = y_cible;
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

zVoxels_coord = z_cible;
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

xVoxels_coord = x_cible;
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

distancesList_TV1 = reshape(distancesList_TV1,[N_XVoxels,N_ZVoxels,Total_T]);

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

distancesList_RV1 = reshape(distancesList_RV1,[N_XVoxels,N_ZVoxels,Total_R]);
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

distancesList_TV2 = reshape(distancesList_TV2,[N_XVoxels,N_YVoxels,Total_T]);

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

distancesList_RV2 = reshape(distancesList_RV2,[N_XVoxels,N_YVoxels,Total_R]);
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
distancesList_TV3 = reshape(distancesList_TV3,[N_YVoxels,N_ZVoxels,Total_T]);
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
distancesList_RV3 = reshape(distancesList_RV3,[N_YVoxels,N_ZVoxels,Total_R]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
waitbar(.55,f,'Calculating Test field');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Champ Test
S = zeros(Nf,Total_T,Total_R,N_cible);
radiation = 3;
fctrue =true;
minr = abs(min(distancesList_R(:,:,1)));

for l = 1:Nf
    for m = 1:Total_T%Number of Transmitter
        for n = 1:Total_R %Number of Receiver
            k = -2i*pi*bande_f(l)/c;
            radiation_pattern = cos(theta).^3;
            radiation_pattern = abs(radiation_pattern/max(radiation_pattern));
            if (bande_f(l) == fc && fctrue == true)
                        fctrue = -fctrue;
                        figure;
                        subplot(2,1,1);
                        plot(rad2deg(theta), (radiation_pattern));
                        xlabel('Angle (degrees)');
                        ylabel('Radiation ');
                        title('Antenna Radiation Pattern Central Frequency');
                        subplot(2,1,2);
                        plot(rad2deg(theta), 10*log10(radiation_pattern));
                        xlabel('Angle (degrees)');
                        ylabel('Radiation ');
                        ylim([-60 0]);
                        title('Antenna Radiation Pattern(dB)');
                        grid on;  % Enable gridlines

                        % Add -3dB reference line
                        hold on;
                        ref_level = -3;  % -3dB reference level
                        plot([rad2deg(min(theta)), rad2deg(max(theta))], [ref_level, ref_level], 'r--', 'LineWidth', 1.5);
                        text(200, -5, '-3dB', 'Color', 'red', 'FontSize', 12);  % Add label
                        hold off;
            end
            g = acos(abs(X_R(n,1,1)-x_cible)/distancesList_R(n,:));
            dg = rad2deg(g(1));
            rad = interp1(rad2deg(theta), radiation_pattern, dg, 'linear');
            
            lambda =c/bande_f(l);
            
            r1 = distancesList_T(m,:);
            r2 = distancesList_R(n,:);
            r = r1+r2;
            newX_R = reshape(X_R,[Total_R,1,3]);
            newX_T = reshape(X_T,[Total_T,1,3]);

            %Transmitter
            angT = acos(abs(newX_T(m,1,1)-x_cible)/distancesList_T(m,:));
            Tangles(m,:) = 90-rad2deg(angT(1));
            radT = interp1(rad2deg(theta), radiation_pattern, Tangles(m,:), 'nearest', 'extrap');

            %Receiver
            angR = acos(abs(newX_R(n,1,1)-x_cible)/distancesList_R(n,:));
            Rangles(n,:) = 90-rad2deg(angR(1));
            radR = interp1(rad2deg(theta), radiation_pattern, Rangles(n,:),  'nearest', 'extrap');
            %S (frequence,which Transmitter,which receiver)
            if radiation == 0
                S(l,m,n,:) = exp(k*r);
            elseif radiation== 1
                S(l,m,n,:) = (lambda/sqrt(4*pi))*exp(k*r)/(4*pi*r1*r2);
            elseif radiation== 2
                S(l,m,n,:) = abs(radR*radT)*exp(k*r);
            else
                S(l,m,n,:) = (lambda/sqrt(4*pi))*exp(k*r)*abs(rad)/(4*pi*r1*r2);
            end
            % S(l,m,n,:) = exp(phase*r)/(4*pi*r);
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
    for n = 1:Total_T % Number of Transmitter
        for o = 1:Total_R % Number of Receiver
            Image_YCTE(k, :, :) = squeeze(Image_YCTE(k, :, :)) + S(k, n, o, :) *exp(((2i * pi * bande_f(k)) / c) * (distancesList_TV1(:, :, n) + distancesList_RV1(:, :, o))); % bande_f(l)
        end
    end
end

TotalImage_YCTE = sum(Image_YCTE,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Image with Z constante 
% The formula of Image is taken from Microwave Imaging in Security-Two
% Decades in Innovation by Sherif S.Ahmed on page 193
Image_ZCTE = zeros(Nf,N_XVoxels,N_YVoxels);

for k = 1:Nf
    for n = 1:Total_T %Number of Transmitter
        for o = 1:Total_R %Number of Receiver)
            Image_ZCTE(k, :, :) = squeeze(Image_ZCTE(k, :, :))+ S(k, n, o, :) * exp(((2i * pi * bande_f(k)) / c) * (distancesList_TV2(:, :, n) + distancesList_RV2(:, :, o))); % bande_f(l)
        end
    end
end

TotalImage_ZCTE = sum(Image_ZCTE,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Image with X constante 
% The formula of Image is taken from Microwave Imaging in Security-Two
% Decades in Innovation by Sherif S.Ahmed on page 193
Image_XCTE = zeros(Nf,N_YVoxels,N_ZVoxels);

for k = 1:Nf
    for n = 1:Total_T %Number of Transmitter
        for o = 1:Total_R %Number of Receiver
            Image_XCTE(k, :, :) = squeeze(Image_XCTE(k, :, :))+ S(k, n, o, :) * exp(((2i * pi * bande_f(k)) / c) * (distancesList_TV3(:, :, n) + distancesList_RV3(:, :, o)));
        end
    end
end
TotalImage_XCTE = sum(Image_XCTE,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
waitbar(.65,f,'Image construction');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define a grid of points for the Antennas, Voxels and Target(s
%For Y cte
figure;
%subplot(1,3,1);
plot3(X_C(:,1), X_C(:,2), X_C(:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
hold on;
plot3(X_T(:,:,1), X_T(:,:,2), X_T(:,:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

plot3(X_R(:,:,1), X_R(:,:,2), X_R(:,:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

%plot3(Voxels1(:,:,1), Voxels1(:,:,2), Voxels1(:,:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
% Add labels and title
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
legend('Target','Transmitter', 'Receiver',  'Location', 'Northeast');
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

% Save data to MAT file
if radiation == 0
    save('data.mat', 'TotalImageReshaped_YCTE', 'TotalImageReshaped_ZCTE', 'TotalImageReshaped_XCTE','BigMax');
else
    save('data_radiation.mat', 'TotalImageReshaped_YCTE', 'TotalImageReshaped_ZCTE', 'TotalImageReshaped_XCTE','BigMax');
end
%% Create a scatter plot Y ct
% Create grids X and Z
[X1,Z1] = meshgrid(ZV, XV);

TotalImageReshaped_YCTE_Normalize =TotalImageReshaped_YCTE/BigMax;
% Find the maximum value and its indices
[maxValueY, maxIndexY] = max(TotalImageReshaped_YCTE_Normalize(:));
[rowY, colY] = ind2sub(size(TotalImageReshaped_YCTE_Normalize), maxIndexY);

% Plotting the 3D image with a heatmap colorbar
figure;
subplot(1,2,1);
surf(X1, Z1,TotalImageReshaped_YCTE_Normalize, 'EdgeColor', 'none');
shading("interp");
xlabel('Z-axis');
ylabel('X-axis');
title('3D Image Plot with Y cte=' + string(yVoxels_coord) +' of target on  [' + string(X_C(:,1)) + ',' + string(X_C(:,2)) + ',' + string(X_C(:,3)) + ']' );

% Adjusting view angle
view(90, -90); % Change the view angle as needed

% Adding a colorbar
colorbar;
clim([0 1]);

dB = 20*log10(TotalImageReshaped_YCTE_Normalize);

subplot(1,2,2);
surf(X1, Z1, dB);
shading("interp");
xlabel('Z-axis');
ylabel('X-axis');
title('3D Image Plot with Y cte=' + string(yVoxels_coord) +' of target on  [' + string(X_C(:,1)) + ',' + string(X_C(:,2)) + ',' + string(X_C(:,3)) + ']' );
colorbar;
clim([-30 0]);
view(90, -90);
% %% Create a scatter plot Z cte
% % Create grids X and Z
% [X2,Y2] = meshgrid(YV, XV);
% TotalImageReshaped_ZCTE_Normalize =TotalImageReshaped_ZCTE/BigMax;
% 
% [maxValueZ, maxIndexZ] = max(TotalImageReshaped_ZCTE_Normalize(:));
% [rowZ,colZ] = ind2sub(size(TotalImageReshaped_ZCTE_Normalize), maxIndexZ);
% 
% % Plotting the 3D image with a heatmap colorbar
% figure;
% subplot(1,2,1);
% surf(X2, Y2, TotalImageReshaped_ZCTE_Normalize, 'EdgeColor', 'none');
% shading("interp");
% xlabel('Y-axis');
% ylabel('X-axis');
% title('3D Image Plot with Z cte = '+ string(zVoxels_coord) + ' of target on  [' + string(X_C(:,1)) + ',' + string(X_C(:,2)) + ',' + string(X_C(:,3)) + ']' );
% 
% % Adjusting view angle
% view(90, -90); % Change the view angle as needed
% 
% % Adding a colorbar
% colorbar;
% clim([0 1]);
% 
% subplot(1,2,2);
% surf(X2, Y2, 20*log10(TotalImageReshaped_ZCTE_Normalize), 'EdgeColor', 'none');
% shading("interp");
% xlabel('Y-axis');
% ylabel('X-axis');
% title('3D Image Plot with Z cte = '+ string(zVoxels_coord) + ' of target on  [' + string(X_C(:,1)) + ',' + string(X_C(:,2)) + ',' + string(X_C(:,3)) + ']' );
% 
% colorbar;
% clim([-30 0]);
% view(90, -90);
% %% Create a scatter plot X cte
% % Create grids Y and Z
% [Y3,Z3] = meshgrid(ZV, YV);
% 
% TotalImageReshaped_XCTE_Normalize =TotalImageReshaped_XCTE/BigMax;
% 
% [maxValueX, maxIndexX] = max(TotalImageReshaped_XCTE_Normalize(:));
% [rowX, colX] = ind2sub(size(TotalImageReshaped_XCTE_Normalize), maxIndexX);
% 
% % Plotting the 3D image with a heatmap colorbar
% figure;
% subplot(1,2,1);
% surf(Y3, Z3, TotalImageReshaped_XCTE_Normalize , 'EdgeColor', 'none');
% shading("interp");
% xlabel('Z-axis');
% ylabel('Y-axis');
% title('3D Image Plot with X cte = ' + string(xVoxels_coord) + ' of target on  [' + string(X_C(:,1)) + ',' + string(X_C(:,2)) + ',' + string(X_C(:,3)) + ']' );
% 
% % Adjusting view angle
% view(90, -90); % Change the view angle as needed
% 
% % Adding a colorbar
% colorbar;
% clim([0 1]);
% 
% subplot(1,2,2);
% surf(Y3, Z3, 20*log10(TotalImageReshaped_XCTE_Normalize ), 'EdgeColor', 'none');
% shading("interp");
% xlabel('Z-axis');
% ylabel('Y-axis');
% title('3D Image Plot with X cte = ' + string(xVoxels_coord) + ' of target on  [' + string(X_C(:,1)) + ',' + string(X_C(:,2)) + ',' + string(X_C(:,3)) + ']' );
% 
% colorbar;
% clim([-30 0]);
% view(90, -90);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
waitbar(1,f,'Finishing');
%% Stop the timer
elapsed_time = toc;
close(f);
% Display the elapsed time
disp(['Elapsed Time: ' num2str(elapsed_time) ' seconds']);