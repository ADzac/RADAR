function ProjetV10IHM(handles,event,radar,ax,ax_XZ,ax_XZ1,ax_YZ,ax_YZ1,ax_XY,ax_XY1)

hradar = guihandles(radar);
a = str2double(get(hradar.N_TX,'String'));
b = str2double(get(hradar.N_TY,'String'));
c = str2double(get(hradar.N_RX,'String'));
d = str2double(get(hradar.N_RY,'String'));
e = str2double(get(hradar.N_cible,'String'));
% Get the selected index from the popup menu
selected_index = get(hradar.Choice, 'Value');

% Get the list of strings for the popup menu
popup_strings = get(hradar.Choice, 'String');

% Get the selected string from the popup menu based on the index
selected_string = popup_strings{selected_index};

% Get the selected index from the popup menu
selected_ind = get(hradar.Choix, 'Value');

% Get the list of strings for the popup menu
popup_str = get(hradar.Choix, 'String');

% Get the selected string from the popup menu based on the index
selected_str = popup_str{selected_ind};

X_cible1 = str2double(get(hradar.X_cible1,'String'));
Y_cible1 = str2double(get(hradar.Y_cible1,'String'));
Z_cible1 = str2double(get(hradar.Z_cible1,'String'));

X_cible2 = str2double(get(hradar.X_cible2,'String'));
Y_cible2 = str2double(get(hradar.Y_cible2,'String'));
Z_cible2 = str2double(get(hradar.Z_cible2,'String'));

X_cible3 = str2double(get(hradar.X_cible3,'String'));
Y_cible3 = str2double(get(hradar.Y_cible3,'String'));
Z_cible3 = str2double(get(hradar.Z_cible3,'String'));

%Plane cuts
X_cut = str2double(get(hradar.X_cut,'String'));
Y_cut = str2double(get(hradar.Y_cut,'String'));
Z_cut = str2double(get(hradar.Z_cut,'String'));

f = str2double(get(hradar.LX,'String'));
g = str2double(get(hradar.LY,'String'));

nf = str2double(get(hradar.Nf,'String'));

minf = str2double(get(hradar.fmin,'String'));
maxf = str2double(get(hradar.fmax,'String'));

minXV = str2double(get(hradar.LXVmin,'String'));
maxXV = str2double(get(hradar.LXVmax,'String'));

minYV = str2double(get(hradar.LYVmin,'String'));
maxYV = str2double(get(hradar.LYVmax,'String'));

minZV = str2double(get(hradar.LZVmin,'String'));
maxZV = str2double(get(hradar.LZVmax,'String'));

%% Error exception

if (~isnan(X_cible2) || ~isnan(Y_cible2) || ~isnan(Z_cible2)) && e<2
    % Display a warning message
    warningMessage = 'Error: Values detected for 2nd target! Please leave it blank or change the target number to 2';
    warndlg(warningMessage, 'Warning');
    % Stop script execution using the error function
    error('Error: Invalid values detected! Script halted.');
end

if (isnan(X_cible2) || isnan(Y_cible2) || isnan(Z_cible2)) && e==2
    % Display a warning message
    warningMessage = 'Error: Please make sure to enter values for each coordinates';
    warndlg(warningMessage, 'Warning');
    % Stop script execution using the error function
    error('Error: Invalid values detected! Script halted.');
end

if (~isnan(X_cible3) || ~isnan(Y_cible3) || ~isnan(Z_cible3)) && e<3
    % Display a warning message
    warningMessage = 'Error: Values detected for 3rd target! Please leave it blank or change the target number to 2';
    warndlg(warningMessage, 'Warning');
    % Stop script execution using the error function
    error('Error: Invalid values detected! Script halted.');
end

if (isnan(X_cible3) || isnan(Y_cible3) || isnan(Z_cible3)) && e==3
    % Display a warning message
    warningMessage = 'Error: Please make sure to enter values for each coordinates';
    warndlg(warningMessage, 'Warning');
    % Stop script execution using the error function
    error('Error: Invalid values detected! Script halted.');
end

if (isnan(f) && ~isnan(g) || ~isnan(f) && isnan(g))
    % Display a warning message
    warningMessage = 'Error: Please make sure to enter values for both LX and LY';
    warndlg(warningMessage, 'Warning');
    % Stop script execution using the error function
    error('Error: Invalid values detected! Script halted.');
end

if ((isnan(minf) && ~isnan(maxf)) || (~isnan(minf) && isnan(maxf)))
    % Display a warning message
    warningMessage = 'Error: Please make sure to enter values for both fmin and fmax';
    warndlg(warningMessage, 'Warning');
    % Stop script execution using the error function
    error('Error: Invalid values detected! Script halted.');
end
%--------------------------------------------------------------------------
% Voxels
if ((isnan(minXV) && ~isnan(maxXV)) || (~isnan(minXV) && isnan(maxXV)))
    % Display a warning message
    warningMessage = 'Error: Please make sure to enter values for both xmin and xmax';
    warndlg(warningMessage, 'Warning');
    % Stop script execution using the error function
    error('Error: Invalid values detected! Script halted.');
end

if ((isnan(minYV) && ~isnan(maxYV)) || (~isnan(minYV) && isnan(maxYV)))
    % Display a warning message
    warningMessage = 'Error: Please make sure to enter values for both ymin and ymax';
    warndlg(warningMessage, 'Warning');
    % Stop script execution using the error function
    error('Error: Invalid values detected! Script halted.');
end

if ((isnan(minZV) && ~isnan(maxZV)) || (~isnan(minZV) && isnan(maxZV)))
    % Display a warning message
    warningMessage = 'Error: Please make sure to enter values for both zmin and zmax';
    warndlg(warningMessage, 'Warning');
    % Stop script execution using the error function
    error('Error: Invalid values detected! Script halted.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~isnan(f) && ~isnan(g)) && (f~=0 && g~=0)
    L = 1;
else
    L= 0;
end
% Calculate
f = waitbar(0,'Please wait...');

% Start the timer
tic;

%% Defining all the parameters
z_coord = 0;
%Structure of T or a square or X
T = selected_string;
Change = selected_str;
%Numbers of antennas
% for square with outer receivers square, change N_Ry for numbers of square
% at each edge, leave N_Rx on one
N_Rx = c; %Receiver on x-axis
N_Ry = d; %Receiver on y-axis
% for T/X structure that encapsulate the receivers,
% value on Y should be the number of transmitter , 
% if want it to be in middle choose 1
%N_Tx can be left as one or zero for middle case
%ex , 6 => 24  or 4*N_Ty transmitter in total
%for square structure, choose number of square at the outer
% each square has 4 transmitters
N_Tx = a;
N_Ty =b; 

Total_T = N_Ty*N_Tx;
Total_R = N_Ry*N_Rx;

%Frequency & Spacing
if ~isnan(nf)
    Nf = nf;  %always odd(impair)
else
    Nf = 101;
end

c = 3e8;
if isnan(minf) || (minf ==0)
    f_min = 0.5e9;
else
    f_min = minf;
end
if isnan(maxf) || (maxf==0)
    f_max = 1.5e9;
else
    f_max = maxf;
end

bande_f = linspace(f_min,f_max,Nf); %If more than 1 frequency
fc = 1e9; %(f_max+f_min)/2;
lc = c/fc;
shannon = lc/2;
k  = 2*pi/lc;

%% Length of the line
% Lx = shannon*N_Rx
Lx = 1.5;
Ly = 1;
space_Rx = Lx/N_Rx;
space_Ry = Ly/N_Ry;

if L == 0
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
end

if L==1
    %Finding the place of the receiving antennas centered on 0
    minRx = -Lx/2;
    maxRx = Lx/2;
    minRy = -Ly/2;
    maxRy = Ly/2;
    if N_Rx == 1
        minRx = 0;
        maxRx = 0;
    end
    %Finding the place of the transmitting antennas centered on 0
    minTx = -Lx/2 -Lx/N_Rx;
    maxTx = Lx/2 + Lx/N_Rx;
    minTy = -Ly/2;
    maxTy = Ly/2;
    if N_Tx == 1
        minTx = 0;
        maxTx = 0;
    end
end

Y_voxels = 0;
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
if T=="T-shaped"
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
if Change == "Normal"
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
else
    maxRX = max(max(X_R(:,:,1)));
    minRX = min(min(X_R(:,:,1)));
    maxRY = 0;
    minRY = 0;
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
    end
end
    X_T =tempt;
    newT = size(X_T);
    Total_T = newT(2)*newT(1);
end
%% X-shaped structure
if T=="X-shaped"
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

%%  Square receiver
if T=="Square X square"
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Placing a Target into the surface z = zc

N_cible = e; %Number of target
%List of target(s) with its coordinates
X_C(1,1,:) = [X_cible1,Y_cible1,Z_cible1]; 
if e == 2
    X_C(1,2,:) = [X_cible2,Y_cible2,Z_cible2];
end
if e == 3
    X_C(1,2,:) = [X_cible2,Y_cible2,Z_cible2];
    X_C(1,3,:) = [X_cible3,Y_cible3,Z_cible3];
end
% x_cible = X_C(1,1,1);
% y_cible = X_C(1,1,2);
% z_cible = X_C(1,1,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
waitbar(.2,f,'Loading your data');
%% Distances for each Receiver with target
nR = size(X_R);
distancesList_R = zeros(nR(1),nR(2),N_cible);
for h = 1:N_cible
    for i = 1:nR(2)
        for j = 1:nR(1)
            distancesList_R(j,i,h) = sqrt((X_R(j,i,1) - X_C(1,h,1)) .^2 + (X_R(j,i,2) - X_C(1,h,2)).^2 + (X_R(j,i,3) - X_C(1,h,3)).^2);
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
            distancesList_T(j,i,h) = sqrt((X_T(j,i,1) - X_C(1,h,1)).^2 + (X_T(j,i,2) - X_C(1,h,2)).^2 + (X_T(j,i,3) - X_C(1,h,3)).^2);
        end
    end
end
distancesList_T = reshape(distancesList_T,[Total_T,N_cible]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Paramaters for the voxels

% L_XVoxels = 6;
if isnan(maxXV) && isnan(minXV)
    XVoxels_min = -3;
    XVoxels_max = 3;
else
    if (minXV == 0 && maxXV~=0) || (minXV ~= 0 && maxXV==0) || (minXV<maxXV)
        XVoxels_min = minXV;
        XVoxels_max = maxXV;
    else
        % Display a warning message
    warningMessage = 'Error: Please make sure both are not 0 and min is less than max';
    warndlg(warningMessage, 'Warning');
    % Stop script execution using the error function
    close(f);
    error('Error: Invalid values detected! Script halted.');
    end
end
N_XVoxels = 51;

% L_ZVoxels = 6;
if isnan(maxZV) && isnan(minZV)
    ZVoxels_min = 1;
    ZVoxels_max = 5;
else
    if (minZV == 0 && maxZV~=0) || (minZV ~= 0 && maxZV==0) || (minZV<maxZV)
        ZVoxels_min = minZV;
        ZVoxels_max = maxZV;
    else
        % Display a warning message
    warningMessage = 'Error: Please make sure both are not 0 and min is less than max';
    warndlg(warningMessage, 'Warning');
    % Stop script execution using the error function
    close(f);
    error('Error: Invalid values detected! Script halted.');
    end
end
N_ZVoxels = 51; 

if isnan(maxYV) && isnan(minYV)
    YVoxels_min = -3;
    YVoxels_max = 3;
else
    if (minYV == 0 && maxYV~=0) || (minYV ~= 0 && maxYV==0)  || (minYV<maxYV)
        YVoxels_min = minYV;
        YVoxels_max = maxXV;
    else
        % Display a warning message
    warningMessage = 'Error: Please make sure both are not 0 and min is less than max';
    warndlg(warningMessage, 'Warning');
    % Stop script execution using the error function
    close(f);
    error('Error: Invalid values detected! Script halted.');
    end
end
N_YVoxels = 51; 

XV = linspace(XVoxels_min, XVoxels_max,N_XVoxels);
ZV = linspace(ZVoxels_min, ZVoxels_max,N_ZVoxels);
YV = linspace(YVoxels_min, YVoxels_max,N_YVoxels);

%% Voxels for X-Z plane with Y cte
Voxels1 = zeros(N_XVoxels,N_ZVoxels,3);

if isnan(Y_cut)
    yVoxels_coord = Y_cible1;
else
    yVoxels_coord = Y_cut;
end

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

if isnan(Z_cut)
    zVoxels_coord = Z_cible1;
else
    zVoxels_coord =Z_cut;
end
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

if isnan(X_cut)
    xVoxels_coord = X_cible1;
else
    xVoxels_coord = X_cut;
end

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

for l = 1:Nf
    for m = 1:Total_T%Number of Transmitter
        for n = 1:Total_R %Number of Receiver
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
    for n = 1:Total_T % Number of Transmitter
        for o = 1:Total_R % Number of Receiver
            Image_YCTE(k, :, :) = squeeze(Image_YCTE(k, :, :)) + S(k, n, o, :) *exp(((2i * pi * bande_f(k)) / c) * (distancesList_TV1(:, :, n) + distancesList_RV1(:, :, o))); % bande_f(l)
        end
    end
end
waitbar(.65,f,'Next Image');
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
%% Define a grid of points for the Antennas, Voxels and Target(s)
%For Y cte

%subplot(1,3,1);
plot3(ax,X_C(1,:,1), X_C(1,:,2), X_C(1,:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
hold(ax,'on');

plot3(ax,X_T(:,:,1), X_T(:,:,2), X_T(:,:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot3(ax,X_R(:,:,1), X_R(:,:,2), X_R(:,:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

%plot3(Voxels1(:,:,1), Voxels1(:,:,2), Voxels1(:,:,3), 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
% Add labels and title
xlabel(ax,'X-axis');
ylabel(ax,'Y-axis');
zlabel(ax,'Z-axis');
legend(ax,'Target','Transmitter', 'Receiver', 'Location', 'Northeast');
%legend('Transmitter', 'Receiver', 'Target', 'Voxels', 'Location', 'Northeast');
view(ax,0,90);
hold(ax,'off');
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

TotalImageReshaped_YCTE_Normalize =TotalImageReshaped_YCTE/BigMax;

[maxValueY, maxIndexY] = max(TotalImageReshaped_YCTE_Normalize(:));
[rowY, colY] = ind2sub(size(TotalImageReshaped_YCTE_Normalize), maxIndexY);

% set(hradar.XZVmaxX,'String',YV(1,rowY));
% set(hradar.XZVmaxZ,'String',ZV(1,colY));
% set(hradar.LevelXZ,'String',maxValueY);

% Plotting the 3D image with a heatmap colorbar
surf(ax_XZ,X1, Z1,TotalImageReshaped_YCTE_Normalize, 'EdgeColor', 'none');
shading(ax_XZ,"interp");
xlabel(ax_XZ,'Z-axis');
ylabel(ax_XZ,'X-axis');
title(ax_XZ,'3D Image Plot with Y cte=' + string(yVoxels_coord) +' of target on  [' + string(X_C(1,:,1)) + ',' + string(X_C(1,:,2)) + ',' + string(X_C(1,:,3)) + ']' );

% Adjusting view angle
view(ax_XZ,90, -90); % Change the view angle as needed

% Adding a colorbar
colorbar(ax_XZ);
caxis(ax_XZ,[0 1]);

dB = 20*log10(TotalImageReshaped_YCTE_Normalize);

surf(ax_XZ1,X1, Z1, dB);
shading(ax_XZ1,"interp");
xlabel(ax_XZ1,'Z-axis');
ylabel(ax_XZ1,'X-axis');
title(ax_XZ1,'3D Image Plot with Y cte=' + string(yVoxels_coord) +' of target on  [' + string(X_C(1,:,1)) + ',' + string(X_C(1,:,2)) + ',' + string(X_C(1,:,3)) + ']' );
colorbar(ax_XZ1);
clim(ax_XZ1,[-30 0]);
view(ax_XZ1,90, -90);
%% Create a scatter plot Z cte
% Create grids X and Z
[X2,Y2] = meshgrid(YV, XV);
TotalImageReshaped_ZCTE_Normalize =TotalImageReshaped_ZCTE/BigMax;

[maxValueZ, maxIndexZ] = max(TotalImageReshaped_ZCTE_Normalize(:));
[rowZ, colZ] = ind2sub(size(TotalImageReshaped_ZCTE_Normalize), maxIndexZ);

% set(hradar.XYVmaxX,'String',XV(1,rowZ));
% set(hradar.XYVmaxY,'String',YV(1,colZ));
% set(hradar.LevelXY,'String',maxValueZ);

% Plotting the 3D image with a heatmap colorbar
surf(ax_XY,X2, Y2, TotalImageReshaped_ZCTE_Normalize, 'EdgeColor', 'none');
shading(ax_XY,"interp");
xlabel(ax_XY,'Y-axis');
ylabel(ax_XY,'X-axis');
title(ax_XY,'3D Image Plot with Z cte = '+ string(zVoxels_coord) + ' of target on  [' + string(X_C(1,:,1)) + ',' + string(X_C(1,:,2)) + ',' + string(X_C(1,:,3)) + ']' );

% Adjusting view angle
view(ax_XY,90, -90); % Change the view angle as needed

% Adding a colorbar
colorbar(ax_XY);
caxis(ax_XY,[0 1]);

surf(ax_XY1,X2, Y2, 20*log10(TotalImageReshaped_ZCTE_Normalize), 'EdgeColor', 'none');
shading(ax_XY1,"interp");
xlabel(ax_XY1,'Y-axis');
ylabel(ax_XY1,'X-axis');
title(ax_XY1,'3D Image Plot with Z cte = '+ string(zVoxels_coord) + ' of target on  [' + string(X_C(1,:,1)) + ',' + string(X_C(1,:,2)) + ',' + string(X_C(1,:,3)) + ']' );

colorbar(ax_XY1);
caxis(ax_XY1,[-30 0]);
view(ax_XY1,90, -90);
%% Create a scatter plot X cte
% Create grids Y and Z
[Y3,Z3] = meshgrid(ZV, YV);

TotalImageReshaped_XCTE_Normalize =TotalImageReshaped_XCTE/BigMax;

[maxValueX, maxIndexX] = max(TotalImageReshaped_XCTE_Normalize(:));
[rowX, colX] = ind2sub(size(TotalImageReshaped_XCTE_Normalize), maxIndexX);

% set(hradar.YZVmaxY,'String',YV(1,rowX));
% set(hradar.YZVmaxZ,'String',ZV(1,colX));
% set(hradar.LevelYZ,'String',maxValueX);
% Plotting the 3D image with a heatmap colorbar
surf(ax_YZ,Y3, Z3, TotalImageReshaped_XCTE_Normalize , 'EdgeColor', 'none');
shading(ax_YZ,"interp");
xlabel(ax_YZ,'Z-axis');
ylabel(ax_YZ,'Y-axis');
title(ax_YZ,'3D Image Plot with X cte = ' + string(xVoxels_coord) + ' of target on  [' + string(X_C(1,:,1)) + ',' + string(X_C(1,:,2)) + ',' + string(X_C(1,:,3)) + ']' );

% Adjusting view angle
view(ax_YZ,90, -90); % Change the view angle as needed

% Adding a colorbar
colorbar(ax_YZ);
clim(ax_YZ,[0 1]);

surf(ax_YZ1,Y3, Z3, 20*log10(TotalImageReshaped_XCTE_Normalize ), 'EdgeColor', 'none');
shading(ax_YZ1,"interp");
xlabel(ax_YZ1,'Z-axis');
ylabel(ax_YZ1,'Y-axis');
title(ax_YZ1,'3D Image Plot with X cte = ' + string(xVoxels_coord) + ' of target on  [' + string(X_C(1,:,1)) + ',' + string(X_C(1,:,2)) + ',' + string(X_C(1,:,3)) + ']' );

colorbar(ax_YZ1);
clim(ax_YZ1,[-30 0]);
view(ax_YZ1,90, -90);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
waitbar(1,f,'Finishing');
%% Stop the timer
elapsed_time = toc;
close(f);
% Display the elapsed time
disp(['Elapsed Time: ' num2str(elapsed_time) ' seconds']);
clear variables;
end