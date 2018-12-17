function [MIMO_Channel,ArrayResponse_TX,ArrayResponse_RX,Alpha] = ChannelGenereationMIMO(Nt,Nr,NumCluster,NumRay,AS)
%%
% Input : 
%         Nt            : Number of transmit antennas 
%         Nr            : Number of receive antennas
%         NumCluster    : Number of clusters
%         NumRay        : Number of rays per cluster
%         AS            : Fixed angular spread at both the transmitter and receiver
% Output : 
%         MIMO_Chan         : Nr x Nt, MIMO channel matrix
%         ArrayResponse_TX  : Nt x (NumCluster x NumRay), transmit steering vectors
%         ArrayResponse_RX  : Nr x (NumCluster x NumRay), receive steering vectors
%         Alpha             : 1 x (NumCluster x NumRay), complex path gain
%
%% Generate random mean AOD, EOD, AOA, EOA of each cluster

% Azimuth   : [-180,180] degree
% Elevation : [0,180] degree

% Transmitter(sectorized)
% Azimuth   : 60 degree wide w.r.t 0 degree
% Elevation : 20 degree wide w.r.t 90 degree
minAOD = -30;
maxAOD = 30; 
minEOD = 80;
maxEOD = 100;
Cluster_AOD = rand(NumCluster,1)*(maxAOD-minAOD)+minAOD; % [-30,30] degree
Cluster_EOD = rand(NumCluster,1)*(maxEOD-minEOD)+minEOD; % [80,100] degree

% Receiver(omni-directional)
Cluster_AOA = (rand(NumCluster,1)-0.5)*360; % [-180,180] degree
Cluster_EOA = rand(NumCluster,1)*180; % [0,180] degree


%% Generate random AOD, EOD, AOA, EOA of rays per cluster (Laplacian distribution)

b = AS/sqrt(2); % Scaling parameter, degree

Randomness = rand(NumRay*NumCluster,1)-0.5;

% Dimension of AOD, EOD, AOA, EOA : (NumRay x NumCluster) x 1 
Ray_AOD = repelem(Cluster_AOD,NumRay,1)-b*sign(Randomness).*log(1-2.*abs(Randomness));
Ray_EOD = repelem(Cluster_EOD,NumRay,1)-b*sign(Randomness).*log(1-2.*abs(Randomness));
Ray_AOA = repelem(Cluster_AOA,NumRay,1)-b*sign(Randomness).*log(1-2.*abs(Randomness));
Ray_EOA = repelem(Cluster_EOA,NumRay,1)-b*sign(Randomness).*log(1-2.*abs(Randomness));

%% Obtain antenna element position vectors (normalized by half of the wavelength)

% Transmitter
Nt_H = sqrt(Nt);
Nt_V = sqrt(Nt);

X_Tx = zeros(1,Nt);
[Y_Tx,Z_Tx] = meshgrid(0:Nt_H-1,0:Nt_V-1);
TxPos = [X_Tx;Y_Tx(:).';Z_Tx(:).']; % 3 x Nt

% Receiver
Nr_H = sqrt(Nr);
Nr_V = sqrt(Nr);

X_Rx = zeros(1,Nr);
[Y_Rx,Z_Rx] = meshgrid(0:Nr_H-1,0:Nr_V-1);
RxPos = [X_Rx;Y_Rx(:).';Z_Rx(:).']; % 3 x Nr

%% Obtain array response vectors at the transmitter and receiver

SphericalUnitVecTx = getSphericalUnitVector(Ray_EOD,Ray_AOD); % 3 x NumRay*NumCluster
SphericalUnitVecRx = getSphericalUnitVector(Ray_EOA,Ray_AOA); % 3 x NumRay*NumCluster

ArrayResponse_TX = (1/sqrt(Nt))*exp(1i*pi*TxPos.'*SphericalUnitVecTx); % Nt x NumRay*NumCluster
ArrayResponse_RX = (1/sqrt(Nr))*exp(1i*pi*RxPos.'*SphericalUnitVecRx); % Nr x NumRay*NumCluster

%% Generate complex path gain

Alpha = sqrt(1/2)*(randn(1,NumRay*NumCluster)+1i*randn(1,NumRay*NumCluster));

%% Generate MIMO channel matrix 

MIMO_Channel = sqrt((Nt*Nr)/(NumRay*NumCluster))*ArrayResponse_RX*diag(Alpha)*ArrayResponse_TX';

end


function SphericalUnitVector = getSphericalUnitVector(theta,phi)
%% 
% Input:
%       theta, phi : M x 1
% Output:
%       SphericalUnitVector: 3 x M
% 
SphericalUnitVector = [(sind(theta).*cosd(phi)).';...
                       (sind(theta).*sind(phi)).';...
                       (cosd(theta)).']; 
end
