%% 
% This is the script to 
%                       1. perform hybrid precoding algorithms for
%                       millimeter wave MIMO systems [1];
%                       2. calculate the spectral efficiency obtained by
%                       optimal precoding, hybrid precoding and beam steering for comparison.
% Basic assumptions:
%                       1. The millimeter wave channel is modelled as a narrowband clustered channnel.  
%                       2. The hybrid precoding algorithm is based on the spatially sparse nature of
%                          millimeter wave propagation.
%                       3. This is the code solely for plotting Fig.3 and 4 of single data stream (Ns = 1)
%                          case in the paper [1] 
%                           
% 
% [1] O. E. Ayach et al., ?Spatially sparse precoding in millimeter wave MIMO systems,? Mar. 2013.
%
%% Clear workplace

clear variables;

%% System parameters

Nt = [64 256]; % Number of transmit antennas
Nr = [16 64]; % Number of receive antennas

Ns = [1]; % Number of data streams
NumRF = [4 6]; % Number of RF chains for precoding and combining 

NumCluster = 8; % Number of clusters
NumRay = 10; % Number of rays per cluster

AS = 7.5; % Angular spread of 7.5 degree

SNR_dB = -40:5:0; % Range of SNR in dB

ITER = 100; % Number of channel generations 

%% Fig.3 

[Fig3_SE_Optimal,Fig3_SE_Hybrid,Fig3_SE_BeamSteering] = SpatiallySparsePrecoding(Nt(1),Nr(1),Ns,NumRF(1),NumCluster,NumRay,AS,SNR_dB,ITER);

figure(); 
l1 = plot(SNR_dB,Fig3_SE_Optimal(1,:),'-s','Color',[0 0.5 0],'LineWidth',2.0,'MarkerSize',8.0);hold on;
l2 = plot(SNR_dB,Fig3_SE_Hybrid(1,:),'-o','Color',[0 0.45 0.74],'LineWidth',2.0,'MarkerSize',8.0);hold on;
l3 = plot(SNR_dB,Fig3_SE_BeamSteering(1,:),'-d','Color',[0.85 0.33 0.1],'LineWidth',2.0,'MarkerSize',8.0);hold on;
grid on;
legend([l1 l2 l3],'Optimal unconstrained precoding','Hybrid precoding and combining','Beam steering','Location','northwest');
title('Fig.3. 64 x 16 mmWave system with 4 RF chains for sparse precoding and MMSE combining (Ns = 1)');
xlabel('SNR (dB)'); ylabel('Spectral efficiency (bits/s/Hz)');

%% Fig.4

[Fig4_SE_Optimal,Fig4_SE_Hybrid,Fig4_SE_BeamSteering] = SpatiallySparsePrecoding(Nt(2),Nr(2),Ns,NumRF(2),NumCluster,NumRay,AS,SNR_dB,ITER);

figure(); 
l1 = plot(SNR_dB,Fig4_SE_Optimal(1,:),'-s','Color',[0 0.5 0],'LineWidth',2.0,'MarkerSize',8.0);hold on;
l2 = plot(SNR_dB,Fig4_SE_Hybrid(1,:),'-o','Color',[0 0.45 0.74],'LineWidth',2.0,'MarkerSize',8.0);hold on;
l3 = plot(SNR_dB,Fig4_SE_BeamSteering(1,:),'-d','Color',[0.85 0.33 0.1],'LineWidth',2.0,'MarkerSize',8.0);hold on;
grid on;
legend([l1 l2 l3],'Optimal unconstrained precoding','Hybrid precoding and combining','Beam steering','Location','northwest');
title('Fig.4. 256 x 64 mmWave system with 6 RF chains for sparse precoding and MMSE combining (Ns = 1)');
xlabel('SNR (dB)'); ylabel('Spectral efficiency (bits/s/Hz)');


