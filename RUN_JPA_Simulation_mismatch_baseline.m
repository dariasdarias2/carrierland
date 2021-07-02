%**************************************************************************
% RUN_JPA_Simulation
%
% This script performs a simulation for the Joint Position and Attitude
% (JPA) estimation. The main idea is to characterize and compare the
% performances of estimators for JPA against estimators which solve in
% parallel the RTK (positioning) and attitude problems.
%  
% CONFIGURATION:
%    Estimators to be used:
%           - useJpa = true/false (default true)
%           - useRtk = true/false (default true)
%           - useAtt = true/false (default true) 
%           - useEkf = true/false (default true)
%    Update rate for prediction: freqImu (default = 1)
%    Number of (GNSS) time samples: nSeconds (default = 100)
%    Number of Monte Carlo iterations: nMonteCarlo (default = 1000)
%    Number of tracked satellites +1: n (default = 7)
%    Number of on-board antennas +1: N (default 3)
%    Separation between on-board antennas: baselineLength (default 1 [m])
%    Initial base-to-rover distance: baseToRoverDistance (default 5000 [m])
%    Probability of occurrence of a cycle slip: probabilityCycleSlip  (default 0.5 / 100 [·])
%    Noise modelling:
%           - sigmaCode: standard deviation for zenith code observations (default 0.30 [m])
%           - sigmaPhase: std for zenith phase obs (default 0.003 [m])
%           - sigmaVel: std for velocity random walk (default 1 m/s/s)
%           - sigmaGyro: std for gyroscope measurements (default sqrt(10)*1E-5 [rad/sqrt(s^3)])
%           - sigmaGyroBias: std for gyroscope bias random walk (default sqrt(10)*1E-7 [rad/sqrt(s^3])
%           - sigmaAcce: (not used at the moment)
%           - sigmaAcceBias: (not used at the moment)
%
% OUTPUTS:
%           - Mean Ambiguity Success Rate (MASR) for JPA, RTK and Attitude (masrJpa, masrRtk, masrAtt)
%           - RMSE for fixed and float position and attitude estimates (rmsePosJpa, rmsePosFixJpa, rmseEulerJpa, ....)
%
% ADDITIONAL NOTES:
%   - Filters which estimate the attitude are based on the
%   Multiplicative/Error State/Indirect Extended Kalman Filter (ESKF) [2], while RTK
%   positioning is solved using EKF.
%   - The order of the elements in the state estimate is as follows:
%   x = [ position, velocity, quaternion, biasGyro, ambiguities ]
%   in case any of the elements do not appear in a filter, the order is
%   still respected (e.g., the state for the ESKF estimating the attitude
%   is xAtt = [ quaternion, biasGyro, ambiguities ]).
%   - The position and velocity are expressed in the ENU frame (inertial
%   frame), and the quaternion expresses the rotation from the local
%   (forward-left-up) to the inertial frames.
%   - Details on the filter configuration and stochastic modeling can be
%   found in [3].
% 
%--------------------------------------------------------------------------
% Date  : 08-FEB-2021                                          
% Authors       : Daniel MEDINA
% Email: daniel.ariasmedina@dlr.de
% Institute of Communications and Navigation, DLR
%--------------------------------------------------------------------------
%
% REFERENCES: 
%  [1] LAMBDA Software Package: Matlab implementation, Version 3.0.
%     Documentation provided with this software package.
%  [2] Sola, Joan. "Quaternion kinematics for the error-state Kalman filter." arXiv preprint arXiv:1711.02508 (2017).
%     IAG General Meeting, Beijing, China
%  [3] Medina, Daniel, et al. "On the Recursive Joint Position and Attitude Determination in Multi-Antenna GNSS Platforms." Remote Sensing 12.12 (2020): 1955.
%
%**************************************************************************
close all;
clear;
dbstop if error;
clc;
rng(3e4)
%% Parameter Configuration
% Method selection
useJpa = 1;
useLCJpa = 1;
useRtk = 0;
useAtt = 1;
useLCAtt = 1;
% Simulation configuration
nMonteCarlo = 10;
nSeconds = 100;
freqImu = 1;
freqGnss = 1; % This CANNOT be changed 
% Body configuration
N = 3; % Number of baselines (number of antennas -1)
n = 7; % Number of DD observations per antenna (number of satellites - 1)
baselineLength = 10; % Distance between antennas mounted of the vehicle (the shorter -> the worse attitude solutions)
baseToRoverDistance = 5000; % Initial distance to the base station (the larger -> the worse positioning solutions)
% Noise of the sensors
sigmaCode = 0.30;
sigmaPhase = 0.003;
sigmaGyro=sqrt(10)*1E-5; % std in (rad/sqrt(s^3))
sigmaBiasGyro=sqrt(10)*1E-7; % std in (rad/sqrt(s^3))
sigmaAcce=sqrt(10)*1E-3; % std in (rad/sqrt(s^3))
sigmaBiasAcce=sqrt(10)*1E-5; % std in (rad/sqrt(s^3))
sigmaAmb = (1E-8).^2; % (std in cycle/sqrt(s))
sigmaVel = 1; % (std in m/s/s)
probabilityCycleSlip = 0.5 / 100; % Probability of cycle slip (independent from each satellite);
% Initial covariances
initPosCov = 10^2;
initVelCov = 1^2;
initAttCov = (10*pi/180)^2;
initBiasGyroCov = (sigmaBiasGyro*100)^2;
initBiasAcceCov = (sigmaBiasAcce*100)^2;
initAmbCov = 5^2; 
% Misc
enablePlot = 1;
solveAmbiguities = 1;


%% Reference simulation
% Definition of the IMU and GNSS times
dtImu = 1/freqImu;
dtGnss = 1/freqGnss;
timeRef = [ 0 : dtImu : nSeconds ].'; % This is the actual reference time
timeGnss = [ 0 : dtGnss: nSeconds ].';
nTime = length( timeRef );
% Reference attitude solution over time
eulerRef = pi * [ (0.1/180)*sin(1/10*timeRef).'; (3/180)*cos(0.5*timeRef).'; 1*cos( 4*2*pi/(nSeconds)*timeRef).' ];
% Initial values for simulation:
posInit = baseToRoverDistance * [ 1; 0; 0 ];
velInit = [ -1; 3; 0 ];
eulerInit = eulerRef(:,1);
quatInit = euler_2_quat( eulerInit );
biasGyroInit = 0.05*pi/180*[ 1; -3; 0.5 ]; % Inertial sensor biases
biasAcceInit = 9.8*1e-2*[ -2; 1; 2 ];
% Location of the base station
posBaseRef = [0;0;0];
% Create (randomly) the vehicle antenna configuration in the body frame
baselineBodyHat = inScript_random_antenna_configuration(N,baselineLength);
baselineBody = [.97, 1.05, 1.02] .* baselineBodyHat;
% Integer ambiguities
M = n*(N+1);
ambInit = randi([ -293 293 ], M, 1);
% Number of unknown parameters (excluding the M ambiguities
p = 3+3+4+3;
% Memory allocation for the reference solution xRef = [ posRef; velRef; quatRef; biasGyroRef; ambRef ] \in \mathbb{R}^{3+3+4+3+M}
posRef = nan( 3, nTime ); 
velRef = nan( 3, nTime );
quatRef = nan( 4, nTime );
biasGyroRef = nan( 3, nTime );
biasAcceRef = nan( 3, nTime );
ambRef = nan( M, nTime );
xRef = nan( p+M, nTime );
cycleSlipRef = zeros( M, nTime );
angularRateRef = nan( 3, nTime );
% Process noise configuration
QVel = diag( [sigmaVel, sigmaVel, sigmaVel/1000].^2 );
Q = blkdiag( QVel, sigmaGyro^2*eye(3), sigmaBiasGyro^2*eye(3), sigmaAmb^2*eye(M) );
QRtk = blkdiag( QVel, sigmaAmb^2*eye(n) );
QAtt = blkdiag( sigmaGyro^2*eye(3), sigmaBiasGyro^2*eye(3), sigmaAmb^2*eye(n*N) );
% Reference initialization
posRef(1:3, 1) = posInit;
velRef(1:3, 1) = velInit;
quatRef(1:4, 1) = quatInit;
eulerRef(1:3, 1) = eulerInit;
biasGyroRef(1:3, 1) = biasGyroInit;
biasAcceRef(1:3, 1) = biasAcceInit;
ambRef(:, 1) = ambInit;
xRef(:,1) = [ posRef(1:3, 1); velRef(1:3, 1); quatRef(1:4, 1); biasGyroRef( 1:3, 1 ); ambRef(:, 1)  ]; 
% Loop for creating the reference trajectory
for iTime = 2:nTime
            velRef( 1:3, iTime ) = velRef(1:3, iTime-1) + dtImu*QVel*randn( 3, 1 ); % Random walk evolution for velocity
            posRef( 1:3, iTime ) = posRef(1:3, iTime-1) + dtImu * velRef(1:3, iTime-1);
%             posRef( 1:3, iTime ) = posRef(1:3, iTime-1) + dtImu * ( velRef(1:3, iTime-1) + velRef(1:3, iTime) )/2;
            quatRef( 1:4,iTime ) = euler_2_quat( eulerRef(:, iTime) );
            biasGyroRef( 1:3, iTime ) = biasGyroRef( 1:3, iTime-1 ) + sigmaBiasGyro*dtImu*randn(3,1); % Random walk evolution for biases
            biasAcceRef( 1:3, iTime ) = biasAcceRef( 1:3, iTime-1 ) + sigmaBiasAcce*dtImu*randn(3,1); % Random walk evolution for biases
            ambRef(:, iTime) = ambRef (:, iTime-1); % Ambiguities do not change over time!           
            % Getting the angular rate from the pre-existing reference attitude (eulerRef)
            deltaQuat = quat_mult(  quat_conj( quatRef( 1:4,iTime-1) ) , quatRef( 1:4,iTime) );
            normAngularRate = acos( deltaQuat(1) ) * 2 / dtImu;
            angularRateRef( 1:3,iTime ) = deltaQuat(2:4) * normAngularRate / sin( normAngularRate*dtImu/2 );
            % Induce cycle slips 
            thereIsCorrection = min( abs( timeGnss - timeRef(iTime) ) ) == 0 ;
            if thereIsCorrection % There are new GNSS observations at this time
                cycleSlipHappened = rand( M, 1 ) < probabilityCycleSlip ; % Create cycle slip occurence to challenge RTK/JPA problem
                if sum( cycleSlipHappened )>0
                    cycleSlipRef( cycleSlipHappened, iTime ) = ones( sum( cycleSlipHappened ),1 );
                    ambRef( cycleSlipHappened, iTime ) = randi( [-293 293], sum( cycleSlipHappened ), 1 ); % If any cycle slip has occurred, create a new ambiguity for that link
                end
            end
        % Updating the overall reference state estimate
        xRef(:,iTime) = [ posRef(1:3, iTime); velRef(1:3, iTime); quatRef(1:4, iTime); biasGyroRef( 1:3, iTime ); ambRef(:, iTime)  ]; 
end % End time recursion


%% GNSS Observations
%  San Fernando IGS Station on 2020/04/03 10:00
pBaseRef = [ 5105515.7657  -555145.6601  3769804.8319 ]'; % Original reference base station location (so that the satellite position make some sense) [Notice that hereinafter we will use posBaseRef]
svPos = [16869458.7696156 -4571120.85679637 19786192.4244003;25825094.4489053 -1143345.84896267 -6665718.14699749;21643438.6067279 -10475877.0319354 -11320937.0081792;18973256.5478105 6898376.19433289 17179884.9146089;13890180.6846935 22328513.4470308 4339505.27713915;17641179.2066494 19377084.5468307 -4620651.24222490;1445746.04037227 -16592974.6249599 20565913.4039738;14347765.7004308 -18312132.8082430 12281607.4295026;10494147.4945647 13154518.8315525 21104041.2333542;3778444.49434751 21414391.3640113 15142134.4168146].';
svElev = [1.25529313212859;0.465399395908134;0.210469879338600;1.09770614106603;0.220919585319286;0.132524707563289;0.367673268496015;0.698695441109834;0.630911754365399;0.153500508281829]; % Notice, the first satellite has the highest elevation
svAzim = [5.83593162385241;3.06170518951831;3.48552529935196;1.23894419986494;1.70757929674904;2.11715311391813;5.44297992037407;4.72703506702490;0.884656095666627;1.08914542364832];
svPos = [ ecef2enu( svPos, pBaseRef ) ]; % Now, the satellite positions are in the ENU frame
n = min( size( svPos, 2 ) -1, n );  
svPos = svPos(:, 1 : n+1 );
svElev  = svElev( 1 : n+1 );
svAzim  = svAzim( 1 : n+1 );
lambda = 0.1904; % Wavelength for GPS L1 [m]
% GNSS Design matrices
G = [ - [ svPos(:,1:n+1) - repmat( pBaseRef, 1, n+1)] ./ vecnorm( svPos(:,1:n+1) - repmat( pBaseRef, 1, n+1), 2, 1 ) ].'; % Geometry matrix in simple RTK problem (without the difference to the pivot satellite!)
A = lambda*eye(n); % A design matrix in RTK problem
DD = [ -ones(n,1), eye(n) ]; % Double-difference matrix
DD_Jpa = ones(N+1)+eye(N+1); % Matrix to create correlations between the master and the base and slave antennas
B = DD*G; % B design matrix in RTK problem
%%%% NOTE: Construction of the covariance matrix of the observations is
%%%% performed by calling [ R_Jpa, R_Rtk, R_att ] = inScript_jpa_noise_model( sigmaPhase, sigmaCode, distance, svElev, n, N, DD );
% % % % Noise covariances 
% % % W = diag( (1 + 10.*exp( -svElev./(10*pi/180) )).^2 );
% % % R_iono = W * (sqrt(2)*0.4*baseDistance/1000/1000)^2; % The ionosphere weighted model is \sigma_I = sqrt(2)*0.4 mm/km and amplified with the elevation dependent weighting matrix
% % % R_phase = W * sigmaPhase^2; % Covariance matrix for code observations (elevation dependent)
% % % R_code = W * sigmaCode^2; % Covariance matrix for phase observations (elevation dependent)
% % % R_Jpa = blkdiag( kron(  DD_Jpa, DD * R_phase * DD.'), kron(  DD_Jpa, DD * R_code * DD.') ); % This is the general JPA covariance matrix w/o the ionospheric residuals due to the baseline length
% % % R_Jpa( 1:n, 1:n ) = R_Jpa( 1:n, 1:n ) + DD*R_iono*DD.';
% % % R_Jpa( n*(N+1)+1:n*(N+1)+n, n*(N+1)+1:n*(N+1)+n ) = R_Jpa( n*(N+1)+1:n*(N+1)+n, n*(N+1)+1:n*(N+1)+n ) + DD*R_iono*DD.';
% % % R_Jpa( 1:n, n*(N+1)+1:n*(N+1)+n ) = DD*R_iono*DD.';
% % % R_Jpa( n*(N+1)+1:n*(N+1)+n, 1:n ) = DD*R_iono*DD.';
% % % R_Rtk = blkdiag( DD * (R_phase+R_iono) * DD.', DD * (R_code+R_iono) * DD.' );
% % % R_att = blkdiag( kron(  1+eye(N), DD * R_phase * DD.'), kron(  1+eye(N), DD * R_code * DD.') );


% [ R_Jpa, R_Rtk, R_att ] = inScript_jpa_noise_model( sigmaPhase, sigmaCode, baseDistance, ones(n+1,1), n, N, DD );
% figure; hold on; grid on; box on; axis equal;
% % clims = [ min(R_Jpa(:)), max(R_Jpa(:)) ];
% % imagesc(R_Jpa(M+1:end,M+1:end), clims)
% aux=R_Jpa(1:M,1:M);
% clims = [ 0, max(aux(:)) ];
% imagesc(R_Jpa(1:M,1:M), clims)
%  a = makeColorMap(0.95*[1 1 1],rgb('gold'),rgb('red'),30);
% colormap(a);
% for i = 1:M+1
%    plot([.5,M+0.5],[i-.5,i-.5],':','LineWidth',0.30,'Color',0.5*[1,1,1]);
%    plot([i-.5,i-.5],[.5,M+0.5],':','LineWidth',0.30,'Color',0.5*[1,1,1]);
% end
% plot([0.5,0.5],[0.5,M+0.5],'k-');
% plot([M+0.5,M+0.5],[0.5,M+0.5],'k-');
% plot([0.5,M+0.5],[0.5,0.5],'k-');
% plot([0.5,M+0.5],[M+0.5,M+0.5],'k-');
% set(gca,'YDir','reverse')
% xlim([0.5,M+0.5]); ylim([0.5,M+0.5]);
% ccc = colorbar;
% ccc.Label.String = '$[\mathbf{Q}_{\Phi}]_{i,i}$ [m]';


%% Variables initialization for the estimators
if useJpa
    errorJpa = nan( p+M, nTime, nMonteCarlo );
    errorPosJpa = nan( 3, nTime, nMonteCarlo );
    errorVelJpa = nan( 3, nTime, nMonteCarlo );
    errorEulerJpa = nan( 3, nTime, nMonteCarlo );
    errorAmbJpa = nan( M, nTime, nMonteCarlo );
    isFixJpa = nan( nTime, nMonteCarlo );
     errorPosFixJpa = nan( 3, nTime, nMonteCarlo );
    errorEulerFixJpa = nan( 3, nTime, nMonteCarlo );
    errorAmbFixJpa = nan( M, nTime, nMonteCarlo );
end % End of JPA-ESKF variable initialization

if useLCJpa
    errorLCJpa = nan( p+M, nTime, nMonteCarlo );
    errorPosLCJpa = nan( 3, nTime, nMonteCarlo );
    errorVelLCJpa = nan( 3, nTime, nMonteCarlo );
    errorEulerLCJpa = nan( 3, nTime, nMonteCarlo );
    errorAmbLCJpa = nan( M, nTime, nMonteCarlo );
    isFixLCJpa = nan( nTime, nMonteCarlo );
     errorPosFixLCJpa = nan( 3, nTime, nMonteCarlo );
    errorEulerFixLCJpa = nan( 3, nTime, nMonteCarlo );
    errorAmbFixLCJpa = nan( M, nTime, nMonteCarlo );
end % End of JPA-ESKF variable initialization

if useRtk
    errorRtk = nan( 3+3+n, nTime, nMonteCarlo );
    errorPosRtk = nan( 3, nTime, nMonteCarlo );
    errorVelRtk = nan( 3, nTime, nMonteCarlo );
    errorAmbRtk = nan( n, nTime, nMonteCarlo );
    isFixRtk = nan( nTime, nMonteCarlo );
    errorPosFixRtk = nan( 3, nTime, nMonteCarlo );
    errorAmbFixRtk = nan( n, nTime, nMonteCarlo );
end % End of RTK-EKF variable initialization

if useAtt
    errorBiasAtt = nan( 3, nTime, nMonteCarlo );
    errorEulerAtt = nan( 3, nTime, nMonteCarlo );
    errorAmbAtt = nan( n*N, nTime, nMonteCarlo );
    isFixAtt = nan( nTime, nMonteCarlo );
    errorEulerFixAtt = nan( 3, nTime, nMonteCarlo );
    errorAmbFixAtt = nan( n*N, nTime, nMonteCarlo );
end % End of Att-ESKF variable initialization

if useLCAtt
    errorBiasLCAtt = nan( 3, nTime, nMonteCarlo );
    errorEulerLCAtt = nan( 3, nTime, nMonteCarlo );
    errorAmbLCAtt = nan( n*N, nTime, nMonteCarlo );
    isFixLCAtt = nan( nTime, nMonteCarlo );
    errorEulerFixLCAtt = nan( 3, nTime, nMonteCarlo );
    errorAmbFixLCAtt = nan( n*N, nTime, nMonteCarlo );
end % End of LCAtt-ESKF variable initialization


%% Monte Carlo Experiment
for iMonteCarlo = 1:nMonteCarlo
    
    % Randomly sampled the initial state estimate around the initial covariance matrix (equal for all evaluated filters):
    P00 = blkdiag( initPosCov*eye(3), initVelCov*eye(3), initAttCov*eye(3), initBiasGyroCov*eye(3), initAmbCov*eye(M) );
    x00 = [ posInit + sqrt(initPosCov)*randn(3,1); ...
        velInit + sqrt(initVelCov)*randn(3,1); ...
        euler_2_quat( eulerInit + sqrt(initAttCov)*randn(3,1) ); ...
        biasGyroInit + sqrt(initBiasGyroCov)*randn(3,1); ...
        ambInit + sqrt(initAmbCov)*randn(M,1); ];
    
    % Filters initialization 
    if useJpa == true % Initialization for JPA - ESKF filter
        xJpa = nan( p+M, nTime );
        fixJpa = nan( p+M, nTime );
        xJpa(:,1) = x00;
        PJpa = P00;
    end
    
    if useLCJpa == true % Initialization for JPA - ESKF filter
        xLCJpa = nan( p+M, nTime );
        fixLCJpa = nan( p+M, nTime );
        xLCJpa(:,1) = x00;
        PLCJpa = P00;
    end
    
    if useRtk == true % Initialization for RTK - EKF filter
        xRtk = nan( 3+3+n, nTime );
        fixRtk = nan( 3+3+n, nTime );
        xRtk(:,1) = [ x00(1:3); x00(4:6); x00(p+1:p+n); ] ;
        PRtk = blkdiag( initPosCov*eye(3), initVelCov*eye(3), initAmbCov*eye(n) );
    end
    
    if useAtt == true % Initialization for Attitude - ESKF filter
        xAtt = nan( 4+3+n*N, nTime );
        fixAtt = nan( 4+3+n*N, nTime );
        xAtt(:,1) = [ x00(7:10); x00(11:13); x00(p+n+1:end); ] ;
        PAtt = blkdiag( initAttCov*eye(3), initBiasGyroCov*eye(3), initAmbCov*eye(n*N) );
    end
    
    if useLCAtt == true % Initialization for Attitude - ESKF filter
        xLCAtt = nan( 4+3+n*N, nTime );
        fixLCAtt = nan( 4+3+n*N, nTime );
        xLCAtt(:,1) = [ x00(7:10); x00(11:13); x00(p+n+1:end); ] ;
        PLCAtt = blkdiag( initAttCov*eye(3), initBiasGyroCov*eye(3), initAmbCov*eye(n*N) );
    end
    
    
    % Main loop
    for iTime = 2:nTime

        % Prediction Step
        omega = angularRateRef( 1:3, iTime ) + biasGyroRef( 1:3, iTime-1 ) + sigmaGyro*randn( 3, 1 );  % Adding the noise and gyroscope bias to the actual angular rate
        
        if useJpa % Prediction for the JPA ESKF
            % Propagate state estimate
            xJpa(:, iTime) = xJpa(:, iTime-1);
            xJpa( 1:3, iTime ) = xJpa( 1:3, iTime ) + dtImu * xJpa( 4:6, iTime-1 );
            omegaHat = omega - xJpa( 11:13, iTime-1 );
            deltaTheta = norm( omegaHat ) * dtImu / 2;
            deltaQuat = [ cos( deltaTheta ); omegaHat/norm(omegaHat) * sin(deltaTheta) ];
            xJpa( 7:10, iTime ) = quat_mult( xJpa( 7:10, iTime-1 ), deltaQuat   );
            % Update covariance
            Fx = [  eye(3), eye(3)*dtImu, zeros(3), zeros(3), zeros(3,M); ...
                zeros(3), eye(3), zeros(3), zeros(3), zeros(3,M);...
                zeros(3), zeros(3), quat_2_DCM(deltaQuat).', -eye(3)*dtImu, zeros(3,M);...
                zeros(3), zeros(3), zeros(3), eye(3), zeros(3,M);...
                zeros(M,3), zeros(M,3), zeros(M,3), zeros(M,3), eye(M) ];
            Fq = [ eye(3)*dtImu^2/2, zeros(3), zeros(3), zeros(3,M); ...
                eye(3)*dtImu, zeros(3), zeros(3), zeros(3,M); ...
                zeros(3), eye(3)*dtImu^2/2, zeros(3), zeros(3,M); ...
                zeros(3), zeros(3), eye(3)*dtImu, zeros(3,M); ...
                zeros(M,3), zeros(M,3), zeros(M,3), eye(M)*dtImu ];
            PJpa = Fx * PJpa * Fx.' + Fq*Q*Fq.';
        end % End of ESKF Prediction
        
        if useLCJpa % Prediction for the JPA ESKF
            % Propagate state estimate
            xLCJpa(:, iTime) = xLCJpa(:, iTime-1);
            xLCJpa( 1:3, iTime ) = xLCJpa( 1:3, iTime ) + dtImu * xLCJpa( 4:6, iTime-1 );
            omegaHat = omega - xLCJpa( 11:13, iTime-1 );
            deltaTheta = norm( omegaHat ) * dtImu / 2;
            deltaQuat = [ cos( deltaTheta ); omegaHat/norm(omegaHat) * sin(deltaTheta) ];
            xLCJpa( 7:10, iTime ) = quat_mult( xLCJpa( 7:10, iTime-1 ), deltaQuat   );
            % Update covariance
            Fx = [  eye(3), eye(3)*dtImu, zeros(3), zeros(3), zeros(3,M); ...
                zeros(3), eye(3), zeros(3), zeros(3), zeros(3,M);...
                zeros(3), zeros(3), quat_2_DCM(deltaQuat).', -eye(3)*dtImu, zeros(3,M);...
                zeros(3), zeros(3), zeros(3), eye(3), zeros(3,M);...
                zeros(M,3), zeros(M,3), zeros(M,3), zeros(M,3), eye(M) ];
            Fq = [ eye(3)*dtImu^2/2, zeros(3), zeros(3), zeros(3,M); ...
                eye(3)*dtImu, zeros(3), zeros(3), zeros(3,M); ...
                zeros(3), eye(3)*dtImu^2/2, zeros(3), zeros(3,M); ...
                zeros(3), zeros(3), eye(3)*dtImu, zeros(3,M); ...
                zeros(M,3), zeros(M,3), zeros(M,3), eye(M)*dtImu ];
            PLCJpa = Fx * PLCJpa * Fx.' + Fq*Q*Fq.';
        end % End of ESKF Prediction
        
         if useRtk % Prediction for RTK (EKF based)
            % Propagate state estimate
            xRtk(:, iTime) = xRtk(:, iTime-1);
            xRtk( 1:3, iTime ) = xRtk( 1:3, iTime ) + dtImu * xRtk( 4:6, iTime-1 );
            % Update covariance
            Fx = [  eye(3), eye(3)*dtImu, zeros(3,n); ...
                zeros(3), eye(3), zeros(3,n);...
                zeros(n,3), zeros(n,3), eye(n) ];
            Fq = [ eye(3)*dtImu^2/2, zeros(3,n); ...
                eye(3)*dtImu, zeros(3,n); ...
                zeros(n,3), eye(n)*dtImu ];
            PRtk = Fx * PRtk * Fx.' + Fq*QRtk*Fq.';
        end % End of RTK Prediction
        
        if useAtt % Prediction for the Attitude ESKF
            % Propagate state estimate
            xAtt(:, iTime) = xAtt(:, iTime-1);
            omegaHat = omega - xAtt( 5:7, iTime-1 );
            deltaTheta = norm( omegaHat ) * dtImu / 2;
            deltaQuat = [ cos( deltaTheta ); omegaHat/norm(omegaHat) * sin(deltaTheta) ];
            xAtt( 1:4, iTime ) = quat_mult( xAtt( 1:4, iTime-1 ), deltaQuat   );
            % Update covariance
            Fx = [  quat_2_DCM(deltaQuat).', -eye(3)*dtImu, zeros(3,n*N);...
                 zeros(3), eye(3), zeros(3,n*N);...
                zeros(n*N,3), zeros(n*N,3), eye(n*N) ];
            Fq = [ eye(3)*dtImu^2/2, zeros(3), zeros(3,n*N); ...
                zeros(3), eye(3)*dtImu, zeros(3,n*N); ...
                zeros(n*N,3), zeros(n*N,3), eye(n*N)*dtImu ];
            PAtt = Fx * PAtt * Fx.' + Fq*QAtt*Fq.';
        end % End of ESKF Prediction
        
        if useLCAtt % Prediction for the Attitude ESKF
            % Propagate state estimate
            xLCAtt(:, iTime) = xLCAtt(:, iTime-1);
            omegaHat = omega - xLCAtt( 5:7, iTime-1 );
            deltaTheta = norm( omegaHat ) * dtImu / 2;
            deltaQuat = [ cos( deltaTheta ); omegaHat/norm(omegaHat) * sin(deltaTheta) ];
            xLCAtt( 1:4, iTime ) = quat_mult( xLCAtt( 1:4, iTime-1 ), deltaQuat   );
            % Update covariance
            Fx = [  quat_2_DCM(deltaQuat).', -eye(3)*dtImu, zeros(3,n*N);...
                 zeros(3), eye(3), zeros(3,n*N);...
                zeros(n*N,3), zeros(n*N,3), eye(n*N) ];
            Fq = [ eye(3)*dtImu^2/2, zeros(3), zeros(3,n*N); ...
                zeros(3), eye(3)*dtImu, zeros(3,n*N); ...
                zeros(n*N,3), zeros(n*N,3), eye(n*N)*dtImu ];
            PLCAtt = Fx * PLCAtt * Fx.' + Fq*QAtt*Fq.';
        end % End of ESKF Prediction
                
                
        % Correction Step
        thereIsCorrection = min( abs( timeGnss - timeRef(iTime) ) ) == 0 ;
        if thereIsCorrection % There are new GNSS observations at this time
            % GNSS Measurements Generator
            baselineRef = nan(3*N,1);
            for iBaseline = 1:3:3*N
                baselineRef( iBaseline:iBaseline+2 ) = quat_rotate( quatRef( :, iTime ), baselineBody( :, ceil(iBaseline/3) ) );
            end
            [ R_Jpa, R_Rtk, R_att ] = inScript_jpa_noise_model( sigmaPhase, sigmaCode, norm( posRef(1:3, iTime)), svElev, n, N, DD );
            noise_Jpa = mgd( 1, 2*M, zeros(2*M,1), R_Jpa ).';
            noise_Rtk = [ noise_Jpa(1:n); noise_Jpa(M+1:M+n) ];
            noise_att = [ noise_Jpa(n+1:M); noise_Jpa(M+n+1:end) ];
            y_Jpa = [ kron( eye(N+1), DD*G ) * [ posRef(1:3, iTime); baselineRef( 1:3*N ) ] + lambda*ambRef( 1:M, iTime ); ...
                          kron( eye(N+1), DD*G ) * [ posRef(1:3, iTime); baselineRef( 1:3*N ) ] ]; 
            y_Jpa = y_Jpa + noise_Jpa;
            y_Rtk = [ y_Jpa(1:n); y_Jpa(M+1:M+n) ] + noise_Rtk;
            y_att = [ y_Jpa(n+1:M); y_Jpa(M+n+1:end) ] + noise_att;
            
            
            if useLCJpa % Correction step for JPA- ESKF
                % Cycle slip repair
                if sum( cycleSlipRef(:,iTime) ) > 0
                    cycleSlipInd = find( cycleSlipRef(:,iTime) == 1 );
                    xLCJpa( p+cycleSlipInd, iTime ) = ( y_Jpa( cycleSlipInd ) - y_Jpa( cycleSlipInd+M ) )/lambda;
                    for ind=1:length( cycleSlipInd )
                        PLCJpa( p-1+cycleSlipInd(ind), :  ) = 0;
                        PLCJpa( :, p-1+cycleSlipInd(ind)  ) = 0;
                        PLCJpa( p-1+cycleSlipInd(ind),  p-1+cycleSlipInd(ind) ) =  10 * R_Jpa( M+cycleSlipInd(ind), M+cycleSlipInd(ind) ); 
                    end
                end
                % Observation model
                baselineHat = nan(3*N,1);
                for iBaseline = 1:3:3*N
                    baselineHat( iBaseline:iBaseline+2 ) = quat_rotate( xLCJpa( 7:10, iTime ), baselineBodyHat( :, ceil(iBaseline/3) ) );
                end
                hx = [ kron( eye(N+1), DD*G ) * [ xLCJpa(1:3, iTime); baselineHat( 1:3*N ) ] + lambda*xLCJpa( p+1:end, iTime ); ...
                          kron( eye(N+1), DD*G ) * [ xLCJpa(1:3, iTime); baselineHat( 1:3*N ) ] ]; 
                % Jacobian for the observation model       
                HQuat = zeros( n, 4 );
                for iBaseline = 1:N
                    Jq = inScript_derivative_Baseline_wrt_Quaternion( xLCJpa( 7:10, iTime ), baselineBodyHat(:,iBaseline) );
                    HQuat = [ HQuat; DD*G*Jq ]; % Getting the Jacobian related to the quaternion
                end
                HQuat = kron( [1;1], HQuat );
                HPos = [ DD*G; zeros( n*N, 3 ); DD*G; zeros( n*N, 3 ) ];
                HVel = zeros( 2*M, 3 );
                HBiasGyro = zeros( 2*M, 3 );
                HAmb = [ lambda * eye(M); zeros(M) ];
                Hx = [ HPos, HVel, HQuat, HBiasGyro, HAmb ];
                % Jacobian for the error state
                HdeltaTheta = 1/2 * quat_productMatrix ( xLCJpa( 7:10, iTime ), 'left' ) * [zeros(1,3); eye(3)];
                HdeltaX = blkdiag( eye(3), eye(3), HdeltaTheta, eye(3),  eye(M) );
                % Jacobian for the total state
                H = Hx * HdeltaX;
                % Covariance matrix update
                S_inv = pinv( H*PLCJpa*H' + R_Jpa );
                K = PLCJpa * H'* S_inv;
                %Constraints
                Delta = zeros(2*(N+1)*n, N);
                for antInd=1:N
                    Delta(n*antInd + 1 : n*(antInd + 1), antInd) = DD*G*baselineHat(3*antInd - 2 : 3*antInd);
                    Delta((N+1 + antInd)*n + 1 : (N+1 + antInd + 1)*n, antInd) = DD*G*baselineHat(3*antInd - 2 : 3*antInd);
                end
                T = zeros(2*(N+1)*n, N);
                Gamma = - K * Delta;
                Psi = Delta' * S_inv * Delta;
                L = K + Gamma * pinv(Psi) * Delta' * S_inv;
                PLCJpa = ( eye(p+M-1) - K*H )*PLCJpa + Gamma * pinv(Psi) * Gamma';
                % Error state update
                errorState = L * ( y_Jpa-hx );
                xLCJpa( 1:6, iTime ) = xLCJpa( 1:6, iTime ) + errorState(1:6);
                xLCJpa( 11:end, iTime ) = xLCJpa( 11:end, iTime ) + errorState(10:end);
                deltaTheta = errorState(7:9);
                deltaQuat = [ cos(norm(deltaTheta)/2); deltaTheta/norm(deltaTheta)*sin(norm(deltaTheta)/2) ];
                xLCJpa(7:10, iTime ) = quat_mult( xLCJpa(7:10, iTime ), deltaQuat );
                
                if solveAmbiguities == true
                    [ aTmp, rTmp ] = LAMBDA( xLCJpa(p+1:end, iTime), PLCJpa( p:end, p:end ) ); % Solving the ILS problem
                    % Fix solution estimation
                    fixLCJpa(:, iTime) = xLCJpa(:, iTime );
                    deltaErrorConditioned = PLCJpa( 1:p-1, p:end )*pinv(PLCJpa( p:end, p:end )) * ( xLCJpa(p+1:end, iTime ) - aTmp(:,1) );
                    fixLCJpa( 1:6, iTime ) = xLCJpa( 1:6, iTime ) - deltaErrorConditioned(1:6);
                    fixLCJpa( 11:13, iTime ) = xLCJpa( 11:13, iTime ) - deltaErrorConditioned(10:12);
                    deltaTheta = deltaErrorConditioned(7:9);
                    deltaQuat = [ cos(norm(deltaTheta)/2); deltaTheta/norm(deltaTheta)*sin(norm(deltaTheta)/2)  ];
                    fixLCJpa( 7:10, iTime ) = quat_mult( xLCJpa(7:10, iTime ), quat_conj(deltaQuat) );
                    fixLCJpa( p+1:end, iTime ) = aTmp(:,1);
                    isFixLCJpa( iTime, iMonteCarlo ) = norm( aTmp(:,1) - ambRef(:, iTime) ) == 0;
                end % End of solving ambiguities for ESKF         
            
            end % End of JPA - ESLCKF correction
            
            
            if useJpa % Correction step for JPA- ESKF
                % Cycle slip repair
                if sum( cycleSlipRef(:,iTime) ) > 0
                    cycleSlipInd = find( cycleSlipRef(:,iTime) == 1 );
                    xJpa( p+cycleSlipInd, iTime ) = ( y_Jpa( cycleSlipInd ) - y_Jpa( cycleSlipInd+M ) )/lambda;
                    for ind=1:length( cycleSlipInd )
                        PJpa( p-1+cycleSlipInd(ind), :  ) = 0;
                        PJpa( :, p-1+cycleSlipInd(ind)  ) = 0;
                        PJpa( p-1+cycleSlipInd(ind),  p-1+cycleSlipInd(ind) ) =  10 * R_Jpa( M+cycleSlipInd(ind), M+cycleSlipInd(ind) ); 
                    end
                end
                % Observation model
                baselineHat = nan(3*N,1);
                for iBaseline = 1:3:3*N
                    baselineHat( iBaseline:iBaseline+2 ) = quat_rotate( xJpa( 7:10, iTime ), baselineBodyHat( :, ceil(iBaseline/3) ) );
                end
                hx = [ kron( eye(N+1), DD*G ) * [ xJpa(1:3, iTime); baselineHat( 1:3*N ) ] + lambda*xJpa( p+1:end, iTime ); ...
                          kron( eye(N+1), DD*G ) * [ xJpa(1:3, iTime); baselineHat( 1:3*N ) ] ]; 
                % Jacobian for the observation model       
                HQuat = zeros( n, 4 );
                for iBaseline = 1:N
                    Jq = inScript_derivative_Baseline_wrt_Quaternion( xJpa( 7:10, iTime ), baselineBodyHat(:,iBaseline) );
                    HQuat = [ HQuat; DD*G*Jq ]; % Getting the Jacobian related to the quaternion
                end
                HQuat = kron( [1;1], HQuat );
                HPos = [ DD*G; zeros( n*N, 3 ); DD*G; zeros( n*N, 3 ) ];
                HVel = zeros( 2*M, 3 );
                HBiasGyro = zeros( 2*M, 3 );
                HAmb = [ lambda * eye(M); zeros(M) ];
                Hx = [ HPos, HVel, HQuat, HBiasGyro, HAmb ];
                % Jacobian for the error state
                HdeltaTheta = 1/2 * quat_productMatrix ( xJpa( 7:10, iTime ), 'left' ) * [zeros(1,3); eye(3)];
                HdeltaX = blkdiag( eye(3), eye(3), HdeltaTheta, eye(3),  eye(M) );
                % Jacobian for the total state
                H = Hx * HdeltaX;
                % Covariance matrix update
                K = PJpa * H'* pinv( H*PJpa*H' + R_Jpa );
                PJpa = ( eye(p+M-1) - K*H )*PJpa;
                % Error state update
                errorState = K * ( y_Jpa-hx );
                xJpa( 1:6, iTime ) = xJpa( 1:6, iTime ) + errorState(1:6);
                xJpa( 11:end, iTime ) = xJpa( 11:end, iTime ) + errorState(10:end);
                deltaTheta = errorState(7:9);
                deltaQuat = [ cos(norm(deltaTheta)/2); deltaTheta/norm(deltaTheta)*sin(norm(deltaTheta)/2) ];
                xJpa(7:10, iTime ) = quat_mult( xJpa(7:10, iTime ), deltaQuat );
                
                if solveAmbiguities == true
                    [ aTmp, rTmp ] = LAMBDA( xJpa(p+1:end, iTime), PJpa( p:end, p:end ) ); % Solving the ILS problem
                    % Fix solution estimation
                    fixJpa(:, iTime) = xJpa(:, iTime );
                    deltaErrorConditioned = PJpa( 1:p-1, p:end )*pinv(PJpa( p:end, p:end )) * ( xJpa(p+1:end, iTime ) - aTmp(:,1) );
                    fixJpa( 1:6, iTime ) = xJpa( 1:6, iTime ) - deltaErrorConditioned(1:6);
                    fixJpa( 11:13, iTime ) = xJpa( 11:13, iTime ) - deltaErrorConditioned(10:12);
                    deltaTheta = deltaErrorConditioned(7:9);
                    deltaQuat = [ cos(norm(deltaTheta)/2); deltaTheta/norm(deltaTheta)*sin(norm(deltaTheta)/2)  ];
                    fixJpa( 7:10, iTime ) = quat_mult( xJpa(7:10, iTime ), quat_conj(deltaQuat) );
                    fixJpa( p+1:end, iTime ) = aTmp(:,1);
                    isFixJpa( iTime, iMonteCarlo ) = norm( aTmp(:,1) - ambRef(:, iTime) ) == 0;
                end % End of solving ambiguities for ESKF         
            
            end % End of JPA - ESKF correction
            
            
            if useRtk % Correction step for RTK-EKF
                % Cycle slip repair
                if sum( cycleSlipRef(1:n,iTime) ) > 0
                    cycleSlipInd = find( cycleSlipRef(1:n,iTime) == 1 );
                    xRtk( 3+3+cycleSlipInd, iTime ) = ( y_Rtk( cycleSlipInd ) - y_Rtk( cycleSlipInd+n ) )/lambda;
                    for ind=1:length( cycleSlipInd )
                        PRtk( 3+3+cycleSlipInd(ind), :  ) = 0;
                        PRtk( :, 3+3+cycleSlipInd(ind)  ) = 0;
                        PRtk( 3+3+cycleSlipInd(ind),  3+3+cycleSlipInd(ind) ) =  10 * R_Rtk( n+cycleSlipInd(ind), n+cycleSlipInd(ind) ); 
                    end
                end
                % Observation model
                hx = [ DD*G*xRtk(1:3, iTime) + lambda*xRtk( 3+3+1:end, iTime );      DD*G*xRtk(1:3, iTime) ]; 
                % Jacobian for the observation model       
                HPos = [ DD*G; DD*G ];
                HVel = zeros( 2*n, 3 );
                HAmb = [ lambda * eye(n); zeros(n) ];
                H = [ HPos, HVel, HAmb ];
                % Covariance matrix update
                K = PRtk * H'* pinv( H*PRtk*H' + R_Rtk );
                PRtk = ( eye(3+3+n) - K*H )*PRtk;
                % Error state update
                 xRtk( :, iTime ) = xRtk( :, iTime ) + K * ( y_Rtk-hx );
                
                if solveAmbiguities == true
                    [ aTmp, rTmp ] = LAMBDA( xRtk(3+3+1:end, iTime), PRtk( 3+3+1:end, 3+3+1:end ) ); % Solving the ILS problem
                    % Fix solution estimation
                    fixRtk(:, iTime) = xRtk(:, iTime );
                    fixRtk( 1:6, iTime ) = xRtk( 1:6, iTime ) - PRtk( 1:3+3, 3+3+1:end )*pinv(PRtk( 3+3+1:end, 3+3+1:end )) * ( xRtk(3+3+1:end, iTime ) - aTmp(:,1) );
                    fixRtk( 3+3+1:end, iTime ) = aTmp(:,1);
                    isFixRtk( iTime, iMonteCarlo ) = norm( aTmp(:,1) - ambRef(1:n, iTime) ) == 0;
                end % End of solving ambiguities for RTK 
                
            end % End of RTK-EKF correction
            
            if useAtt % Correction step for Attitude - ESKF
                % Cycle slip repair
                if sum( cycleSlipRef(:,iTime) ) > 0
                    cycleSlipInd = find( cycleSlipRef(n+1:end,iTime) == 1 );
                    xAtt( 4+3+cycleSlipInd, iTime ) = ( y_att( cycleSlipInd ) - y_att( cycleSlipInd+n*N ) )/lambda;
                    for ind=1:length( cycleSlipInd )
                        PAtt( 4+3-1+cycleSlipInd(ind), :  ) = 0;
                        PAtt( :, 4+3-1+cycleSlipInd(ind)  ) = 0;
                        PAtt( 4+3-1+cycleSlipInd(ind),  4+3-1+cycleSlipInd(ind) ) =  10 * R_att( n*N+cycleSlipInd(ind), n*N+cycleSlipInd(ind) ); 
                    end
                end
                % Observation model
                baselineHat = nan(3*N,1);
                for iBaseline = 1:3:3*N
                    baselineHat( iBaseline:iBaseline+2 ) = quat_rotate( xAtt( 1:4, iTime ), baselineBody( :, ceil(iBaseline/3) ) );
                end
                hx = [ kron( eye(N), DD*G ) * baselineHat( 1:3*N ) + lambda*xAtt( 4+3+1:end, iTime ); ...
                          kron( eye(N), DD*G ) * baselineHat( 1:3*N ) ]; 
                % Jacobian for the observation model       
                HQuat = [];
                for iBaseline = 1:N
                    Jq = inScript_derivative_Baseline_wrt_Quaternion( xAtt( 1:4, iTime ), baselineBody(:,iBaseline) );
                    HQuat = [ HQuat; DD*G*Jq ]; % Getting the Jacobian related to the quaternion
                end
                HQuat = kron( [1;1], HQuat );
                HBiasGyro = zeros( 2*n*N, 3 );
                HAmb = [ lambda * eye(n*N); zeros(n*N) ];
                Hx = [ HQuat, HBiasGyro, HAmb ];
                % Jacobian for the error state
                HdeltaTheta = 1/2 * quat_productMatrix ( xAtt( 1:4, iTime ), 'left' ) * [zeros(1,3); eye(3)];
                HdeltaX = blkdiag( HdeltaTheta, eye(3),  eye(n*N) );
                % Jacobian for the total state
                H = Hx * HdeltaX;
                % Covariance matrix update
                K = PAtt * H'* pinv( H*PAtt*H' + R_att);
                PAtt = ( eye(4+3+n*N-1) - K*H )*PAtt;
                % Error state update
                errorState = K * ( y_att-hx );
                xAtt( 5:end, iTime ) = xAtt( 5:end, iTime ) + errorState(4:end);
                deltaTheta = errorState(1:3);
                deltaQuat = [ cos(norm(deltaTheta)/2); deltaTheta/norm(deltaTheta)*sin(norm(deltaTheta)/2) ];
                xAtt(1:4, iTime ) = quat_mult( xAtt(1:4, iTime ), deltaQuat );
                
                if solveAmbiguities == true
                    [ aTmp, rTmp ] = LAMBDA( xAtt(4+3+1:end, iTime), PAtt( 4+3:end, 4+3:end ) ); % Solving the ILS problem
                    % Fix solution estimation
                    fixAtt(:, iTime) = xAtt(:, iTime );
                    deltaErrorConditioned = PAtt( 1:4+3-1, 4+3:end )*pinv(PAtt( 4+3:end, 4+3:end )) * ( xAtt(4+3+1:end, iTime ) - aTmp(:,1) );
                    fixAtt( 5:7, iTime ) = xAtt( 5:7, iTime ) - deltaErrorConditioned(4:6);
                    deltaTheta = deltaErrorConditioned(1:3);
                    deltaQuat = [ cos(norm(deltaTheta)/2); deltaTheta/norm(deltaTheta)*sin(norm(deltaTheta)/2)  ];
                    fixAtt( 1:4, iTime ) = quat_mult( xAtt(1:4, iTime ), quat_conj(deltaQuat) );
                    fixAtt( 4+3+1:end, iTime ) = aTmp(:,1);
                    isFixAtt( iTime, iMonteCarlo ) = norm( aTmp(:,1) - ambRef(n+1:end, iTime) ) == 0;
                end % End of solving ambiguities for Attitude - ESKF         
                
            end % End of Attitude - ESKF correction
            
 
            if useLCAtt % Correction step for Attitude - ESKF
                % Cycle slip repair
                if sum( cycleSlipRef(:,iTime) ) > 0
                    cycleSlipInd = find( cycleSlipRef(n+1:end,iTime) == 1 );
                    xLCAtt( 4+3+cycleSlipInd, iTime ) = ( y_att( cycleSlipInd ) - y_att( cycleSlipInd+n*N ) )/lambda;
                    for ind=1:length( cycleSlipInd )
                        PLCAtt( 4+3-1+cycleSlipInd(ind), :  ) = 0;
                        PLCAtt( :, 4+3-1+cycleSlipInd(ind)  ) = 0;
                        PLCAtt( 4+3-1+cycleSlipInd(ind),  4+3-1+cycleSlipInd(ind) ) =  10 * R_att( n*N+cycleSlipInd(ind), n*N+cycleSlipInd(ind) ); 
                    end
                end
                % Observation model
                baselineHat = nan(3*N,1);
                for iBaseline = 1:3:3*N
                    baselineHat( iBaseline:iBaseline+2 ) = quat_rotate( xLCAtt( 1:4, iTime ), baselineBodyHat( :, ceil(iBaseline/3) ) );
                end
                hx = [ kron( eye(N), DD*G ) * baselineHat( 1:3*N ) + lambda*xLCAtt( 4+3+1:end, iTime ); ...
                          kron( eye(N), DD*G ) * baselineHat( 1:3*N ) ]; 
                % Jacobian for the observation model       
                HQuat = [];
                for iBaseline = 1:N
                    Jq = inScript_derivative_Baseline_wrt_Quaternion( xLCAtt( 1:4, iTime ), baselineBodyHat(:,iBaseline) );
                    HQuat = [ HQuat; DD*G*Jq ]; % Getting the Jacobian related to the quaternion
                end
                HQuat = kron( [1;1], HQuat );
                HBiasGyro = zeros( 2*n*N, 3 );
                HAmb = [ lambda * eye(n*N); zeros(n*N) ];
                Hx = [ HQuat, HBiasGyro, HAmb ];
                % Jacobian for the error state
                HdeltaTheta = 1/2 * quat_productMatrix ( xLCAtt( 1:4, iTime ), 'left' ) * [zeros(1,3); eye(3)];
                HdeltaX = blkdiag( HdeltaTheta, eye(3),  eye(n*N) );
                % Jacobian for the total state
                H = Hx * HdeltaX;
                % Covariance matrix update
                S_inv = pinv( H*PLCAtt*H' + R_att );
                K = PLCAtt * H'* S_inv;
                %Constraints
                Delta = zeros(2*N*n, N);
                for antInd=1:N
                    Delta(n*(antInd - 1) + 1 : n*antInd, antInd) = DD*G*baselineHat(3*antInd - 2 : 3*antInd);
                    Delta((N + antInd - 1)*n + 1 : (N + antInd)*n, antInd) = DD*G*baselineHat(3*antInd - 2 : 3*antInd);
                end
                T = zeros(2*(N+1)*n, N);
                Gamma = - K * Delta;
                Psi = Delta' * S_inv * Delta;
                L = K + Gamma * pinv(Psi) * Delta' * S_inv;
                PLCAtt = ( eye(4+3+n*N-1) - K*H )*PLCAtt + Gamma * pinv(Psi) * Gamma';
                
                % Error state update
                errorState = L * ( y_att-hx );
                xLCAtt( 5:end, iTime ) = xLCAtt( 5:end, iTime ) + errorState(4:end);
                deltaTheta = errorState(1:3);
                deltaQuat = [ cos(norm(deltaTheta)/2); deltaTheta/norm(deltaTheta)*sin(norm(deltaTheta)/2) ];
                xLCAtt(1:4, iTime ) = quat_mult( xLCAtt(1:4, iTime ), deltaQuat );
                
                if solveAmbiguities == true
                    [ aTmp, rTmp ] = LAMBDA( xLCAtt(4+3+1:end, iTime), PLCAtt( 4+3:end, 4+3:end ) ); % Solving the ILS problem
                    % Fix solution estimation
                    fixLCAtt(:, iTime) = xLCAtt(:, iTime );
                    deltaErrorConditioned = PLCAtt( 1:4+3-1, 4+3:end )*pinv(PLCAtt( 4+3:end, 4+3:end )) * ( xLCAtt(4+3+1:end, iTime ) - aTmp(:,1) );
                    fixLCAtt( 5:7, iTime ) = xLCAtt( 5:7, iTime ) - deltaErrorConditioned(4:6);
                    deltaTheta = deltaErrorConditioned(1:3);
                    deltaQuat = [ cos(norm(deltaTheta)/2); deltaTheta/norm(deltaTheta)*sin(norm(deltaTheta)/2)  ];
                    fixLCAtt( 1:4, iTime ) = quat_mult( xLCAtt(1:4, iTime ), quat_conj(deltaQuat) );
                    fixLCAtt( 4+3+1:end, iTime ) = aTmp(:,1);
                    isFixLCAtt( iTime, iMonteCarlo ) = norm( aTmp(:,1) - ambRef(n+1:end, iTime) ) == 0;
                end % End of solving ambiguities for Attitude - LCESKF         
                
            end % End of Attitude - ESKF correction
            
        end % End of correction step
                
            
                        
            
        
    end % End of time loop

    % Saving results
    if useJpa % Results for JPA- ESKF
        errorJpa( 1:p+M, :, iMonteCarlo ) = xJpa - xRef;
        errorPosJpa( 1:3, :, iMonteCarlo ) = xJpa( 1:3, : ) - posRef;
        errorVelJpa( 1:3, :, iMonteCarlo ) = xJpa( 4:6, : ) - velRef;
        errorAmbJpa( 1:M, :, iMonteCarlo ) = xJpa( p+1:end, : ) - ambRef;
        for aux=1:nTime,  errorEulerJpa( 1:3, aux, iMonteCarlo ) = inScript_minimal_angular_error( quat_2_euler( xJpa( 7:10, aux ) ), eulerRef(:, aux) ); end
        errorPosFixJpa( 1:3, :, iMonteCarlo ) = fixJpa( 1:3, : ) - posRef;
        errorAmbFixJpa( 1:M, :, iMonteCarlo ) = fixJpa( p+1:end, : ) - ambRef;
        for aux=1:nTime,  errorEulerFixJpa( 1:3, aux, iMonteCarlo ) = inScript_minimal_angular_error( quat_2_euler( fixJpa( 7:10, aux ) ), eulerRef(:, aux) ); end
    end % End of results for JPA- ESKF
    
    if useLCJpa % Results for JPA- ESKF
        errorLCJpa( 1:p+M, :, iMonteCarlo ) = xLCJpa - xRef;
        errorPosLCJpa( 1:3, :, iMonteCarlo ) = xLCJpa( 1:3, : ) - posRef;
        errorVelLCJpa( 1:3, :, iMonteCarlo ) = xLCJpa( 4:6, : ) - velRef;
        errorAmbLCJpa( 1:M, :, iMonteCarlo ) = xLCJpa( p+1:end, : ) - ambRef;
        for aux=1:nTime,  errorEulerLCJpa( 1:3, aux, iMonteCarlo ) = inScript_minimal_angular_error( quat_2_euler( xLCJpa( 7:10, aux ) ), eulerRef(:, aux) ); end
        errorPosFixLCJpa( 1:3, :, iMonteCarlo ) = fixLCJpa( 1:3, : ) - posRef;
        errorAmbFixLCJpa( 1:M, :, iMonteCarlo ) = fixLCJpa( p+1:end, : ) - ambRef;
        for aux=1:nTime,  errorEulerFixLCJpa( 1:3, aux, iMonteCarlo ) = inScript_minimal_angular_error( quat_2_euler( fixLCJpa( 7:10, aux ) ), eulerRef(:, aux) ); end
    end % End of results for JPA- LCESKF
    
     if useRtk % Results for RTK- EKF
        errorPosRtk( 1:3, :, iMonteCarlo ) = xRtk( 1:3, : ) - posRef;
        errorVelRtk( 1:3, :, iMonteCarlo ) = xRtk( 4:6, : ) - velRef;
        errorAmbRtk( 1:n, :, iMonteCarlo ) = xRtk( 3+3+1:end, : ) - ambRef(1:n, :);
        errorPosFixRtk( 1:3, :, iMonteCarlo ) = fixRtk( 1:3, : ) - posRef;
        errorAmbFixRtk( 1:n, :, iMonteCarlo ) = fixRtk( 3+3+1:end, : ) - ambRef(1:n, :);
    end % Enf of results for RTK- EKF
    
    if useAtt % Results for Att- ESKF
        errorBiasAtt( 1:3, :, iMonteCarlo ) = xAtt( 5:7, : ) - biasGyroRef;
        errorAmbAtt( 1:n*N, :, iMonteCarlo ) = xAtt( 4+3+1:end, : ) - ambRef(n+1:end,:);
        for aux=1:nTime,  errorEulerAtt( 1:3, aux, iMonteCarlo ) = inScript_minimal_angular_error( quat_2_euler( xAtt( 1:4, aux ) ), eulerRef(:, aux) ); end
        errorAmbFixAtt( 1:n*N, :, iMonteCarlo ) = fixAtt( 4+3+1:end, : ) - ambRef(n+1:end,:);
        for aux=1:nTime,  errorEulerFixAtt( 1:3, aux, iMonteCarlo ) = inScript_minimal_angular_error( quat_2_euler( fixAtt( 1:4, aux ) ), eulerRef(:, aux) ); end
    end % End of results for Att- ESKF
    
    if useLCAtt % Results for Att- LCESKF
        errorBiasLCAtt( 1:3, :, iMonteCarlo ) = xLCAtt( 5:7, : ) - biasGyroRef;
        errorAmbLCAtt( 1:n*N, :, iMonteCarlo ) = xLCAtt( 4+3+1:end, : ) - ambRef(n+1:end,:);
        for aux=1:nTime,  errorEulerLCAtt( 1:3, aux, iMonteCarlo ) = inScript_minimal_angular_error( quat_2_euler( xLCAtt( 1:4, aux ) ), eulerRef(:, aux) ); end
        errorAmbFixLCAtt( 1:n*N, :, iMonteCarlo ) = fixLCAtt( 4+3+1:end, : ) - ambRef(n+1:end,:);
        for aux=1:nTime,  errorEulerFixLCAtt( 1:3, aux, iMonteCarlo ) = inScript_minimal_angular_error( quat_2_euler( fixLCAtt( 1:4, aux ) ), eulerRef(:, aux) ); end
    end % End of results for Att- LCESKF
    
    
end % End of Monte Carlo loop


%% RMSE Estimation
% Initialization of variables for RMSE analysis
if useJpa, rmsePosJpa=nan(nTime,1); rmsePosFixJpa=nan(nTime,1); rmseEulerJpa=nan(nTime,1); rmseEulerFixJpa=nan(nTime,1);  rmseBiasJpa=nan(nTime,1); masrJpa=nan(nTime, 1); end
if useLCJpa, rmsePosLCJpa=nan(nTime,1); rmsePosFixLCJpa=nan(nTime,1); rmseEulerLCJpa=nan(nTime,1); rmseEulerFixLCJpa=nan(nTime,1);  rmseBiasLCJpa=nan(nTime,1); masrLCJpa=nan(nTime, 1); end
if useRtk, rmsePosRtk=nan(nTime,1); rmsePosFixRtk=nan(nTime,1); masrRtk=nan(nTime, 1); end
if useAtt, rmseEulerAtt=nan(nTime,1); rmseEulerFixAtt=nan(nTime,1);  rmseBiasAtt=nan(nTime,1);  masrAtt=nan(nTime, 1); end
if useLCAtt, rmseEulerLCAtt=nan(nTime,1); rmseEulerFixLCAtt=nan(nTime,1);  rmseBiasLCAtt=nan(nTime,1);  masrLCAtt=nan(nTime, 1); end

for iTime = 1:nTime
    
    if useJpa
        tmp = vecnorm( reshape( errorPosJpa (1:3, iTime, : ), nMonteCarlo, 3 ) , 2, 2 ).^2;
        rmsePosJpa( iTime ) = sqrt(  mean( tmp(:), 'omitnan' )  );
        tmp = vecnorm( reshape( errorPosFixJpa(1:3, iTime, : ), nMonteCarlo, 3 ), 2,  2 ).^2;
        rmsePosFixJpa( iTime ) = sqrt(  mean( tmp(:), 'omitnan' )  );
        tmp = vecnorm( reshape( errorEulerJpa (1:3, iTime, : ), nMonteCarlo, 3 ), 2, 2 ).^2;
        rmseEulerJpa( iTime ) = sqrt(  mean( tmp(:), 'omitnan' )  );
        tmp = vecnorm( reshape( errorEulerFixJpa(1:3, iTime, : ), nMonteCarlo, 3 ) , 2, 2 ).^2;
        rmseEulerFixJpa( iTime ) = sqrt(  mean( tmp(:), 'omitnan' )  );
        tmp = vecnorm( reshape( errorJpa (11:13, iTime, : ), nMonteCarlo, 3 ) , 2, 2 ).^2;
        rmseBiasJpa( iTime ) = sqrt(  mean( tmp(:) )  );
        masrJpa( iTime ) = mean( isFixJpa(iTime, :), 'omitnan' );
    end
    
    if useLCJpa
        tmp = vecnorm( reshape( errorPosLCJpa (1:3, iTime, : ), nMonteCarlo, 3 ) , 2, 2 ).^2;
        rmsePosLCJpa( iTime ) = sqrt(  mean( tmp(:), 'omitnan' )  );
        tmp = vecnorm( reshape( errorPosFixLCJpa(1:3, iTime, : ), nMonteCarlo, 3 ), 2,  2 ).^2;
        rmsePosFixLCJpa( iTime ) = sqrt(  mean( tmp(:), 'omitnan' )  );
        tmp = vecnorm( reshape( errorEulerLCJpa (1:3, iTime, : ), nMonteCarlo, 3 ), 2, 2 ).^2;
        rmseEulerLCJpa( iTime ) = sqrt(  mean( tmp(:), 'omitnan' )  );
        tmp = vecnorm( reshape( errorEulerFixLCJpa(1:3, iTime, : ), nMonteCarlo, 3 ) , 2, 2 ).^2;
        rmseEulerFixLCJpa( iTime ) = sqrt(  mean( tmp(:), 'omitnan' )  );
        tmp = vecnorm( reshape( errorLCJpa (11:13, iTime, : ), nMonteCarlo, 3 ) , 2, 2 ).^2;
        rmseBiasLCJpa( iTime ) = sqrt(  mean( tmp(:) )  );
        masrLCJpa( iTime ) = mean( isFixLCJpa(iTime, :), 'omitnan' );
    end
    
    if useRtk
        tmp = vecnorm( reshape( errorPosRtk (1:3, iTime, : ), nMonteCarlo, 3 ) , 2, 2 ).^2;
        rmsePosRtk( iTime ) = sqrt(  mean( tmp(:), 'omitnan' )  );
        tmp = vecnorm( reshape( errorPosFixRtk(1:3, iTime, : ), nMonteCarlo, 3 ) , 2, 2 ).^2;
        rmsePosFixRtk( iTime ) = sqrt(  mean( tmp(:), 'omitnan' )  );
        masrRtk( iTime ) = mean( isFixRtk(iTime, :), 'omitnan' );
    end
    
    if useAtt
        tmp = vecnorm( reshape( errorEulerAtt (1:3, iTime, : ), nMonteCarlo, 3 ) , 2, 2 ).^2;
        rmseEulerAtt( iTime ) = sqrt(  mean( tmp(:), 'omitnan' )  );
        tmp = vecnorm( reshape( errorEulerFixAtt(1:3, iTime, : ), nMonteCarlo, 3 ) , 2, 2 ).^2;
        rmseEulerFixAtt( iTime ) = sqrt(  mean( tmp(:), 'omitnan' )  );
        tmp = vecnorm( reshape( errorBiasAtt (1:3, iTime, : ), nMonteCarlo, 3 ) , 2, 2 ).^2;
        rmseBiasAtt( iTime ) = sqrt(  mean( tmp(:), 'omitnan' )  );
        masrAtt( iTime ) = mean( isFixAtt(iTime, :), 'omitnan' );
    end

    if useLCAtt
        tmp = vecnorm( reshape( errorEulerLCAtt (1:3, iTime, : ), nMonteCarlo, 3 ) , 2, 2 ).^2;
        rmseEulerLCAtt( iTime ) = sqrt(  mean( tmp(:), 'omitnan' )  );
        tmp = vecnorm( reshape( errorEulerFixLCAtt(1:3, iTime, : ), nMonteCarlo, 3 ) , 2, 2 ).^2;
        rmseEulerFixLCAtt( iTime ) = sqrt(  mean( tmp(:), 'omitnan' )  );
        tmp = vecnorm( reshape( errorBiasLCAtt (1:3, iTime, : ), nMonteCarlo, 3 ) , 2, 2 ).^2;
        rmseBiasLCAtt( iTime ) = sqrt(  mean( tmp(:), 'omitnan' )  );
        masrLCAtt( iTime ) = mean( isFixLCAtt(iTime, :), 'omitnan' );
    end
    
end


%% Visualization
set(groot,'defaultFigureColor','w')
set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',13)

if enablePlot == true
    
    % Illustrating here the complete reference solution
    figure; hold on; grid on; box on; axis equal; % Reference trajectory
    plot( posRef(1,:).', posRef(2,:).', '-', 'LineWidth', 1.5 );
    plot( posBaseRef(1), posBaseRef(2), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r' );
    xlabel('East (m)'); ylabel('North (m)');
    figure; hold on; grid on; box on;  % Reference velocity over time
    plot( velRef(1,:).', '-', 'LineWidth', 1.5 ); plot( velRef(2,:).', '-', 'LineWidth', 1.5 ); plot( velRef(3,:).', '-', 'LineWidth', 1.5 );
    xlabel('time (s)'); ylabel('speed (m/s)'); legend('$v_e$','$v_n$','$v_u$' );
    figure; hold on; grid on; box on;  % Gyroscope biases over time
    plot( biasGyroRef(1,:).', '-', 'LineWidth', 1.5 ); plot( biasGyroRef(2,:).', '-', 'LineWidth', 1.5 ); plot( biasGyroRef(3,:).', '-', 'LineWidth', 1.5 );
    xlabel('time (s)'); ylabel('bias (rad/s)'); legend('$b_{\omega,x}$','$v_{\omega,y}$','$v_{\omega,z}$' );
    figure; hold on; grid on; box on; % Ambiguities over time
    plot( ambRef.', '-' );
    xlabel('time (s)'); ylabel('ambiguities (cycles)');
    figure; hold on; grid on; box on;  % Euler angles (orientation) over time
    plot( 180/pi*eulerRef(1,:).', '-', 'LineWidth', 1.5 ); plot( 180/pi*eulerRef(2,:).', '-', 'LineWidth', 1.5 ); plot( 180/pi*eulerRef(3,:).', '-', 'LineWidth', 1.5 );
    xlabel('time (s)'); ylabel('attitude (deg)'); legend('$\theta_{x}$ (roll)','$\theta_{x}$ (pitch)','$\theta_{z}$ (yaw)' );
    figure; hold on; grid on; box on; %  Angular rate over time
    plot( 180/pi*angularRateRef(1,:).', '-', 'LineWidth', 1.5 ); plot( 180/pi*angularRateRef(2,:).', '-', 'LineWidth', 1.5 ); plot( 180/pi*angularRateRef(3,:).', '-', 'LineWidth', 1.5 );
    xlabel('time (s)'); ylabel('angular rate (deg/s)'); legend('$\omega_{x}$','$\omega_{x}$','$\omega_{z}$' );
    
    % Plot reference with error ellipses
    if useJpa
        figure; hold on; grid on; box on; axis equal; xlabel('East [m]'); ylabel('North [m]');
        plot( posRef(1,:).', posRef(2,:).', 'k-', 'LineWidth', 1.5 );
        plot( posBaseRef(1), posBaseRef(2), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r' );
        confidenceLevel = 90; % 9973 % 68 % 90 %% 99
        for iTime=2:20:nTime
            plot( posRef(1,iTime), posRef(2,iTime), 'kx','MarkerSize',12,'LineWidth',2);
            errorXY = reshape( errorPosJpa( 1:2, iTime, : ), nMonteCarlo, 2  );
            ellipse90 = errorEllipseFromData( errorXY, confidenceLevel );
            h1=plot( posRef(1,iTime)+ellipse90(:,1), posRef(2,iTime)+ellipse90(:,2), '-', 'LineWidth', 1.50, 'Color',rgb('DarkBlue'));
            errorXY = reshape( errorPosFixJpa( 1:2, iTime, : ), nMonteCarlo, 2  );
            ellipse90 = errorEllipseFromData( errorXY, confidenceLevel );
            h2=plot( posRef(1,iTime)+ellipse90(:,1), posRef(2,iTime)+ellipse90(:,2), '-', 'LineWidth', 1.50, 'Color',rgb('DarkOliveGreen'));
        end
    end
    
    % Plot reference with error ellipses with constraints
    if useLCJpa
        figure; hold on; grid on; box on; axis equal; xlabel('East [m]'); ylabel('North [m]');
        plot( posRef(1,:).', posRef(2,:).', 'k-', 'LineWidth', 1.5 );
        plot( posBaseRef(1), posBaseRef(2), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r' );
        confidenceLevel = 90; % 9973 % 68 % 90 %% 99
        for iTime=2:20:nTime
            plot( posRef(1,iTime), posRef(2,iTime), 'kx','MarkerSize',12,'LineWidth',2);
            errorXY = reshape( errorPosLCJpa( 1:2, iTime, : ), nMonteCarlo, 2  );
            ellipse90 = errorEllipseFromData( errorXY, confidenceLevel );
            h1=plot( posRef(1,iTime)+ellipse90(:,1), posRef(2,iTime)+ellipse90(:,2), '-', 'LineWidth', 1.50, 'Color',rgb('DarkBlue'));
            errorXY = reshape( errorPosFixLCJpa( 1:2, iTime, : ), nMonteCarlo, 2  );
            ellipse90 = errorEllipseFromData( errorXY, confidenceLevel );
            h2=plot( posRef(1,iTime)+ellipse90(:,1), posRef(2,iTime)+ellipse90(:,2), '-', 'LineWidth', 1.50, 'Color',rgb('DarkOliveGreen'));
        end
    end
    
    % Estimators Position RMSE
    figure; hold on; grid on; box on; xlim([1 nTime]); leg={};
    if useJpa,  ax=gca; ax.ColorOrderIndex=1; leg = [leg; 'JPA - Float'];  plot( rmsePosJpa, '-', 'LineWidth', 1.5 ); end
    if useJpa, ax=gca; ax.ColorOrderIndex=1;  leg = [leg; 'JPA - Fix'];  plot( rmsePosFixJpa, '--', 'LineWidth', 1.5 ); end
    if useLCJpa,  ax=gca; ax.ColorOrderIndex=4; leg = [leg; 'LCJPA - Float'];  plot( rmsePosLCJpa, '-', 'LineWidth', 1.5 ); end
    if useLCJpa, ax=gca; ax.ColorOrderIndex=4;  leg = [leg; 'LCJPA - Fix'];  plot( rmsePosFixLCJpa, '--', 'LineWidth', 1.5 ); end
    if useRtk,  ax=gca; ax.ColorOrderIndex=2; leg = [leg; 'RTK - Float'];  plot( rmsePosRtk, '-', 'LineWidth', 1.5); end
    if useRtk,  ax=gca; ax.ColorOrderIndex=2; leg = [leg; 'RTK - Fix'];  plot( rmsePosFixRtk, '--', 'LineWidth', 1.5 ); end
    legg = legend( leg, 'FontSize', 16 );
    xlabel('time (s)'); ylabel('RMSE$_{\mathbf{p} }$ (m)');
    set(gca,'YScale','log');
    
    % Estimators Attitude RMSE
    figure; hold on; grid on; box on; xlim([1 nTime]); leg={};
    if useJpa,   leg = [leg; 'JPA - Float']; ax=gca; ax.ColorOrderIndex=1; plot( rmseEulerJpa*180/pi, '-', 'LineWidth', 1.5 ); end
    if useJpa,   leg = [leg; 'JPA - Fix']; ax=gca; ax.ColorOrderIndex=1; plot( rmseEulerFixJpa*180/pi, '--' ); end
    if useLCJpa,   leg = [leg; 'LCJPA - Float']; ax=gca; ax.ColorOrderIndex=5; plot( rmseEulerLCJpa*180/pi, '-', 'LineWidth', 1.5 ); end
    if useLCJpa,   leg = [leg; 'LCJPA - Fix']; ax=gca; ax.ColorOrderIndex=5; plot( rmseEulerFixLCJpa*180/pi, '--' ); end
    if useAtt,   leg = [leg; 'Att - Float']; ax=gca; ax.ColorOrderIndex=3;  plot( rmseEulerAtt*180/pi, '-', 'LineWidth', 1.5 ); end
    if useAtt,   leg = [leg; 'Att - Fix']; ax=gca; ax.ColorOrderIndex=3;  plot( rmseEulerFixAtt*180/pi, '--', 'LineWidth', 1.5 ); end
    if useLCAtt,   leg = [leg; 'LCAtt - Float']; ax=gca; ax.ColorOrderIndex=4;  plot( rmseEulerLCAtt*180/pi, '-', 'LineWidth', 1.5 ); end
    if useLCAtt,   leg = [leg; 'LCAtt - Fix']; ax=gca; ax.ColorOrderIndex=4;  plot( rmseEulerFixLCAtt*180/pi, '--', 'LineWidth', 1.5 ); end
    legg = legend( leg, 'FontSize', 16 );
    ylabel('$\| \bar{\theta} \|_{rmse}$ ($^\circ$)')
    xlabel('time (s)'); ylabel('RMSE$_{\theta }$ (deg)');
    set(gca,'YScale','log');
    
    % Bias errors against covariance
    
    % Mean ambiguity success rate (MASR)
    figure; hold on; grid on; box on; xlim([2 nTime]); leg={};
    if useJpa,   leg = [leg; 'JPA'];  plot( mean( isFixJpa, 2, 'omitnan' )*100, 'LineWidth', 2 ); end
    if useRtk,   leg = [leg; 'RTK'];  plot( mean( isFixRtk, 2, 'omitnan' )*100, 'LineWidth', 2 ); end
    if useAtt,   leg = [leg; 'Att'];  plot( mean( isFixAtt, 2, 'omitnan' )*100, 'LineWidth', 2 ); end
    if useLCAtt,   leg = [leg; 'LCAtt'];  plot( mean( isFixLCAtt, 2, 'omitnan' )*100, 'LineWidth', 2 ); end
    if useLCJpa,   leg = [leg; 'LCJPA'];  plot( mean( isFixLCJpa, 2, 'omitnan' )*100, 'LineWidth', 2 ); end
    legg = legend( leg, 'FontSize', 13 );
    xlabel('time (s)'); ylabel('MASR (\%)');
    
end

%% Additional functions
function Jq = inScript_derivative_Baseline_wrt_Quaternion( q, b )
    Jq =  2 * [ q(1)*b + cross(q(2:4), b), q(2:4)'*b*eye(3) + q(2:4)*b' - b*q(2:4)' - q(1)*skewMatrix(b)  ];
end

function XiMatrix = inScript_get_Xi_from_quaternion( q )
    XiMatrix = [ -q(2:4).'; q(1)*eye(3) + skewMatrix(q(2:4))  ];
end

function baselineBody = inScript_random_antenna_configuration(Nbaselines,baselineLength)
    if Nbaselines<2
        error('Not enough number of antennas!')
    end
    baselineBody = randn(3, Nbaselines);
    baselineBody = baselineLength * baselineBody./vecnorm(baselineBody, 2, 1);
    dot_products = baselineBody' * baselineBody;
    count = 0;
    while any( abs(dot_products(~eye(Nbaselines))) > 0.9 )
        baselineBody = randn(3, Nbaselines);
        baselineBody = baselineLength * baselineBody./vecnorm(baselineBody, 2, 1);
        dot_products = baselineBody' * baselineBody;
        count = count + 1;
        if count > 3
            break
        end
    end
end

function [ R_Jpa, R_Rtk, R_att ] = inScript_jpa_noise_model( sigmaPhase, sigmaCode, distance, svElev, n, N, DD )
    DD_Jpa = ones(N+1)+eye(N+1); 
    W = diag( (1 + 10.*exp( -svElev./(10*pi/180) )).^2 );
    R_iono = W * (sqrt(2)*0.4*distance/1000/1000)^2; % The ionosphere weighted model is \sigma_I = sqrt(2)*0.4 mm/km and amplified with the elevation dependent weighting matrix
    R_phase = W * sigmaPhase^2;
    R_code = W * sigmaCode^2;
    R_Jpa = blkdiag( kron(  DD_Jpa, DD * R_phase * DD.'), kron(  DD_Jpa, DD * R_code * DD.') );
    R_Jpa( 1:n, 1:n ) = R_Jpa( 1:n, 1:n ) + DD*R_iono*DD.';
    R_Jpa( n*(N+1)+1:n*(N+1)+n, n*(N+1)+1:n*(N+1)+n ) = R_Jpa( n*(N+1)+1:n*(N+1)+n, n*(N+1)+1:n*(N+1)+n ) + DD*R_iono*DD.';
    R_Jpa( 1:n, n*(N+1)+1:n*(N+1)+n ) = DD*R_iono*DD.';
    R_Jpa( n*(N+1)+1:n*(N+1)+n, 1:n ) = DD*R_iono*DD.';
    R_Rtk = blkdiag( DD * (R_phase+R_iono) * DD.', DD * (R_code+R_iono) * DD.' );
    R_att = blkdiag( kron(  1+eye(N), DD * R_phase * DD.'), kron(  1+eye(N), DD * R_code * DD.') );
end

function qNew = inScript_quat_integrate_omega( qOld, omega, dt )
    deltaTheta = norm( omega ) * dt / 2;
    deltaQuat = [ cos( deltaTheta ); gyroCorr/norm(gyroCorr) * sin(deltaTheta) ];
    qNew = quat_mult( qOld, deltaQuat   );
end

function errorAngular = inScript_minimal_angular_error( eulerHat, eulerRef )
    aux = eulerHat - eulerRef;
    aux = sign(aux).* [ min( 2*pi - abs(aux(1)), abs( aux(1) ) );  min( 2*pi - abs(aux(2)), abs( aux(2) ) ); min( 2*pi - abs(aux(3)), abs( aux(3) ) ) ]; % To find the shortest angular error
    errorAngular = aux;
end