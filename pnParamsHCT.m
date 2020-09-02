function [params] = pnParamsHCT(~)
% Requires V2Struct
% Generates parameters structure for pnDrift
% For use with pindrift

if nargin>0
   
else 
    %Physical constants
    kB = 8.6173324e-5;    % Boltzmann constant [eV K^-1]
    T = 277;              % Temperature [K]
    epp0 = 552434;        % e^2 eV^-1 cm^-1 -Checked (02-11-15)
    q = 1;                % in e
    e = 1.61917e-19;      % Charge of an electron in Coulombs for current calculations

    % General Parameters
    Int =3;             % Bias Light intensity
    tmax = 5e-5;        % Time
    pulseon =0;         % Switch pulse on for TPV
    pulselen = 2e-10;    % Transient pulse length
    tstart=1e-10;
    pulseint =5;        % Transient pulse intensity- for BL and TM models, 100 mW/cm2 assumed
    Vapp =0;            % Applied bias
    Vtransient=0.01;
    wAC=1e3;
    fastrec = 0;        % Can be used to accelerate finding initial conditions
    BC =4;                  % Boundary Conditions. Must be set to one for first solution % BC=5 Impedance measurement
    figson =1;              % Toggle figures on/off 
    side = 1;                % illumination side 1 = EE, 2 = SE
    calcJ =1;              % Calculates Currents- slows down solving calcJ = 1, calculates DD currents at every position, calcJ = 2, calculates DD at boundary.
    mobset =1;             % Switch on/off electron hole mobility- MUST BE SET TO ZERO FOR INITIAL SOLUTION
    symm=0;
    equilibrium=0;
    discretetrap=0;
    % OM = Optical Model 
    % 0 = Uniform Generation 
    % 1 = Beer-Lamber (Requires pre calculation using Igor code & gen profile in base workspace)
    % 2 = Transfer Matrix (Stanford)
    OM  =0;
    % Parameters for time point outputs
    tmesh_type = 2;  %%% 1 means linear spaced time, 2 is logaritmic
    t0 = tmax/1000;
    tpoints = 1000;
    AbsTol=1e-6;
    RelTol=1e-3;
    % Material Properties
   
    [layer,layers]=readlayersCT('layerdata.xlsx');
    
    
    % SRH parameters
    se = 1e-15;             % [cm^2] Electron capture cross-section
    v = 1e9;                % [cm/s] Carrier group velocity. Estimates from: http://www.tf.uni-kiel.de/matwis/amat/semi_en/kap_2/advanced/t2_3_1.html


    % Langevin trapping parameters
    cn = 6e-7;
    cp = 6e-7;
    Tn = 1e-3;

    % Define geometry of the system (m=0 1D, m=1 cylindrical polar coordinates,
    % m=2 spherical polar coordinates).
    m = 0; 

   
    varlist = who('*')'; %%%Basically orders the variable names in a cell array and the asterisc takes them all
    varstr = strjoin(varlist, ','); %%%The previous array now is a string and it is separated by commas

    % Pack parameters in to structure 'params'
    varcell = who('*')';                    % Store variables names in cell array
    varcell = ['fieldnames', varcell];      % adhere to syntax for v2struct

    params = v2struct(varcell);
    toz=hstruct(params);
    params=toz.data;
end 
end