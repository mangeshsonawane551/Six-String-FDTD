%-------------------------------------------------------------------------%
%_*Program Description*_  
%This program designs a script that calls  stiff string FDTD function for 
% six different strings, in order to build a short music note piece of a 
%six-string EADGBE-tuned PHOSPHOR BRONZE STRING guitar. The output of each
%string is added in the output vector 'piece_music'. The output is generated in 
%mono.
%-------------------------------------------------------------------------%
clear all
dbstop if error

%%%%% options
opts.plot_on = false;
opts.useforloop = true;
opts.add_stiffness = true;
opts.input_type = 'plucked';
opts.output_type = 'displacement';
opts.bctype = 'clamped';
%opts.bctype = 'simply_supported';

%-------------------------------------------------------------------------%
% E
%-------------------------------------------------------------------------%

%%%%% physical string parameters 
phys_param.T = 119.73e1;                 % tension (N)
phys_param.r = 0.0003429;                % string radius (m)
phys_param.rho = 8890;                   % density (kg/m^3)
phys_param.T60 = 2;                      % T60 (s)
phys_param.L = 1;                        % length (m)
phys_param.E = 1.2e11;                   % Young's modulus (Pa)

%%%%% simulation parameters 
sim_param.SR = 44100*1;                % sample rate (Hz)
sim_param.Tf = 1;                      % duration of simulation (s)
sim_param.xi = 0.1;                    % coordinate of excitation (normalised, 0-1)
sim_param.famp = 1;                    % peak amplitude of excitation (N)
sim_param.dur = 0.009;                 % duration of excitation (s)
sim_param.exc_st = 0.09;               % start time of excitation (s)
sim_param.xo = 0.9;                    % coordinate of output (normalised, 0-1)

y1 = string_fdtd_s1889125(opts,phys_param,sim_param);
%-------------------------------------------------------------------------%
% A
%-------------------------------------------------------------------------%
%%%%% physical string parameters 
phys_param.T = 104.14e1;               % tension (N)
phys_param.r = 0.00023;                % string radius (m)
phys_param.rho = 8890;                 % density (kg/m^3)
phys_param.T60 = 2;                    % T60 (s)
phys_param.L = 1;                      % length (m)
phys_param.E = 1.2e11;                 % Young's modulus (Pa)

%%%%% simulation parameters 
sim_param.SR = 44100*1;                % sample rate (Hz)
sim_param.Tf = 1;                      % duration of simulation (s)
sim_param.xi = 0.2;                    % coordinate of excitation (normalised, 0-1)
sim_param.famp = 1;                    % peak amplitude of excitation (N)
sim_param.dur = 0.009;                 % duration of excitation (s)
sim_param.exc_st = 0.01;               % start time of excitation (s)
sim_param.xo = 0.9;                    % coordinate of output (normalised, 0-1)

y2 = string_fdtd_s1889125(opts,phys_param,sim_param);

%-------------------------------------------------------------------------%
% D
%-------------------------------------------------------------------------%
%%%%% physical string parameters 
phys_param.T = 82.47e1;                 % tension (N)
phys_param.r = 0.00015;                 % string radius (m)
phys_param.rho = 8890;                  % density (kg/m^3)
phys_param.T60 = 2;                     % T60 (s)
phys_param.L = 1;                       % length (m)
phys_param.E = 1.2e11;                  % Young's modulus (Pa)

%%%%% simulation parameters 
sim_param.SR = 44100*1;                % sample rate (Hz)
sim_param.Tf = 1;                      % duration of simulation (s)
sim_param.xi = 0.1;                    % coordinate of excitation (normalised, 0-1)
sim_param.famp = 1;                    % peak amplitude of excitation (N)
sim_param.dur = 0.007;                 % duration of excitation (s)
sim_param.exc_st = 0.01;               % start time of excitation (s)
sim_param.xo = 0.9;                    % coordinate of output (normalised, 0-1)

y3 = string_fdtd_s1889125(opts,phys_param,sim_param);

%-------------------------------------------------------------------------%
% G
%-------------------------------------------------------------------------%
%%%%% physical string parameters 
phys_param.T = 82.87e1;                % tension (N)
phys_param.r = 0.000114;               % string radius (m)
phys_param.rho = 8890;                 % density (kg/m^3)
phys_param.T60 = 2;                    % T60 (s)
phys_param.L = 1;                      % length (m)
phys_param.E = 1.2e11;                 % Young's modulus (Pa)

%%%%% simulation parameters 
sim_param.SR = 44100*1;                % sample rate (Hz)
sim_param.Tf = 1;                      % duration of simulation (s)
sim_param.xi = 0.1;                    % coordinate of excitation (normalised, 0-1)
sim_param.famp = 1;                    % peak amplitude of excitation (N)
sim_param.dur = 0.008;                 % duration of excitation (s)
sim_param.exc_st = 0.01;               % start time of excitation (s)
sim_param.xo = 0.9;                    % coordinate of output (normalised, 0-1)

y4 = string_fdtd_s1889125(opts,phys_param,sim_param);

%-------------------------------------------------------------------------%
% B
%-------------------------------------------------------------------------%
%%%%% physical string parameters 
phys_param.T = 79.33e1;                % tension (N)
phys_param.r = 0.00018;                % string radius (m)
phys_param.rho = 8890;                 % density (kg/m^3)
phys_param.T60 = 2;                    % T60 (s)
phys_param.L = 1;                      % length (m)
phys_param.E = 1.2e11;                 % Young's modulus (Pa)

%%%%% simulation parameters 
sim_param.SR = 44100*1;                % sample rate (Hz)
sim_param.Tf = 1;                      % duration of simulation (s)
sim_param.xi = 0.1;                    % coordinate of excitation (normalised, 0-1)
sim_param.famp = 1;                    % peak amplitude of excitation (N)
sim_param.dur = 0.005;                 % duration of excitation (s)
sim_param.exc_st = 0.01;               % start time of excitation (s)
sim_param.xo = 0.9;                    % coordinate of output (normalised, 0-1)

y5 = string_fdtd_s1889125(opts,phys_param,sim_param);

%-------------------------------------------------------------------------%
% E
%-------------------------------------------------------------------------%
%%%%% physical string parameters 
phys_param.T = 72.179e1;               % tension (N)
phys_param.r = 0.000125;               % string radius (m)
phys_param.rho = 8890;                 % density (kg/m^3)
phys_param.T60 = 2;                    % T60 (s)
phys_param.L = 1;                      % length (m)
phys_param.E = 1.2e11;                 % Young's modulus (Pa)

%%%%% simulation parameters 
sim_param.SR = 44100*1;                % sample rate (Hz)
sim_param.Tf = 1;                      % duration of simulation (s)
sim_param.xi = 0.1;                    % coordinate of excitation (normalised, 0-1)
sim_param.famp = 1;                    % peak amplitude of excitation (N)
sim_param.dur = 0.008;                 % duration of excitation (s)
sim_param.exc_st = 0.01;               % start time of excitation (s)
sim_param.xo = 0.9;                    % coordinate of output (normalised, 0-1)

y6 = string_fdtd_s1889125(opts,phys_param,sim_param);

%-------------------------------------------------------------------------%
% MONO Output
%-------------------------------------------------------------------------%

%G major chord
Gmaj=y3+y4+y5;

%D major
Dmaj=y3+y2;

%Em chord
Em=y4+y5+y6;

%Cmajor chord
Cmaj= y6+y4;

%Music Piece vector
piece_music = cat(1,y3,y5,Gmaj,y3,y2,Dmaj,y4,y5,y6,Em,Dmaj,Cmaj,Gmaj,Dmaj);

%Normalised output
piece_music = piece_music/max(abs(piece_music));

soundsc(piece_music,sim_param.SR);
