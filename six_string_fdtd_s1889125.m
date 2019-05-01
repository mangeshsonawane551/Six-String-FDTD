%-------------------------------------------------------------------------%
% _*Program Description*_ 
%This program designs a script that calls   
% stiff string FDTD function for six different strings, in order to
%six-string EADGBE-tuned steel-string guitar. The output of each string is
% build a simplified model of a  added in the output vector 'total'. 
%The output is generated in mono as well as in stereo. 
%Panning has also been implemented for the output from left to right 
%with ratio of -1 to 1.
%-------------------------------------------------------------------------%
clear all
dbstop if error

%%%%% options
opts.plot_on = false;
opts.useforloop = true;
opts.add_stiffness = true;
opts.input_type = 'plucked';
opts.output_type = 'displacement';
%opts.bctype = 'clamped';
opts.bctype = 'simply_supported';

%-------------------------------------------------------------------------%
% E
%-------------------------------------------------------------------------%

%%%%% physical string parameters 
phys_param.T = 85.592e1;               % tension (N)
phys_param.r = 0.0005969;              % string radius (m)
phys_param.rho = 7850;                 % density (kg/m^3)
phys_param.T60 = 2;                    % T60 (s)
phys_param.L = 1;                      % length (m)
phys_param.E = 2e11;                   % Young's modulus (Pa)

%%%%% simulation parameters 
sim_param.SR = 44100*1;                % sample rate (Hz)
sim_param.Tf = 1;                      % duration of simulation (s)
sim_param.xi = 0.1;                    % coordinate of excitation (normalised, 0-1)
sim_param.famp = 1;                    % peak amplitude of excitation (N)
sim_param.dur = 0.001;                 % duration of excitation (s)
sim_param.exc_st = 0.01;               % start time of excitation (s)
sim_param.xo = 0.9;                    % coordinate of output (normalised, 0-1)

y1 = string_fdtd_s1889125(opts,phys_param,sim_param);

% Stereo generation
y1L=y1;
y1R=y1*0;
y1Stereo=[y1L y1R];


%-------------------------------------------------------------------------%
% A
%-------------------------------------------------------------------------%
%%%%% physical string parameters 
phys_param.T = 99.635e1;               % tension (N)
phys_param.r = 0.0004826;              % string radius (m)
phys_param.rho = 7850;                 % density (kg/m^3)
phys_param.T60 = 2;                    % T60 (s)
phys_param.L = 1;                      % length (m)
phys_param.E = 2e11;                   % Young's modulus (Pa)

%%%%% simulation parameters 
sim_param.SR = 44100*1;                % sample rate (Hz)
sim_param.Tf = 1;                      % duration of simulation (s)
sim_param.xi = 0.1;                    % coordinate of excitation (normalised, 0-1)
sim_param.famp = 1;                    % peak amplitude of excitation (N)
sim_param.dur = 0.001;                 % duration of excitation (s)
sim_param.exc_st = 0.01;               % start time of excitation (s)
sim_param.xo = 0.9;                    % coordinate of output (normalised, 0-1)

y2 = string_fdtd_s1889125(opts,phys_param,sim_param);

% Stereo generation
y2L=y2 * sqrt(2)/2;
y2R=y2 * sqrt(2)/2;
y2Stereo=[y2L y2R];


%-------------------------------------------------------------------------%
% D
%-------------------------------------------------------------------------%
%%%%% physical string parameters 
phys_param.T = 94.14e1;                % tension (N)
phys_param.r = 0.0003556;              % string radius (m)
phys_param.rho = 7850;                 % density (kg/m^3)
phys_param.T60 = 2;                    % T60 (s)
phys_param.L = 1;                      % length (m)
phys_param.E = 2e11;                   % Young's modulus (Pa)

%%%%% simulation parameters 
sim_param.SR = 44100*1;                % sample rate (Hz)
sim_param.Tf = 1;                      % duration of simulation (s)
sim_param.xi = 0.1;                    % coordinate of excitation (normalised, 0-1)
sim_param.famp = 1;                    % peak amplitude of excitation (N)
sim_param.dur = 0.001;                 % duration of excitation (s)
sim_param.exc_st = 0.01;               % start time of excitation (s)
sim_param.xo = 0.9;                    % coordinate of output (normalised, 0-1)

y3 = string_fdtd_s1889125(opts,phys_param,sim_param);

% Stereo generation
y3L=y3*0;
y3R=y3;
y3Stereo=[y3L y3R];


%-------------------------------------------------------------------------%
% G
%-------------------------------------------------------------------------%
%%%%% physical string parameters 
phys_param.T = 122.28e1;               % tension (N)
phys_param.r = 0.0002921;              % string radius (m)
phys_param.rho = 7850;                 % density (kg/m^3)
phys_param.T60 = 2;                    % T60 (s)
phys_param.L = 1;                      % length (m)
phys_param.E = 2e11;                   % Young's modulus (Pa)

%%%%% simulation parameters 
sim_param.SR = 44100*1;                % sample rate (Hz)
sim_param.Tf = 1;                      % duration of simulation (s)
sim_param.xi = 0.1;                    % coordinate of excitation (normalised, 0-1)
sim_param.famp = 1;                    % peak amplitude of excitation (N)
sim_param.dur = 0.001;                 % duration of excitation (s)
sim_param.exc_st = 0.01;               % start time of excitation (s)
sim_param.xo = 0.9;                    % coordinate of output (normalised, 0-1)

y4 = string_fdtd_s1889125(opts,phys_param,sim_param);

% Stereo generation
y4L=y4 * sqrt(2)/2;
y4R=y4 * sqrt(2)/2;
y4Stereo=[y4L y4R];


%-------------------------------------------------------------------------%
% B
%-------------------------------------------------------------------------%
%%%%% physical string parameters 
phys_param.T = 79.33e1;                % tension (N)
phys_param.r = 0.00018;                % string radius (m)
phys_param.rho = 7850;                 % density (kg/m^3)
phys_param.T60 = 2;                    % T60 (s)
phys_param.L = 1;                      % length (m)
phys_param.E = 2e11;                   % Young's modulus (Pa)

%%%%% simulation parameters 
sim_param.SR = 44100*1;                % sample rate (Hz)
sim_param.Tf = 1;                      % duration of simulation (s)
sim_param.xi = 0.1;                    % coordinate of excitation (normalised, 0-1)
sim_param.famp = 1;                    % peak amplitude of excitation (N)
sim_param.dur = 0.001;                 % duration of excitation (s)
sim_param.exc_st = 0.01;               % start time of excitation (s)
sim_param.xo = 0.9;                    % coordinate of output (normalised, 0-1)

y5 = string_fdtd_s1889125(opts,phys_param,sim_param);

% Stereo generation
y5L=y5;
y5R=y5*0;
y5Stereo=[y5L y5R];



%-------------------------------------------------------------------------%
% E
%-------------------------------------------------------------------------%
%%%%% physical string parameters 
phys_param.T = 87.279e1;               % tension (N)
phys_param.r = 0.0001395;              % string radius (m)
phys_param.rho = 7850;                 % density (kg/m^3)
phys_param.T60 = 2;                    % T60 (s)
phys_param.L = 1;                      % length (m)
phys_param.E = 2e11;                   % Young's modulus (Pa)

%%%%% simulation parameters 
sim_param.SR = 44100*1;                % sample rate (Hz)
sim_param.Tf = 1;                      % duration of simulation (s)
sim_param.xi = 0.1;                    % coordinate of excitation (normalised, 0-1)
sim_param.famp = 1;                    % peak amplitude of excitation (N)
sim_param.dur = 0.001;                 % duration of excitation (s)
sim_param.exc_st = 0.01;               % start time of excitation (s)
sim_param.xo = 0.9;                    % coordinate of output (normalised, 0-1)

y6 = string_fdtd_s1889125(opts,phys_param,sim_param);

% Stereo generation
y6L=y6 * sqrt(2)/2;
y6R=y6 * sqrt(2)/2;
y6Stereo=[y6L y6R];


%-------------------------------------------------------------------------%
% MONO Output (Uncomment to use)
%-------------------------------------------------------------------------%

%%Pause time is 0.5s  between two notes
% tp = 0.5;
% pause_time = zeros(floor(tp*44100),1);
%%Output sound
% total = cat(1,y1,pause_time,y2,pause_time,y3,pause_time,y4,pause_time,y5...
%     ,pause_time,y6);
%%Normalised vector
%total= total/max(abs(total));
% soundsc(total,44100);
%-------------------------------------------------------------------------%
% Stereo Output with panning effect
%-------------------------------------------------------------------------%
%%Pause time is 0.5s  between two notes
pause_time = zeros(floor(0.5*44100),2);
%%Output sound vector in stereo
total = cat(1,y1Stereo,pause_time,y2Stereo,pause_time,y3Stereo,pause_time,y4Stereo,pause_time,y5Stereo...
    ,pause_time,y6Stereo);
%%Normalised vector
total= total/max(abs(total));
soundsc(total,44100);



%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%







