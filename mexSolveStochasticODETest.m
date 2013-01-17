%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Model Parameters, Program Parameters, Initial Conditions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Set up default values in structures/arrays %%%

%%% ModelParamaters contains rate constants for the ODEs %%%

handles.ModelParameters.gamma_V     = 510;
handles.ModelParameters.gamma_VA    = 619.2;
handles.ModelParameters.gamma_VH    = 1.02;
handles.ModelParameters.alpha_V     = 1.7;
handles.ModelParameters.a_V1        = 100;
handles.ModelParameters.a_V2        = 23000;
handles.ModelParameters.b_HD        = 4;
handles.ModelParameters.a_R         = 1;
handles.ModelParameters.gamma_HV    = 0.34;
handles.ModelParameters.b_HF        = 0.01;
handles.ModelParameters.b_IE        = 0.066;
handles.ModelParameters.a_I         = 1.5;
handles.ModelParameters.b_MD        = 1;
handles.ModelParameters.b_MV        = 0.0037;
handles.ModelParameters.a_M         = 1;
handles.ModelParameters.b_F         = 250000;
handles.ModelParameters.c_F         = 2000;
handles.ModelParameters.b_FH        = 17;
handles.ModelParameters.a_F         = 8;
handles.ModelParameters.b_EM        = 8.3;
handles.ModelParameters.b_EI        = 2.72;
handles.ModelParameters.a_E         = 0.4;
handles.ModelParameters.b_PM        = 11.5;
handles.ModelParameters.a_P         = 0.4;
handles.ModelParameters.b_A         = 0.043;
handles.ModelParameters.gamma_AV    = 146.2;
handles.ModelParameters.a_A         = 0.043;
handles.ModelParameters.r           = 3e-5;

%%% ProgramParameters contains information such as start and end times %%%

handles.ProgramParameters.t_start   = 0;        % Start time
handles.ProgramParameters.t_end     = 30;       % End time


%%% InitialConditions contains the initial values of the variables %%%

handles.InitialConditions(1)    = 0.01;         % V(0)
handles.InitialConditions(2)    = 1;            % H(0)
handles.InitialConditions(3)    = 0;            % I(0)
handles.InitialConditions(4)    = 0;            % M(0)
handles.InitialConditions(5)    = 0;            % F(0)
handles.InitialConditions(6)    = 0;            % R(0)
handles.InitialConditions(7)    = 1;            % E(0)
handles.InitialConditions(8)    = 1;            % P(0)
handles.InitialConditions(9)    = 1;            % A(0)
handles.InitialConditions(10)   = 0.1;          % S(0)

%%% Timestep size
options.stepSize                = 0.01;         % Timestep size

%%% NoiseParameters contains the initial values of the variables %%%

handles.NoiseParameters(1,1)    = 1;            % V(0)
handles.NoiseParameters(2,1)    = 0;            % H(0)
handles.NoiseParameters(3,1)    = 0;            % I(0)
handles.NoiseParameters(4,1)    = 0;            % M(0)
handles.NoiseParameters(5,1)    = 0;            % F(0)
handles.NoiseParameters(6,1)    = 0;          % R(0)
handles.NoiseParameters(7,1)    = 0;            % E(0)
handles.NoiseParameters(8,1)    = 0;            % P(0)
handles.NoiseParameters(9,1)    = 0;            % A(0)
handles.NoiseParameters(10,1)   = 0;            % S(0)

handles.NoiseParameters(1,2)    = 0.5;            % V(0)
handles.NoiseParameters(2,2)    = 0;            % H(0)
handles.NoiseParameters(3,2)    = 0;            % I(0)
handles.NoiseParameters(4,2)    = 0;            % M(0)
handles.NoiseParameters(5,2)    = 0;            % F(0)
handles.NoiseParameters(6,2)    = 0;            % R(0)
handles.NoiseParameters(7,2)    = 0;            % E(0)
handles.NoiseParameters(8,2)    = 0;            % P(0)
handles.NoiseParameters(9,2)    = 0;            % A(0)
handles.NoiseParameters(10,2)   = 0;            % S(0)

handles.NoiseParameters(1,3)    = 1;            % V(0)
handles.NoiseParameters(2,3)    = 0;            % H(0)
handles.NoiseParameters(3,3)    = 0;            % I(0)
handles.NoiseParameters(4,3)    = 0;            % M(0)
handles.NoiseParameters(5,3)    = 0;            % F(0)
handles.NoiseParameters(6,3)    = 0;            % R(0)
handles.NoiseParameters(7,3)    = 0;            % E(0)
handles.NoiseParameters(8,3)    = 0;            % P(0)
handles.NoiseParameters(9,3)    = 0;            % A(0)
handles.NoiseParameters(10,3)   = 0;            % S(0)


profile on

profile clear

[times, Solution] = mexSolveStochasticODE(handles.ModelParameters,...
    [handles.ProgramParameters.t_start,handles.ProgramParameters.t_end],...
    handles.InitialConditions, options, handles.NoiseParameters);

profile off