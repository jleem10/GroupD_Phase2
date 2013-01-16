function [ Y, times ] = sderatesolver(ModelParameters,ProgramParameters,InitialConditions,NoiseParameters)
%SDEratesolver a stochastic differential equation solver that uses the
%Euler-Maruyama method to solve the influenza model.
%   Takes as inputs struct ModelParameters, struct ProgramParameters (these
%   define the time limits and steps), struct InitialConditions and struct
%   NoiseParameters. NoiseParameters is configured as struct with two
%   elements, VNoise and SNoise, describing the SD of the Gaussian Noise
%   added to the production of V and S, respectively.

times=ProgramParameters.t_start:ProgramParameters.t_step:ProgramParameters.t_end;

Y=zeros(10,length(times));
Y(:,1)=InitialConditions;

for t=2:length(times)
    
    Y(:,t) = Y(:,t-1) + ProgramParameters.t_step * DifferentialEquations(Y(:,t-1),ModelParameters,NoiseParameters);
    
    Y(:,t)=Y(:,t).*(0.^(Y(:,t)<0));
    
end


end

function [ dYdt] = DifferentialEquations(Y, ModelParameters,NoiseParameters)
%
% Differential Equations Governing the Model Populations
% 
% Input Arguments: Time, t - type DOUBLE
%				   Model Variables, Y, - type VECTOR (as defined in SolveODE)
%				   Model Parameters, ModelParameters - type STRUCT (as defined in SolveODE)
%                  Noise Parameters - Struct
%
% Output Arguments: VECTOR dYdt containing point differentials
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Assign Data to Local Variables %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

V = Y(1);               
H = Y(2);
I = Y(3);
M = Y(4);
F = Y(5);
R = Y(6);
E = Y(7);
P = Y(8);
A = Y(9);
S = Y(10);
D = 1 - H - R - I;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Local Differentials %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dYdt(1,1)  = (ModelParameters.gamma_V * I)*(1+NoiseParameters.VNoise*randn) - (ModelParameters.gamma_VA * S * A * V) - (ModelParameters.gamma_VH * H * V) - (ModelParameters.alpha_V * V) - (ModelParameters.a_V1 * V)/(1 + ModelParameters.a_V2 * V);
dYdt(2,1)  = (ModelParameters.b_HD * D)*(H + R) + (ModelParameters.a_R * R) - (ModelParameters.gamma_HV * V * H) - (ModelParameters.b_HF * F * H);
dYdt(3,1)  = (ModelParameters.gamma_HV * V * H) - (ModelParameters.b_IE * I * E) - (ModelParameters.a_I * I);
dYdt(4,1)  = (ModelParameters.b_MD * D + ModelParameters.b_MV * V)*(1 - M) - (ModelParameters.a_M * M);
dYdt(5,1)  = (ModelParameters.b_F * M) + (ModelParameters.c_F * I) - (ModelParameters.b_FH * H * F) - (ModelParameters.a_F * F);
dYdt(6,1)  = (ModelParameters.b_HF * F * H) - (ModelParameters.a_R * R);
dYdt(7,1)  = (ModelParameters.b_EM * M * E) - (ModelParameters.b_EI * I * E) + (ModelParameters.a_E)*(1 - E);
dYdt(8,1)  = (ModelParameters.b_PM * M * P) + (ModelParameters.a_P)*(1 - P);
dYdt(9,1)  = (ModelParameters.b_A * P) - (ModelParameters.gamma_AV * S * A * V) - (ModelParameters.a_A * A);
dYdt(10,1) = (ModelParameters.r*P)*(1 - S)*(1+NoiseParameters.SNoise*randn);

end

