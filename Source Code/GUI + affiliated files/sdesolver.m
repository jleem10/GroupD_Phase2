function [ times, Y ] = sdesolver(ModelParameters,ProgramParameters,InitialConditions,NoiseParameters)
%SDEsolver a stochastic differential equation solver that uses the
%   Euler-Maruyama method to solve the influenza model.
%   Takes as inputs struct ModelParameters, struct ProgramParameters (these
%   define the time limits and steps), struct InitialConditions and struct
%   NoiseParameters. NoiseParameters in configured as a 10 by 3 matrix,
%   each row being a variable. The first column returns 1 if the variable
%   needs noise, zero otherwise. The second column gives A and the third B,
%   such that the noise for each variable x scales as A*x^B. Outputs
%   are [Solution, times].

times=ProgramParameters.t_start:ProgramParameters.t_step:ProgramParameters.t_end;

Y=zeros(12,length(times));
Y(:,1)=InitialConditions;

for t=2:length(times)
    
    Y(:,t) = Y(:,t-1) + ProgramParameters.t_step * DifferentialEquations(Y(:,t-1),ModelParameters)...
        + sqrt(ProgramParameters.t_step).*NoiseParameters(:,1).*((Y(:,t-1)).^(NoiseParameters(:,2))).*randn;
    
    Y(:,t)=Y(:,t).*(0.^(Y(:,t)<0));
    
end


end

function [ dYdt] = DifferentialEquations(Y, ModelParameters)
%
% Differential Equations Governing the Model Populations
% 
% Input Arguments: Time, t - type DOUBLE
%				   Model Variables, Y, - type VECTOR (as defined in SolveODE)
%				   Model Parameters, ModelParameters - type STRUCT (as defined in SolveODE)
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
TV = Y(11); % Change in helper T cell %
Bi = Y(12); % Change in B cell %
D = 1 - H - R - I;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Local Differentials %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (ModelParameters.t_peak == 0)
    c_1 = 0; c_2 = 0; c_3 = 0;
else
    c_3 = 1/ModelParameters.t_peak;
    c_2 = (c_3^2*ModelParameters.t_half^2-4*c_3*ModelParameters.t_half+1)/ModelParameters.t_half;
    c_1 = ModelParameters.d_max*(2*c_3+c_2);
end;
% add this c_1*t./(1+c_2*t+c_3^2*t.^2)*

dYdt(1,1)  = (ModelParameters.gamma_V * I)...
    - (c_1*t./(1+c_2*t+c_3^2*t.^2)*ModelParameters.gamma_V * I)...  %-CHANGE!-%
    - (ModelParameters.gamma_VA * S * A * V)...
    - (ModelParameters.gamma_VH * H * V)...
    - (ModelParameters.alpha_V * V)...
    - (ModelParameters.a_V1 * V)/(1 + ModelParameters.a_V2 * V);

%dYdt(1,1)  = (ModelParameters.gamma_V * I) - (ModelParameters.gamma_VA * S * A * V) - (ModelParameters.gamma_VH * H * V) - (ModelParameters.alpha_V * V) - (ModelParameters.a_V1 * V)/(1 + ModelParameters.a_V2 * V);
dYdt(2,1)  = (ModelParameters.b_HD * D)*(H + R) + (ModelParameters.a_R * R) - (ModelParameters.gamma_HV * V * H) - (ModelParameters.b_HF * F * H);
dYdt(3,1)  = (ModelParameters.gamma_HV * V * H) - (ModelParameters.b_IE * I * E) - (ModelParameters.a_I * I);
dYdt(4,1)  = (ModelParameters.b_MD * D + ModelParameters.b_MV * V)*(1 - M) - (ModelParameters.a_M * M);
dYdt(5,1)  = (ModelParameters.b_F * M) + (ModelParameters.c_F * I) - (ModelParameters.b_FH * H * F) - (ModelParameters.a_F * F);
dYdt(6,1)  = (ModelParameters.b_HF * F * H) - (ModelParameters.a_R * R);
dYdt(7,1)  = (ModelParameters.b_EM * M * E) - (ModelParameters.b_EI * I * E) + (ModelParameters.a_E)*(1 - E);
%dYdt(8,1)  = (ModelParameters.b_PM * M * P) + (ModelParameters.a_P)*(1 - P);
dYdt(8,1) = (ModelParameters.b_PB * Bi * P) + (ModelParameters.a_P)*(1-P); % Revised equation for accounting B cells %
dYdt(9,1)  = (ModelParameters.b_A * P) - (ModelParameters.gamma_AV * S * A * V) - (ModelParameters.a_A * A);
dYdt(10,1) = (ModelParameters.r*P)*(1 - S);
dYdt(11,1) = (ModelParameters.kPT - ModelParameters.kDT) * TV * M;
dYdt(12,1) = (-ModelParameters.kDB * Bi) + (ModelParameters.kPB) * ((TV * Bi) / (TV + Bi));

end

