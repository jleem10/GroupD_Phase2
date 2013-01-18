%INFLUENZAA_GUI creates a new copy of the GUI influenzaA_GUI or raises the
% existing copy. From the GUI, the user can set the value of various
% parameters and run one of two ODE solvers to solve the system of ODEs
% representing the progression of an influenza A infection.
%
%Initialisation: Default values for each of the initial conditions and
% parameters can be set.
%
%Parameter Selection: With nsweeps = 1, the user can produce graphs showing
% disease progression using the default values set in the initialisation
% tab. With nsweeps > 1, the user can choose linear or logarithmic spacing
% of the selected parameter between a minimum and maximum value. This
% parameter is then swept over, producing nsweeps series of data. Minimum
% and maximum values of time are also set here.
%
%IMPORTANT: When nsweeps is set to 1, the value of that parameter is
% selected by the menu in the Initialisation section - the min and max
% values in the Parameter Selection tab have no effect.
%
%Plotting: The user can choose which variables to plot from the sets of
% data produced. Limits of the x and y axes can be changed to improve the
% appearance of the plots.
%
%Save: The save buttons allow the user to save the data as a text file or
% to save the current plot. The text files are best opened in an internet
% browser to ensure that all line breaks display properly. A variety of
% file types are available when saving the current plot.
%
%MEX solver: An alternative ODE solver has been implemented in C and
% inserted as a MEX function into the main code. The function has been
% compiled to support a 64-bit Windows opertaing system. Recompilation will
% be necessary to support its use on other operating systems.
%
% -------------------------------------------------------------------------
%
%   A Program to implement and reproduce the results of Hancioglu et al.'s Paper
%   'A Dynamical Model of Human Immune Response to Influenza A Virus Infection'
%   Copyright (C) 2013  Paul Taylor et al.
% 
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

function varargout = influenzaA_GUI(varargin)
%
% INFLUENZAA_GUI MATLAB code for influenzaA_GUI.fig
%      INFLUENZAA_GUI, by itself, creates a new INFLUENZAA_GUI or raises the existing
%      singleton*.
%
%      H = INFLUENZAA_GUI returns the handle to a new INFLUENZAA_GUI or the handle to
%      the existing singleton*.
%
%      INFLUENZAA_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INFLUENZAA_GUI.M with the given input arguments.
%
%      INFLUENZAA_GUI('Property','Value',...) creates a new INFLUENZAA_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before influenzaA_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to influenzaA_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%

% Edit the above text to modify the response to help influenzaA_GUI

% Last Modified by GUIDE v2.5 16-Jan-2013 16:36:08

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Begin initialization code - DO NOT EDIT %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @influenzaA_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @influenzaA_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);

if nargin && ischar(varargin{1})

    gui_State.gui_Callback = str2func(varargin{1});

end

if nargout

    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});

else

    gui_mainfcn(gui_State, varargin{:});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% End initialization code - DO NOT EDIT %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GUI-Opening Function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes just before influenzaA_GUI is made visible.
function influenzaA_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
%
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to influenzaA_GUI (see VARARGIN)
%

%%% Choose default command line output for influenzaA_GUI %%%

handles.output = hObject;

%%% General %%%

handles.success = 0;        % Check is output can be/has been produced
                            % Only plot if true
handles.line = 1;           % 'LineWidth' for graphs in points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialise variables for GUI %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

handles.initial_parameter = 1;      % Number of popupmenu selection for initialisation
handles.initialise_variable = 1;
handles.initialise_value = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameter selection tab %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% No extra handles variables required here
% All inputs are stored in structures listed below for easy input to the ODE
% functions

%%%%%%%%%%%%%%%%%%%%
%%% Plotting tab %%%
%%%%%%%%%%%%%%%%%%%%

% These values apply to axis1 (may introduce more axes later)

handles.axis1_xaxis = 12;   % Time from popupmenu
handles.axis1_yaxis = 1;    % V from popupmenu
handles.axis1_log = [0 0];  % true/false values for plotting [x y] on log scale
handles.xlim = [0 1];
handles.ylim = [0 1];

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
handles.ModelParameters.d_max       = 4;  %----CHANGE----%
handles.ModelParameters.t_peak      = 0.15;  %----CHANGE----%
handles.ModelParameters.t_half      = 0.4;  %----CHANGE----%
handles.ModelParameters.t_on        = 5.1;  %----CHANGE----%

%%% ProgramParameters contains information such as start and end times %%%

handles.ProgramParameters.t_start   = 0;        % Start time
handles.ProgramParameters.t_end     = 15;       % End time
handles.ProgramParameters.min       = 0.01;     % V(0) min value
handles.ProgramParameters.max       = 0.01;     % V(0) max value (same as default)
handles.ProgramParameters.nsweeps   = 1;        % Number of different values in parameter sweep
handles.ProgramParameters.sweepvar = 1;         % V(0) from popupmenu
handles.ProgramParameters.logsweep = 0;         % log spacing of sweeps
handles.ProgramParameters.mex = 0;              % option for using mex ODE solver instead of ODE15s

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

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes influenzaA_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GUI Output Function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Outputs from this function are returned to the command line.
function varargout = influenzaA_GUI_OutputFcn(hObject, eventdata, handles) 
%
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Callback functions %%%%%%%%%
% Actions following button press etc. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The first use of each type of uiobject (button, menu etc. still has the
% comments from GUIDE. Repeated copies of this information have been deleted
% for clarity.

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialisation tab %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on selection change in popupmenu_initialisation_parameter.
function popupmenu_initialisation_parameter_Callback(hObject, eventdata, handles)
%
% hObject    handle to popupmenu_initialisation_parameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_initialisation_parameter contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_initialisation_parameter

handles.initial_parameter = get(hObject, 'Value');

guidata(hObject, handles);  % Update handles structure
                            % Do this at end of every callback which changes handles

current_value = getInitialValue(handles,handles.initial_parameter);

set(handles.edit_initialisation_value,'String',num2str(current_value));

guidata(hObject, handles);
    
% --- Executes on value change in edit_initialisation_value.
function edit_initialisation_value_Callback(hObject, eventdata, handles)
%
% hObject    handle to edit_initialisation_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%

% Hints: get(hObject,'String') returns contents of edit_initialisation_value as text
%        str2double(get(hObject,'String')) returns contents of edit_initialisation_value as a double

handles.new_value = str2double(get(hObject,'String'));
guidata(hObject, handles); 

handles = setInitialValue(handles);

guidata(hObject, handles);    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameter Selection tab %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on selection change in popupmenu_parameter_sweepvar.
function popupmenu_parameter_sweepvar_Callback(hObject, eventdata, handles)
handles.ProgramParameters.sweepvar = get(hObject, 'Value');
current_value = getInitialValue(handles,handles.ProgramParameters.sweepvar);
set(handles.edit_parameter_min,'String',num2str(current_value));
set(handles.edit_parameter_max,'String',num2str(current_value));
guidata(hObject, handles);

function edit_parameter_min_Callback(hObject, eventdata, handles)
tmp = str2double(get(hObject,'String'));
if (tmp>=0 && isnan(tmp) ~=1)
    handles.ProgramParameters.min = tmp;
else 
	set(hObject,'String',num2str(handles.ProgramParameters.min));
end
guidata(hObject, handles);

function edit_parameter_max_Callback(hObject, eventdata, handles)
tmp = str2double(get(hObject,'String'));
if (tmp>=0 && isnan(tmp) ~=1)
    handles.ProgramParameters.max = tmp;
else 
	set(hObject,'String',num2str(handles.ProgramParameters.min));
end
guidata(hObject, handles);

function edit_parameter_nsweeps_Callback(hObject, eventdata, handles)
tmp = round(str2double(get(hObject,'String')));
if (tmp>=1 && isnan(tmp) ~=1)
    handles.ProgramParameters.nsweeps = tmp;
end
set(hObject,'String',num2str(handles.ProgramParameters.nsweeps)); % Show rounded number
guidata(hObject, handles);

function edit_parameter_tmin_Callback(hObject, eventdata, handles)
tmp = str2double(get(hObject,'String'));
if (tmp>=0 && isnan(tmp) ~=1)
    handles.ProgramParameters.t_start = tmp;
end
guidata(hObject, handles);

function edit_parameter_tmax_Callback(hObject, eventdata, handles)
tmp = str2double(get(hObject,'String'));
if (tmp>0 && isnan(tmp) ~=1)
    handles.ProgramParameters.t_end = tmp;
end
guidata(hObject, handles);

function checkbox_parameter_log_sweep_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_parameter_log_sweep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_parameter_log_sweep
handles.ProgramParameters.logsweep = get(hObject, 'Value');
guidata(hObject, handles);

% --- Executes on button press in pushbutton_parameter_run.
function pushbutton_parameter_run_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_parameter_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.success = prepareRun(handles);    % Perform checks to ensure
                                            % valid inputs (success = 1 is all okay)
if ((handles.success == 1) && (handles.ProgramParameters.nsweeps == 1))
    handles.result = cell(1,1);
    handles.result{1} = SolveODE(handles.InitialConditions,...
        handles.ProgramParameters,handles.ModelParameters);
    handles.SweepVariableName = 'Only one set of results - no sweep';
    handles.SweepValues = 1;
elseif ((handles.success == 1) && (handles.ProgramParameters.nsweeps ~= 1))
    handles.result = cell(handles.ProgramParameters.nsweeps);
    [handles.result handles.legend_entries handles.SweepVariableName, handles.SweepValues] = parameterSweep(handles.InitialConditions,...
        handles.ProgramParameters,handles.ModelParameters);
else fprintf(1,'There is a problem with the values sent to Solve ODE\n');
    fprintf(1,'Check for input errors such as min > max\n\n');
end
handles.result = checkZeros(handles.result,handles.ProgramParameters.nsweeps);
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Plotting tab %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on selection change in popupmenu_plot_xaxis.
function popupmenu_plot_xaxis_Callback(hObject, eventdata, handles)
handles.axis1_xaxis = get(hObject, 'Value');
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_plot_yaxis.
function popupmenu_plot_yaxis_Callback(hObject, eventdata, handles)
handles.axis1_yaxis = get(hObject, 'Value');
guidata(hObject, handles);

% --- Executes on button press in checkbox_plot_logx.
function checkbox_plot_logx_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot_logx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox_plot_logx
handles.axis1_log(1) = get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in checkbox_plot_logy.
function checkbox_plot_logy_Callback(hObject, eventdata, handles)
handles.axis1_log(2) = get(hObject,'Value');
guidata(hObject, handles);

function edit_plot_minx_Callback(hObject, eventdata, handles)
handles.xlim(1) = str2double(get(hObject,'String'));
xlim(handles.xlim);
guidata(hObject, handles);

function edit_plot_maxx_Callback(hObject, eventdata, handles)
handles.xlim(2) = str2double(get(hObject,'String'));
xlim(handles.xlim);
guidata(hObject, handles);

function edit_plot_miny_Callback(hObject, eventdata, handles)
handles.ylim(1) = str2double(get(hObject,'String'));
ylim(handles.ylim);
guidata(hObject, handles);

function edit_plot_maxy_Callback(hObject, eventdata, handles)
handles.ylim(2) = str2double(get(hObject,'String'));
ylim(handles.ylim);
guidata(hObject, handles);

% --- Executes on button press in pushbutton_plot_plot.
function pushbutton_plot_plot_Callback(hObject, eventdata, handles)
if (handles.success==1)
    plotSweep(handles);
    x_limits = xlim;
    y_limits = ylim;
    set(handles.edit_plot_minx,'String',sprintf('%.3g',x_limits(1)));
    set(handles.edit_plot_maxx,'String',sprintf('%.3g',x_limits(2)));
    set(handles.edit_plot_miny,'String',sprintf('%.3g',y_limits(1)));
    set(handles.edit_plot_maxy,'String',sprintf('%3.g',y_limits(2)));
else
    fprintf(1,'Can''t Plot Yet - No Results\n');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Save buttons %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pushbutton_save_Callback(hObject, eventdata, handles)
FileName = uiputfile('*.txt');
SaveToFile (handles.result, FileName, handles.SweepVariableName, handles.SweepValues);

function pushbutton_save_figure_Callback(hObject, eventdata, handles)
saveFigure(handles);

function checkbox_mex_Callback(hObject, eventdata, handles)
handles.ProgramParameters.mex = get(hObject,'Value');
guidata(hObject,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Create functions %%%%%%%%%%
%%%%%%%%%%%% Made by GUIDE %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Not in use here (apart from colour change), but do not delete!

% --- Executes during object creation, after setting all properties.
function popupmenu_initialisation_parameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_initialisation_parameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_initialisation_value_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu_parameter_sweepvar_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_parameter_min_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_parameter_max_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_parameter_nsweeps_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_parameter_tmin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_parameter_tmax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu_plot_xaxis_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu_plot_yaxis_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_plot_minx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_plot_maxx_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_plot_miny_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_plot_maxy_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on popupmenu_initialisation_parameter and none of its controls.
function popupmenu_initialisation_parameter_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_initialisation_parameter (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
