%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    CHRONIC AORTIC REGURGITATION (2D)                    %
%                 RECONSTRUCT GITHUB REDUCED-ORDER MODELS                 %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Giuseppe Di Labbio                                                      %
% Department of Mechanical Engineering                                    %
% Polytechnique Montréal, Montréal, Québec, Canada                        %
%                                                                         %
% Last Update: December 16th, 2019 by Giuseppe Di Labbio                  %
%                                                                         %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Copyright (C) 2019 Giuseppe Di Labbio                                   %
%                                                                         %
% This program is free software: you can redistribute it and/or modify it %
% under the terms of the GNU General Public License as published by the   %
% Free Software Foundation, either version 3 of the License, or (at your  %
% option) any later version.                                              %
%                                                                         %
% This program is distributed in the hope that it will be useful, but     %
% WITHOUT ANY WARRANTY; without even the implied warranty of              %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        %
% General Public License for more details.                                %
%                                                                         %
% You should have received a copy of the GNU General Public License along %
% with this program. If not, see <https://www.gnu.org/licenses/>.         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% CLEAR & CLOSE %%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% clear all;
% close all;
% cla;
% clf;
clc;

%% %%%%%%%%%%%%%%%%%%%%%%%%%% INPUT VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Specify the path to the reduced-order models.
path  = 'Y:\Documents\GitHub\ROMs--LV_flow_with_AR\';

% Specify the case of interest. Choose from 'Healthy', 'Mild', 'Moderate1',
% 'Moderate2', 'Severe1', 'Severe1DMD', 'Severe2' and 'Severe2DMD'.
state = 'Severe2DMD';

% Refine data to original resolution ('cubic', 'linear', 'none').
refine = 'none';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

DMDcases = {'Healthy'; 'Mild'; 'Moderate1'; 'Moderate2'; 'Severe1DMD';  ...
            'Severe2DMD'};
DMDnums  = [49; 41; 61; 47; 133; 119];
PODcases = {'Severe1'; 'Severe2'};
PODnums  = [84; 71];
t        = (((0:342).' - 15)/400)/(6/7);

% Read the X data.
fileID = fopen([path 'X.dat'], 'rt');
format = repmat('%f', 1, 1);
X      = cell2mat(textscan(fileID, format));
fclose(fileID);

% Read the Y data.
fileID = fopen([path 'Y.dat'], 'rt');
format = repmat('%f', 1, 1);
Y      = cell2mat(textscan(fileID, format));
fclose(fileID);

% Read the mask data.
if max(strcmp(state, {'Severe1DMD'; 'Severe2DMD'}))
    fileID = fopen([path state(1:7) '_Mask.dat'], 'rt');
else
    fileID = fopen([path state '_Mask.dat'], 'rt');
end
format = repmat('%f', 1, 343);
mask   = cell2mat(textscan(fileID, format));
mask   = reshape(mask, 75, 82, 343);
fclose(fileID);

[isDMD, ind] = max(strcmp(state, DMDcases));
if isDMD
    % Read the real part of the spatial mode.
    fileID  = fopen([path state '_Modes_Re.dat'], 'rt');
    format  = repmat('%f', 1, DMDnums(ind));
    ModesRe = cell2mat(textscan(fileID, format));
    fclose(fileID);
    
    % Read the imaginary part of the spatial mode.
    fileID  = fopen([path state '_Modes_Im.dat'], 'rt');
    format  = repmat('%f', 1, DMDnums(ind));
    ModesIm = cell2mat(textscan(fileID, format));
    fclose(fileID);
    
    % Read the real part of the temporal dynamics.
    fileID = fopen([path state '_Dynamics_Re.dat'], 'rt');
    format = repmat('%f', 1, DMDnums(ind));
    TdynRe = cell2mat(textscan(fileID, format));
    fclose(fileID);
    
    % Read the imaginary part of the temporal dynamics.
    fileID = fopen([path state '_Dynamics_Im.dat'], 'rt');
    format = repmat('%f', 1, DMDnums(ind));
    TdynIm = cell2mat(textscan(fileID, format));
    fclose(fileID);
    
    % Combine the spatial modes into a complex variable.
    Modes = ModesRe + 1i*ModesIm;
    clear ModesRe ModesIm;
    
    % Combine the temporal dynamics into a complex variable.
    Tdyn = (TdynRe + 1i*TdynIm).';
    clear TdynRe TdynIm;
    
    % Reconstruct the flow.
    Flow = real(Modes*Tdyn);
    U    = mask.*reshape(Flow(1:6150,:), 75, 82, 343);
    V    = mask.*reshape(Flow(6151:end,:), 75, 82, 343);
    clear Flow Modes Tdyn;
end
clear DMDcases DMDnums ind isDMD;

[isPOD, ind] = max(strcmp(state, PODcases));
if isPOD
    % Read the spatial mode.
    fileID = fopen([path state '_Modes.dat'], 'rt');
    format = repmat('%f', 1, PODnums(ind));
    Modes  = cell2mat(textscan(fileID, format));
    fclose(fileID);
    
    % Read the temporal dynamics.
    fileID = fopen([path state '_Dynamics.dat'], 'rt');
    format = repmat('%f', 1, PODnums(ind));
    Tdyn   = cell2mat(textscan(fileID, format));
    fclose(fileID);
    
    % Reconstruct the flow.
    Flow = Modes*(Tdyn.');
    U    = mask.*reshape(Flow(1:6150,:), 75, 82, 343);
    V    = mask.*reshape(Flow(6151:end,:), 75, 82, 343);
    clear Flow Modes Tdyn;
end
clear PODcases PODnums ind isPOD

% Refine the domain to its original size.
if max(strcmp(refine, {'cubic'; 'linear'}))
    % Refine the X and Y coordinates.
    Xr = interp1(1:2:163, X, 1:163).';
    Yr = interp1(1:2:149, Y, 1:149).';
    [Xg, Yg] = meshgrid(Xr, Yr);
    
    % Refine the velocity components.
    Ur = zeros(length(Yr), length(Xr), length(t)); Vr = Ur;
    for k = 1:length(t)
        Ur(:,:,k) = interp2(X, Y, U(:,:,k), Xg, Yg, refine);
        Vr(:,:,k) = interp2(X, Y, V(:,:,k), Xg, Yg, refine);
    end
    X = Xr; Y = Yr; U = Ur; V = Vr; mask = logical(U);
end
clear k refine Ur Vr Xg Xr Yg Yr;

% Clear remaining variables.
clear ans fileID format path state;

%% %%%%%%%%%%%%%%%%%%%%%%%%% SUPPRESS MESSAGES %%%%%%%%%%%%%%%%%%%%%%%%% %%

%#ok<*CLALL>
% Line(s) 41
% Message(s)
% * Using CLEAR ALL usually decreases code performance and is often
%   unnecessary.
% Reason(s)
% * The 'clear all' command is here used to simply leave the workspace
%   uncrowded.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%                                                                         %
% Line(s) N/A                                                             %
% * N/A.                                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
