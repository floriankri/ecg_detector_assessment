function R_DECO()
% R-DECO: An open-source Matlab based graphical user interface for the detection and correction of R-peaks
%
% R_DECO() Starts the R-peak detection GUI tool
%
% Author(s):    Jonathan Moeyersons       (Jonathan.Moeyersons@esat.kuleuven.be)
%               Sabine Van Huffel         (Sabine.Vanhuffel@esat.kuleuven.be)
%               Carolina Varon            (Carolina.Varon@esat.kuleuven.be)
%
% Version History:
% - 06/05/2019   JM      Initial version
%
% Copyright (c) 2019,  Jonathan Moeyersons, KULeuven-ESAT-STADIUS 
%
% This software is made available for non commercial research purposes only 
% under the GNU General Public License. However, notwithstanding any 
% provision of the GNU General Public License, this software may not be 
% used for commercial purposes without explicit written permission after 
% contacting jonathan.moeyersons@esat.kuleuven.be
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.


% Determine where your m-file's folder is.
filename = mfilename('fullpath');

% Get only the folder
[folder,~,~] = fileparts(filename);

% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

% Check version
if verLessThan('matlab', '9.0')
    warndlg('You are using an old matlab version, the GUI might not be fully functional')
end

% Run R-peak detection GUI
R_peak_detection