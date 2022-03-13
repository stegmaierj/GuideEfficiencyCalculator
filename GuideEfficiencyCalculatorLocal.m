%%
% GuideEfficiencyCalculator.
% Copyright (C) 2017 Christelle Etard, Swarnima Joshi, Johannes Stegmaier, Ralf Mikut, and Uwe Strähle
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the Liceense at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
% Please refer to the documentation for more information about the software
% as well as for installation instructions.
%
% If you use this application for your work, please cite the repository and one
% of the following publications:
%
% Etard, C., Joshi, S., Stegmaier, J., Mikut, R., & Strähle, U. "Tracking
% of indels by decomposition is a simple and effective method to assess
% efficiency of guide RNAs in zebrafish". Zebrafish, 14(6), 586-588. 2017.
%
%%

%% add third party tools
addpath('toolbox/');

%% get the input files
disp('Please select one or more *.scf files for processing. Folders are skipped but you can select files from multiple folders.');
inputFiles = uipickfiles('FilterSpec', '*.scf');

%% 
if (~iscell(inputFiles))
    disp('No input files selected, aborting the processing. Please select *.scf files for processing!');
end

%% process all input files sequentially
for i=1:length(inputFiles)

    if (isfolder(inputFiles{i}))
        disp(['Skipping folder ' inputFiles{i} '. Please select individual files!']);
        continue;
    end

    %% perform preview
    disp(['Generating preview for ' inputFiles{i} ' ...']);
    [efficiency, stdDev, stdErr, figureHandle] = GuideEfficiencyCalculator(inputFiles{i}, -1, -1, true);

    %% maximize preview figure for PAM and radius estimation
    set(figureHandle, 'units', 'normalized', 'OuterPosition', [0 0 1 1]);
    figure(figureHandle);
    title(strrep(inputFiles{i}, '_', '-'));

    %% bring figure to focus
    figure(figureHandle);

    %% ask the user for the PAM location and the radius to use
    prompt = {'Enter PAM Location:','Enter Radius:'};
    dlgtitle = 'Parameter Specification';
    dims = [1 35];
    definput = {'-1','50'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    close(figureHandle);

    %% perform processing with the user-checked parameters
    disp(['Generating result images for ' inputFiles{i} ' with L=' answer{1} ' and R=' answer{2} '...']);
    GuideEfficiencyCalculator(inputFiles{i}, str2double(answer{1}), str2double(answer{2}), false);

    [folder, file, ext] = fileparts(inputFiles{i});
    disp(['Results for ' inputFiles{i} ' saved to ' folder filesep 'results' filesep]);
end