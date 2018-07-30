function [Model] = ReadModel_adapted(File)

% Read and parse the parameter file
% Input: 
%   - File  : Path to the parameter file
% Output:
%   - Model     : Structure containing the voxel related parameters
%   - Seq       : Structure containing the sequence related parameters
%
% Code by N.Pannetier, C.Debacker, INSERM, Grenoble, 2012

Model = struct;

%% read lines
fid = fopen(File,'rt');
if fid > 0
    C = textscan(fid, '%s', 'Delimiter',''); C = C{1};
else
    fprintf('Cannot open Parameter file\n');
end
fclose(fid);

%% Read Model strucutre
for a=1:numel(C)
   eval(C{a}); 
end
