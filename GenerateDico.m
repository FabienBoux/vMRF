function [] = GenerateDico(FileConfig, FileSeq, FileParams, PathOut)

if nargin > 4
    error('GenerateDico:TooManyInputs');
end
switch nargin
    case 1
        error('GenerateDico:NotEnoughInputs')
    case 2
        PathOut = pwd();
end

addpath(genpath(fullfile(pwd(), 'functions')))


%% Simulations
tic
[X, Y, Param_names, AltModel] = GenerateSignals(FileConfig, FileSeq, FileParams);
Exec_time = toc;

filename = [datestr(now,'yyyy-mm-dd-HH:MM') '_' num2str(size(X,2)) '-samples' '_' num2str(size(Y,2)) '-parameters' '_' num2str(size(X,1)) '-signals' '.mat'];
save(fullfile(PathOut, filename), 'X','Y', 'Param_names', 'AltModel', 'Exec_time')






