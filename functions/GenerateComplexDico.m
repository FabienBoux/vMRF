function [] = GenerateComplexDico(FileConfig, FileSeq, FileParams, PathOut)
addpath(genpath(fullfile(pwd(), 'functions')))
% First simulation
fprintf(['Simulation for sequence : ' FileSeq{1} '\n']);
tic
[X, Y, Param_names, AltModel] = GenerateSignals(FileConfig, FileSeq{1}, FileParams);
Exec_time = toc;
fprintf(['\t Completed in ' num2str(Exec_time) ' s\n']);
filenames{1} = [datestr(now,'yyyy-mm-dd-HH:MM') '_' num2str(size(X,2)) '-samples' '_' num2str(size(Y,2)) '-parameters' '_' num2str(size(X,1)) '-signals' '.mat'];
save(fullfile(PathOut, filenames{1}), 'X','Y', 'Param_names', 'AltModel', 'Exec_time')
% Save parameter configuration
NewFileParams = [strtok(FileParams,'.') '_copy.txt'];
WriteModel(fullfile(PathOut, filenames{1}), NewFileParams, 'AltModel')
% Other simulations
for s = 2:length(FileSeq)
    clear X Y Param_names AltModel Exec_time
    
    fprintf(['Simulation for sequence : ' FileSeq{s} '\n']);
    
    tic
    [X, Y, Param_names, AltModel] = GenerateSignals(FileConfig, FileSeq{s}, NewFileParams);
    Exec_time = toc;
    
    fprintf(['\t Completed in ' num2str(Exec_time) ' s\n']);
    
    filenames{s} = [datestr(now,'yyyy-mm-dd-HH:MM-ss') '_' num2str(size(X,2)) '-samples' '_' num2str(size(Y,2)) '-parameters' '_' num2str(size(X,1)) '-signals' '.mat'];
    save(fullfile(PathOut, filenames{s}), 'X','Y', 'Param_names', 'AltModel', 'Exec_time')
end
Xcomplex = [];
for s = 1:length(filenames)
    load(fullfile(PathOut, filenames{s}))
    
    Xcomplex = [Xcomplex X];
    if s > 1
        if any(Ycomplex ~= Y), error('Parameters are not equals'); end
    end
    Ycomplex = Y;   
end
X = Xcomplex;
Y = Ycomplex;
filename = [datestr(now,'yyyy-mm-dd-HH:MM-ss') '_complex_dico' '.mat'];
if exist(PathOut,'dir')
    save(fullfile(PathOut, filename), 'X','Y', 'Param_names', 'AltModel', 'FileSeq')
else
    save(fullfile('.', filename), 'X','Y', 'Param_names', 'AltModel', 'FileSeq')
end