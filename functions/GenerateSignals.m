function [X, Y, Param_names, AltModel] = GenerateSignals(DicoFile, SeqFile, ParamFile)

addpath(genpath(fullfile(pwd, 'tools/mrvox')));

% Read both SeqFile and ModelFile
Model   = ReadModel(DicoFile);
Seq     = SelectSeq(SeqFile,1);

% Read alternative model
AltModel = ReadModel_adapted(ParamFile);

field_names = fieldnamesr(AltModel);

Ylen = 0;
for f = 1:length(field_names)
    tmp = eval(['AltModel.' field_names{f} ';']);
    if length(tmp) == 1
        eval(['Model.' field_names{f} ' = ' num2str(tmp) ';']);
    else
        Ylen              = Ylen + 1;
        Params{Ylen}      = tmp;
        Param_names{Ylen} = field_names{f};
    end
end

for f = 2:Ylen
    if length(Params{f-1}) ~= length(Params{f})
        error('Error parameter array sizes')
    else
        Y(:,f) = Params{f};
    end
end
Y(:,1) = Params{1};

parfor_progress(size(Y,1));
parfor sim = 1:size(Y,1)
    % Simulate Signal    
    [Sa(sim,:), Sphi(sim,:)]  = VoxelSim2D_do_one_adapted(Model,Seq,Param_names,Y(sim,:));
    parfor_progress;
end
parfor_progress(0);

X = Sa .* exp(1i * Sphi);







