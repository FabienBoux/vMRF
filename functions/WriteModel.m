function [] = WriteModel(FileDico, FileOut, ModelVarName)

if ~exist(ModelVarName, 'var')
    ModelVarName = 'AltModel';
end

load(FileDico);
Model = eval(ModelVarName);

fields = fieldnamesr(Model, 'prefix');

for f = 1:length(fields)
    lines{f} = [fields{f} ' = [' num2str(eval(fields{f})) '];'];
end

fid = fopen(FileOut,'w');
for l = 1:length(lines)
    fprintf(fid, '%s\n', lines{l});
end
fclose(fid);