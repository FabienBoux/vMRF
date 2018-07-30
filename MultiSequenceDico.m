
sequence_files = {'files/myconf/seqpar_GESFIDE.txt',...
    'files/myconf/seqpar_MGE.txt',...
    'files/myconf/seqpar_MSME.txt'};

path_out = './files/multidico';

GenerateComplexDico('files/myconf/voxpar_dico.txt', sequence_files, 'files/myconf/params_rand.txt', path_out)
GenerateComplexDico('files/myconf/voxpar_dico.txt', sequence_files, 'files/myconf/params_grid.txt', path_out)