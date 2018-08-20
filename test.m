

% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_grid   	= 'dictionaries/multidico/2018-07-25-13:41_complex_dico.mat';
file_regression	= 'dictionaries/multidico/2018-07-25-11:22_complex_dico.mat';
seq_names       = {'GEFIDSE','MGE','MSME'};

cstr.Sigma      = 'd';
cstr.Gammat     = 'd';
cstr.Gammaw     = '';
Lw              = 0;
K               = 20;

coord           = [1 2 3 4];
seq_sizes       = [32 32 30];
signal_test     = 1000;
nb_tests        = 100;

snr_train       = [0 inf];
snr_test        = [0 inf];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Load dico
load(file_regression)
Xregr = X;
Yregr = Y(:,coord);

load(file_grid)
Xdico = X;
Ydico = Y(:,coord);
clear X Y

X           = [Xregr; Xdico];
Y           = [Yregr; Ydico];


%% Convert seq_sizes to logical vectors
clear mat nmat
vect    = ones(1,length(seq_sizes));
mat     = vect;
for i = 1:length(seq_sizes)-1
    vect(1:i)   = 0;
    mat         = [mat; unique(perms(vect), 'rows')];
end
for i = 1:size(mat,1), nmat(i,:) = repelem(mat(i,:), seq_sizes); end
%nmat   	= logical(repmat(nmat, 1,2));
nmat   	= logical([nmat zeros(size(nmat))]);



for seq = 1:length(seq_sizes)
    
    fprintf(['Iter ' num2str(seq) '\n'])
    
    % Select trainning data and train model
    Xtrain              = X(:,nmat(seq,:));
%     Xtrain              = AddNoise(Xtrain, snr_train(1)+(snr_train(2)-snr_train(1))*rand(size(Xtrain,1),1), 0);
    Ytrain              = Y;
    [theta,~]           = EstimateInverseFunction(Ytrain, Xtrain, K, Lw, 100, cstr, 0);
    
    subplot(3,3,seq)
    
    hold on
    for c = 1:length(coord)
        plot(squeeze(abs(theta.A(:,c))) ./ max(abs(squeeze(theta.A(:,c)))), 'linewidth',1.5)
    end
    legend({'BVf', 'VSI', 'T2', 'ADC'})
end