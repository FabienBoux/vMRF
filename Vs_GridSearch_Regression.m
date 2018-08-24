
file_grid       = 'dictionaries/compare_methods/23-Aug-2018_Regular.mat';
file_trainning  = 'dictionaries/compare_methods/23-Aug-2018_RandquasiMC.mat';
file_test       = 'dictionaries/compare_methods/24-Aug-2018_Test.mat';

snr_levels      = [10 20 50 100];
generation      = 0;
coord           = 1:3;


%%
if generation == 1
    voxpar_file     = 'dictionaries/compare_methods/conf/voxpar_dico.txt';
    sequence_file   = 'files/myconf/seqpar_GESFIDE.txt';
    
    [X,Y,~,~] = GenerateSignals(voxpar_file, sequence_file, 'dictionaries/compare_methods/conf/params_grid.txt');
    Y = Y(1:end/2,coord); X = abs([X(1:end/2,:) X(end/2+1:end,:)]);
    save(file_grid,'X','Y')
    
    [X,Y,~,~] = GenerateSignals(voxpar_file, sequence_file, 'dictionaries/compare_methods/conf/params_quasiMC.txt');
    Y = Y(1:end/2,coord); X = abs([X(1:end/2,:) X(end/2+1:end,:)]);
    save(file_trainning,'X','Y')
    
    [X,Y,~,~] = GenerateSignals(voxpar_file, sequence_file, 'dictionaries/compare_methods/conf/params_test.txt');
    Y = Y(1:end/2,coord); X = abs([X(1:end/2,:) X(end/2+1:end,:)]);
    save(file_test,'X','Y')
end


% Load file generated (or pre-generated if 'generation'=0)
load(file_test)
Xtest = [X(:,end/2+1:end) X(:,1:end/2)];
Ytest = Y;

load(file_trainning)
Xtrain = [X(:,end/2+1:end) X(:,1:end/2)];
Ytrain = Y;

load(file_grid)
Xgrid = [X(:,end/2+1:end) X(:,1:end/2)];
Ygrid = Y;
clear X Y


for snr = 1:length(snr_levels)
    
    Xtrain_noisy = AddNoise(Xtrain, snr_levels(snr), 0);
    
    Xtest_noisy = AddNoise(Xtest, snr_levels(snr), 0);
    real_snr(snr,:) = max(Xtest_noisy,[],2) ./ std(Xtest_noisy - Xtest,[],2);
    
    [Rslt{snr}.Estimation, Rslt{snr}.Parameters, Rslt{snr}.Errors] = EstimateParametersBothMethod(Xtrain_noisy, Ytrain, Xgrid, Ygrid, Xtest_noisy, Ytest);
    
end


%%

% figure
% for snr = 1:length(snr_levels)
%     err_r{snr} = (Rslt{snr}.Estimation.Regression.Y - Rslt{snr}.Estimation.Real.Y);
%     err_g{snr} = (Rslt{snr}.Estimation.GridSearch.Y - Rslt{snr}.Estimation.Real.Y);
%     
%     subplot(2,length(snr_levels),snr)
%     hist([err_g{snr}(:,1) err_r{snr}(:,1)], -.02:.002:.02)
%     xlim([-.02 .02])
%     legend('Grid Search','Regression')
%     title(['Snr = ' num2str(mean(real_snr(snr),2),3)])
%     
%     subplot(2,length(snr_levels),length(snr_levels)+snr)
%     hist([err_g{snr}(:,2) err_r{snr}(:,2)], (-.5:0.05:.5)*1e-5)
%     xlim([-.5 .5] *1e-5)
%     legend('Grid Search','Regression')
% end


figure
for snr = 1:length(snr_levels)
    err_r{snr} = (Rslt{snr}.Estimation.Regression.Y - Rslt{snr}.Estimation.Real.Y) ./ Rslt{snr}.Estimation.Real.Y;
    err_g{snr} = (Rslt{snr}.Estimation.GridSearch.Y - Rslt{snr}.Estimation.Real.Y) ./ Rslt{snr}.Estimation.Real.Y;
    
    err2_r{snr} = (Rslt{snr}.Estimation.Regression.Y - Rslt{snr}.Estimation.Real.Y).^2 ./ Rslt{snr}.Estimation.Real.Y.^2;
    err2_g{snr} = (Rslt{snr}.Estimation.GridSearch.Y - Rslt{snr}.Estimation.Real.Y).^2 ./ Rslt{snr}.Estimation.Real.Y.^2;
    
    
    subplot(2,length(snr_levels),snr)
    hist([mean(err_g{snr},2) mean(err_r{snr},2)], -.8:.05:.8)
    xlim([-.8 .8])
    legend('Grid Search','Regression')
    title(['Snr = ' num2str(mean(real_snr(snr),2),3)])
    if snr == 1, ylabel('nAE'); end 
    
    
    subplot(2,length(snr_levels),4+snr)
    hist([mean(err2_g{snr},2) mean(err2_r{snr},2)], -.005:.005:.25)
    xlim([-.005 .255])
    legend('Grid Search','Regression')
    title(['Snr = ' num2str(mean(real_snr(snr),2),3)])
    if snr == 1, ylabel('nQE'); end 
end
















