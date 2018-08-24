
%% Generate dicos

voxpar_file  	= 'dictionaries/compare_grids/conf/voxpar_dico.txt';
sequence_file   = 'files/myconf/seqpar_GESFIDE.txt';
params_files    = {'dictionaries/compare_grids/conf/params_grid.txt',...
                   'dictionaries/compare_grids/conf/params_rand.txt',...
                   'dictionaries/compare_grids/conf/params_quasiMC.txt',...
                   'dictionaries/compare_grids/conf/params_test.txt',...
                  };
file_names      = {'Regular','RandUnif','RandquasiMC', 'Test'};
method_name     = {'Regular','RandUnif','RandquasiMC'};

path_out        = 'dictionaries/compare_grids';

%%

for f = 4:length(params_files)
    [X, Y, ~, AltModel] = GenerateSignals(voxpar_file, sequence_file, params_files{f});
    Y = Y(1:end/2,1:2); X = abs([X(1:end/2,:) X(end/2+1:end,:)]);
    save([path_out '/' date '_' method_name{f} '.mat'], 'X','Y','AltModel')
end
clear X Y AltModel


%%

dico_files  = {'dictionaries/compare_grids/23-Aug-2018_Regular.mat',...
               'dictionaries/compare_grids/23-Aug-2018_RandUnif.mat',...
               'dictionaries/compare_grids/23-Aug-2018_RandquasiMC.mat',...
              };
          
test_file   = 'dictionaries/compare_grids/24-Aug-2018_Test.mat';
          
cstr.Sigma 	= 'd';
cstr.Gammat	= 'd';
cstr.Gammaw	= '';
K        	= 20;
Lw          = 0;

load(test_file)
Xtest = X;
Ytest = Y;
clear X Y


figure

for f = 1:length(dico_files)

    load(dico_files{f})

    subplot(4,3,[f 3+f]); plot(Y(:,1), Y(:,2), '.', 'MarkerSize',10)
    xlim([0 0.16]); ylim([0 16*1e-6])
    xlabel('BVf'); ylabel('VSI')

    theta       = EstimateInverseFunction(Y, X, K, 1, 150, cstr, 0);
    
    for rep = 1:100
        ind         = randperm(size(Xtest,1),50);
        Ypredict    = EstimateParametersFromModel(Xtest(ind,:), theta, 0);
        
        for i = 1:size(Y,2)
                Ypredict(Ypredict(:,i) > max(Y(:,i)),i) = nan;
                Ypredict(Ypredict(:,i) < min(Y(:,i)),i) = nan;
        end

        [Rmse(f,:,rep), Nrmse(f,:,rep), Mae(f,:,rep)] = EvaluateEstimation(Ytest(ind,:), Ypredict, Y);
    end
    title([method_name{f} ' - Average NRMSE = ' num2str(mean(mean(Nrmse(f,:,:),3),2),2)])
end


subplot(4,3,[7 10])
hist(squeeze(mean(Nrmse,2))', 0.02:0.02:0.5)
legend(method_name)
xlim([0.02 0.5])
title('Average NRMSE histogram')

subplot(4,3,[8 9])
hist(squeeze(Mae(:,1,:))', .0005:(.007-.0005)/30:.007)
xlim([0 .007])
title('MAE histogram on BVf estimations')

subplot(4,3,[11 12])
hist(squeeze(Mae(:,2,:))', (0:1.4/30:1.4)*1e-6)
xlim([0 1.4]*1e-6)
title('MAE histogram on VSI estimations')





