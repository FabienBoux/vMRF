
%% Generate dicos

nb_repetitions  = 200;

voxpar_file  	= 'dictionaries/compare_grids/conf/voxpar_dico.txt';
sequence_file   = 'files/myconf/seqpar_GESFIDE.txt';
params_files    = {'dictionaries/compare_grids/conf/params_grid.txt',...
                   'dictionaries/compare_grids/conf/params_rand.txt',...
                   'dictionaries/compare_grids/conf/params_quasiMC.txt',...
                   'dictionaries/compare_grids/conf/params_test.txt',...
                  };
method_name     = {'Regular','RandUnif','RandquasiMC'};


%%

cstr.Sigma 	= 'd';
cstr.Gammat	= 'd';
cstr.Gammaw	= '';
K        	= 20;
Lw          = 0;


% Generate test signals
[X, Y, ~, ~] = GenerateSignals(voxpar_file, sequence_file, params_files{4});
Ytest   = Y(1:end/2,1:2);
Xtest   = abs([X(1:end/2,:) X(end/2+1:end,:)]);


% Generate grid one time
[X, Y, ~, ~] = GenerateSignals(voxpar_file, sequence_file, params_files{1});
Ygrid   = Y(1:end/2,1:2);
Xgrid   = abs([X(1:end/2,:) X(end/2+1:end,:)]);


%%
for rep = 1:nb_repetitions
    
    fprintf([num2str(rep) 'th repetition \n'])
    
    [X, Y, ~, ~] = GenerateSignals(voxpar_file, sequence_file, params_files{2});
    Yrandunif   = Y(1:end/2,1:2);
    Xrandunif   = abs([X(1:end/2,:) X(end/2+1:end,:)]);
    
    [X, Y, ~, ~] = GenerateSignals(voxpar_file, sequence_file, params_files{3});
    YrandqMC    = Y(1:end/2,1:2);
    XrandqMC    = abs([X(1:end/2,:) X(end/2+1:end,:)]);

    clear X Y 
    
    % Regular grid
    Ypred_grid  = EstimateParametersFromGrid(Xtest, Xgrid, Ygrid, 0);
    
    [Rmse{1}(rep,:), Nrmse{1}(rep,:), Mae{1}(rep,:)] = EvaluateEstimation(Ytest, Ypred_grid, Ygrid);
    
    % Random Uniform
    theta       = EstimateInverseFunction(Yrandunif, Xrandunif, K, Lw, 200, cstr, 0);
    Ypred_randunif = EstimateParametersFromModel(Xtest, theta, 0);

    for i = 1:size(Ygrid,2)
            Ypred_randunif(Ypred_randunif(:,i) > max(Ygrid(:,i)),i) = nan;
            Ypred_randunif(Ypred_randunif(:,i) < min(Ygrid(:,i)),i) = nan;
    end
        
    [Rmse{2}(rep,:), Nrmse{2}(rep,:), Mae{2}(rep,:)] = EvaluateEstimation(Ytest, Ypred_randunif, Yrandunif);
    
    % Random quasiMC
    theta       = EstimateInverseFunction(YrandqMC, XrandqMC, K, Lw, 200, cstr, 0);
    Ypred_randqMC = EstimateParametersFromModel(Xtest, theta, 0);

    for i = 1:size(Ygrid,2)
            Ypred_randqMC(Ypred_randqMC(:,i) > max(Ygrid(:,i)),i) = nan;
            Ypred_randqMC(Ypred_randqMC(:,i) < min(Ygrid(:,i)),i) = nan;
    end
        
    [Rmse{3}(rep,:), Nrmse{3}(rep,:), Mae{3}(rep,:)] = EvaluateEstimation(Ytest, Ypred_randqMC, YrandqMC);   
end


%%

figure
subplot(2,2,[1 3])
hist([mean(Nrmse{1},2) mean(Nrmse{2},2) mean(Nrmse{3},2)], 30)
legend(method_name)
title('Average NRMSE')

subplot(2,2,2)
hist([mean(Mae{1}(:,1),2) mean(Mae{2}(:,1),2) mean(Mae{3}(:,1),2)], 30)
legend(method_name)
title('Average MAE - BVf')

subplot(2,2,4)
hist([mean(Mae{1}(:,2),2) mean(Mae{2}(:,2),2) mean(Mae{3}(:,2),2)], 30)
legend(method_name)
title('Average MAE - VSI')



