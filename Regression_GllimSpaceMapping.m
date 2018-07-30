% 23/01/18 
% Le but de ce script est d'étudier le paramètre d'erreur r renvoyé par
% gllim au moment de l'estimation de theta



% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_dico       = 'files/dico_random_220000signals.mat';

signal_number   = 40000;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(file_dico)

K           = 2:13;
Lt          = 0;

rand_perm   = randperm(length(X), signal_number);

Xtrain      = abs(X(rand_perm, :));
Ytrain      = Y(rand_perm, :);

colors      = [	0           0.4470      0.7410;
            	0.8500      0.3250      0.0980;
              	0.9290      0.6940      0.1250;
                0.4940      0.1840      0.5560;
                0.4660      0.6740      0.1880;
                0.3010      0.7450      0.9330;
                0.6350      0.0780      0.1840;
                0           0           1.0000;
                0           0.5000   	0;
                1.0000         0    	0;
                0           0.7500      0.7500;
                0.7500      0           0.7500;
                0.7500      0.7500      0;
                0.2500      0.2500      0.2500;
                ];
colors      = [colors; colors];

% Ytrain = [Ytrain zeros(length(Ytrain), 1)];
% for i = 1:length(Xtrain)
%     Ytrain(i,3)   = rand* 100;
%     Xtrain(i,:)   = AddNoise(Xtrain(i,:), Ytrain(i,3), 0);
% end

for elt = 1:length(K)
    cstr.Sigma = 'd';
    [theta, r]  = EstimateInverseFunction(Ytrain, Xtrain, K(elt), Lt, 30, cstr, 0);

    [~, mm] = max(r');

    subplot(3,4,elt)
    hold on
    for i = 1:size(r,2)
        plot(Ytrain(mm == i, 1), Ytrain(mm == i, 2), '.', 'Color',colors(i,:))
    end
%     for i = 1:size(r,2)
%         text(mean(Ytrain(mm == i, 1)), mean(Ytrain(mm == i, 2)), num2str(theta.Sigma(1,1,i),3))
%     end
    %legend(string(squeeze(theta.Sigma(1,1,:))))
    xlim([0.01 0.25]); ylim([1e-6 1e-5])
    title('K = ' + string(K(elt)))
end






