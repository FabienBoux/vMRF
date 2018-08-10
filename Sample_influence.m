% SAMPLE_INFLUENCE
%
% We are investigated on the influence of the sample number
%
% Note:
%
% Fabien Boux - 08/2018

addpath(genpath(fullfile(pwd(), 'functions')))
addpath(genpath(fullfile(pwd(), 'tools')))


%%% Parameters
%sequence
file        = 'dictionaries/multidico/2018-07-25-11:22_complex_dico.mat';
coord       = [1 2];
snr         = [50 200];

%regression
K           = 15;
Lw          = 0;
cstr.Sigma  = 'd';
cstr.Gammat = 'd';
cstr.Gammaw = '';


%%% Load
load(file)
N       = size(X,2);
L       = length(coord);
X   	= AddNoise(X, snr(1)+(snr(2)-snr(1))*rand(size(X,1),1), 0);


%%% Initialization
M0      = 10;
prob    = .99;
tol   	= [0.01 1e-6];
limit   = 100;


%%% Loop
O       = 100;
acc     = inf(O,L);
iter    = 1;
M       = zeros(1,limit);
ind     = randperm(N,N);

x       = X(randperm(size(X,1),O),:);
while all(mean(acc) > tol) || (iter < limit)
    
    %learning
    if iter == 1
        M(iter) = M0;
    else
        M(iter) = ceil(M(iter-1) * lambda);
    end
    fprintf('Iteration %3d - M = %d \n', iter, M(iter))
    
    %laerning
    [theta, r]  = EstimateInverseFunction(Y(:,coord), X(:,ind(1:M(iter))), K, Lw, 100, cstr, 0);
    
    x_samples   = [0.001:0.001:0.3; (0.01:0.04:12) *1e-6];
    x1          = 0:.0005:.3;
    x2          = (0:.05:13) *1e-6;
    [X1,X2]     = meshgrid(x1,x2);
    
    for obs = 1:O
        %estimation step
        [x_dens, psi, theta]    = gllim_inverse_dens_modified(x(obs,ind(1:M(iter)))', theta, ones(size(x(obs,ind(M(iter)))')), x_samples, 0);
        [x_predict, alpha]      = gllim_inverse_map(x(obs,ind(1:M(iter)))',theta,0);
        
        %accuracy computing
        D       = zeros(size(X1));
        Dtot    = D;
        for k = 1:size(theta.A_s,3) %should be K but it doesn't work and i don't know why ??????
            D       = mvnpdf([X1(:) X2(:)],(theta.A_s(:,:,k)*x(obs,ind(1:M(iter)))'+theta.b_s(:,k))',theta.Sigma_s(:,:,k));
            D       = reshape(D,length(x2),length(x1));
            Dtot    = Dtot + alpha(k)*D;
        end
        D       = D ./ sum(sum(D));
        acc(obs,:) = ComputeAccuracy({x1, x2}, D, prob, 0);
    end
    
    Acc(iter,1) = mean(acc(:,1));
    Acc(iter,2) = mean(acc(:,2));
    
    fprintf('\t %.4f / %f \n', Acc(iter,1), Acc(iter,2))
    
    %update
    lambda  = 1.2; %TODO
    iter    = iter + 1;
end









