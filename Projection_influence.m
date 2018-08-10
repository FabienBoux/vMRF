% PROJECTION_INFLUENCE
%
% We are investigating on the projection of each subspace. Some
% interrested plots here
%
% Note: 
%
% Fabien Boux - 08/2018

addpath(genpath(fullfile(pwd(), 'functions')))
addpath(genpath(fullfile(pwd(), 'tools')))

%% Learning step

% 3 sequences
% file_rand   = 'dictionaries/multidico/2018-07-25-11:22_complex_dico.mat';
% file_grid   = 'dictionaries/multidico/2018-07-25-13:41_complex_dico.mat';
% only MGEFIDSE sequence
file_rand   = 'dictionaries/multidico/2018-07-25-09:35_32-samples_5-parameters_80000-signals.mat';
file_grid   = 'dictionaries/multidico/2018-07-25-11:58_32-samples_5-parameters_80000-signals.mat';
coord       = [1 2 3 4];

K           = 4;
Lw          = 0;

cstr.Sigma  = 'd';
cstr.Gammat = 'd';
cstr.Gammaw = '';

load(file_rand)
Xr          = X;
Yr          = Y;
load(file_grid)
Xg          = X;
Yg          = Y;

X           = [Xr; Xg];
Y           = [Yr; Yg];
% snr         = [50 200];
% X           = AddNoise(X, snr(1)+(snr(2)-snr(1))*rand(size(X,1),1), 0);
[theta, r]	= EstimateInverseFunction(Y(:,coord), X, K, Lw, 100, cstr, 0);


%% First plot

colors  = [	0           0.4470      0.7410;            	0.8500      0.3250      0.0980;              	0.9290      0.6940      0.1250;                0.4940      0.1840      0.5560;                0.4660      0.6740      0.1880;                0.3010      0.7450      0.9330;                0.6350      0.0780      0.1840;                0           0           1.0000;                0           0.5000   	0;                1.0000         0    	0;                0           0.7500      0.7500;                0.7500      0           0.7500;                0.7500      0.7500      0;                0.2500      0.2500      0.2500;                ];

figure
subplot(2*K,2,1:2:2*K)
hold on
[~, mm] = max(r');
for k = 1:K
    plot(Y(mm == k, 1), Y(mm == k, 2), '.', 'Color',colors(k,:))
    leg{k}      = ['Subspace ' num2str(k)];
    
    meank(k,:)  = [mean(Y(mm == k, 1)) mean(Y(mm == k, 2))];
    maxk(k,:)   = [max(Y(mm == k, 1))  max(Y(mm == k, 2))];
    mink(k,:)   = [min(Y(mm == k, 1))  min(Y(mm == k, 2))];
end
xlim([0.01 0.25]); ylim([1e-6 1e-5])
legend(leg)
title('Space mapping')


nb_signal   = 10;
for k = 1:K
    subplot(2*K,2,2*k)
    hold on
    
    loc = find(mm == k);
    plot(X(randperm(length(loc), nb_signal),:)', 'Color',colors(k,:))
    title(['Signals from subspace ' num2str(k)])
end


%% Second plot

figure
subplot(2*K,2,1:2:2*K)
hold on
[~, mm] = max(r');
for k = 1:K
    plot(Y(mm == k, 1), Y(mm == k, 2), '.', 'Color',colors(k,:))
    leg{k}      = ['Subspace ' num2str(k)];
    
    meank(k,:)  = [mean(Y(mm == k, 1)) mean(Y(mm == k, 2))];
    stdk(k,:)   = [std(Y(mm == k, 1))  std(Y(mm == k, 2))];
    
    stdk_signal(k,:)   = std(X(mm == k, :));
end
xlim([0.01 0.2]); ylim([1e-6 1e-5])
legend(leg)
title('Space mapping')

for k = 1:K    
    subplot(2*K,2,2*k)
    
    plot(squeeze(abs(theta.A(:,1,k))) ./ max(abs(squeeze(theta.A(:,1,k)))), 'linewidth',1.5)
    hold on
    plot(squeeze(abs(theta.A(:,2,k))) ./ max(abs(squeeze(theta.A(:,2,k)))), 'linewidth',1.5)
    title(['Subspace ' num2str(k)])
    legend({'BVf', 'VSI'})
end



vect    = unique(Yg(:,3));
t2      = vect(ceil(length(vect)/4));
vect    = unique(Yg(:,4));
adc     = vect(ceil(length(vect)/2));

for k = 1:K
    v       = unique(Yg(:,2));
    [~,l]   = min(abs(v - meank(k,2)));
    vsi     = v(l);
    vect    = (Yg(:,2) == vsi) & (Yg(:,3) == t2) & (Yg(:,4) == adc);
    vect(vect) = (Yg(vect,1) <= meank(k,1)+1.5*stdk(k,1)) & (Yg(vect,1) >= meank(k,1)-1.5*stdk(k,1));
%     vect(vect) = (Yg(vect,1) <= maxk(k,1)) & (Yg(vect,1) >= mink(k,1));
    
    subplot(2*K,2,2*K + 2*k)
    plot(Xg(vect,:)', 'Color',colors(k,:))
    v = Yg(vect,1);
    title(['VSI = ' num2str(vsi,2) '   -   BVf \in [' num2str([v(1) v(end)],4) ']'])
    
    subplot(2*K,2,1:2:2*K)
    hold on
    plot(Yg(vect,1), Yg(vect,2), 'k.', 'markersize',12)
end


for k = 1:K
    v       = unique(Yg(:,1));
    [~,l]   = min(abs(v - meank(k,1)));
    bvf     = v(l);
    vect    = (Yg(:,1) == bvf) & (Yg(:,3) == t2) & (Yg(:,4) == adc);
    vect(vect) = (Yg(vect,2) <= meank(k,2)+1.5*stdk(k,2)) & (Yg(vect,2) >= meank(k,2)-1.5*stdk(k,2));
%     vect(vect) = (Yg(vect,2) <= maxk(k,2)) & (Yg(vect,2) >= mink(k,2));
    
    subplot(2*K,2,2*K -1+ 2*k)
    plot(Xg(vect,:)', 'Color',colors(k,:))
    v = Yg(vect,2);
    title(['BVf = ' num2str(bvf,3) '   -   VSI \in [' num2str([v(1) v(end)],2) ']'])
    
    subplot(2*K,2,1:2:2*K)
    hold on
    plot(Yg(vect,1), Yg(vect,2), 'k.', 'markersize',12)
end


%% Third plot

for k = 1:K
    subplot(K,2,2*k-1)
    hold on
    plot(squeeze(abs(theta.A(:,1,k))) ./ max(abs(squeeze(theta.A(:,1,k)))), 'k', 'linewidth',1.5)
%     plot(stdk_signal(k,:) ./ max(stdk_signal(k,:)), 'm', 'linewidth',1.5)
    subplot(K,2,2*k)
    hold on
    plot(squeeze(abs(theta.A(:,2,k))) ./ max(abs(squeeze(theta.A(:,2,k)))), 'k', 'linewidth',1.5)
%     plot(stdk_signal(k,:) ./ max(stdk_signal(k,:)), 'm', 'linewidth',1.5)
end




vect    = unique(Yg(:,3));
t2      = vect(ceil(length(vect)/4));
vect    = unique(Yg(:,4));
adc     = vect(ceil(length(vect)/2));

for k = 1:K
    v       = unique(Yg(:,2));
    [~,l]   = min(abs(v - meank(k,2)));
    vsi     = v(l);
    vect    = (Yg(:,2) == vsi) & (Yg(:,3) == t2) & (Yg(:,4) == adc);
    vect(vect) = (Yg(vect,1) <= meank(k,1)+1.5*stdk(k,1)) & (Yg(vect,1) >= meank(k,1)-1.5*stdk(k,1));
    
    subplot(K,2,2*k-1)
    plot(Xg(vect,:)', 'Color',colors(k,:))
    v = Yg(vect,1);
    title(['VSI = ' num2str(vsi,2) '   -   BVf \in [' num2str([v(1) v(end)],4) ']'])
    xlim([0 length(X(1,:))])
end


for k = 1:K
    v       = unique(Yg(:,1));
    [~,l]   = min(abs(v - meank(k,1)));
    bvf     = v(l);
    vect    = (Yg(:,1) == bvf) & (Yg(:,3) == t2) & (Yg(:,4) == adc);
    vect(vect) = (Yg(vect,2) <= meank(k,2)+1.5*stdk(k,2)) & (Yg(vect,2) >= meank(k,2)-1.5*stdk(k,2));
    
    subplot(K,2,2*k)
    plot(Xg(vect,:)', 'Color',colors(k,:))
    v = Yg(vect,2);
    title(['BVf = ' num2str(bvf,3) '   -   VSI \in [' num2str([v(1) v(end)],2) ']'])
    xlim([0 length(X(1,:))])
end



