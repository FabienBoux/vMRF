
addpath(genpath(fullfile(pwd(), 'functions')))
addpath(genpath(fullfile(pwd(), 'tools')))


%% Learning step

% complete 3 sequences
% file        = 'dictionaries/multidico/2018-07-25-11:22_complex_dico.mat';
% only GESFIDSE sequence
file        = 'dictionaries/multidico/2018-07-25-09:35_32-samples_5-parameters_80000-signals.mat';
% file        = 'dictionaries/mydicos/dico_random_220000signals_4parameters.mat';
coord       = [1 2];

v = randperm(30,16);
GESFIDSE    = zeros(1,32); %GESFIDSE(11:end)  = 1;
% MSME        = zeros(1,32); MSME(v)      = 1;
% MGE         = zeros(1,30); MGE(v(1:end-1)) = 1;
% remove_samples = logical(repmat([GESFIDSE MSME MGE],1,2));
remove_samples = logical(repmat(GESFIDSE,1,2));

K           = 6;
Lw          = 0;

cstr.Sigma  = 'd';
cstr.Gammat = 'd';
cstr.Gammaw = '';

load(file)
X           = X(:,~remove_samples);
snr         = [50 200];
X           = AddNoise(X, snr(1)+(snr(2)-snr(1))*rand(size(X,1),1), 0);
[theta, r]	= EstimateInverseFunction(Y(:,[1 2]), X, K, Lw, 100, cstr, 0);


%% Choosing test signal

nb_test = 15000;
ind     = randi([1 length(X)], 1,nb_test);
y       = X(ind,:)';
x       = Y(ind,:);

[x_predict,alpha]   = gllim_inverse_map(y,theta,0);
x_predict           = x_predict';

% Estimating model error
errorS = zeros(size(theta.Sigma(:,:,1)));
for k = 1:K
    errorS = errorS + alpha(k) * theta.Sigma(:,:,k);
end

errorG = zeros(nb_test, length(theta.Gamma(:,:,1)));
for k = 1:K
    for c = 1:length(coord)
        errorG(:,c) = errorG(:,c) + alpha(:,k) * theta.Gamma(c,c,k);
    end
end
quad_error = (x(:,1:length(coord)) - x_predict(:,1:length(coord))).^2;


[~,loc] = max(alpha,[],2);
for k = 1:K
    error_expected{k}   = errorG(loc==k,1:length(coord));
    error_real{k}       = abs(x(loc==k,1:length(coord)) - x_predict(loc==k,1:length(coord)));
    [Rmse{k}, Nrmse{k}, Mae{k}]	= EvaluateEstimation(x(loc==k,1:length(coord)),  x_predict(loc==k,:), x);
end
figure
for k = 1:K
    subplot(5,4,2*k-1)
    hist(error_expected{k}(:,1),(0.1:0.01:5)*1e-3)
    title(num2str(mean(error_expected{k}(:,1)),2))
    subplot(5,4,2*k)
    hist(error_real{k}(:,1),0:0.00005:0.005)
    title(num2str(mean(error_real{k}(:,1)),2))
    xlim([0 0.005])
end



%% Interpretation - test 1

colors  = [	0           0.4470      0.7410;            	0.8500      0.3250      0.0980;              	0.9290      0.6940      0.1250;                0.4940      0.1840      0.5560;                0.4660      0.6740      0.1880;                0.3010      0.7450      0.9330;                0.6350      0.0780      0.1840;                0           0           1.0000;                0           0.5000   	0;                1.0000         0    	0;                0           0.7500      0.7500;                0.7500      0           0.7500;                0.7500      0.7500      0;                0.2500      0.2500      0.2500;                ];
colors  = repmat(colors, 3,1);

% plotting: 1.space mapping, 2.gamma by k-space, 3.real estimation error
figure
subplot(221); hold on
[~,loc] = max(alpha,[],2);
for k = 1:K
    plot(x(loc == k,1), x(loc == k,2), '.')
end
legend(num2str((1:K)'))

subplot(222); hold on
for k = 1:K
    plot(squeeze(theta.Gamma(1,1,k)),squeeze(theta.Gamma(2,2,k)),'o','Color',colors(k,:)); 
end
legend(num2str((1:K)'))

subplot(223); hold on
for k = 1:K
    plot(error_real{k}(:,1), error_real{k}(:,2), '.', 'Color',colors(k,:))
end
legend(num2str((1:K)'))

xl      = xlim;
yl      = ylim;
subplot(224)
hold on
for k = 1:K
    xtmp    = error_real{k}(:,1);
    ytmp    = error_real{k}(:,2);
    mu      = [mean(xtmp) mean(ytmp)];
    sigma   = 5*cov(xtmp,ytmp);
    
    x1      = xl(1):.001:xl(2);
    x2      = yl(1):.00000001:yl(2);
    [X1,X2] = meshgrid(x1,x2);
    F       = mvnpdf([X1(:) X2(:)],mu,sigma);
    F       = reshape(F,length(x2),length(x1));
    
    %contour(x1,x2,F,3, 'Color',colors(k,:),'LineWidth',0.2);
    contour(x1,x2,F,1, 'Color',colors(k,:),'LineWidth',2.5);
    xlim(xl)
    ylim(yl)
end
%xlabel('x'); ylabel('y');
%line([0 0 1 1 0],[1 0 0 1 1],'linestyle','--','color','k');
% axis square;
% grid on;

for k =1:K
	meank(k,:)  = mean(error_real{k});
    gammak(k,:) = diag(theta.Gamma(:,:,k))';
end


%% Interpretation - test 2

figure
% plotting Sigma by k-space
subplot(3,3,7); hold on
[~,loc] = max(alpha,[],2);
for k = 1:K
    plot(x(loc == k,1), x(loc == k,2), '.')
end
xlim([min(x(:,1)) max(x(:,1))])
ylim([min(x(:,2)) max(x(:,2))])
legend(num2str((1:K)'))

subplot(3,3,[1:3])
hold on
for k = 1:K
    tmp     = y(:,loc ==k);
    plot(tmp(:,randperm(size(tmp,2), 5)), 'LineWidth', 1.5, 'Color', colors(k,:))
end
clear tmp
%legend(['Typical signal: ' num2str((1:K)')])


subplot(3,3,[4:6])
hold on
for k = 1:K
    meank(k,:)  = mean(error_real{k});
    nrmsek(k)   = mean(Nrmse{k});
    plot(diag(theta.Sigma(:,:,k)), 'LineWidth', 2)
end
legend([num2str((1:K)') repelem(' - NRMSE = ',K,1) num2str(nrmsek')])
title(['Mean NRMSE = ' num2str(mean(nrmsek))])


%%

xbvf        = 0.001:0.001:0.3;
xvsi        = (0.01:0.04:12) *1e-6;
x_samples	= [xbvf; xvsi];

ind = randi(1000,1,4);
ind = [586 287 888 747 366]; %621 880 747

figure
for i = 1:length(ind)
    yt = X(ind(i),:);
    xt = Y(ind(i),:);
    
    [x_dens,psi,theta]  = gllim_inverse_dens_modified(yt',theta,ones(size(yt)),x_samples,0);
    [xt_predict,alpha]  = gllim_inverse_map(yt',theta,0);
    x1      = 0:.001:.3;
    x2      = (0:.001:13) *1e-6;
    [X1,X2] = meshgrid(x1,x2);
    Ftot    = zeros(size(X1));
    Ftot_w  = zeros(size(X1));
    for k =1:K
        F       = mvnpdf([X1(:) X2(:)],(theta.A_s(:,:,k)*yt'+theta.b_s(:,k))',theta.Sigma_s(:,:,k));
        F       = reshape(F,length(x2),length(x1));
        set(gca,'YDir','normal')
%         Ftot    = Ftot + F; 
        Ftot_w  = Ftot_w + alpha(k)* F;
    end
%     subplot(2,length(ind),2*i-1)
%     imagesc(x1,x2,Ftot)
%     set(gca,'YDir','normal')
    Ftot_w = Ftot_w ./ sum(sum(Ftot_w));
    
    subplot(2,length(ind),i)
    imagesc(x1,x2,Ftot_w)
    hold on
    plot(xt(1),xt(2), 'gx')
    plot(xt_predict(1), xt_predict(2), 'rx')
    set(gca,'YDir','normal')
    
    subplot(2,length(ind),length(ind)+i)
    contour(x1,x2,Ftot_w,[1e-5 1e-9 1e-13 1e-17]);
    hold on
    plot(xt(1),xt(2), 'gx')
    plot(xt_predict(1), xt_predict(2), 'rx')
end



%%
% figure
% plot(quad_error)
% hold on; plot(errorG)

disp(mean(quad_error))
disp(mean(errorG))

% % Generating distribution and estimation
% [X1,X2] = meshgrid(xbvf,xvsi);
% F       = zeros(size(X1));
% for k = 1:K
%     F = F + alpha(k) * reshape(mvnpdf([X1(:) X2(:)],psi.mu(:,k)',psi.S(:,:,k)),length(xvsi),length(xbvf));
% end
% 
% 
% % Plotting
% figure
% subplot(121)
% surf(xbvf,xvsi,F);
% title(['Predict = ' num2str(x_predict(1),3) ' / ' num2str(x_predict(2),3)...
%     '\newline Real = ' num2str(Y(ind,1),3) ' / ' num2str(Y(ind,2),3)])
% 
% subplot(122)
% imagesc(xbvf,xvsi,F)
% hold on
% plot(x_predict(1), x_predict(2),'rx','LineWidth', 2)
% plot(Y(ind,1), Y(ind,2),'gx','LineWidth', 2)
% legend({'Prediction','Real value'})
% set(gca,'YDir','normal')





