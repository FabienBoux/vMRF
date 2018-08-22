% CONFIDENCE_INDEX_V3
%
% The objective of this script is to idenbtify if it's possible to predict
% the quality of the estimations.
%
% Note:
%
% Fabien Boux - 08/2018


addpath(genpath(fullfile(pwd(), 'functions')))
addpath(genpath(fullfile(pwd(), 'tools')))


%% Learning step

% complete 3 sequences
file        = 'dictionaries/multidico/18-08-21/2018-08-22-09:53_32-samples_4-parameters_36000-signals.mat';
% only GESFIDSE sequence
% file        = 'dictionaries/multidico/2018-07-25-09:35_32-samples_5-parameters_80000-signals.mat';

coord       = [1 2 3];
snr         = [100 500];

K           = 20;
Lw          = 0;
cstr.Sigma  = 'd';
cstr.Gammat = 'd';
cstr.Gammaw = '';


load(file)
X           = AddNoise(X, snr(1)+(snr(2)-snr(1))*rand(size(X,1),1), 0);
[theta, r]	= EstimateInverseFunction(Y(:,coord), X, K, Lw, 100, cstr, 0);



%% Test 1

xbvf        = 0.001:0.001:0.3;
xvsi        = (0.01:0.04:12) *1e-6;
xt2         = (1:1:300) *1e-3;
%xadc        = (503:3:1400) *1e-12;
% x_samples	= [xbvf; xvsi; xt2; xadc];
x_samples	= [xbvf; xvsi; xt2];

ind         = randi(1000,1,1); %278

for i = 1:length(ind)    
    yt = X(ind(i),:);
    xt = Y(ind(i),:);
    
    [x_dens,psi,theta]  = gllim_inverse_dens_modified(yt',theta,ones(size(yt')),x_samples,0);
    [xt_predict,alpha]  = gllim_inverse_map(yt',theta,0);
    x1      = 0:.0001:.3;
    x2      = (0:.005:13) *1e-6;
    [X1,X2] = meshgrid(x1,x2);
    Ftot    = zeros(size(X1));
    for k = 1:K
        F       = mvnpdf([X1(:) X2(:)],(theta.A_s(1:2,:,k)*yt'+theta.b_s(1:2,k))',theta.Sigma_s(1:2,1:2,k)');
        F       = reshape(F,length(x2),length(x1));
        Ftot  = Ftot + alpha(k)* F;
    end
    Ftot = Ftot ./ sum(sum(Ftot));
    
%     subplot(221)
%     imagesc(x1,x2,Ftot)
%     hold on
%     plot(xt(1),xt(2), 'gx')
%     plot(xt_predict(1), xt_predict(2), 'rx')
%     set(gca,'YDir','normal')
    
    figure
    
    subplot(221)
    contour(x1,x2,Ftot,[1e-5 1e-9 1e-13 1e-17]);
    hold on
    plot(xt(1),xt(2), 'gx')
    plot(xt_predict(1), xt_predict(2), 'rx')
    
    xl = xlim;
    yl = ylim;
    
    [~,l1] = min(abs(x1 - xt_predict(1)));    
    [~,l2] = min(abs(x2 - xt_predict(2)));
    
    subplot(223)
    plot(x1, Ftot(l2,:), 'LineWidth', 1.5)
    hold on
    line('XData', [xt_predict(1) xt_predict(1)], 'YData', ylim, 'Color','r')
    line('XData', [xt(1) xt(1)], 'YData', ylim, 'Color','g')
    xlim(xl)
    
    subplot(222)
    plot(Ftot(:,l1), x2, 'LineWidth', 1.5)
    hold on
    line('YData', [xt_predict(2) xt_predict(2)], 'XData', xlim, 'Color','r')
    line('YData', [xt(2) xt(2)], 'XData', xlim, 'Color','g')
    ylim(yl)
end


%% Test 1 bis

xbvf        = 0.001:0.001:0.3;
xvsi        = (0.01:0.04:12) *1e-6;
xt2         = (1:1:300) *1e-3;
%xadc        = (503:3:1400) *1e-12;
% x_samples	= [xbvf; xvsi; xt2; xadc];
x_samples	= [xbvf; xvsi; xt2];

ind         = randi(1000,1,10);
% ind         = [446 160 134 983 1000 674];

for i = 1:length(ind)    
    yt = X(ind(i),:);
    xt = Y(ind(i),:);
    
    [x_dens,psi,theta]  = gllim_inverse_dens_modified(yt',theta,ones(size(yt')),x_samples,0);
    [xt_predict,alpha]  = gllim_inverse_map(yt',theta,0);
    x1      = 0:.001:.3;
    x2      = (0:.01:13) *1e-6;
    [X1,X2] = meshgrid(x1,x2);
    Ftot    = zeros(size(X1));
    for k = 1:K
        F       = mvnpdf([X1(:) X2(:)],(theta.A_s(1:2,:,k)*yt'+theta.b_s(1:2,k))',theta.Sigma_s(1:2,1:2,k)');
        F       = reshape(F,length(x2),length(x1));
        Ftot  = Ftot + alpha(k)* F;
    end 
    Ftot = Ftot ./ max(max(Ftot));
    
    
    figure    
    ha(1) = subplot(221);
    contour(x1,x2,Ftot,[.9 .5 .1 .001 .00001]);
    hold on
    plot(xt(1),xt(2), 'gx', 'LineWidth', 2)
    plot(xt_predict(1), xt_predict(2), 'rx', 'LineWidth', 2)
    
    
    xl = xlim;
    yl = ylim;
    
    ha(3) = subplot(223);
    plot(x1, sum(Ftot,1), 'LineWidth', 2)
    hold on
    line('XData', [xt_predict(1) xt_predict(1)], 'YData', ylim, 'Color','r', 'LineWidth', 1.5)
    line('XData', [xt(1) xt(1)], 'YData', ylim, 'Color','g', 'LineWidth', 1.5)
    xlim(xl)
    
    
    ha(2) = subplot(222);
    plot(sum(Ftot,2), x2, 'LineWidth', 2)
    hold on
    line('YData', [xt_predict(2) xt_predict(2)], 'XData', xlim, 'Color','r', 'LineWidth', 1.5)
    line('YData', [xt(2) xt(2)], 'XData', xlim, 'Color','g', 'LineWidth', 1.5)
    ylim(yl)
    
    
    linkaxes([ha(1), ha(2)], 'y')
    linkaxes([ha(1), ha(3)], 'x')
end
