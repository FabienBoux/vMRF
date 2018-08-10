
addpath(genpath(fullfile(pwd(), 'functions')))
addpath(genpath(fullfile(pwd(), 'tools')))


%% Learning step

% complete 3 sequences
file        = 'dictionaries/multidico/2018-07-25-11:22_complex_dico.mat';
% only GESFIDSE sequence
% file        = 'dictionaries/multidico/2018-07-25-09:35_32-samples_5-parameters_80000-signals.mat';

coord       = [1 2];
snr         = [20 200];

K           = 20;
Lw          = 0;
cstr.Sigma  = 'd';
cstr.Gammat = 'd';
cstr.Gammaw = '';


load(file)
X           = AddNoise(X, snr(1)+(snr(2)-snr(1))*rand(size(X,1),1), 0);
[theta, r]	= EstimateInverseFunction(Y(:,coord), X, K, Lw, 100, cstr, 0);



%%
xbvf        = 0:0.001:0.3;
xvsi        = (0:0.04:12) *1e-6;
% xt2         = 0.02:0.001:0.32;
x_samples	= [xbvf; xvsi];
% x_samples	= [xbvf; xvsi; xt2];

ind         = randperm(size(X,1), 40000);

temp = zeros(size(coord));
parfor i = 1:length(ind)    
    t = theta;
    
    [~,~,t]  = gllim_inverse_dens_modified(X(ind(i),:)',t,ones(size(X(ind(i),:)')),x_samples,0);
    [xt_predict{i},alpha{i}]  = gllim_inverse_map(X(ind(i),:)',t,0);
    
    mae(i,:)        = Y(ind(i),coord) - xt_predict{i}';
    
    tmp = temp;
    for k = 1:K
        tmp = tmp + alpha{i}(k) * diag(t.Sigma_s(:,:,k))';
    end
    confidence(i,:) = tmp;
end


%% First investigation

ind     = 2;
bands   = 100;
xbins   = (0:.001:.5) *1e-5;

[counts,centers] = hist(confidence(:,ind),bands);
dw          = unique(diff(centers));
centers     = centers + dw(1)/2;

% figure
% hold on
clear m mm
for i = 2:length(centers)-1
    mm(i) = i;
    if counts(i) <= 20
        m(i) = nan;
        s(i) = nan;
        fprintf([num2str(i) ' skipped \n'])
    else
        d = mae((confidence(:,ind) <= centers(i)) & (confidence(:,ind) > centers(i-1)),ind);
        m(i) = mean(abs(d));
        s(i) = std(d);
    end
end
m(1) = nan;
s(1) = nan;
% legend(num2str((1:bands)'))

figure
subplot(121)
plot(mm(~isnan(m)), m(~isnan(m)), 'o')
mdl = fitlm(mm(~isnan(m)), m(~isnan(m)));
p = mdl.Coefficients.Estimate;
hold on;
plot(mm(~isnan(m)), mm(~isnan(m))*p(2) + p(1))
title(['Rsquared = ' num2str(mdl.Rsquared.Ordinary,2)])

subplot(122)
plot(mm(~isnan(s)), s(~isnan(s)), 'o')
mdl = fitlm(mm(~isnan(s)), s(~isnan(s)));
p = mdl.Coefficients.Estimate;
hold on;
plot(mm(~isnan(s)), mm(~isnan(s))*p(2) + p(1))
title(['Rsquared = ' num2str(mdl.Rsquared.Ordinary,2)])



%% Second investigation

ind     = 1;
bands   = 10;
xbins   = (0:.001:.5) *1e-5;

[counts,centers] = hist(confidence(:,ind),bands);
dw          = unique(diff(centers));
centers     = centers + dw(1)/2;

clear class
c = 1;
for i = 2:length(centers)-1
    mm(i) = i;
    if counts(i) <= 50
        m(i) = nan;
        s(i) = nan;
        fprintf([num2str(i) ' skipped \n'])
    else
        d = mae((confidence(:,ind) <= centers(i)) & (confidence(:,ind) > centers(i-1)),ind);
        m(i) = mean(d);
        s(i) = std(d);
        class_conf{c} = confidence((confidence(:,ind) <= centers(i)) & (confidence(:,ind) > centers(i-1)),ind);
        class{c} = d;
        c = c +1;
    end
end
m(1) = nan;
s(1) = nan;


figure
for i = 1:length(class)
    subplot(4,4,2*i-1)
    hist(abs(class{i}),0:.001:.1)
    xlim([0 .1])
    
    subplot(4,4,2*i)
    hist(class_conf{i},500)
end







