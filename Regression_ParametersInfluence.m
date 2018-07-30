% 22/01/18 
% Le but de ce script est d'imager l'influence des paramètres à partir de
% l'apprentissage. Pour ce faire on fera un apprentissage sur BVf et VSI 
%
% 14/05/18
% le file_dico est modifié pour essayer un dictionnaire composé de signaux
% MGEFIDSE pre et post injection d'AC
% sequence_file


% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file_dico = 'files/dico_random_220000signals.mat';

% sequence_file   = 'files/config/seqpar_GESFIDE.txt';
% params_file     = 'files/voxpar_dico_2params.txt';

% 14/05/18 update
file_dico = 'files/06-Apr-2018_64-samples_4-parameters_185367-signals.mat';
file_grid = 'files/04-Apr-2018_64-samples_4-parameters_185367-signals.mat';

sequence_file   = 'files/config/seqpar_GESFIDE.txt';
params_file     = 'files/myconf/voxpar_dico.txt';

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load(file_dico)

coord   = 1:3;
Y       = Y(:,coord);
X       = abs(X);


%% Model learning
K           = 20;
lw          = 0;
cstr.Sigma  = '';

[theta, r]  = EstimateInverseFunction(Y, X, K, lw, 50, cstr, 0);


%% Parameter influence

grid = load(file_grid);

grid.Y  = grid.Y(:,coord);
grid.X  = abs(grid.X);

% BVf
vect = find(all([grid.Y(:,2) == 5e-6, grid.Y(:,3) == 1e-9]'));

Y_bvf     	= grid.Y(vect,1)';
Y_vsi       = 5e-6 *ones(size(Y_bvf));
Y_adc       = 1e-9 *ones(size(Y_vsi));

[X_bvf, ~] = gllim_forward_map([Y_bvf; Y_vsi; Y_adc], theta, 1);


figure
subplot(222); mesh(repmat(1:size(X_bvf,1), size(X_bvf,2), 1)', repmat(Y_bvf, size(X_bvf,1), 1), abs(X_bvf))
xlabel('Samples'); ylabel('BVf'); zlabel('Signal intensity'); title('Map from learnt model')

subplot(221); mesh(repmat(1:size(grid.X,2), length(vect), 1)', repmat(grid.Y(vect,1)', size(grid.X,2), 1), abs(grid.X(vect,:))')
xlabel('Samples'); ylabel('BVf'); zlabel('Signal intensity'); title('Map from simulation tool')

subplot(2,2,3); mesh(repmat(1:size(X_bvf,1), size(X_bvf,2), 1)', repmat(Y_bvf, size(X_bvf,1), 1), ((grid.X(vect,:)' - abs(X_bvf)).^2))
xlabel('Samples'); ylabel('BVf'); zlabel('Squared error'); title('Map error')
subplot(2,2,4); mesh(repmat(1:size(X_bvf,1), size(X_bvf,2), 1)', repmat(Y_bvf, size(X_bvf,1), 1), abs((grid.X(vect,:)' - abs(X_bvf))))
xlabel('Samples'); ylabel('BVf'); zlabel('Absolute error'); title('Map error')

colormap(jet)


% VSI
vect = find(all([grid.Y(:,1) == 0.1, grid.Y(:,3) == 1e-9]'));

Y_vsi       = grid.Y(vect,2)';
Y_bvf       = 0.1 *ones(size(Y_vsi));
Y_adc       = 1e-9 *ones(size(Y_vsi));

[X_vsi, alpha] = gllim_forward_map([Y_bvf; Y_vsi; Y_adc], theta, 0);


figure
subplot(222); mesh(repmat(1:size(X_vsi,1), size(X_vsi,2), 1)', repmat(Y_vsi, size(X_vsi,1), 1), abs(X_vsi))
xlabel('Samples'); ylabel('VSI'); zlabel('Signal intensity'); title('Map from learnt model')

subplot(221); mesh(repmat(1:size(grid.X,2), length(vect), 1)', repmat(grid.Y(vect,2)', size(grid.X,2), 1), abs(grid.X(vect,:))')
xlabel('Samples'); ylabel('VSI'); zlabel('Signal intensity'); title('Map from simulation tool')

subplot(2,2,3); mesh(repmat(1:size(X_vsi,1), size(X_vsi,2), 1)', repmat(Y_vsi, size(X_vsi,1), 1), ((grid.X(vect,:)' - abs(X_vsi)).^2))
xlabel('Samples'); ylabel('BVf'); zlabel('Squared error'); title('Map error')
subplot(2,2,4); mesh(repmat(1:size(X_vsi,1), size(X_vsi,2), 1)', repmat(Y_vsi, size(X_vsi,1), 1), abs((grid.X(vect,:)' - abs(X_vsi))))
xlabel('Samples'); ylabel('BVf'); zlabel('Absolute error'); title('Map error')

colormap(jet)

% ADC
vect = find(all([grid.Y(:,1) == 0.1, grid.Y(:,2) == 5e-6]'));

Y_adc       = grid.Y(vect,3)';
Y_bvf       = 0.1 *ones(size(Y_adc));
Y_vsi       = 5e-6 *ones(size(Y_bvf));

[X_adc, alpha] = gllim_forward_map([Y_bvf; Y_vsi; Y_adc], theta, 0);

figure
subplot(222); mesh(repmat(1:size(X_adc,1), size(X_adc,2), 1)', repmat(Y_adc, size(X_adc,1), 1), abs(X_adc))
xlabel('Samples'); ylabel('ADC'); zlabel('Signal intensity'); title('Map from learnt model')

subplot(221); mesh(repmat(1:size(grid.X,2), length(vect), 1)', repmat(grid.Y(vect,3)', size(grid.X,2), 1), abs(grid.X(vect,:))')
xlabel('Samples'); ylabel('ADC'); zlabel('Signal intensity'); title('Map from simulation tool')

subplot(2,2,3); mesh(repmat(1:size(X_adc,1), size(X_adc,2), 1)', repmat(Y_adc, size(X_adc,1), 1), ((grid.X(vect,:)' - abs(X_adc)).^2))
xlabel('Samples'); ylabel('BVf'); zlabel('Squared error'); title('Map error')
subplot(2,2,4); mesh(repmat(1:size(X_adc,1), size(X_adc,2), 1)', repmat(Y_adc, size(X_adc,1), 1), abs((grid.X(vect,:)' - abs(X_adc))))
xlabel('Samples'); ylabel('BVf'); zlabel('Absolute error'); title('Map error')

colormap(jet)


%% Parameter influence (variante)

grid = load(file_grid);

grid.X = abs(grid.X);

% BVf
vect = find(all([grid.Y(:,1) == 0.1, grid.Y(:,2) == 5e-6]'));

Y_adc     	= grid.Y(vect,3)';
Y_bvf       = 0.1 *ones(size(Y_adc));

[X_adc, ~] = gllim_forward_map([Y_bvf; Y_adc], theta, 1);


figure
subplot(222); mesh(repmat(1:size(X_adc,1), size(X_adc,2), 1)', repmat(Y_adc, size(X_adc,1), 1), abs(X_adc))
xlabel('Samples'); ylabel('ADC'); zlabel('Signal intensity'); title('Map from learnt model')

subplot(221); mesh(repmat(1:size(grid.X,2), length(vect), 1)', repmat(grid.Y(vect,3)', size(grid.X,2), 1), abs(grid.X(vect,:))')
xlabel('Samples'); ylabel('ADC'); zlabel('Signal intensity'); title('Map from simulation tool')

subplot(2,2, 3:4); mesh(repmat(1:size(X_adc,1), size(X_adc,2), 1)', repmat(Y_adc, size(X_adc,1), 1), 100 * abs(grid.X(vect,:)' - abs(X_adc)) ./ grid.X(vect,:)')
xlabel('Samples'); ylabel('ADC'); zlabel('Relative error (%)'); title('Map error')
zlim([0 120])
colormap(jet)


