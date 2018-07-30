

grid_bvf        = 0.01:0.01:0.25;
grid_vsi        = 1e-6:5*1e-7:1e-5;
[Xgrid, Ygrid]  = meshgrid(grid_bvf, grid_vsi);

rand_bvf        = 0.01 + 0.24   * rand(1, size(Xgrid,1) * size(Xgrid,2));
rand_vsi        = 1e-6 + 9*1e-6 * rand(1, size(Xgrid,1) * size(Xgrid,2));

figure
subplot(121); plot(reshape(Xgrid,1,[]), reshape(Ygrid,1,[]), '.', 'LineWidth', 2)
xlabel('BVf'); ylabel('VSI (in m)')
title('Dictionary')
subplot(122); plot(rand_bvf, rand_vsi, '.', 'LineWidth', 2)
xlabel('BVf'); ylabel('VSI (in m)')
title('Learning dataset')