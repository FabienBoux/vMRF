function [acc] = ComputeAccuracy(v, D, e, disp)

if nargin > 4
    error('ComputeAccuracy: too much input arguments')
end
if ~exist('disp','var')
    disp = 0;
end

N       = length(size(D));
s       = size(D);

% if N ~= length(tol)
%     error('ComputeAccuracy: Input sizes incorrect')
% end

d       = reshape(D, 1,[]);
[ds, ind] = sort(d, 'descend');
p       = cumsum(ds);
n       = sum(~(p >= e));

dd      = zeros(size(d));
dd(ind(1:n)) = 1;

if disp == 1
    imagesc(v{1}, v{2}, reshape(dd, s))
    set(gca,'YDir','normal')
end

acc     = zeros(1,N);
for i = 1:N
    loc         = find(max(reshape(dd,s),[],i) == 1);
    acc(i)      = v{i}(max(loc)) - v{i}(min(loc));
end