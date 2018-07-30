
clear

files.path      = '/home/bouxf/Bureau/fingerprint/IRM/Image_Analyses_data';
files.animal    = 'C6-K89';

files.pre       = 'MGEFIDSEpre-60ms-S.mat';
files.post      = 'MGEFIDSEpost-60ms-S.mat';
files.mask      = 'VOI-Atlas_Brain.mat';

files.bvf_ref   = 'BVf-S.mat';
files.bvf_fingerpt = 'BVf-Num-dense_se.mat';
files.vsi_ref   = 'VSI-S.mat';
files.vsi_fingerpt = 'R-Num-dense_se.mat';

list_files      = {'pre', 'post', 'mask', 'bvf_ref', 'vsi_ref', 'bvf_fingerpt', 'vsi_fingerpt'};

folder_files    = dir(files.path);
for i = 1:length(list_files)
    for j = 1:length(folder_files)
        if (contains(folder_files(j).name, getfield(files, list_files{i})) && contains(folder_files(j).name, files.animal))
            files = setfield(files, list_files{i}, [files.path '/' folder_files(j).name]);
        end
    end
end
clear i j

pre  = load(files.pre);
post = load(files.post);
mask = load(files.mask);

temp = zeros(size(mask.uvascroi(1).value));
for i = 1:5 %length(mask.uvascroi)
    temp = temp + mask.uvascroi(i).value;
end
mask = (temp ~= 0);
clear temp

dico = post.uvascim.image.reco.data ./ pre.uvascim.image.reco.data;
clear pre post

for i = 1:size(dico,3)
    for j = 1:size(dico,4)
        dico(:, :, i, j) = dico(:, :, i, j) .* mask;
    end
end
clear i j mask

for i = 1:size(dico,3)
    temp(i,:) = reshape(squeeze(dico(:,:,i,:)), 1,[]);
end
retained = (sum(temp) ~= 0);
X = temp(:,retained)';
clear dico temp i

% Ref values
load(files.bvf_ref)
Y(:,1)  = reshape(squeeze(uvascim.image.reco.data), 1,[])';
load(files.vsi_ref)
Y(:,2)  = reshape(squeeze(uvascim.image.reco.data), 1,[])';
Y       = Y(retained,:);

% Fingerprint values
load(files.bvf_fingerpt)
Yf(:,1) = reshape(squeeze(uvascim.image.reco.data), 1,[])';
load(files.vsi_fingerpt)
Yf(:,2) = reshape(squeeze(uvascim.image.reco.data), 1,[])';
Yf      = Yf(retained,:);
clear uvascim retained


% Save
save(['dico_' files.animal '_S.mat'], 'X','Y','Yf')
