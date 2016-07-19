function [] = mozek(mri,E_mri_size,el)
% nakresli obrazky mozku pacienta s elektrodami
% starsi skript nahrazen funkci mozek2

figure('Name','Mozek');
% if nargin >= 3
%     n = em{el,1};
%     x = em{el,3};
%     y = em{el,2};
%     z = em{el,4};
% else
%     x = round(size(mri,1)/2);
%     y  = round(size(mri,2)/2);
%     z  = round(size(mri,3)/2);
% end

x=55; y=256;z=256;
mri_vox_dim = [0.500000000000000,0.419999986886978,0.419999986886978];
%[x,y,z] = ind2sub(size(E_mri_size),find(E_mri_size==el));

%y x z
%coronal section = zpredu dozadu
subplot(1,3,1);
section(squeeze(mri(:,:,round(z))),mri_vox_dim([1 2])); %256
hold on;
plot(y,x,'or');

%sagital section = zleva doprava
subplot(1,3,2);
section(squeeze(mri(:,round(y),:)),mri_vox_dim([1 3])); %250
hold on;
plot(z,x,'or');

%horizontal section = shoral dolu
subplot(1,3,3);
section(squeeze(mri(round(x),:,:)),mri_vox_dim([2 3])); %55 Y;   x y z 
hold on;
plot(z,y,'or');


% if nargin >= 2
%     figure('Name','Elektrody');
%     subplot(1,3,1);
%     hold on;
%     for j = 1:size(em,1)
%         plot(em{j,3},em{j,2},'o');
%     end
%     subplot(1,3,2);
%     hold on;
%     for j = 1:size(em,1)
%         plot(em{j,3},em{j,4},'o');
%     end
%     subplot(1,3,3);
%     hold on;
%     for j = 1:size(em,1)
%         plot(em{j,2},em{j,4},'o');
%     end
% end




end
function []=section(mri0,vox_dim)
%jeden obrazek jednoho rezu mozkem. m je matice 2D
    mri0 = 500-mri0; %aby bila hmota byla bila
    imagesc([0 vox_dim(2)*size(mri0,2)],[0 vox_dim(1)*size(mri0,1)],mri0);
    colormap(gray);
    axis square;
end