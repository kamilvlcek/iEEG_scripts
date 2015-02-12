function [ ] = mozek2( mri,implanted_electrode,mri_vox_dim,els )
%MOZEK2 Summary of this function goes here
%   els jsou jmena elektrod, napriklad els = {'C2','E5','E3','P3'};


%figure('Name','elektrody');
% sp = 1;
% for xx = 2:4;
%     for yy = 2:4;
% %         subplot(3,3,sp);
%         sp = sp + 1;
%         for ii = 1:size(implanted_electrode,1)
% %           plot(implanted_electrode{ii,xx},implanted_electrode{ii,yy},'or');
% %           hold on;
%         end
%     end
% end




%spocitam jmena a pozice elektrod
if ~exist('els','var') , els = {0}; end %7.11.2014
electrodes = cell(1,numel(els)); %tohle bude struktura
for ee = 1:numel(els) %electrodes je cell array
    EE = els{ee};
    if ischar(EE)
        electrodes{ee}.el_name = EE;
        electrodes{ee}.el = findelbyname(implanted_electrode,EE);
    elseif EE > 0
         electrodes{ee}.el_name = implanted_electrode{EE,1};
         electrodes{ee}.el = EE;
    else
        electrodes{ee}.el_name = 'zadna elektroda';
        electrodes{ee}.el = 0;
    end
    
    if ee==1 %pozice v mozku podle prvni elektrody
        el = electrodes{ee}.el;
        el_name = electrodes{ee}.el_name;
        if el == 0
            x=size(mri,1)/2; y=size(mri,2)/2; z=size(mri,3)/2;
        else
            x = implanted_electrode{el,3}/mri_vox_dim(2);
            y = implanted_electrode{el,2}/mri_vox_dim(1);
            z = implanted_electrode{el,4}/mri_vox_dim(3);
             
        end
    end
end

figure('Name',['Mozek2 ' el_name '=' num2str(el)]);
%coronal section = zpredu dozadu
subplot(1,3,1);
section(squeeze(mri(:,:,round(z))),mri_vox_dim([1 2])); %256
hold on;

for ii = 1:size(implanted_electrode,1)
    %plot(implanted_electrode{ii,2},implanted_electrode{ii,3},'ob');
    hold on;
end
for ee = 1:size(implanted_electrode,1)
    xy =  [implanted_electrode{ ee,2} implanted_electrode{ee,3}];
    %plot(xy(1),xy(2), 'o','MarkerSize',1.5);
    circle(1,xy(1),xy(2),'b',0); 
end
    
for ee = numel(electrodes):-1:1
    if electrodes{ee}.el>0
        xy =  [implanted_electrode{ electrodes{ee}.el,2} implanted_electrode{ electrodes{ee}.el,3}];
        circle(3,xy(1),xy(2),iff(ee==1,'r','c'),1); 
        text(xy(1)+10,xy(2)+10, electrodes{ee}.el_name);
    end
end

 
%sagital section = zleva doprava
subplot(1,3,2);
section(squeeze(mri(:,round(y),:)),mri_vox_dim([1 3])); %250
hold on;

for ii = 1:size(implanted_electrode,1)
    %plot(implanted_electrode{ii,4},implanted_electrode{ii,3},'ob');
    hold on;
end
for ee = 1:size(implanted_electrode,1)
    xy =  [implanted_electrode{ ee,4} implanted_electrode{ee,3}];
    circle(1,xy(1),xy(2),'b',0); 
end
for ee = numel(electrodes):-1:1
    if electrodes{ee}.el>0
        xy =  [implanted_electrode{ electrodes{ee}.el,4} implanted_electrode{ electrodes{ee}.el,3}];
        circle(3,xy(1),xy(2),iff(ee==1,'r','c'),1); 
        text(xy(1)+10,xy(2)+10, electrodes{ee}.el_name);
    end
end

%horizontal section = shoral dolu
subplot(1,3,3);
section(squeeze(mri(round(x),:,:)),mri_vox_dim([2 3])); %55 Y;   x y z 
hold on;
for ii = 1:size(implanted_electrode,1)
    %plot(implanted_electrode{ii,4},implanted_electrode{ii,2},'ob');
    %circle(3,implanted_electrode{ii,4},implanted_electrode{ii,2},iff(ii==el,'r','b'),iff(ii==el,1,0));
    hold on;
end
for ee = 1:size(implanted_electrode,1)
    xy =  [implanted_electrode{ ee,4} implanted_electrode{ee,2}];
    circle(1,xy(1),xy(2),'b',0); 
end
for ee = numel(electrodes):-1:1
    if electrodes{ee}.el>0
        xy =  [implanted_electrode{ electrodes{ee}.el,4} implanted_electrode{ electrodes{ee}.el,2}];
        circle(3,xy(1),xy(2),iff(ee==1,'r','c'),1); 
        text(xy(1)+10,xy(2)+10, electrodes{ee}.el_name);
    end
end

tightfig;
%http://blogs.mathworks.com/pick/2012/12/21/figure-margins-subplot-spacings-and-more/#2
end

function []=section(mri0,vox_dim)
%jeden obrazek jednoho rezu mozkem. m je matice 2D
    mri0 = 500-mri0; %aby bila hmota byla bila
    imagesc([0 size(mri0,2)*vox_dim(2)],[0 size(mri0,1)*vox_dim(1)],mri0);
   colormap(gray);
   axis equal;
   %axis xy;
end

function [el]=findelbyname(implanted_electrode,elname)
    el = 0;
    for ee = 1:size(implanted_electrode,1)
        if strcmp(implanted_electrode{ee,1},elname)
            el = ee;
            break;
        end
    end
end

