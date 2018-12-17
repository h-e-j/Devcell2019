%Heath Johnson 12/17/18
%Calculates the punctate myosin intensity as a function of angular 
%posistion around the embryo

%%Set up directories and varibles
Filepath='D:\Program Files\MATLAB\R2013a\work\'; %Directory where the m-files are located
addpath('D:\Program Files\MATLAB\R2013a\work\')
savnam='ratioredsubhkred';
basenames={'myo'};
columnW=[177];
FolderNames={'myo squeeze confocal'};
bstep=31;
angles=361;
fileFolder = ['D:\Program Files\MATLAB\R2013a\work\' sprintf('%s',FolderNames{1}) '\1\'];

cd(fileFolder)
load psudeowhite %load colormap
%%
for uu=1:length(columnW)
    HSLog=zeros(columnW(uu),angles);
    basename=basenames{uu};
    for w=1:columnW(uu)
        
        if  w<10
            LSn=[basename '000'  num2str(w)];
            
        elseif w<100
            LSn=[basename '00'  num2str(w)];
        else
            LSn=[basename '0'  num2str(w)];
        end
        LS=imread([LSn '.tif']);
        noback=LS>700;
        [a b]=size(LS);
        noback=imclose(noback,strel('disk',2));
        noback=bwareaopen(noback,1000);
        noback=imfill(noback,'holes');
        rbroE= imopen(noback,strel('Disk',100));
        rbroB=noback;
        lrbro=length(noback);
        rbroB(:,9*lrbro/10:lrbro)=0;
        rbroE(1:lrbro/10)=0;
        bestbro=rbroB+rbroE;
        noback=bestbro;
        LSb=LS.*uint16(ones(size(LS))-noback);
        figure(1)
        imshow(noback)
        truesize
        LSbb=double(mean(nonzeros(LSb)));
        myo=uint16(noback).*LS;
        myo=double(myo);
        myohole=imerode(noback,strel('disk',50));
        countmask=noback-myohole;
        countmask=countmask>0;
        clear STATS AllCenter;
        labelcell=bwlabel(noback);
        STATS=regionprops(labelcell,'Area','Centroid');
        OGSize=cell2mat({STATS.Area});
        [OGSortSize index]=sort(OGSize,'descend');
        AllCenter=cell2mat({STATS.Centroid});
        Center(1)=AllCenter(2*index(1)-1);
        Center(2)=AllCenter(2*index(1));
        dcc=double(countmask);
        [SpotYb,SpotXb]=find(countmask);%find indicies of all peripheral pixels
        HSalphasav=zeros(1,length(SpotXb));
        %% Find orientation of peripheral pixels
        HSalphasav=zeros(1,length(SpotXb));
        for i=1:length(SpotXb) %for each pixel
            HSalpha=atan2((SpotYb(i))-Center(2),(SpotXb(i)-Center(1)))*180/pi+181;
            HSalphasav(i)=fix(HSalpha);
            HSLog(w,fix(HSalpha))=HSLog(w,fix(HSalpha))+myo(SpotYb(i),SpotXb(i));
            if fix(HSalpha)==30
                dcc(SpotYb(i),SpotXb(i))= 2;
            end
            
        end
        normcount=histc(HSalphasav,1:361);
        for gg=1:length(normcount)
            if normcount(gg)~=0
                HSLog(w,gg)=HSLog(w,gg)/normcount(gg);
            end
        end
        
        
        if mod(w,20)==0
            figure(2)
            imagesc(HSLog(80:columnW(1),:))
            pause(0.01)
        end
    end
end
%% Plot results
figure(2)
binsiz=60;
sublog=HSLog(1:columnW(1),:);
sublogm=sublog;
imagesc(sublogm)

colormap(A)

caxis([860 1300])
set(gca,'TickDir','out')
set(gcf, 'Color', 'w')
set(gca,'Fontsize',16)

set(gca,'YTick',[1  31 61 91 121]);
set(gca,'YTickLabel',{'0','30','60','90','120'});
set(gca,'XTick',[1 91 181 271 361]);
set(gca,'XTickLabel',{'anterior','dorsal','posterior','ventral','anterior'});






