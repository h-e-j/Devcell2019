%Heath Johnson 12-17-18
%Rotate embryos and calculates the fluorescent intensity of FISH expression in a contour
%around the embryo for one channel and the number of RNA puncta around the
%contour in a second channel of the same image


Filepath='D:\Program Files\MATLAB\R2013a\work\'; %Directory where the m-files are located
addpath('D:\Program Files\MATLAB\R2013a\work\')
savnam='ratioredsubhkredc';
cluster=3; %number of k-means clusters
basenames={'dark' '2hr','40min','15min','30min'};
columnW=[76 79 75 81 71]; %number of images in each set
medsiz=4; %median filtering size
FolderNames={'2017-5-19 tll FISH duration'};
loadangles=1; % Load embryo rotation angles
numpoints=100; %number of bins
bstep=11; %pixel depth of contours
numpointshk=300; %bins for second channel
fileFolder = ['D:\Program Files\MATLAB\R2013a\work\' sprintf('%s',FolderNames{1})];
%%
cd(fileFolder)
load psudeowhite
for uu=1:length(columnW)
    
    basename=basenames{uu};
    totalcells=zeros(columnW(uu),1);
    totalred= totalcells;
    redtotalcelloverlap= totalcells;
    totalredthick= totalcells;
    plog=zeros(columnW(uu),360);
    plogs=zeros(columnW(uu),360);
    if loadangles==1
        load(['rotatebrosflip' sprintf('%s',basename)])
    else
        angle=zeros(columnW(uu),1);
        flipr=angle;
    end
    
    Nhk=zeros(1,numpointshk);
    Ntotmat=zeros(columnW(uu),numpoints);
    Ntothkmat=zeros(columnW(uu),numpointshk);
    
    Ntot=zeros(1,numpoints);
    Ntothk=zeros(1,numpointshk);
    for w=1:columnW(uu)
        w
        Ntots=zeros(1,numpoints);
        Ntothks=zeros(1,numpointshk);
        
        if  w<10
            dapin=[basename '_t00'  num2str(w) '_c001'];
            hkbn=[basename '_t00'  num2str(w) '_c003'];
            
            redn=[basename '_t00'  num2str(w) '_c002'];
            
            mistn=[basename '_t00'  num2str(w) '_c004'];
        else
            dapin=[basename '_t0'  num2str(w) '_c001'];
            hkbn=[basename '_t0'  num2str(w) '_c003'];
            
            redn=[basename '_t0'  num2str(w) '_c002'];
            mistn=[basename '_t0'  num2str(w) '_c004'];
        end
        
        dapi=imread([dapin '.tif']);
        mist=imread([mistn '.tif']);
        hkb=imread([hkbn '.tif']);
        red=imread([redn '.tif']);
        B = medfilt2(mist, [medsiz medsiz]);
        outlie=mist-B;
        
        
        [v1 , c1] = kmeanssub(dapi, cluster+3,floor(numel(dapi)/5));
        noback=c1>1;
        noback=imclose(noback,strel('disk',2));
        noback=bwareaopen(noback,20);
        noback=imfill(noback,'holes');
        
        hkbb=hkb.*uint16(ones(size(hkb))-noback);
        redb=red.*uint16(ones(size(hkb))-noback);
        hkbb=double(mean(nonzeros(hkbb)));
        redb=double(mean(nonzeros(redb)));
        
        
        hkb=double(hkb)-hkbb;
        wholebro=noback;
        lbro=bwlabel(wholebro);
        STATS=regionprops(wholebro,'Area','Centroid');
        OGSize=cell2mat({STATS.Area});
        [OGSortSize index]=sort(OGSize,'descend');
        wholebro=zeros(size(wholebro));
        wholebro(find(lbro==index(1)))=1;
        hkb=(hkb-hkbb);
        red=double(red)-redb;
        hkb=hkb/mean(nonzeros(double(wholebro).*red));
        wholebroring=wholebro-imerode(wholebro,strel('disk',11));
        wholebro=wholebroring;
        kmeanssubbreak(mistG, cluster+2,floor(numel(dapi)/10));
        [v3 , c3, breakbin3] = kmeanssubbreak(outlie, cluster+2,floor(numel(dapi)/10));
        redpunct=c3>3;
        redpunct=redpunct-bwareaopen(redpunct,4);
        
        if loadangles~=1
            
            figure(22)
            imagesc(hkb)
            figure(222)
            imagesc(mist)
            figure(2)
            imagesc(dapi)
            colormap('gray')
            truesize([1024 1024])
            h = imline;
            pos = h.getPosition;
            px=pos(:,1);
            py=pos(:,2);
            angle(w)=atan2((py(2)-py(1)),(px(2)-px(1)))*180/pi;
            
            
        end
        rbro=imrotate(wholebro,angle(w));
        
        hkb=imrotate(hkb,angle(w));
        if loadangles~=1
            figure(2)
            imshow(rbro)
            truesize
            h = imline;
            pos = h.getPosition;
            px=pos(:,1);
            py=pos(:,2);
            if py>256
                flipr(w)=1;
            else
                flipr(w)=0;
                
            end
        end
        if flipr(w)==1
            
            rbro=flip(rbro);
            hkb=flip(hkb);
        end
        rbro=flip(rbro,2);
        
        rbroE= imopen(imfill(rbro,'holes'),strel('Disk',50));
        rbroB=rbro;
        lrbro=length(rbro);
        rbroB(:,lrbro/2:lrbro)=0;
        rbroE(1:lrbro/2)=0;
        bestbro=rbroB+rbroE;
        bestbro=imfill(bestbro,'holes');
        rbro=bestbro-imerode(bestbro,strel('disk',11));
        
        hkb=flip(hkb,2);
        STATS=regionprops(rbro,'Area','Centroid');
        
        
        rmist=imrotate(redpunct,angle(w));
        
        if flipr(w)==1
            rmist=flip(rmist);
        end
        rmist=flip(rmist,2);
        rmist= rmist&rbro;
        [yy xx]=find(rmist);
        stepbro=imfill(rbro,'holes');
        for bs=1:bstep
            B = bwboundaries(stepbro,'noholes');
            
            BB=(B{1});
            
            Xout=BB(:,2);
            Yout=BB(:,1);
            hkints=zeros(length(length(Yout)),1);
            for hh=1:length(Yout)
                stepbro(BB(hh,1),BB(hh,2))=0;
                hkints(hh)=hkb(BB(hh,1),BB(hh,2));
            end
            
            ploglist=zeros(length(BB(:,1)),1);
            
            for j=1:length(xx)
                iix=find(Xout==xx(j));
                iiy=find(Yout==yy(j));
                ol=intersect(iix,iiy);
                if ~isempty(ol)
                    ii=ol;
                else
                    ii=0;
                end
                ploglist(j)=ii(1);
                
            end
            
            [N,edges] = histcounts(ploglist,linspace(1,length(Xout),numpoints+1));
            
            hkbin=linspace(1,length(Xout),numpointshk+1);
            for zz=1:length(hkbin)-1
                Nhk(zz)=sum(hkints(floor(hkbin(zz):hkbin(zz+1))));
            end
            
            Ntot=Ntot+N;
            Ntothk=Ntothk+Nhk;
            
            Ntots=Ntots+N;
            Ntothks=Ntothks+Nhk;
            
            if bs==bstep
                figure(3)
                dapim=imrotate(dapi,angle(w));
                if flipr(w)==1
                    dapim=flip(dapim);
                end
                dapim=flip(dapim,2);
                imagesc(dapim)
                colormap('gray')
                [yr xr]=find(imfill(rbro,'holes')-stepbro);
                
                hold on
                plot(xr,yr,'.b')
                
                plot(xx,yy,'.r')
                hold off
                truesize([1024 1024])
                figure(66)
                bar(Ntothk)
                pause(0.001)
                
            end
            
        end
        Ntotmat(w,:)=Ntots;
        Ntothkmat(w,:)=Ntothks;
    end
    if loadangles ~=1
        save(['rotatebrosflip' basename],'angle','flipr')
    end
    figure
    bar(Ntot/columnW(uu))
    xlabel('P-A')
    ylabel('Counts per embryo')
    
    if uu==1
        dark=Ntot/columnW(uu);
        darkhk=Ntothk/columnW(uu);
        darkmat=Ntotmat;
        darkhkmat=Ntothkmat;
        save(['dark' sprintf('%s',savnam)],'dark','darkhk','darkmat','darkhkmat')
    elseif uu==2
        hr2opto=Ntot/columnW(uu);
        hr2optohk=Ntothk/columnW(uu);
        hr2optomat=Ntotmat;
        hr2optohkmat=Ntothkmat;
        save(['2hr' sprintf('%s',savnam)],'hr2opto','hr2optohk','hr2optomat','hr2optohkmat')
    elseif uu==3
        min40opto=Ntot/columnW(uu);
        min40optohk=Ntothk/columnW(uu);
        min40optomat=Ntotmat;
        min40optohkmat=Ntothkmat;
        save(['min40' sprintf('%s',savnam)],'min40opto','min40optohk','min40optomat','min40optohkmat' )
    elseif uu==4
        min15opto=Ntot/columnW(uu);
        min15optohk=Ntothk/columnW(uu);
        min15optomat=Ntotmat;
        min15optohkmat=Ntothkmat;
        save(['min15' sprintf('%s',savnam)],'min15opto','min15optohk','min15optomat','min15optohkmat')
    elseif uu==5
        min30opto=Ntot/columnW(uu);
        min30optohk=Ntothk/columnW(uu);
        min30optomat=Ntotmat;
        min30optohkmat=Ntothkmat;
        save(['min30' sprintf('%s',savnam)],'min30opto','min30optohk','min30optomat','min30optohkmat')
    end
    
end

%%Plot results
figure(56)
plot(1:300,hr2optohk/mean(darkhk),'LineWidth',3,'Color','b')
hold on
plot(1:300,min40optohk/mean(darkhk),'LineWidth',3,'Color','g')
plot(1:300,min30optohk/mean(darkhk),'LineWidth',3,'Color','m')
plot(1:300,min15optohk/mean(darkhk),'LineWidth',3,'Color','r')
plot(1:300,darkhk/mean(darkhk),'LineWidth',3,'Color','k')
hold off
ylabel('normalized tll expression')
xlim([1 300])
ylim([0.5 3.5])
set(gcf,'color','white')
xticks([1 0.25*300 150 0.75*300 300])
xticklabels({'anterior','dorsal','posterior','ventral','anterior'})
set(gca,'Fontsize',18)
box('off')
