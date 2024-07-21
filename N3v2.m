 clear all

%% Input parameters
VideoName='2-1.mp4';
txtname='N2_0.1.txt';
figure1='N2_0.1_f1.avi';
figure2='N2_0.1_f2.avi';
figure3='N2_0.1_f3.avi';
Ncell=2; % no. of unit cells
Rmin=5;
Rmax=15;
area=[9,38,716-9,679-38]; % x-coord, y-coord, width, height
%% 
%% Program
v1 = VideoWriter(figure1);
open(v1);
v2 = VideoWriter(figure2);
open(v2);
v3 = VideoWriter(figure3);
open(v3);
txt = fopen(txtname,'w');
videoObj=VideoReader(VideoName);
videoFrames=get(videoObj,'NumberOfFrames');
fprintf(txt,'%8s  %8s  %8s  %8s \n','e','k','ka/e','k/e');

Nnode=(Ncell+1)^2; % no. of nodes
for i=1:floor(videoFrames/30):videoFrames
    RGB=read(videoObj,i);
    figure(1)
    
    RGB=imcrop(RGB,area);
    imshow(RGB);
    
    Grey0 = rgb2gray(RGB);
    %     Grey=imadjust(Grey,[0.1 0.6],[],0.55);
    %     imshow(Grey)
    %     [centers, radii] = imfindcircles(Grey,[8 20],'ObjectPolarity','dark','Method','TwoStage');
    %     n=length(radii)
    %     viscircles(centers, min(radii)*ones(size(radii)),'EdgeColor','y');
    
    %% Find at least N circles
    sensitivity=0.85;
    adjst=0.3;
    n=0;
    while n<Nnode
        %         sensitivity=sensitivity-0.02;
        %         if sensitivity<0.6
        %             break
        %         end
        
        adjst=adjst+0.1;
        Grey=imadjust(Grey0,[0.1 0.8],[],adjst);
        imshow(Grey)
        
        [centers, radii] = imfindcircles(Grey,[Rmin Rmax],'ObjectPolarity','dark','Method','TwoStage');
        n=length(radii);
        viscircles(centers, min(radii)*ones(size(radii)),'EdgeColor','y');
        pause(0.1)
    end
    
    % Delete extra circles
    while n~=Nnode
        distance=zeros(n,n);
        for j=1:n
            distance(j,:)=sort(vecnorm(transpose(centers(j,:)-centers)));
        end
        minDist=sum(distance(:,1:3),2);
        minDistSort=sort(minDist);
        minDist1=minDistSort(1);
        minDist2=minDistSort(2);
        if or(centers(minDist==minDist2,1)==max(centers(:,1)),centers(minDist==minDist2,2)==min(centers(:,2)))
           centers(minDist==minDist2,:)=[];
           radii(minDist==minDist2,:)=[];
        else
            centers(minDist==minDist1,:)=[];
            radii(minDist==minDist1,:)=[];
        end
%         minDistSort=sort(minDist);
%         stdDist=minDistSort(floor(length(minDist2)/2));
%         deviation=abs(minDist-stdDist);
%         centers(deviation==max(deviation),:)=[];
%         radii(deviation==max(deviation))=[];
        n=length(radii);
    end
    viscircles(centers, min(radii)*ones(size(radii)),'EdgeColor','r');
    frame = getframe;
    writeVideo(v1,frame)
    
    % Find the points at corner, and calibrate the square
    if i==1
        distance=zeros(Nnode,Nnode);
        corners=zeros(4,2);
        for j=1:n
            distance(j,:)=sort(vecnorm(transpose(centers(j,:)-centers)));
        end
        maxDist=sum(distance,2);
        for j=1:4
            corners(j,:)=centers(maxDist==max(maxDist),:);
            maxDist(maxDist==max(maxDist),:)=0;
        end
        coordSum=sum(corners,2);
        A=corners(coordSum==min(coordSum),:);
        C=corners(coordSum==max(coordSum),:);
        corners(coordSum==min(coordSum) | coordSum==max(coordSum),:)=[];
        B=corners(corners(:,1)==min(corners(:,1)),:);
        corners(corners(:,1)==min(corners(:,1)),:)=[];
        D=corners;
        x1=A(1); x2=B(1); x3=C(1); x4=D(1);
        y1=A(2); y2=B(2); y3=C(2); y4=D(2);
        Ly=91.675; % length of the sample, mm 
        Lx=88.785;
        X1=Lx; X2=0; X3=0; X4=Lx;
        Y1=0; Y2=0; Y3=Ly; Y4=Ly;
        T=[ x1 y1 1 0 0 0 -X1*x1 -X1*y1;
            0 0 0 x1 y1 1 -Y1*x1 -Y1*y1;
            x2 y2 1 0 0 0 -X2*x2 -X2*y2;
            0 0 0 x2 y2 1 -Y2*x2 -Y2*y2;
            x3 y3 1 0 0 0 -X3*x3 -X3*y3;
            0 0 0 x3 y3 1 -Y3*x3 -Y3*y3;
            x4 y4 1 0 0 0 -X4*x4 -X4*y4;
            0 0 0 x4 y4 1 -Y4*x4 -Y4*y4];
        X=[X1 Y1 X2 Y2 X3 Y3 X4 Y4]';
        fa=T\X;
        
        a=fa(1);b=fa(2);c=fa(3);
        d=fa(4);e=fa(5);f=fa(6);
        g=fa(7);h=fa(8);
        rot=[d e f;
            a b c;
            g h 1];
        
        %                 pix1=rot*[x1 y1 1]'/(g*x1+h*y1+1)  % verify
    end
    
    % Get the real coordinates, mm
    Centers=zeros(size(centers));
    for j=1:length(radii)
        xi=centers(j,1); yi=centers(j,2);
        pix=rot*[xi yi 1]'/(g*xi+h*yi+1);
        Centers(j,:)=pix(1:2);
    end
    
    figure(2)
    plot(Centers(:,1),Centers(:,2),'ro')
    set(gca,'XTick',0:90/Ncell:Lx)
    set(gca,'YTick',0:90/Ncell:Ly)
    pax = gca;
    pax.FontWeight="bold";
    pax.FontSize = 16;
    grid on
    axis([-15,105,-5,95])
    title('Deformed nodes')
    frame = getframe(figure(2));
    writeVideo(v2,frame)
    
    
    % Curve fitting
    Centers=sortrows(Centers,1);
    Centerline=Centers(floor(Ncell/2)*(Ncell+1)+1:floor(Ncell/2+1)*(Ncell+1),:);
    Centerline=sortrows(Centerline,2);
    Centerline2=Centerline;
    for j=1:length(Centerline)
        Centerline2(j,1)=Centerline(j,1)-((Centerline(end,1)-Centerline(1,1))/(Centerline(end,2)-Centerline(1,2))*(Centerline(j,2)-Centerline(1,2)));
    end
    Centerline=Centerline2;
    if i==1
        L0=max(Centerline(:,2))-min(Centerline(:,2));
        Centerline0=Centerline;
    end
    
    x=0:Lx/Ncell:Lx;
    y=Centerline(:,1)-Centerline0(:,1);
    P= polyfit(x, y', 2);
    k=2*P(1);
    
    xi=linspace(min(x),max(x),20);
    yi= polyval(P, xi);
    
    f3=figure(3);
    
    plot(x,y,'r*',xi,yi,'k','LineWidth',1.2);
    pax = gca;
    pax.FontWeight="bold";
    pax.FontSize = 16;
    axis([0,90,-15,15])
    view(90,-90)
    % axis equal
    legend('Experiment','Curve fitting','Location','northeast')
    if abs(e-0.03)<0.001
        times=1
        x'
        y
        pause()
    end
    % Calculate ka/e
    e=(L0-(max(Centerline(:,2))-min(Centerline(:,2))))/L0;
    kaOVERe=k/e*Ly/Ncell;
    kOVERe=k/e;
    fprintf(txt,'%f  %f  %f  %f \n',e,k,kaOVERe,kOVERe);
    figure(3)
    title('Horizontal displacements')
    %     title(['e=-',num2str(e*100),'%, ka/e=',num2str(kaOVERe)])
    
    
%     f3.Position=[235 89 560 420];
    frame = getframe(figure(3));
    writeVideo(v3,frame)

end
close(v1)
close(v2)
close(v3)
fclose(txt);
