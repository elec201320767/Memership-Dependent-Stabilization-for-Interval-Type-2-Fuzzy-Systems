%% Plot membership funciotns
%% 
function PlotPolytopes(LowerMF,UpperMF,mfsum_min,mfsum_max)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig.1 : Box + standard simplex + upper/lower simplex
% Fig.2 : membership distribution
% Fig.3 : upper membership distribution
% Fig.4 : lower membershup disrtibution
% Vertex Colors:    r- inter sections btw upper simplex and Box
%                   g- inter sections btw standard simplex and Box
%                   b- inter sections btw lower simplex and Box
%                   k- Box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vertices = Type2FuzzyVertices([LowerMF; UpperMF]);% Find vertex set for theta
VerticesLower = Type2FuzzyVerticesBound([LowerMF; UpperMF],mfsum_min);% Find vertex set for thetaU
VerticesUpper = Type2FuzzyVerticesBound([LowerMF; UpperMF],mfsum_max);% Find vertex set for thetaL
ExtremeRect = ExtremeBox([LowerMF; UpperMF]);% Find vertex set for extrema

VerticesL = [Vertices VerticesLower];
% VerticesL = [VerticesL ExtremeRect(2,:)' ExtremeRect(5,:)'];
VerticesU = [Vertices VerticesUpper]; 
SizeVertices = size(Vertices); SizeVerticesLower = size(VerticesLower); SizeVerticesUpper = size(VerticesUpper); SizeBox = size(ExtremeRect);
SizeVerticesL = size(VerticesL); SizeVerticesU = size(VerticesU); 

vertex=[1 0 0; 0 1 0; 0 0 1];
face=[1 2 3];
MembershipColor = [0.5 0.5 0.5]; BoundColor = [0.7 0.7 0.7];
%% Fig.1 (Vertex + Box) Plot
figure(1) 
MrkSz = 3;
for i = 1 : SizeVertices(2)
    plot3(Vertices(1,i),Vertices(2,i),Vertices(3,i),'o','MarkerEdgeColor','g','MarkerSize',MrkSz);
    hold on;  
end
for i = 1 : SizeVerticesUpper(2)
    plot3(VerticesUpper(1,i),VerticesUpper(2,i),VerticesUpper(3,i),'o','MarkerEdgeColor','r','MarkerSize',MrkSz);
    hold on;  
end    
for i = 1 : SizeVerticesLower(2)
    plot3(VerticesLower(1,i),VerticesLower(2,i),VerticesLower(3,i),'o','MarkerEdgeColor','b','MarkerSize',MrkSz);
    hold on;   
end    
for i = 1 : SizeBox(1)
    plot3(ExtremeRect(i,1),ExtremeRect(i,2),ExtremeRect(i,3),'o','MarkerEdgeColor','k','MarkerSize',MrkSz);
    hold on;  
end

Edge(:,:,1) = [ExtremeRect(1,:);ExtremeRect(2,:)];  Edge(:,:,2) = [ExtremeRect(1,:);ExtremeRect(3,:)];  Edge(:,:,3) = [ExtremeRect(1,:);ExtremeRect(5,:)];
Edge(:,:,4) = [ExtremeRect(4,:);ExtremeRect(2,:)];  Edge(:,:,5) = [ExtremeRect(4,:);ExtremeRect(3,:)];  Edge(:,:,6) = [ExtremeRect(4,:);ExtremeRect(8,:)];
Edge(:,:,7) = [ExtremeRect(6,:);ExtremeRect(2,:)];  Edge(:,:,8) = [ExtremeRect(6,:);ExtremeRect(5,:)];  Edge(:,:,9) = [ExtremeRect(6,:);ExtremeRect(8,:)];
Edge(:,:,10) = [ExtremeRect(7,:);ExtremeRect(3,:)]; Edge(:,:,11) = [ExtremeRect(7,:);ExtremeRect(5,:)]; Edge(:,:,12) = [ExtremeRect(7,:);ExtremeRect(8,:)];

for i = 1 : length(Edge)
    plot3(Edge(:,1,i),Edge(:,2,i),Edge(:,3,i),'k');
    hold on;  
end
axis([0 1.5 0 1.5 0 1.5])
vertex=[1 0 0; 0 1 0; 0 0 1];
face=[1 2 3];
patch('Faces',face,'Vertices',vertex,'Facecolor',MembershipColor,'FaceAlpha',.5);
hold on
vertex=[mfsum_max 0 0; 0 mfsum_max 0; 0 0 mfsum_max];
patch('Faces',face,'Vertices',vertex,'Facecolor',BoundColor,'FaceAlpha',.5);
vertex=[mfsum_min 0 0; 0 mfsum_min 0; 0 0 mfsum_min];
patch('Faces',face,'Vertices',vertex,'Facecolor',BoundColor,'FaceAlpha',.5);
view([1.5 1.2 1.2]);
face=[1 2 4 6 5 3];
% patch('Faces',face,'Vertices',Vertices','Facecolor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','b');
grid on
set(gcf, 'Position',  [300, 300, 800 , 600])
% set(gcf, 'Position',  [100, 100, 400 , 300])
xlabel('$h_1$','Fontsize',15,'Interpreter','latex')
ylabel('$h_2$','Fontsize',15,'Interpreter','latex')
zlabel('$h_3$','Fontsize',15,'Interpreter','latex')

%% Fig.2 : Distribution for MF
figure(2) 
for i = 1 : SizeVertices(2)
    plot3(Vertices(1,i),Vertices(2,i),Vertices(3,i),'o','MarkerEdgeColor','g','MarkerSize',MrkSz);
    hold on;  
end
patch('Faces',face,'Vertices',Vertices','Facecolor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','k');
ax_view = [0.1 0.7 0.1 0.4 0.1 0.7];
axis(ax_view)
view([1.5 1.2 1.2]);
face=[1 2 4 6 5 3];
grid on
% set(gcf, 'Position',  [100, 100, 400 , 300])
set(gcf, 'Position',  [300, 300, 800 , 600])
xlabel('$h_1$','Interpreter','latex')
ylabel('$h_2$','Interpreter','latex')
zlabel('$h_3$','Interpreter','latex')

%% Fig.3 : Distribution for UMF
figure(3)
for i = 1 : SizeVerticesU(2)
    plot3(VerticesU(1,i),VerticesU(2,i),VerticesU(3,i),'o','MarkerEdgeColor','r','MarkerSize',MrkSz);
    hold on;  
end   

patch('Faces',face,'Vertices',Vertices','Facecolor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','k');
hold on
patch('Faces',face,'Vertices',VerticesUpper','Facecolor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','k');
hold on
face = [1 2 3 4];
v = [VerticesU(:,3) VerticesU(:,9) VerticesU(:,11) VerticesU(:,5)]';
patch('Faces',face,'Vertices',v,'Facecolor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','k');
hold on
v = [VerticesU(:,11) VerticesU(:,5) VerticesU(:,6) VerticesU(:,12)]';
patch('Faces',face,'Vertices',v,'Facecolor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','k');
hold on
v = [VerticesU(:,6) VerticesU(:,12) VerticesU(:,10) VerticesU(:,4)]';
patch('Faces',face,'Vertices',v,'Facecolor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','k');
hold on
v = [VerticesU(:,10) VerticesU(:,4) VerticesU(:,2) VerticesU(:,8)]';
patch('Faces',face,'Vertices',v,'Facecolor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','k');
hold on
v = [VerticesU(:,8) VerticesU(:,2) VerticesU(:,1) VerticesU(:,7)]';
patch('Faces',face,'Vertices',v,'Facecolor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','k');
hold on
v = [VerticesU(:,1) VerticesU(:,7) VerticesU(:,9) VerticesU(:,3)]';
patch('Faces',face,'Vertices',v,'Facecolor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','k');
hold on
axis(ax_view)
view([1.5 1.2 1.2]);
face=[1 2 4 6 5 3];
grid on
% set(gcf, 'Position',  [100, 100, 400 , 300])
set(gcf, 'Position',  [300, 300, 800 , 600])
xlabel('$h_1$','Interpreter','latex')
ylabel('$h_2$','Interpreter','latex')
zlabel('$h_3$','Interpreter','latex')
%% Fig.4 : Distribution for UMF
figure(4)
for i = 1 : SizeVerticesL(2)
    plot3(VerticesL(1,i),VerticesL(2,i),VerticesL(3,i),'o','MarkerEdgeColor','b','MarkerSize',MrkSz);
    hold on;  
end   

patch('Faces',face,'Vertices',Vertices','Facecolor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','k');
hold on
face = [1 2 3 4];
patch('Faces',face,'Vertices',VerticesLower','Facecolor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','k');
hold on
v = [VerticesL(:,3) VerticesL(:,7) VerticesL(:,10) VerticesL(:,5)]';
patch('Faces',face,'Vertices',v,'Facecolor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','k');
hold on 
v = [VerticesL(:,7) VerticesL(:,1) VerticesL(:,2) VerticesL(:,8)]';
patch('Faces',face,'Vertices',v,'Facecolor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','k');
hold on 
v = [VerticesL(:,2) VerticesL(:,8) VerticesL(:,9) VerticesL(:,4)]';
patch('Faces',face,'Vertices',v,'Facecolor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','k');
hold on 
v = [VerticesL(:,9) VerticesL(:,4) VerticesL(:,6) VerticesL(:,10)]';
patch('Faces',face,'Vertices',v,'Facecolor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','k');
hold on 
face = [1 2 3];
v = [VerticesL(:,11) VerticesL(:,3) VerticesL(:,7)]';
patch('Faces',face,'Vertices',v,'Facecolor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','k');
hold on 
v = [VerticesL(:,11) VerticesL(:,1) VerticesL(:,7)]';
patch('Faces',face,'Vertices',v,'Facecolor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','k');
hold on 
v = [VerticesL(:,11) VerticesL(:,1) VerticesL(:,3)]';
patch('Faces',face,'Vertices',v,'Facecolor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','k');
hold on 
v = [VerticesL(:,10) VerticesL(:,5) VerticesL(:,12)]';
% v = [VerticesL(:,10) VerticesL(:,5) VerticesL(:,12) VerticesL(:,6)]';
patch('Faces',face,'Vertices',v,'Facecolor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','k');
hold on 
v = [VerticesL(:,5) VerticesL(:,12) VerticesL(:,6)]';
patch('Faces',face,'Vertices',v,'Facecolor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','k');
hold on 
v = [VerticesL(:,10) VerticesL(:,12) VerticesL(:,6)]';
patch('Faces',face,'Vertices',v,'Facecolor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','k');
hold on 

axis(ax_view)
view([1.5 1.2 1.2]); 
grid on
% set(gcf, 'Position',  [100, 100, 400 , 300])
set(gcf, 'Position',  [300, 300, 800 , 600])
xlabel('$h_1$','Interpreter','latex')
ylabel('$h_2$','Interpreter','latex')
zlabel('$h_3$','Interpreter','latex')
end
 

