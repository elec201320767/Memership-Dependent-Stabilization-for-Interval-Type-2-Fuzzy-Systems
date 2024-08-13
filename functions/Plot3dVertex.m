%% Plot membership funciotns
%% 
function [VerticesL,VerticesU] = Plot3dVertex(LowerMF,UpperMF,mfsum_min,mfsum_max)
Vertices = Type2FuzzyVertices([LowerMF; UpperMF]);% Find vertex set for theta
 
VerticesLower = Type2FuzzyVerticesBound([LowerMF; UpperMF],mfsum_min);% Find vertex set for thetaU
VerticesUpper = Type2FuzzyVerticesBound([LowerMF; UpperMF],mfsum_max);% Find vertex set for thetaL
ExtremeRect = ExtremeBox([LowerMF; UpperMF]);% Find vertex set for extrema

VerticesL = [Vertices VerticesLower];
% VerticesL = [VerticesL ExtremeRect(2,:)' ExtremeRect(5,:)'];
VerticesU = [Vertices VerticesUpper]; 
SizeVertices = size(Vertices); SizeVerticesL = size(VerticesL); SizeVerticesU = size(VerticesUpper); SizeBox = size(ExtremeRect);

   
%% Vertex plot
figure()  
for i = 1 : SizeVertices(2)
    plot3(Vertices(1,i),Vertices(2,i),Vertices(3,i),'o','MarkerEdgeColor','b');
    hold on;  
end
for i = 1 : SizeVerticesU(2)
    plot3(VerticesUpper(1,i),VerticesUpper(2,i),VerticesUpper(3,i),'o','MarkerEdgeColor','r');
    hold on;  
end   
for i = 1 : SizeVerticesL(2)
    plot3(VerticesL(1,i),VerticesL(2,i),VerticesL(3,i),'o','MarkerEdgeColor','g');
    hold on;   
end  
axis([0 1.5 0 1.5 0 1.5])
vertex=[1 0 0; 0 1 0; 0 0 1];
face=[1 2 3];
MembershipColor = [0.5 0.5 0.5]; BoundColor = [0.7 0.7 0.7];
patch('Faces',face,'Vertices',vertex,'Facecolor',MembershipColor,'FaceAlpha',.5);
hold on
vertex=[mfsum_max 0 0; 0 mfsum_max 0; 0 0 mfsum_max];
patch('Faces',face,'Vertices',vertex,'Facecolor',BoundColor,'FaceAlpha',.5);
vertex=[mfsum_min 0 0; 0 mfsum_min 0; 0 0 mfsum_min];
patch('Faces',face,'Vertices',vertex,'Facecolor',BoundColor,'FaceAlpha',.5);
view([1.5 1.5 1.5]);
face=[1 2 4 6 5 3];
patch('Faces',face,'Vertices',Vertices','Facecolor',[0.5 0.5 0.5],'FaceAlpha',.5);
grid on
set(gcf, 'Position',  [300, 300, 800 , 600])
xlabel('$\theta_1(t)$','Interpreter','latex')
ylabel('$\theta_2(t)$','Interpreter','latex')
zlabel('$\theta_3(t)$','Interpreter','latex')
% legend('$\dot\theta_1(x_1(t))$ Boundary','$\dot\theta_2(x_1(t))$ Boundary','$\dot\theta_3(x_1(t))$ Boundary','Interpreter','Latex')
%% (Vertex + Box) Plot
figure() 
for i = 1 : SizeVertices(2)
    plot3(Vertices(1,i),Vertices(2,i),Vertices(3,i),'o','MarkerEdgeColor','b');
    hold on;  
end
for i = 1 : SizeVerticesU(2)
    plot3(VerticesUpper(1,i),VerticesUpper(2,i),VerticesUpper(3,i),'o','MarkerEdgeColor','r');
    hold on;  
end    
for i = 1 : SizeVerticesL(2)
    plot3(VerticesL(1,i),VerticesL(2,i),VerticesL(3,i),'o','MarkerEdgeColor','g');
    hold on;   
end    
for i = 1 : SizeBox(1)
    plot3(ExtremeRect(i,1),ExtremeRect(i,2),ExtremeRect(i,3),'o','MarkerEdgeColor','k');
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
patch('Faces',face,'Vertices',Vertices','Facecolor',[0.5 0.5 0.5],'FaceAlpha',.5,'EdgeColor','b');
grid on
set(gcf, 'Position',  [300, 300, 800 , 600])
xlabel('$h_1(x(t))$','Interpreter','latex')
ylabel('$h_2(x(t))$','Interpreter','latex')
zlabel('$h_3(x(t))$','Interpreter','latex')


end
 

