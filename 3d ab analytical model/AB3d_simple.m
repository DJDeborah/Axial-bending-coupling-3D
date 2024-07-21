%% 3D axial bending coupling
clear all;
close all;
% parameters
global L epsilon nelx nely nelz D macroWidth macroHeight macroLength nn IEN coords
L = 20e-3;
a = L/2; % beam length
theta = 30;
epsilon = -0.1;
DeformFactor = 1;

% homogenized material bulk - output
nelx = 8;
nely = 8;
nelz = 8;
macroLength = L;
macroWidth = L; 
macroHeight = L; 

[coords,IEN] = Cubemesh();

Latticeplot(coords,IEN);

[D] = Q_ab3D(a,theta);

% B matrix
    % gauss integration informations
    gp = 0.7746;
    xiRow = [-gp,0,gp];
    etaRow = xiRow;
    zetaRow = xiRow;
    
    w = [0.5556 0.8889 0.5556];

    nn = size(coords,1);
    ne = size(IEN,1);
    nen = size(IEN,2);  % number of element node, 20   
    ndof = nn*6;    % number of degrees of freedom per node, 6

    
    K = sparse(ndof,ndof);
    LM = zeros(ne,nen*6);
    
    for p = 1:ne
        for t = 1:nen
            for r = 1:6
                ind = (t-1)*6 + r;
                LM(p,ind) = 6*IEN(p,t) - 6 + r ;                           % create the LM matrix
            end 
        end
    end    

    for e = 1:ne
        coord = coords(IEN(e,:),:);
        
        ke = sparse(nen*6,nen*6);
        for i = 1:3
            for j = 1:3
                for k = 1:3
                    xi = xiRow(i);
                    eta = etaRow(j);
                    zeta = zetaRow(k);
                    [B, J_det] = Bmatrix3D(xi,eta,zeta,coord);
                    ke = ke + w(i)*w(j)*w(k)*B'*D*B*J_det;
                end
            end
        end
        K(LM(e,:),LM(e,:)) = K(LM(e,:),LM(e,:)) + ke; 

    end
% Boundary condition
BCtype = 'stretchZfree';
[u,free,essential,~]=BoundaryCondition3D(coords,BCtype);

% Solution
K_ff = K(free,free);
F_f = zeros(length(free),1);
u(free) = K_ff \ (F_f-K(free,essential)*u(essential));
StrainEnergyCont = u'*K*u/2;

%% Plot
% CubVect = [1,3,8,6,13,15,20,18,2,5,7,4,14,17,19,16,9,10,12,11];
CubDraw = [1,9,2,10,3,11,4,12,1,...
         1,17,5,16,8,20,4,...
         4,20,8,15,7,19,3,...
         3,19,7,14,6,18,2,...
         2,18,6,13,5,17,1];
figure()

for e = 1:ne
    coord = coords(IEN(e,:),:);

    % Plot the element outline and displacement
        coordD = zeros(size(coord,1),3);
        nodes = IEN(e,:);
        for temp = 1:size(coord,1)
            coordD(temp,3)=coord(temp,3)+DeformFactor*u(6*nodes(temp)-3);
            coordD(temp,2)=coord(temp,2)+DeformFactor*u(6*nodes(temp)-4);
            coordD(temp,1)=coord(temp,1)+DeformFactor*u(6*nodes(temp)-5);
        end
           
%         coordI = coord;
        % use the connectivity Cubdraw to draw the 3d lattice
        coordDraw = coord(CubDraw,:);
        coordDDraw = coordD(CubDraw,:);
        plot3(coordDraw(:,1),coordDraw(:,2),coordDraw(:,3),'-g');
        hold on;
        plot3(coordDDraw(:,1),coordDDraw(:,2),coordDDraw(:,3), '-b');
        xlabel('x')
        ylabel('y')
        zlabel('z')
end
% view([-90 0 0])
% % Set up 3D plot to record
% daspect([1,1,1]);axis tight;
% 
% % Set up recording parameters (optional), and record
% OptionZ.FrameRate=15;OptionZ.Duration=50;OptionZ.Periodic=true;
% CaptureFigVid([-20,0;-110,0;-120,80;-200,90;-330,90], 'ab',OptionZ)
% z= 1;
