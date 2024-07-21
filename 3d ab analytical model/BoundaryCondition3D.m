function [u,free,essential,F] = BoundaryCondition3D(Pos,Btype)
% Pos: All the global position of nodes 'coords'
% Btype: Boundary condition type
% False: isLattice
global L epsilon
% Variables 
nn = size(Pos,1);           % number of total nodes
macroLength = max(Pos(:,1));
macroWidth = max(Pos(:,2));
macroHeight = max(Pos(:,3));

xfixed = zeros(nn,1);
yfixed = zeros(nn,1);
zfixed = zeros(nn,1);
alphafixed = zeros(nn,1);
betafixed = zeros(nn,1);
gammafixed = zeros(nn,1);
F = zeros(nn*3,1);

xmoves = zeros(nn,1);       xMag = 0;      
ymoves = zeros(nn,1);       yMag = 0;      
zmoves = zeros(nn,1);       zMag = 0;      
alphamoves = zeros(nn,1);   alphaMag = 0;      
betamoves = zeros(nn,1);    betaMag = 0;      
gammamoves = zeros(nn,1);   gammaMag = 0;      

uRot = zeros(nn,3);
tol = 10^-5*L;

% Difinition of different types of boundary conditions
switch Btype
    case 'stretchZ'
    % fix y? how to define a overconstraint?
%         bottomnodes = abs(Pos(:,3))<tol;
%         rot0 = 
    case 'stretchZfree'
    % bottom plane - fix z
    % bottom-beam lay towards the x direction - fix y
    % bottom-beam lay towards the y direction - fix x
        bottomnodes = abs(Pos(:,3)) <= tol;
        topnodes = abs(Pos(:,3)-macroHeight) <= tol;
        zfixed = bottomnodes;
        yfixed = (abs(Pos(:,1)) < tol & abs(Pos(:,3)) < tol);
        xfixed = (abs(Pos(:,2)) < tol & abs(Pos(:,3)) < tol);
        rot0 = xfixed;
        gammafixed = rot0;
%         betafixed = rot0;
        zmoves = topnodes;

        zMag = epsilon*L;
    case 'stretchYfree'
    % bottom plane fix y
        bottomnodes = abs(Pos(:,2)) <= tol;
        topnodes = abs(Pos(:,2)-macroWidth) <= tol;
        zfixed = (abs(Pos(:,1)) < tol & abs(Pos(:,2)) < tol);
        yfixed = bottomnodes;
        xfixed = (abs(Pos(:,2)) < tol & abs(Pos(:,3)) < tol);
        ymoves = topnodes;
        rot0 = xfixed;
        gammafixed = rot0;
        yMag = epsilon*L;
end

% Assign value to boundary condition
essential = zeros(nn*6,1);
essential(1:6:end,1) = or(xfixed,xmoves);
essential(2:6:end,1) = or(yfixed,ymoves);
essential(3:6:end,1) = or(zfixed,zmoves);
essential(4:6:end,1) = or(alphafixed,alphamoves);
essential(5:6:end,1) = or(betafixed,betamoves);
essential(6:6:end,1) = or(gammafixed,gammamoves);

free = find(1-essential);
essential = find(essential);

% nodes with displacement
u = zeros(nn*6,1);
u(1:6:end,1) = xmoves*xMag+uRot(:,1);
u(2:6:end,1) = ymoves*yMag+uRot(:,2);
u(3:6:end,1) = zmoves*zMag+uRot(:,3);
u(4:6:end,1) = alphamoves*alphaMag;
u(5:6:end,1) = betamoves*betaMag;
u(6:6:end,1) = gammamoves*gammaMag;





end