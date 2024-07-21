  function [B, J_det] = Bmatrix3D(xi,eta,zeta,coord)
%% Shape function of 20 nodes
N = zeros(1,20);
% Shape function derivative
dN = zeros(3,20);
    % eight nodes at the corner
    Loc8 = [-1, -1, -1;
            1, -1, -1;
            1,  1, -1;
           -1,  1, -1;
           -1, -1,  1;
            1, -1,  1;
            1,  1,  1;
           -1,  1,  1];
    for i = 1:8
        cx = Loc8(i,1);
        ce = Loc8(i,2); 
        cz = Loc8(i,3);

        N(i) = 1/8*(1+cx*xi)*(1+ce*eta)*(1+cz*zeta)...
            *(cx*xi+ce*eta+cz*zeta-2);

        dN(1,i) = 1/8*cx*(1+ce*eta)*(1+cz*zeta)*(cx*xi+ce*eta+cz*zeta-2)...
                    +1/8*(1+cx*xi)*(1+ce*eta)*(1+cz*zeta)*cx;
        dN(2,i) = 1/8*(1+cx*xi)*ce*(1+cz*zeta)*(cx*xi+ce*eta+cz*zeta-2)...
                    +1/8*(1+cx*xi)*(1+ce*eta)*(1+cz*zeta)*ce;
        dN(3,i) = 1/8*(1+cx*xi)*(1+ce*eta)*cz*(cx*xi+ce*eta+cz*zeta-2)...
                    +1/8*(1+cx*xi)*(1+ce*eta)*(1+cz*zeta)*cz;
    end

    im_local = [9,11,15,13, 10,12,16,14, 17,18,19,20];
    
    % 4 mid nodes towards xi 
    Loc4 = [0, -1, -1;
            0, 1, -1;
            0, 1, 1;
            0, -1, 1];
    for j = 1:4
%         cx = Loc4(j,1);
        ce = Loc4(j,2);
        cz = Loc4(j,3);

        N(im_local(j)) = 1/4*(1-xi^2)*(1+ce*eta)*(1+cz*zeta);

        dN(1,im_local(j)) = 1/4*(-2*xi)*(1+ce*eta)*(1+cz*zeta);
        dN(2,im_local(j)) = 1/4*(1-xi^2)*(ce)*(1+cz*zeta);
        dN(3,im_local(j)) = 1/4*(1-xi^2)*(1+ce*eta)*cz;
    end

    % 4 mid nodes towards eta
    Loc4 = [1, 0 ,-1;
           -1, 0, -1;
           -1, 0, 1;
           1, 0, 1];
    for j = 1:4
        cx = Loc4(j,1);
%         ce = Loc4(j,2);
        cz = Loc4(j,3);

        N(im_local(j+4)) = 1/4*(1-eta^2)*(1+cx*xi)*(1+cz*zeta);

        dN(1,im_local(j+4)) = 1/4*(1-eta^2)*cx*(1+cz*zeta);
        dN(2,im_local(j+4)) = 1/4*(-2*eta)*(1+cx*xi)*(1+cz*zeta);
        dN(3,im_local(j+4)) = 1/4*(1-eta^2)*(1+cx*xi)*cz;
    end

    % 4 mid nodes towards zeta
    Loc4 = [-1, -1, 0;
            1, -1, 0;
            1, 1, 0;
            -1, 1,0];
    for j = 1:4
        cx = Loc4(j,1);
        ce = Loc4(j,2);
%         cz = Loc4(j,3);
        
        N(im_local(j+8)) = 1/4*(1-zeta^2)*(1+cx*xi)*(1+ce*eta);

        dN(1,im_local(j+8)) = 1/4*(1-zeta^2)*cx*(1+ce*eta);
        dN(2,im_local(j+8)) = 1/4*(1-zeta^2)*(1+cx*xi)*ce;
        dN(3,im_local(j+8)) = 1/4*(-2*zeta)*(1+cx*xi)*(1+ce*eta);
    end

%% B matrix
% Jacobian matrix
    J = dN*coord;
    J_det = det(J);

% B matrix for an element
    Be = J\dN;

% B for micropolar
    B = sparse(12,20*6);

    B(1,1:6:120) = Be(1,:); % dN/dx
    B(2,2:6:120) = Be(2,:); % dN/dy
    B(3,3:6:120) = Be(3,:); % dN/dz

    B(4,2:6:120) = Be(3,:);
    B(4,4:6:120) = -N;

    B(5,3:6:120) = Be(1,:);
    B(5,5:6:120) = N;

    B(6,2:6:120) = Be(1,:);
    B(6,6:6:120) = -N;

    B(7,3:6:120) = Be(2,:);
    B(7,4:6:120) = N;

    B(8,1:6:120) = Be(3,:);
    B(8,5:6:120) = -N;

    B(9,1:6:120) = Be(2,:);
    B(9,6:6:120) = N;

% k rotation gradient
    
    B(10,4:6:120) = Be(1,:);
    B(11,5:6:120) = Be(2,:);
    B(12,6:6:120) = Be(3,:);
    B(13,5:6:120) = Be(3,:);
    B(14,4:6:120) = Be(3,:);
    B(15,4:6:120) = Be(2,:);
    B(16,6:6:120) = Be(2,:);
    B(17,6:6:120) = Be(1,:);
    B(18,5:6:120) = Be(1,:);    

end

