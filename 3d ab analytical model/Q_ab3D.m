function [Q] = Q_ab3D(Lx,theta_verti)
% theta_verti = 30;
% theta_hori = 0;
% theta_trans = 0;
% Q matrix of axial bending 3D

nen = 2;    % number of element node
ndof = 6;
n = 2^5;     % number of discretized ligaments
ien = zeros(n,nen);
K_curLx = zeros((n+1)*ndof);
K_curLy = zeros((n+1)*ndof);
K_curLz = zeros((n+1)*ndof);
    for k = 1:n
        ien(k,1:2) = [k,k+1];
    end
im = zeros(n,nen*ndof);
    for p = 1:n
        for q = 1:nen
            for w = 1:ndof*nen
                im(p,w) = (ien(p,1)-1)*ndof + w;                               % The location matrix for single beam
            end 
        end
    end
 
% Material properties
E=9.3;                                                                     % Young's modulus (MPa)
nu=0.3;
G=E/2/(1+nu);                                                                     % Shear modulus
% Geometric parameters Lxyz txyz hxyz 
% start with the same L t h /crosssection
% Lx=1 to verify Q matrix
% Lx = 20e-3;
% Lx = Lx/1000;
t = Lx*0.1;
h = 2*t;
Iz = h*t^3/12; 
Iy = t*h^3/12;
A = t*h; 

%% Stiffness matrix K_beamx for straight beam
% We first start with beams in xyz are of the same legth
L = Lx/n;
K_elz = K_el(Iz,Iy,A,L);

%% Single beam element
% General rotation consists of three angle a(alpha-x),b(beta-y),g(gamma-z)
% with in this case an undulation angle added

% Vertical beam z, all of the beams in x,y,z are designed to have the same
% configuration first
% Undulation matrix
un_verti = zeros(n,1);        %amm z alpha beam in z direction, orientation of alpha
% here we consider achiral
for j=1:n     % achiral
    un_verti(j)=theta_verti*(1-1/n-(j-1)*2/n);
end

un_verti = un_verti*pi/180;
T = eye(12);
for i = 1:n
    % 1 undulation degree 30 degree 'T'
        t_3by3 = [cos(un_verti(i)) sin(un_verti(i)) 0;
            -sin(un_verti(i)) cos(un_verti(i)) 0;
            0                 0                1];
        T(1:3,1:3) = t_3by3;
        T(4:6,4:6) = t_3by3;
        T(7:9,7:9) = t_3by3;
        T(10:12,10:12) =  t_3by3;
        K_ieLz = T'*K_elz*T;
        K_curLz(im(i,:),im(i,:)) = K_curLz(im(i,:),im(i,:)) + K_ieLz;
end

K1(:,1:6)=K_curLz(:,1:6); K1(:,7:12)=K_curLz(:,6*n+1:6*n+6); K1(:,13:6*n+6)=K_curLz(:,7:6*n);
K2(1:6,:)=K1(1:6,:);      K2(7:12,:)=K1(6*n+1:6*n+6,:);      K2(13:6*n+6,:)=K1(7:6*n,:);
Kc_z = K2(1:12,1:12) - K2(1:12,13:end)*inv(K2(13:end,13:end))*K2(13:end,1:12);


% Horizontal beam y

for i = 1:n
        K_curLy(im(i,:),im(i,:)) = K_curLy(im(i,:),im(i,:)) + K_elz;
end

K1(:,1:6)=K_curLy(:,1:6); K1(:,7:12)=K_curLy(:,6*n+1:6*n+6); K1(:,13:6*n+6)=K_curLy(:,7:6*n);
K2(1:6,:)=K1(1:6,:);      K2(7:12,:)=K1(6*n+1:6*n+6,:);      K2(13:6*n+6,:)=K1(7:6*n,:);
Kc_y = K2(1:12,1:12) - K2(1:12,13:end)*inv(K2(13:end,13:end))*K2(13:end,1:12);

Kc_x = Kc_y;
% ###################################################################### %
%% RVE beam assembly
% Vertical beams z
ien = [1,7;
       6,1];
K_O = zeros(7*6);
T = eye(12);
    % transfor to global system        
        % expression for local coordinates
        ez=[0,0,1;
            0,-1,0;
            1,0,0];
        rotz = ez;
for i = 1:size(ien,1)
    % 1 undulation degree 30 degree 'T'
        t_3by3 = rotz;
        T(1:3,1:3) = t_3by3;
        T(4:6,4:6) = t_3by3;
        T(7:9,7:9) = t_3by3;
        T(10:12,10:12) = t_3by3;
        Kz= T'*Kc_z*T;
%         Kz(abs(Kz)<1e-15)=0;
        K_O(6*(ien(1,i)-1)+1:6*(ien(1,i)-1)+6, 6*(ien(1,i)-1)+1:6*(ien(1,i)-1)+6) = ...
        K_O(6*(ien(1,i)-1)+1:6*(ien(1,i)-1)+6, 6*(ien(1,i)-1)+1:6*(ien(1,i)-1)+6) + Kz(1:6,1:6);
        K_O(6*(ien(1,i)-1)+1:6*(ien(1,i)-1)+6, 6*(ien(2,i)-1)+1:6*(ien(2,i)-1)+6) = ...
        K_O(6*(ien(1,i)-1)+1:6*(ien(1,i)-1)+6, 6*(ien(2,i)-1)+1:6*(ien(2,i)-1)+6) + Kz(1:6,7:12);
        K_O(6*(ien(2,i)-1)+1:6*(ien(2,i)-1)+6, 6*(ien(1,i)-1)+1:6*(ien(1,i)-1)+6) = ...
        K_O(6*(ien(2,i)-1)+1:6*(ien(2,i)-1)+6, 6*(ien(1,i)-1)+1:6*(ien(1,i)-1)+6) + Kz(7:12,1:6);
        K_O(6*(ien(2,i)-1)+1:6*(ien(2,i)-1)+6, 6*(ien(2,i)-1)+1:6*(ien(2,i)-1)+6) = ...
        K_O(6*(ien(2,i)-1)+1:6*(ien(2,i)-1)+6, 6*(ien(2,i)-1)+1:6*(ien(2,i)-1)+6) + Kz(7:12,7:12);
end

% Horizontal beam y
        % expression for local coordinates
        ey=[0,1,0;
            0,0,1;
            1,0,0];
ien = [1,2;
       4,1];
T = eye(12);
for i = 1:size(ien,1)
    % 1 undulation degree 30 degree 'T'
        t_3by3 = ey;
        T(1:3,1:3) = t_3by3;
        T(4:6,4:6) = t_3by3;
        T(7:9,7:9) = t_3by3;
        T(10:12,10:12) = t_3by3;
        Ky = T'*Kc_y*T;
%         Ky(abs(Ky)<1e-15)=0;
        K_O(6*(ien(1,i)-1)+1:6*(ien(1,i)-1)+6, 6*(ien(1,i)-1)+1:6*(ien(1,i)-1)+6) = ...
        K_O(6*(ien(1,i)-1)+1:6*(ien(1,i)-1)+6, 6*(ien(1,i)-1)+1:6*(ien(1,i)-1)+6) + Ky(1:6,1:6);
        K_O(6*(ien(1,i)-1)+1:6*(ien(1,i)-1)+6, 6*(ien(2,i)-1)+1:6*(ien(2,i)-1)+6) = ...
        K_O(6*(ien(1,i)-1)+1:6*(ien(1,i)-1)+6, 6*(ien(2,i)-1)+1:6*(ien(2,i)-1)+6) + Ky(1:6,7:12);
        K_O(6*(ien(2,i)-1)+1:6*(ien(2,i)-1)+6, 6*(ien(1,i)-1)+1:6*(ien(1,i)-1)+6) = ...
        K_O(6*(ien(2,i)-1)+1:6*(ien(2,i)-1)+6, 6*(ien(1,i)-1)+1:6*(ien(1,i)-1)+6) + Ky(7:12,1:6);
        K_O(6*(ien(2,i)-1)+1:6*(ien(2,i)-1)+6, 6*(ien(2,i)-1)+1:6*(ien(2,i)-1)+6) = ...
        K_O(6*(ien(2,i)-1)+1:6*(ien(2,i)-1)+6, 6*(ien(2,i)-1)+1:6*(ien(2,i)-1)+6) + Ky(7:12,7:12);
end

% transverse beam x
        % expression for local coordinates
        ex=[1,0,0;
            0,0,1;
            0,-1,0];
ien = [1,5;
       3,1];
T = eye(12);
for i = 1:size(ien,1)
    % 1 undulation degree 30 degree 'T'
        t_3by3 = ex;
        T(1:3,1:3) = t_3by3;
        T(4:6,4:6) = t_3by3;
        T(7:9,7:9) = t_3by3;
        T(10:12,10:12) = t_3by3;
        Kx= T'*Kc_x*T;
%         Kx(abs(Kx)<1e-15)=0;
        K_O(6*(ien(1,i)-1)+1:6*(ien(1,i)-1)+6, 6*(ien(1,i)-1)+1:6*(ien(1,i)-1)+6) = ...
        K_O(6*(ien(1,i)-1)+1:6*(ien(1,i)-1)+6, 6*(ien(1,i)-1)+1:6*(ien(1,i)-1)+6) + Kx(1:6,1:6);
        K_O(6*(ien(1,i)-1)+1:6*(ien(1,i)-1)+6, 6*(ien(2,i)-1)+1:6*(ien(2,i)-1)+6) = ...
        K_O(6*(ien(1,i)-1)+1:6*(ien(1,i)-1)+6, 6*(ien(2,i)-1)+1:6*(ien(2,i)-1)+6) + Kx(1:6,7:12);
        K_O(6*(ien(2,i)-1)+1:6*(ien(2,i)-1)+6, 6*(ien(1,i)-1)+1:6*(ien(1,i)-1)+6) = ...
        K_O(6*(ien(2,i)-1)+1:6*(ien(2,i)-1)+6, 6*(ien(1,i)-1)+1:6*(ien(1,i)-1)+6) + Kx(7:12,1:6);
        K_O(6*(ien(2,i)-1)+1:6*(ien(2,i)-1)+6, 6*(ien(2,i)-1)+1:6*(ien(2,i)-1)+6) = ...
        K_O(6*(ien(2,i)-1)+1:6*(ien(2,i)-1)+6, 6*(ien(2,i)-1)+1:6*(ien(2,i)-1)+6) + Kx(7:12,7:12);
end
% Assembled stiffness
    K_O(abs(K_O)<1e-15)=0;
    K=K_O(7:end,7:end)-K_O(7:end,1:6)/(K_O(1:6,1:6))*K_O(1:6,7:end);

%% Taylor expansion
V = 4*Lx^3; %(*3? decide the coefficient later)
% directional angle matrix ori, different from Euler angle
% for beam 1~6 ori[alpha,beta,gamma] 
ori = [0, pi/2, pi/2;
       pi/2, pi/2, pi;
       pi, pi/2, pi/2;
       pi/2, pi/2, 0;
       pi/2, pi, pi/2;
       pi/2, 0, pi/2];

syms e11 e22 e33 e12 e21 e13 e31 e23 e32 real
syms k11 k22 k33 k12 k21 k13 k31 k23 k32 real
syms a real
Nodedisp = sym('nd', [6 6]);
E12 = (e12+e21)/2;
E21 = E12;
E13 = (e13+e31)/2;
E31 = E13;
E23 = (e23+e32)/2;
E32 = E23;
% A12 = (e12-e21)/2;
% A21 = -A12;
% A13 = (e13-e32)/2;
% A31 = -A13;
% A23 = (e23-e32)/2;
% A32 = -A23;
phi_3 = (e21-e12)/2;
phi_2 = (e13-e31)/2;
phi_1 = (e32-e23)/2;
PHI = [0,0,0,phi_1,phi_2,phi_3]';
% from beam 1-6
NodeVector = [0,a,0;
              -a,0,0;
              0,-a,0;
              a,0,0;
              0,0,-a;
              0,0,a]';
% Gradient of kinematic relationship
KineGrad = [e11, E21, E31;
            E12, e22, E32;
            E13, E23, e33;
            k11, k12, k13;
            k21, k22, k23;
            k31, k32, k33];

    for b = 1:6
        NodeDisp(:,b) = KineGrad*(NodeVector(:,b))*Lx+PHI;
    end
    NodeDisp = subs(NodeDisp,a,1);
    
    % Calculate the strain energy of RVE
    d = reshape(NodeDisp,36,1);
    w = d'*K*d/(2*V);
    
    strain = [e11 e22 e33 e23 e13 e12 e32 e31 e21 ...
        k11 k22 k33 k23 k13 k12 k32 k31 k21];
    for m = 1:length(strain)
        for n = 1:length(strain)
            Q(m,n) = diff(w,strain(m),strain(n));
        end
    end
    
    Q = double(Q);
    Q(abs(Q)<1e-5) = 0;


    inv_Q = inv(Q);
    Q(abs(Q)<1e-4) = 0;    
    v = -inv_Q(1,2)/inv_Q(1,1);
    E_effect = 1/inv_Q(3,3);
    G_effect = 0.5/(inv_Q(3,3)+inv_Q(3,4));
    eta1=inv_Q(3,14)/inv_Q(3,3);
end