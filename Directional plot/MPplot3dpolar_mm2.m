% 3d plot
clear all; close all;
clc; 
%% Input 3
% Division / Elements of the angle theta, phi, psi
n = 500; 
m = n/2;
o = 100;
% 3 rotation dofs in 3d 
theta = linspace(0,2*pi,n);
phi = linspace(0,pi,m);
% psi = linspace(0,2*pi,o);

%% Q matrix 3d
% Vertical curved beam
C1 = [2840000000	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	2840000000	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	2382327204	-0.005316491	0	0	-0.006038021	0	0	0	0	0	0	-690212.3196	0	0	0	0;
0	0	-0.005316491	19855563.5	0	0	8635110.759	0	0	0	0	0	0	0	-2429.258193	0	0	2429.258193;
0	0	0	0	51176919	0.000113024	0	15665577.7	-0.000113023	0	0	22641.47552	0	0	0	0	0	0;
0	0	0	0	0.000113024	73507338.53	0	0.00010532	40092661.47	0	0	0	-16894.83436	0	0	-653.1088588	0	0;
0	0	-0.006038021	8635110.759	0	0	19845168.33	0	0	0	0	0	0	0	2427.007613	0	0	-2427.007613;
0	0	0	0	15665577.7	0.00010532	0	24552066.99	-0.00010532	0	0	6930.698459	0	0	0	0	0	0;
0	0	0	0	-0.000113023	40092661.47	0	-0.00010532	73507338.53	0	0	0	16894.83436	0	0	653.1088596	0	0;
0	0	0	0	0	0	0	0	0	450	0	0	0	0	0	0	0	0;
0	0	0	0	0	0	0	0	0	0	450	0	0	0	0	0	0	0;
0	0	0	0	22641.47552	0	0	6930.698459	0	0	0	462.4100012	0	0	0	0	0	0;
0	0	0	0	0	-16894.83436	0	0	16894.83436	0	0	0	3627.320976	0	0	104.029644	0	0;
0	0	-690212.3196	0	0	0	0	0	0	0	0	0	0	1136.91999	0	0	27.07108491	0;
0	0	0	-2429.258193	0	0	2427.007613	0	0	0	0	0	0	0	937.7112473	0	0	8.955419364;
0	0	0	0	0	-653.1088588	0	0	653.1088596	0	0	0	104.029644	0	0	3681.457404	0	0;
0	0	0	0	0	0	0	0	0	0	0	0	0	27.07108491	0	0	3678.482393	0;
0	0	0	2429.258193	0	0	-2427.007613	0	0	0	0	0	0	0	8.955419364	0	0	937.7112473]; 

% 3 Curved beams
C2 = [2840000000	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	2840000000	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	1614870649	-0.003909795	0	0	-0.003540604	0	0	0	0	0	0	-927215.344	0	0	0	0;
0	0	-0.003909795	20035610.04	0	0	8730288.8	0	0	0	0	0	0	0	-4729.187425	0	0	4729.187425;
0	0	0	0	50663686.87	0	0	15971593.31	0	0	0	44626.04188	0	0	0	0	0	0;
0	0	0	0	0	73843360.01	0	0	39756639.99	0	0	0	-32283.8987	0	0	-1267.242395	0	0;
0	0	-0.003540604	8730288.8	0	0	19994045.16	0	0	0	0	0	0	0	4711.800209	0	0	-4711.800209;
0	0	0	0	15971593.31	0	0	24682087.3	0	0	0	14068.24169	0	0	0	0	0	0;
0	0	0	0	0	39756639.99	0	0	73843360.01	0	0	0	32283.8987	0	0	1267.242396	0	0;
0	0	0	0	0	0	0	0	0	450	0	0	0	0	0	0	0	0;
0	0	0	0	0	0	0	0	0	0	450	0	0	0	0	0	0	0;
0	0	0	0	44626.04188	0	0	14068.24169	0	0	0	498.8556841	0	0	0	0	0	0;
0	0	0	0	0	-32283.8987	0	0	32283.8987	0	0	0	3472.34176	0	0	100.6966594	0	0;
0	0	-927215.344	0	0	0	0	0	0	0	0	0	0	1460.463008	0	0	27.14013555	0;
0	0	0	-4729.187425	0	0	4711.800209	0	0	0	0	0	0	0	931.9642997	0	0	14.70236698;
0	0	0	0	0	-1267.242395	0	0	1267.242396	0	0	0	100.6966594	0	0	3681.388551	0	0;
0	0	0	0	0	0	0	0	0	0	0	0	0	27.14013555	0	0	3678.501231	0;
0	0	0	4729.187425	0	0	-4711.800209	0	0	0	0	0	0	0	14.70236698	0	0	931.9642997;];

C3=[2840000000	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	2840000000	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0;
0	0	1063340071	0.000396958	0	0	0.000378361	0	0	0	0	0	0	-901921.9939	0	0	0	0;
0	0	0.000396958	20332531.8	0	0	8902842.311	0	0	0	0	0	0	0	-6791.376231	0	0	6791.376231;
0	0	0	0	49890871.67	0.000160984	0	16455740	-0.000160984	0	0	65398.33434	0	0	0	0	0	0;
0	0	0	0	0.000160984	74347792.67	0	0.000147701	39252207.33	0	0	0	-44950.99964	0	0	-1810.402935	0	0;
0	0	0.000378361	8902842.311	0	0	20239311.09	0	0	0	0	0	0	0	6735.985673	0	0	-6735.985674;
0	0	0	0	16455740	0.000147701	0	24909072.38	-0.000147701	0	0	21570.63908	0	0	0	0	0	0;
0	0	0	0	-0.000160984	39252207.33	0	-0.000147701	74347792.67	0	0	0	44950.99964	0	0	1810.402935	0	0;
0	0	0	0	0	0	0	0	0	450	0	0	0	0	0	0	0	0;
0	0	0	0	0	0	0	0	0	0	450	0	0	0	0	0	0	0;
0	0	0	0	65398.33434	0	0	21570.63908	0	0	0	557.1004519	0	0	0	0	0	0;
0	0	0	0	0	-44950.99964	0	0	44950.99964	0	0	0	3236.72395	0	0	95.55243299	0	0;
0	0	-901921.9939	0	0	0	0	0	0	0	0	0	0	1678.286196	0	0	27.25056511	0;
0	0	0	-6791.376231	0	0	6735.985673	0	0	0	0	0	0	0	923.1156711	0	0	23.55099559;
0	0	0	0	0	-1810.402935	0	0	1810.402935	0	0	0	95.55243299	0	0	3681.284275	0	0;
0	0	0	0	0	0	0	0	0	0	0	0	0	27.25056511	0	0	3678.533415	0;
0	0	0	6791.376231	0	0	-6735.985674	0	0	0	0	0	0	0	23.55099559	0	0	923.1156711];
dimap = [0.403922000	0.000000000	0.121569000
	0.698039000	0.094118000	0.168627000
	0.839216000	0.376471000	0.301961000
	0.956863000	0.647059000	0.509804000
	0.992157000	0.858824000	0.780392000
	0.968627000	0.968627000	0.968627000
	0.819608000	0.898039000	0.941176000
	0.572549000	0.772549000	0.870588000
	0.262745000	0.576471000	0.764706000
	0.129412000	0.400000000	0.674510000
	0.019608000	0.188235000	0.380392000];
N = 30;
newmap = interp1(1:size(dimap,1),dimap,linspace(1,size(dimap,1),N),'linear');


%% Initialization
E = zeros (m,n); 
% Shear modulus and Poisson's ratio 
G = zeros (m,n); 
niu = zeros(m,n);
% Normal vector 
nn = zeros (3,1); 
% coupling coefficient
ab = zeros(m,n);

% Initialize of transformation matrices
d_v = zeros(18,1);
d_vs = zeros(18,1);
n_v = zeros(18,1);
T_v = zeros(18,1);
n_vs = zeros(18,1);
t_v = zeros(18,1);
m_v = zeros(18,1);
d_ab = zeros(18,1);
% x_E = zeros (m,n); 
% y_E = zeros (m,n); 
% z_E = zeros (m,n); 
% x_K = zeros (m,n); 
% y_K = zeros (m,n); 
% z_K = zeros (m,n); 
% x_niu = zeros (m,n); 
% y_niu = zeros (m,n); 
% z_niu = zeros (m,n); 
%% Normal Moduli 
S1 = inv(C1);
S2 = inv(C2);
S3 = inv(C3);
U_V = zeros(18,1);
U_V(1:3) = 1;
for i = 1:n
    for j = 1:m
    % position vector r 
    d = [sin(phi(j))*cos(theta(i));
        sin(phi(j))*sin(theta(i));
        cos(phi(j))];
    % d dyac d 
    D = d*d';
    Ds = 0.5*(D+D');
    Da = 0.5*(D-D');
    % In Voigt without consider coefficient first
    d_v(1:9) = [D(1,1);
                D(2,2);
                D(3,3);
                D(2,3);
                D(1,3);
                D(1,2);
                D(3,2);
                D(3,1);
                D(2,1);];
    d_vs(1:9) = [D(1,1);
                D(2,2);
                D(3,3);
                Ds(2,3);
                Ds(1,3);
                Ds(1,2);
                Da(3,2);
                Da(3,1);
                Da(2,1);];
    E1(j,i) = 1/(d_vs'*(S1*d_vs)); 
    E2(j,i) = 1/(d_vs'*(S2*d_vs)); 
    E3(j,i) = 1/(d_vs'*(S3*d_vs)); 

    % shear vector n
    psi = pi;
%     nn = -[sin(theta(i))*sin(psi)-cos(phi(j))*cos(theta(i))*cos(psi);
%         -cos(theta(i))*sin(psi)-cos(phi(j))*sin(theta(i))*cos(psi);
%         sin(phi(j))*cos(psi)];
    nn = [cos(phi(j))*cos(theta(i)); cos(phi(j))*sin(theta(i));-sin(phi(j))];
    % The third direction:
    tt = cross(d, nn);
    T = tt*tt';
    T_v(1:9) = [T(1,1);
           T(2,2);
           T(3,3);
           T(2,3);
           T(1,3);
           T(1,2);
           T(3,2);
           T(3,1);
           T(2,1);];
    N = nn*nn';
    Ns = 0.5*(N+N');
    Na = 0.5*(N-N'); 
    n_v(1:9) = [N(1,1);
                N(2,2);
                N(3,3);
                N(2,3);
                N(1,3);
                N(1,2);
                N(3,2);
                N(3,1);
                N(2,1);];
    n_vs(1:9) = [N(1,1);
                N(2,2);
                N(3,3);
                Ns(2,3);
                Ns(1,3);
                Ns(1,2);
                Na(3,2);
                Na(3,1);
                Na(2,1);];
    niu1(j,i) = -E1(j,i)*d_vs'*(S1*n_vs);
    niu2(j,i) = -E2(j,i)*d_vs'*(S2*n_vs);    
    niu3(j,i) = -E3(j,i)*d_vs'*(S3*n_vs);

    M = sqrt(2)/2*(d*nn'+nn*d');
    m_v(1:9) = [M(1,1);
                M(2,2);
                M(3,3);
                M(2,3);
                M(1,3);
                M(1,2);
                M(3,2);
                M(3,1);
                M(2,1);];
    G1(j,i) = 1/(2* m_v'*(S1* m_v)); 
    G2(j,i) = 1/(2* m_v'*(S2* m_v));
    G3(j,i) = 1/(2* m_v'*(S3* m_v)); 
    K1(j,i) = 1/3/(d_v'*(S1*d_v)+d_v'*(S1*n_v)+d_v'*(S1*T_v)); 
    K2(j,i) = 1/3/(d_v'*(S2*d_v)+d_v'*(S2*n_v)+d_v'*(S2*T_v));
    K3(j,i) = 1/3/(d_v'*(S3*d_v)+d_v'*(S3*n_v)+d_v'*(S3*T_v));
%
%     t = cross(d,nn);
%     T = t*t';
%     t_v(1:9) = [T(1,1);
%                 T(2,2);
%                 T(3,3);
%                 T(2,3);
%                 T(1,3);
%                 T(1,2);
%                 T(3,2);
%                 T(3,1);
%                 T(2,1);];
%     d_ab(10:18) = [0,0,0,0, ,0,0,0,0]';
    d_ab(10:18) = [0;
                0;
                0;
                1;
                1;
                1;
                1;
                1;
                1];
    ab1(j,i) = E1(j,i)*d_v'*(S1*d_ab); 
    ab2(j,i) = E2(j,i)*d_v'*(S2*d_ab); 
    ab3(j,i) = E3(j,i)*d_v'*(S3*d_ab); 
    end
end
%     
%     nG = [sin(theta(i));-cos(theta(i))];
%     NG = nG*nG';
%     n_v = [NG(1,1);
%            NG(2,2);
%            NG(1,2);
%            NG(2,1);
%            0;
%            0;];
%     niu(i) = -E(i)*d_v'*(S*n_v);
%     ab(i) = E(i)*n_v'*(S*d_ab);

% transform to xy
for i = 1:n
    for j = 1:m
        x_E1(j,i) = E1(j,i)*sin(phi(j))* cos(theta(i));
        y_E1(j,i) = E1(j,i)*sin(phi(j))* sin(theta(i));
        z_E1(j,i) = E1(j,i)*cos(phi(j));
        x_E2(j,i) = E2(j,i)*sin(phi(j))* cos(theta(i));
        y_E2(j,i) = E2(j,i)*sin(phi(j))* sin(theta(i));
        z_E2(j,i) = E2(j,i)*cos(phi(j));
        x_E3(j,i) = E3(j,i)*sin(phi(j))* cos(theta(i));
        y_E3(j,i) = E3(j,i)*sin(phi(j))* sin(theta(i));
        z_E3(j,i) = E3(j,i)*cos(phi(j));

        x_K1(j,i) = K1(j,i)*sin(phi(j))* cos(theta(i));
        y_K1(j,i) = K1(j,i)*sin(phi(j))* sin(theta(i));
        z_K1(j,i) = K1(j,i)*cos(phi(j));
        x_K2(j,i) = K2(j,i)*sin(phi(j))* cos(theta(i));
        y_K2(j,i) = K2(j,i)*sin(phi(j))* sin(theta(i));
        z_K2(j,i) = K2(j,i)*cos(phi(j));
        x_K3(j,i) = K3(j,i)*sin(phi(j))* cos(theta(i));
        y_K3(j,i) = K3(j,i)*sin(phi(j))* sin(theta(i));
        z_K3(j,i) = K3(j,i)*cos(phi(j));

        x_niu1(j,i) = niu1(j,i)*sin(phi(j))* cos(theta(i));
        y_niu1(j,i) = niu1(j,i)*sin(phi(j))* sin(theta(i));
        z_niu1(j,i) = niu1(j,i)*cos(phi(j));
        x_niu2(j,i) = niu2(j,i)*sin(phi(j))* cos(theta(i));
        y_niu2(j,i) = niu2(j,i)*sin(phi(j))* sin(theta(i));
        z_niu2(j,i) = niu2(j,i)*cos(phi(j));
        x_niu3(j,i) = niu3(j,i)*sin(phi(j))* cos(theta(i));
        y_niu3(j,i) = niu3(j,i)*sin(phi(j))* sin(theta(i));
        z_niu3(j,i) = niu3(j,i)*cos(phi(j));        
    
        x_G1(j,i) = G1(j,i)*sin(phi(j))* cos(theta(i));
        y_G1(j,i) = G1(j,i)*sin(phi(j))* sin(theta(i));
        z_G1(j,i) = G1(j,i)*cos(phi(j));
        x_G2(j,i) = G2(j,i)*sin(phi(j))* cos(theta(i));
        y_G2(j,i) = G2(j,i)*sin(phi(j))* sin(theta(i));
        z_G2(j,i) = G2(j,i)*cos(phi(j));
        x_G3(j,i) = G3(j,i)*sin(phi(j))* cos(theta(i));
        y_G3(j,i) = G3(j,i)*sin(phi(j))* sin(theta(i));
        z_G3(j,i) = G3(j,i)*cos(phi(j));

        x_ab1(j,i) = ab1(j,i)*sin(phi(j))* cos(theta(i));
        y_ab1(j,i) = ab1(j,i)*sin(phi(j))* sin(theta(i));
        z_ab1(j,i) = ab1(j,i)*cos(phi(j));
        x_ab2(j,i) = ab2(j,i)*sin(phi(j))* cos(theta(i));
        y_ab2(j,i) = ab2(j,i)*sin(phi(j))* sin(theta(i));
        z_ab2(j,i) = ab2(j,i)*cos(phi(j));
        x_ab3(j,i) = ab3(j,i)*sin(phi(j))* cos(theta(i));
        y_ab3(j,i) = ab3(j,i)*sin(phi(j))* sin(theta(i));
        z_ab3(j,i) = ab3(j,i)*cos(phi(j));
    end
end
%     figure(1)
% plot(x_E,y_E)
% figure(2)
% plot(x_niu,y_niu)

%% Ploting 
% Young's modulus
figure(1)
s11=surf(x_E1/10e6,y_E1/10e6,z_E1/10e6,E1/10e6);
s11.FaceAlpha = 0.2;
hold on
plot3(x_E1(:,[125 375])/10e6,y_E1(:,[125 375])/10e6,z_E1(:,[125 375])/10e6,':','Color',[0.019608000,0.188235000,0.380392000],'LineWidth',1.5)
s12=surf(x_E2/10e6,y_E2/10e6,z_E2/10e6,E2/10e6);
plot3(x_E2(:,[125 375])/10e6,y_E2(:,[125 375])/10e6,z_E2(:,[125 375])/10e6,'--','Color',[0.403922000,0.000000000,0.121569000],'LineWidth',1.5)
s12.FaceAlpha = 0.4;
s13=surf(x_E3/10e6,y_E3/10e6,z_E3/10e6,E3/10e6);
plot3(x_E3(:,[125 375])/10e6,y_E3(:,[125 375])/10e6,z_E3(:,[125 375])/10e6,'-k','LineWidth',2.5)
s13.FaceAlpha = 0.8;
colormap(newmap);
    pax = gca;
    pax.FontWeight="bold";
    pax.FontSize = 16;
title('Youngs modulus E','FontSize',18,'FontWeight','bold')
xlabel('X_1','FontSize',16,'FontWeight','bold','FontName','Times New Roman','FontAngle','italic')
ylabel('X_2','FontSize',16,'FontWeight','bold','FontName','Times New Roman','FontAngle','italic')
zlabel('X_3','FontSize',16,'FontWeight','bold','FontName','Times New Roman','FontAngle','italic')
view([2,0.5,0.8])
box off;grid off;axis equal
shading flat
c = colorbar;
c.Label.String = 'MPa'; 
% c.Label.Position = [0 0 0];
c.Label.FontSize = 16;
c.Label.FontWeight = 'bold';
pax = gca;
pax.FontSize = 16;

% xlim([-60,60])
% ylim([-60,60])
% zlim([-60,60])


% Bulk modulus
figure(2)
s21 = surf(x_K1/10e6,y_K1/10e6,z_K1/10e6,K1/10e6);
s21.FaceAlpha = 0.2;
hold on
plot3(x_K1(:,[125 375])/10e6,y_K1(:,[125 375])/10e6,z_K1(:,[125 375])/10e6,':','Color',[0.019608000,0.188235000,0.380392000],'LineWidth',1.8)
s22 = surf(x_K2/10e6,y_K2/10e6,z_K2/10e6,K2/10e6);
plot3(x_K2(:,[125 375])/10e6,y_K2(:,[125 375])/10e6,z_K2(:,[125 375])/10e6,'--','Color',[0.403922000,0.000000000,0.121569000],'LineWidth',1.5)
s22.FaceAlpha = 0.4;
s23=surf(x_K3/10e6,y_K3/10e6,z_K3/10e6,K3/10e6);
plot3(x_K3(:,[125 375])/10e6,y_K3(:,[125 375])/10e6,z_K3(:,[125 375])/10e6,'-k','LineWidth',2)
s23.FaceAlpha = 0.8;
colormap(newmap);
    pax = gca;
    pax.FontWeight="bold";
    pax.FontSize = 16;
title('Bulk modulus B','FontSize',18,'FontWeight','bold')
xlabel('X_1','FontSize',16,'FontWeight','bold','FontName','Times New Roman','FontAngle','italic')
ylabel('X_2','FontSize',16,'FontWeight','bold','FontName','Times New Roman','FontAngle','italic')
zlabel('X_3','FontSize',16,'FontWeight','bold','FontName','Times New Roman','FontAngle','italic')
view([2,0.5,0.8])
box off;grid off;axis equal
shading flat
c = colorbar;
c.Label.String = 'MPa'; 
% c.Ticks =[18.9 19 19.1 19.2 19.3];
% c.TickLabels= [18.9 19 19.1 19.2  19.3];
c.Label.FontSize = 16;
c.Label.FontWeight = 'bold';
pax = gca;
pax.FontSize = 16;
% xlim([-110,110])
% ylim([-110,110])
% zlim([-110,110])

% Poisson's ratio
figure(3)
s31 = surf(x_niu1,y_niu1,z_niu1,niu1);
s31.FaceAlpha = 0.2;
hold on
plot3(x_niu1(:,[125 375]),y_niu1(:,[125 375]),z_niu1(:,[125 375]),':','Color',[0.019608000,0.188235000,0.380392000],'LineWidth',1.5)
s32 = surf(x_niu2,y_niu2,z_niu2,niu2);
plot3(x_niu2(:,[125 375]),y_niu2(:,[125 375]),z_niu2(:,[125 375]),'--','Color',[0.403922000,0.000000000,0.121569000],'LineWidth',1.5)
s32.FaceAlpha = 0.4;
s33 = surf(x_niu3,y_niu3,z_niu3,niu3);
plot3(x_niu3(:,[125 375]),y_niu3(:,[125 375]),z_niu3(:,[125 375]),'-k','LineWidth',2)
s33.FaceAlpha = 0.8;
colormap(newmap);
    pax = gca;
    pax.FontWeight="bold";
    pax.FontSize = 16;
title('Poissons ratio \nu','FontSize',18,'FontWeight','bold')
xlabel('X_1','FontSize',16,'FontWeight','bold','FontName','Times New Roman','FontAngle','italic')
ylabel('X_2','FontSize',16,'FontWeight','bold','FontName','Times New Roman','FontAngle','italic')
zlabel('X_3','FontSize',16,'FontWeight','bold','FontName','Times New Roman','FontAngle','italic')
view([2,0.5,0.8])
box off;grid off;axis equal
shading flat
c = colorbar;
% c.Label.String = 'MPa'; 
% c.Ticks= [0.2 0.4 0.6 0.8];
% c.TickLabels= [0.2 0.4 0.6 0.8];
c.Label.FontSize = 16;
c.Label.FontWeight = 'bold';
pax = gca;
pax.FontSize = 16;

% Shear modulus
figure(4)
s41 = surf(x_G1/10e6,y_G1/10e6,z_G1/10e6,G1/10e6);
s41.FaceAlpha = 0.2;
hold on
plot3(x_G1(:,[125 375])/10e6,y_G1(:,[125 375])/10e6,z_G1(:,[125 375])/10e6,':','Color',[0.019608000,0.188235000,0.380392000],'LineWidth',1.5)
s42 = surf(x_G2/10e6,y_G2/10e6,z_G2/10e6,G2/10e6);
plot3(x_G2(:,[125 375])/10e6,y_G2(:,[125 375])/10e6,z_G2(:,[125 375])/10e6,'--','Color',[0.403922000,0.000000000,0.121569000],'LineWidth',1.5)
s42.FaceAlpha = 0.4;
s43 = surf(x_G3/10e6,y_G3/10e6,z_G3/10e6,G3/10e6);
plot3(x_G3(:,[125 375])/10e6,y_G3(:,[125 375])/10e6,z_G3(:,[125 375])/10e6,'-k','LineWidth',2)
s43.FaceAlpha = 0.8;
colormap(newmap);
    pax = gca;
    pax.FontWeight="bold";
    pax.FontSize = 16;
title('Shear modulus G','FontSize',18,'FontWeight','bold')
xlabel('X_1','FontSize',16,'FontWeight','bold','FontName','Times New Roman','FontAngle','italic')
ylabel('X_2','FontSize',16,'FontWeight','bold','FontName','Times New Roman','FontAngle','italic')
zlabel('X_3','FontSize',16,'FontWeight','bold','FontName','Times New Roman','FontAngle','italic')
view([2,0.5,0.8])
box off;grid off;axis equal
shading flat
c = colorbar;
c.Label.String = 'MPa'; 
% c.Label.Position = [0 0 0];
c.Label.FontSize = 16;
c.Label.FontWeight = 'bold';
pax = gca;
pax.FontSize = 16;
% xlim([-45,45])
% ylim([-45,45])
% zlim([-45,45])

% Axial bending 
figure(5)
s51 = surf(x_ab1,y_ab1,z_ab1,ab1);
s51.FaceAlpha = 0.2;
hold on
plot3(x_ab1(:,[125 375]),y_ab1(:,[125 375]),z_ab1(:,[125 375]),':','Color',[0.019608000,0.188235000,0.380392000],'LineWidth',1.8)
s52 = surf(x_ab2,y_ab2,z_ab2,ab2);
plot3(x_ab2(:,[125 375]),y_ab2(:,[125 375]),z_ab2(:,[125 375]),'--','Color',[0.403922000,0.000000000,0.121569000],'LineWidth',1.5)
s52.FaceAlpha = 0.35;
s53 = surf(x_ab3,y_ab3,z_ab3,ab3);
plot3(x_ab3(:,[125 375]),y_ab3(:,[125 375]),z_ab3(:,[125 375]),'-k','LineWidth',1.2)
s53.FaceAlpha = 0.5;
colormap(newmap);
    pax = gca;
    pax.FontWeight="bold";
    pax.FontSize = 16;
title('Axial bending coupling coeffcient \zeta','FontSize',18,'FontWeight','bold')
xlabel('X_1','FontSize',16,'FontWeight','bold','FontName','Times New Roman','FontAngle','italic')
ylabel('X_2','FontSize',16,'FontWeight','bold','FontName','Times New Roman','FontAngle','italic')
zlabel('X_3','FontSize',16,'FontWeight','bold','FontName','Times New Roman','FontAngle','italic')
view([2,0.5,0.8])
box off;grid off;
shading flat

c = colorbar;
c.Label.Position = [2 0 0];
c.Label.FontSize = 16;
c.Label.FontWeight = 'bold';

pax = gca;
pax.FontSize = 16;
pax.XTick=[-300 0 300];axis equal
pax.XTickLabel={'-300','0','300'};
pax.YTick=[-400 0 400];
pax.YTickLabel={'-400','0','400'};
% xlim([-540,540])
% ylim([-540,540])
% zlim([-540,540])

