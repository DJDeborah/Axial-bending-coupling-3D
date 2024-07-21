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

% dimapb = [	0.819608000	0.898039000	0.941176000
% 	0.572549000	0.772549000	0.870588000
% 	0.262745000	0.576471000	0.764706000
% 	0.129412000	0.400000000	0.674510000
% 	      
% 	0.019608000	0.188235000	0.380392000];
% N = 30;
% newmap1 = interp1(1:size(dimapb,1),dimapb,linspace(1,size(dimapb,1),N),'linear');

%% Q matrix 3d
% Vertical curved beam
C1 = [2382327204	0	0	0	0.006187112	0	0	0.005418094	0	0	0	0	0	0	0	0	0	690212.3196;
0	2382327204	0	0.006264163	0	0	0.005090521	0	0	0	0	0	0	0	-690212.3196	0	0	0;
0	0	2382327204	-0.005552363	0	0	-0.005802323	0	0	0	0	0	0	-690212.3196	0	0	0	0;
0	0.006264163	-0.005552363	19915343.99	0	0	8655616.722	0	0	0	0	0	0	2508.560605	-2420.403701	0	5018.524333	2420.40374;
0.006187112	0	0	0	51392223.94	108998.1227	0	15668799.11	-108998.1218	48.3452145	48.3452154	22719.27849	-17385.4698	0	0	-9678.664411	0	0;
0	0	0	0	108998.1227	73140198.85	0	-27249.20269	40067979.93	-17644.67529	32313.58431	48.25340249	-16726.93657	0	0	-742.9835342	0	0;
0	0.005090521	-0.005802323	8655616.722	0	0	19915343.99	0	0	0	0	0	0	-2508.560683	2420.403782	0	-5018.524336	-2420.40374;
0.005418094	0	0	0	15668799.11	-27249.20269	0	24599548.19	27249.20421	-12.08615919	-12.08615939	6904.54024	4346.315235	0	0	2419.637071	0	0;
0	0	0	0	-108998.1218	40067979.93	0	27249.20421	73140198.85	-32313.58431	17644.67529	-48.25340103	16726.93657	0	0	742.9836333	0	0;
0	0	0	0	48.3452145	-17644.67529	0	-12.08615919	-32313.58431	466.6693697	-7.770040916	0.021402397	-7.419094299	0	0	-0.32954422	0	0;
0	0	0	0	48.3452154	32313.58431	0	-12.08615939	17644.67529	-7.770040916	466.6693697	0.021402397	-7.419094293	0	0	-0.329544222	0	0;
0	0	0	0	22719.27849	48.25340249	0	6904.54024	-48.25340103	0.021402397	0.021402397	462.4367326	-7.696536784	0	0	-4.284738801	0	0;
0	0	0	0	-17385.4698	-16726.93657	0	4346.315235	16726.93657	-7.419094299	-7.419094293	-7.696536784	3604.235712	0	0	118.5077172	0	0;
0	0	-690212.3196	2508.560605	0	0	-2508.560683	0	0	0	0	0	0	1134.95131	0.713373496	0	30.23297655	-0.713373496;
0	-690212.3196	0	-2420.403701	0	0	2420.403782	0	0	0	0	0	0	0.713373496	1134.724686	0	1.427145982	8.969334593;
0	0	0	0	-9678.664411	-742.9835342	0	2419.637071	742.9836333	-0.32954422	-0.329544222	-4.284738801	118.5077172	0	0	3614.73706	0	0;
0	0	0	5018.524333	0	0	-5018.524336	0	0	0	0	0	0	30.23297655	1.427145982	0	3609.245494	-1.427145983;
690212.3196	0	0	2420.40374	0	0	-2420.40374	0	0	0	0	0	0	-0.713373496	8.969334593	0	-1.427145983	1134.724686]; 

% 3 Curved beams
C2=[1614870649	0	0	0	0.003580151	0	0	0.003937209	0	0	0	0	0	0	0	0	0	927215.344;
0	1614870649	0	0.003536619	0	0	0.003913963	0	0	0	0	0	0	0	-927215.344	0	0	0;
0	0	1614870649	-0.003907189	0	0	-0.003543391	0	0	0	0	0	0	-927215.344	0	0	0	0;
0	0.003536619	-0.003907189	20271832.32	0	-0.000659794	8818514.981	0	0	0	0	0	0	4907.081411	-4659.408577	0	9337.140612	4659.40857;
0.003580151	0	0	0	51490103.85	407698.8764	0	15995196.45	-407698.8736	362.8082451	362.8082492	45215.16471	-32588.28145	0	0	-18566.65452	0	0;
0	0	0	-0.000659794	407698.8764	72404514.24	0.000659794	-101945.796	39664275.05	-34291.42714	63426.7266	360.064925	-31046.3649	0	0	-1899.752435	0	0;
0	0.003913963	-0.003543391	8818514.981	0	0.000659794	20271832.32	0	0	0	0	0	0	-4907.081396	4659.408562	0	-9337.140611	-4659.40857;
0.003937209	0	0	0	15995196.45	-101945.796	0	24870758.46	101945.7988	-90.72082165	-90.7208226	13867.34442	8148.755175	0	0	4642.623551	0	0;
0	0	0	0	-407698.8736	39664275.05	0	101945.7988	72404514.24	-63426.7266	34291.42714	-360.0649184	31046.3649	0	0	1899.752508	0	0;
0	0	0	0	362.8082451	-34291.42714	0	-90.72082165	-63426.7266	515.1139927	-29.6389281	0.32041914	-27.62793312	0	0	-1.690575832	0	0;
0	0	0	0	362.8082492	63426.7266	0	-90.7208226	34291.42714	-29.6389281	515.1139927	0.320419144	-27.6279331	0	0	-1.690575845	0	0;
0	0	0	0	45215.16471	360.064925	0	13867.34442	-360.0649184	0.32041914	0.320419144	499.2543098	-28.78079316	0	0	-16.39739874	0	0;
0	0	0	0	-32588.28145	-31046.3649	0	8148.755175	31046.3649	-27.62793312	-27.6279331	-28.78079316	3388.628225	0	0	151.851457	0	0;
0	0	-927215.344	4907.081411	0	0	-4907.081396	0	0	0	0	0	0	1452.783028	2.775502796	0	38.99832624	-2.775502796;
0	-927215.344	0	-4659.408577	0	0	4659.408562	0	0	0	0	0	0	2.775502796	1452.49697	0	5.281196245	14.77585599;
0	0	0	0	-18566.65452	-1899.752435	0	4642.623551	1899.752508	-1.690575832	-1.690575845	-16.39739874	151.851457	0	0	3426.948976	0	0;
0	0	0	9337.140612	0	0	-9337.140611	0	0	0	0	0	0	38.99832624	5.281196245	0	3414.639613	-5.281196246;
927215.344	0	0	4659.40857	0	0	-4659.40857	0	0	0	0	0	0	-2.775502796	14.77585599	0	-5.281196246	1452.49697;];

C3=[1063340071	0	0	0	-0.000334375	0	0	-0.000365457	0	0	0	0	0	0	0	0	0	901921.9939;
0	1063340071	0	-0.000336631	0	0	-0.000438227	0	0	0	0	0	0	0	-901921.9939	0	0	0;
0	0	1063340071	0.000347845	0	0	0.000427015	0	0	0	0	0	0	-901921.9939	0	0	0	0;
0	-0.000336631	0.000347845	20854491.05	0	0	9123589.307	0	-0.000123997	0	0	0	0	7083.241335	-6561.583622	0	12463.57877	6561.583621;
-0.000334375	0	0	0	51628518.39	820795.3741	0	16546846.49	-820795.3748	1101.286564	1101.286564	67212.17073	-43966.94756	0	0	-25851.7024	0	0;
0	0	0	0	820795.3741	71218415.31	0	-205487.3895	39069475.11	-49089.19736	92224.42569	1082.699617	-41277.63868	0	0	-3547.343257	0	0;
0	-0.000438227	0.000427015	9123589.307	0	0	20854491.05	0	0.000123996	0	0	0	0	-7083.241332	6561.58362	0	-12463.57877	-6561.583621;
-0.000365457	0	0	0	16546846.49	-205487.3895	0	25329597.23	205487.3894	-275.7087918	-275.7087915	20936.43144	11007.19322	0	0	6472.013623	0	0;
0	0	0	-0.000123997	-820795.3748	39069475.11	0.000123996	205487.3894	71218415.31	-92224.42569	49089.19736	-1082.699618	41277.63868	0	0	3547.343262	0	0;
0	0	0	0	1101.286564	-49089.19736	0	-275.7087918	-92224.42569	590.8461353	-61.59575755	1.452691595	-55.38348569	0	0	-4.759580271	0	0;
0	0	0	0	1101.286564	92224.42569	0	-275.7087915	49089.19736	-61.59575755	590.8461353	1.452691594	-55.38348569	0	0	-4.759580266	0	0;
0	0	0	0	67212.17073	1082.699617	0	20936.43144	-1082.699618	1.452691595	1.452691594	558.8924035	-57.99618126	0	0	-34.10061651	0	0;
0	0	0	0	-43966.94756	-41277.63868	0	11007.19322	41277.63868	-55.38348569	-55.38348569	-57.99618126	3075.315114	0	0	190.0179508	0	0;
0	0	-901921.9939	7083.241335	0	0	-7083.241332	0	0	0	0	0	0	1661.701306	5.955468793	0	51.31576262	-5.955468793;
0	-901921.9939	0	-6561.583622	0	0	6561.58362	0	0	0	0	0	0	5.955468793	1661.396297	0	10.47916496	23.75728928;
0	0	0	0	-25851.7024	-3547.343257	0	6472.013623	3547.343262	-4.759580271	-4.759580266	-34.10061651	190.0179508	0	0	3149.834584	0	0;
0	0	0	12463.57877	0	0	-12463.57877	0	0	0	0	0	0	51.31576262	10.47916496	0	3128.402302	-10.47916496;
901921.9939	0	0	6561.583621	0	0	-6561.583621	0	0	0	0	0	0	-5.955468793	23.75728928	0	-10.47916496	1661.396297;
];
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
colormap(newmap);
s13.FaceAlpha = 0.8;
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
s33.FaceAlpha = 0.7;
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
c.Label.String = 'MPa'; 
c.Label.Position = [2 0 0];
c.Label.FontSize = 16;
c.Label.FontWeight = 'bold';

pax = gca;
pax.FontSize = 16;
pax.XTick=[-300 0 300];
pax.XTickLabel={'-300','0','300'};
pax.YTick=[-400 0 400];
pax.YTickLabel={'-400','0','400'};
% xlim([-540,540])
% ylim([-540,540])
% zlim([-540,540])

