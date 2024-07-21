function [K] = K_el(Iz,Iy,A,L)
E=9.3e6;G=E/2/(1+0.3); 
J = Iz+Iy;
K = [A*E/L  0           0           0       0           0          -A*E/L   0           0            0       0          0;
         0      12*E*Iz/L^3 0           0       0           6*E*Iz/L^2  0      -12*E*Iz/L^3  0            0       0          6*E*Iz/L^2;
         0      0           12*E*Iy/L^3 0       -6*E*Iy/L^2 0           0       0          -12*E*Iy/L^3  0      -6*E*Iy/L^2 0;
         0      0           0           G*J/L   0           0           0       0           0           -G*J/L   0          0;
         0      0          -6*E*Iy/L^2  0       4*E*Iy/L    0           0       0           6*E*Iy/L^2   0       2*E*Iy/L   0;
         0      6*E*Iz/L^2  0           0       0           4*E*Iz/L    0      -6*E*Iz/L^2  0            0       0          2*E*Iz/L;
        -A*E/L  0           0           0       0           0           A*E/L   0           0            0       0          0;
         0     -12*E*Iz/L^3 0           0       0          -6*E*Iz/L^2  0       12*E*Iz/L^3 0            0       0         -6*E*Iz/L^2;
         0      0          -12*E*Iy/L^3 0       6*E*Iy/L^2  0           0       0           12*E*Iy/L^3  0       6*E*Iy/L^2 0;
         0      0           0          -G*J/L   0           0           0       0           0            G*J/L   0          0; 
         0      0          -6*E*Iy/L^2  0       2*E*Iy/L    0           0       0           6*E*Iy/L^2   0       4*E*Iy/L   0;
         0      6*E*Iz/L^2  0           0       0           2*E*Iz/L    0      -6*E*Iz/L^2  0            0       0          4*E*Iz/L;] ;
end

