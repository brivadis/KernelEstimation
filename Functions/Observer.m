function zp = Observer(z, t, Omega, Nx, dx, W11, W12, W21, W22, alpha, beta1, beta2, c)
% Champ de vecteur du système et de l'observateur
    g1 = cont1(t, c, Nx, dx, Omega, z);
    g2 = cont2(t, c, Nx, dx, Omega, z);
    zp(1:Nx) = - z(1:Nx) + c*sin(t*(dx+Omega')) + W11*S(z(1:Nx)) + W12*S(z((Nx+1):2*Nx));
    zp((Nx+1):2*Nx) = - z((Nx+1):2*Nx) + c*sin(sqrt(2)*t*(dx+Omega')) + W21*S(z(1:Nx))+ W22*S(z((Nx+1):2*Nx));
    zp((2*Nx+1):3*Nx) = - z((2*Nx+1):3*Nx) + c*sin(t*(dx+Omega')) + W11*S(z((2*Nx+1):3*Nx)) + W12*S(z((Nx+1):2*Nx));
    zp((3*Nx+1):4*Nx) = - z((Nx+1):2*Nx) + c*sin(sqrt(2)*t*(dx+Omega')) + reshape(z((4*Nx+1):(4*Nx+Nx^2)), [Nx, Nx])*g1 + reshape(z((4*Nx+Nx^2+1):(4*Nx+2*Nx^2)), [Nx, Nx])*g2 - alpha*(z((3*Nx+1):4*Nx)-z((Nx+1):2*Nx));
    W21p = -beta1*(z((3*Nx+1):4*Nx)-z((Nx+1):2*Nx))*g1';
    W22p = -beta2*(z((3*Nx+1):4*Nx)-z((Nx+1):2*Nx))*g2';
    zp((4*Nx+1):(4*Nx+Nx^2)) = reshape(W21p, [1, Nx^2]);
    zp((4*Nx+Nx^2+1):(4*Nx+2*Nx^2)) = reshape(W22p, [1, Nx^2]);
end