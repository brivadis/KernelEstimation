function zp = Observer(z, t, Z, Omega, Nx, dx, W11, W12, W21, W22, alpha, beta1, beta2, c, feedback, noz2)
% Champ de vecteur du système et de l'observateur


    zr = Z(:,1);

    z1 = z(1:Nx);
    z2 = z((Nx+1):2*Nx);
    zhat1 = z((2*Nx+1):3*Nx);
    zhat2 = z((3*Nx+1):4*Nx);
    What11 = reshape(z((4*Nx+1):(4*Nx+Nx^2)), [Nx, Nx]);
    What12 = reshape(z((4*Nx+Nx^2+1):(4*Nx+2*Nx^2)), [Nx, Nx]);
    zr1 = zr(1:Nx);
    zr2 = zr((Nx+1):2*Nx);
    zhatr2 = z((3*Nx+1):4*Nx);

    u1 = c*sin(100*t*(dx+Omega'));%100*
    u2 = c*sin(sqrt(2)*100*t*(dx+Omega'));%100*

    if noz2
        W12 = 0;
        What12 = 0;
    end

    if feedback
        u = cont(z1, zr1, zr2, What11, What12, alpha);
    else
        u = 0;
    end

    z1p = -z1 + u1 + W11*S(zr1) + W12*S(zr2) + u;
    z2p = -z2 + u2 + W21*S(zr1) + W22*S(zr2);
    zhat1p = -alpha*(zhat1-z1) - z1 + u1 + What11*S(zr1) + What12*S(zhatr2) + u;
    zhat2p = -zhat2 + u2 + W21*S(zr1) + W22*S(zhatr2);
    What11p = -beta1*(zhat1-z1)*S(zr1)';
    What12p = -beta2*(zhat1-z1)*S(zhatr2)';

    zp = [z1p', z2p', zhat1p', zhat2p', reshape(What11p, [1, Nx^2]), reshape(What12p, [1, Nx^2])];
end