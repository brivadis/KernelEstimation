function W = gauss(Omega, mu, sigma, dx, Nx)
    % Noyau gaussien sur le tore
W = exp(-sigma*min(abs(Omega-Omega'), Omega(end)-Omega(1)-abs(Omega-Omega')).^2);
W = mu/norm(W)*W;
end