%% Main
clear all;

addpath('./Functions/');

% Type de système
feedback = true; % true => input = stabilizing feedback law
v_input = false; % true => input += c*sin(100*lambda*t*r)
if feedback && ~v_input
    c = 0;
else
    c = 1000;
end
% c = 1000;%1000 or 0 % input = c*sin(100*lambda*t*r)
noz2 = true; % true => n2 = 0 (full sate measurment)
perts = [linspace(0, 100, 10)];%2 % Input perturbation
movie = false;

d = 0.1; % Delay
dperts = d;%linspace(0.5, 1, 10);%0.099; % = d; Delay perturbation

n_pert = size(perts, 2);
n_d = size(dperts, 2);
limsup = zeros(1, max(n_pert,n_d));
limsup_mean = zeros(1, max(n_pert,n_d));

for i = 1:n_d
dpert = dperts(i);

for k = 1:n_pert
inputpert = perts(k);

% Discr??tisation 
Nx = 20;
Omega = linspace(0, 1, Nx);

if Nx>1
    dx = Omega(2) - Omega(1);
else
    dx=1;
end

% Param??tres syst??me
sigma = 60;
W11 = gauss(Omega, 2, sigma, dx, Nx);
W12 = gauss(Omega, 2, sigma, dx, Nx);
W21 = gauss(Omega, -2, sigma, dx, Nx);
W22 = gauss(Omega, 0.1, sigma, dx, Nx);

% Param??tres observateur
alpha = 100;%100
beta1 = 100;%100
beta2 = 100;%100

% Temps et CI
tmax = 10;%10
% z0 = [ones(2*Nx, 1); zeros(2*Nx+2*Nx^2, 1)];
z0 = [ones(3*Nx, 1); zeros(Nx+2*Nx^2, 1)];

if noz2
    W12 = 0;
    W21 = 0;
    W22 = 0;
    z0((3*Nx+1):4*Nx) = z0((Nx+1):2*Nx);
    z0((4*Nx+Nx^2+1):(4*Nx+2*Nx^2)) = zeros(Nx^2, 1);
end



% W22hat0 = exp(-sigma*min(abs(Omega+Omega'-1), Omega(end)-Omega(1)-abs(Omega+Omega'-1)).^2)*dx;
% W22hat0 = 2/norm(W22hat0)*W22hat0;
% z0((4*Nx+Nx^2+1):(4*Nx+2*Nx^2)) = reshape(W22hat0, [1, Nx^2]);

% R??solution syst??me et observateur
f = @(t, z, Z) Observer(z, t, Z, Omega, Nx, dx, W11, W12, W21, W22, alpha, beta1, beta2, c, feedback, noz2, inputpert)';

if d~=0 || dpert ~= 0
    options=ddeset('RelTol',1e-2);
    if d ~= dpert
        sol = dde23(f,[d, dpert],z0,[0 tmax],options);
    else
        sol = dde23(f,d,z0,[0 tmax],options);
    end
    T = sol.x';
    z = sol.y';
else
    [T, z] = ode45(@(t, z) f(t, z, z), [0 tmax], z0);
end
Nt = size(T, 1);

% g1 = cont1(T, c, Nx, dx, Omega, z')';
% g2 = cont2(T, c, Nx, dx, Omega, z')';

%% Plot
close all;
% 
% if Nx>1
% 
% figure
% surf(Omega, T, z(:, 1:Nx))
% title('z1')
% figure
% surf(Omega, T, z(:, (2*Nx+1):3*Nx))
% title('z1hat')
% figure
% surf(Omega, T, z(:, (2*Nx+1):3*Nx) - z(:, 1:Nx))
% title('z1tilde')
% 
% figure
% surf(Omega, T, z(:, (Nx+1):2*Nx))
% title('z2')
% figure
% surf(Omega, T, z(:, (3*Nx+1):4*Nx))
% title('z2hat')
% figure
% surf(Omega, T, z(:, (3*Nx+1):4*Nx) - z(:, (Nx+1):2*Nx))
% title('z2tilde')
% 
% else
% 
% figure
% plot(T, z(:, 1))
% hold on
% plot(T, z(:, 3))
% title('z1, z1hat')
% 
% figure
% plot(T, z(:, 3) - z(:, 1))
% title('z1tilde')
% 
% figure
% plot(T, z(:, 2))
% hold on
% plot(T, z(:, 4))
% title('z2, z2hat')
% 
% figure
% plot(T, z(:, 4) - z(:, 2))
% title('z2tilde')
% 
% end
% 
% figure
% plot(T, sqrt(sum(z(:, (2*Nx+1):3*Nx) - z(:, 1:Nx), 2).^2))
% hold on
% plot(T, sqrt(sum(z(:, (3*Nx+1):4*Nx) - z(:, Nx+1:2*Nx), 2).^2))
% title('norm(ztilde)')
% 
% figure
% plot(T, z(:, (4*Nx+1):(4*Nx+Nx^2)) - reshape(W21, [1, Nx^2]))
% title('W21tilde')
% figure
% plot(T, z(:, (4*Nx+Nx^2+1):(4*Nx+2*Nx^2)) - reshape(W22, [1, Nx^2]))
% title('W22tilde')

ErrW11 = sqrt(sum((z(:, 4*Nx+1:4*Nx+Nx^2) - reshape(W11, [1, Nx^2])).^2, 2));

if not(noz2)
    ErrW12 = sqrt(sum((z(:, 4*Nx+Nx^2+1:4*Nx+2*Nx^2) - reshape(W12, [1, Nx^2])).^2, 2));
end

% define figure properties
opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'img/';
opts.width      = 8;
opts.height     = 6;
opts.fontType   = 'Times';
opts.fontSize   = 9;

% create new figure
lw = 2.5;
fig = figure;
plot(T, ErrW11, '-', 'Linewidth',lw)
if not(noz2)
    hold on
    plot(T, ErrW12, '--', 'Linewidth',lw)
end
hold on
plot(T, sqrt(sum((z(:, (2*Nx+1):3*Nx) - z(:, 1:Nx)).^2, 2)), '-.', 'Linewidth',lw)
if not(noz2)
    hold on
    plot(T, sqrt(sum((z(:, (3*Nx+1):4*Nx) - z(:, Nx+1:2*Nx)).^2, 2)), ':', 'Linewidth',lw)
end
% title('Evolution of the error between the state and the observer based on $T^d_\lambda$', 'interpreter', 'latex')
xlabel('Time $t$', 'interpreter', 'latex', 'FontSize', 12)
ylabel('Estimation error', 'interpreter', 'latex', 'FontSize', 12)
if noz2
    legend({'$\|\tilde w(t)\|_{L^2}$', '$\|\tilde z(t)\|_{L^2}$'}, 'interpreter', 'latex', 'FontSize', 14)
else
    legend({'$\|\tilde w_{11}(t)\|_{L^2}$', '$\|\tilde w_{12}(t)\|_{L^2}$', '$\|\tilde z_1(t)\|_{L^2}$' '$\|\tilde z_2(t)\|_{L^2}$'}, 'interpreter', 'latex', 'FontSize', 14)
end
ylim([0, 15])

% scaling
fig.Units               = 'centimeters';
fig.Position(3)         = opts.width;
fig.Position(4)         = opts.height;

% set text properties
set(fig.Children, ...
    'FontName',     'Times', ...
    'FontSize',     9);

% remove unnecessary white space
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))

print('Images/Observer_error', '-depsc')

if feedback

    % create new figure
    fig2 = figure;
    normz = sqrt(sum(z(:, 1:Nx).^2, 2));
    plot(T, normz, '-', 'Linewidth',lw)
    if not(noz2)
        hold on
        plot(T, sqrt(sum(z(:, Nx+1:2*Nx).^2, 2)), '--', 'Linewidth',lw)
    end
    % title('Evolution of the error between the state and the observer based on $T^d_\lambda$', 'interpreter', 'latex')
    xlabel('Time $t$', 'interpreter', 'latex', 'FontSize', 12)
    ylabel('Norm of the state', 'interpreter', 'latex', 'FontSize', 12)
    if not(noz2)
        legend({'$\|z_{1}(t)\|_{L^2}$', '$\|z_{2}(t)\|_{L^2}$'}, 'interpreter', 'latex', 'FontSize', 14)
    else
        legend({'$\|z(t)\|_{L^2}$'}, 'interpreter', 'latex', 'FontSize', 14)
    end
    % ylim([0, 15])

    limsup(max(k,i)) = max(normz(floor(Nt/2):end));
    limsup_mean(max(k,i)) = mean(normz(floor(Nt/2):end));
    
    % scaling
    fig2.Units               = 'centimeters';
    fig2.Position(3)         = opts.width;
    fig2.Position(4)         = opts.height;
    
    % set text properties
    set(fig2.Children, ...
        'FontName',     'Times', ...
        'FontSize',     9);
    
    % remove unnecessary white space
    set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
    
    print('Images/State', '-depsc')
    
    % create new figure
    U = zeros(Nt, Nx);
    dt = T(2) - T(1);
    jr = max(1, floor(dpert/dt));
    for j = 1:Nt
        U(j, :) = c*sin(100*T(j)*(dx+Omega')) + cont(z(j, 1:Nx)', z(max(j-jr, j), 1:Nx)', z(max(j-jr, j), (Nx+1):2*Nx)', reshape(z(j, (4*Nx+1):(4*Nx+Nx^2))', [Nx, Nx]), (not(noz2))*reshape(z(j, (4*Nx+Nx^2+1):(4*Nx+2*Nx^2))', [Nx, Nx]), alpha);
%         U(i, :) = - reshape(z(i, (4*Nx+1):(4*Nx+Nx^2)), [Nx, Nx])*S(z(i, 1:Nx)');
    end
    
    fig3 = figure;
    % [X, Y] = meshgrid(Omega, T);
    % surf(X, Y, U, 'Linewidth',lw)
%     semilogy(T, sqrt(sum(U(:, 1:Nx), 2).^2), '-', 'Linewidth',lw)
    plot(T, sqrt(sum(U(:, 1:Nx).^2, 2)), '-', 'Linewidth',lw)
    xlabel('Time $t$', 'interpreter', 'latex', 'FontSize', 12)
    ylabel('Input norm', 'interpreter', 'latex', 'FontSize', 12)
    legend({'$\|u_1(t)\|_{L^2}$'}, 'interpreter', 'latex', 'FontSize', 14)
    
    % ylim([0, 15])
    
    % scaling
    fig3.Units               = 'centimeters';
    fig3.Position(3)         = opts.width;
    fig3.Position(4)         = opts.height;
    
    % set text properties
    set(fig3.Children, ...
        'FontName',     'Times', ...
        'FontSize',     9);
    
    % remove unnecessary white space
    set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
    
    print('Images/Input', '-depsc')

end

% figure
% surf(Omega, Omega, W22)
% figure
% surf(Omega, Omega, reshape(z(1, (4*Nx+Nx^2+1):(4*Nx+2*Nx^2)), [Nx, Nx]))
% figure
% surf(Omega, Omega, reshape(z(10000, (4*Nx+Nx^2+1):(4*Nx+2*Nx^2)), [Nx, Nx]))
% figure
% surf(Omega, Omega, reshape(z(100000, (4*Nx+Nx^2+1):(4*Nx+2*Nx^2)), [Nx, Nx]))
% 
% figure
% imagesc(W22)
% figure
% imagesc(reshape(z(1, (4*Nx+Nx^2+1):(4*Nx+2*Nx^2)), [Nx, Nx]))
% figure
% imagesc(reshape(z(10000, (4*Nx+Nx^2+1):(4*Nx+2*Nx^2)), [Nx, Nx]))
% figure
% imagesc(reshape(z(100000, (4*Nx+Nx^2+1):(4*Nx+2*Nx^2)), [Nx, Nx]))
% 
% drawnow

if movie

    nb_img = 100;
    pas = floor(Nt/nb_img);
    cpt_img = 1;
    fig4 = figure;
    for i = 1:pas:Nt
    %     imagesc(reshape(z(i, (4*Nx+Nx^2+1):(4*Nx+2*Nx^2)), [Nx, Nx]))
        fig4;
        surf(Omega, Omega, reshape(z(i, (4*Nx+1):(4*Nx+Nx^2)), [Nx, Nx]))
        xlabel('$r$', 'interpreter', 'latex', 'FontSize', 16)
        ylabel('$r''$', 'interpreter', 'latex', 'FontSize', 16)
        if not(noz2)
            zlabel('$\hat{w}_{11}(t, r, r'')$', 'interpreter', 'latex', 'FontSize', 16)
        else
            zlabel('$\hat{w}(t, r, r'')$', 'interpreter', 'latex', 'FontSize', 16)
        end
        axis([0 1 0 1 -0.5 1])
        drawnow
        print(['Images/Fig_' num2str(cpt_img) '_t' num2str(floor(T(i))) ',' num2str(floor(100*(T(i) - floor(T(i)))))], '-dpng')
        cpt_img = cpt_img+1;
    end

end




end
end

%% Plot final
% create new figure
    fig5 = figure;
    plot(perts, limsup, '*', 'Linewidth',lw)
    hold on
    plot(perts, limsup_mean, 'o', 'Linewidth',lw)
    % title('Evolution of the error between the state and the observer based on $T^d_\lambda$', 'interpreter', 'latex')
    xlabel('Perturbation amplitude', 'interpreter', 'latex', 'FontSize', 12)
    ylabel('Steady state value', 'interpreter', 'latex', 'FontSize', 12)
    legend({'$\max_{5\leq t\leq 10} \|z(t)\|_{L^2}$', '$\frac15\int_{t=5}^{t=10} \|z(t)\|_{L^2} \mathrm{d}t$'}, 'interpreter', 'latex', 'FontSize', 14)
    % ylim([0, 15])
    
    % scaling
    fig5.Units               = 'centimeters';
    fig5.Position(3)         = opts.width;
    fig5.Position(4)         = opts.height;
    
    % set text properties
    set(fig5.Children, ...
        'FontName',     'Times', ...
        'FontSize',     9);
    
    % remove unnecessary white space
    set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
    
    print('Images/Perts', '-depsc')
