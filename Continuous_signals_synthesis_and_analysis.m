clear ; set(0,'DefaultFigureWindowStyle','docked')
%% Data and analysis

data = 'Green function_Thin Isotropic Homegeneous';

switch data
    case 'plaque Guilherme'
        load('C:\Users\Francois_Fabre\Documents\Thèse ANR\ManipFrancoisFabre\plaque\plaque_guilherme.mat')

        S = reshape(frf, [], 9, 5) ;
        S = permute(S, [2 3 1]) ;

        N_freq = 1000 ;
        w = w(1:round(size(S,3)/(N_freq-1)):size(S,3)) ;
        S = S(:,:, 1:round(size(S,3)/(N_freq-1)):size(S,3)) ;

        steps = [4.3 6.1]*1e-2 ;
        Lx = size(S, 2); Ly = size(S, 1);
        x0 = 0 ; y0 = 0 ;
        xgrid = x0 + (0:Lx-1)*steps(1) ; ygrid = y0 + (0:Ly-1)*steps(2) ;
        [Xgrid, Ygrid] = meshgrid(xgrid, ygrid) ;

        decayThreshold = inf ;

    case 'random components'

        N_comp = 10;
        mod_Re_k_synth = 0.99*(rand(N_comp,1)*2-1)*pi ;
        theta_Re_k_synth = 2*pi./(1 + randi(20, N_comp, 1)) ;
        Re_k_synth = mod_Re_k_synth .* [cos(theta_Re_k_synth) sin(theta_Re_k_synth)] ;

        mod_Im_k_synth = 0.029*(rand(N_comp,1)*2-1)*pi ;
        theta_Im_k_synth = 2*pi./(1 + randi(20, N_comp, 1)) ;
        Im_k_synth = mod_Im_k_synth .* [cos(theta_Im_k_synth) sin(theta_Im_k_synth)] ;

        k_synth = Re_k_synth + 1j * Im_k_synth ;

        Amp = ones(1,N_comp) ;

        steps = [1 1] ;
        Lx = 31; Ly = 30;
        x0 = 0 ; y0 = 0 ;
        xgrid = x0 + (0:Lx-1) * steps(1) ; ygrid = y0 + (0:Ly-1) * steps(2);
        [Ygrid, Xgrid] = meshgrid(ygrid, xgrid) ;

        w = 1 ;
        S = Amp * exp(1j*k_synth * [Xgrid(:), Ygrid(:)].') ;
        S = reshape(S, size(Xgrid)) ;

        if Lx==1
            S = repmat(S, 2, 1);
        elseif Ly==1
            S = repmat(S, 1, 2);
        end
        
        decayThreshold = inf ;

    case 'Green function_Thin Isotropic Homegeneous'
        
        E = 210e9; % Young's modulus
        Nu = 0.3 ; % Poisson ratio
        rho = 7850 ; % mass density
        h = 1e-3 ; % thickness
        D = 1/12 * h^3 * E / (1-Nu^2) ; % bending rigidity

        F0 = 1; % load magnitude

        % natural frequency
        w_fcn = @(f) 2*pi*f ;

        % Distance from the source
        R = @(X, Y) sqrt(X.^2 + Y.^2) ;

        % Hankel funtion and modified bessel function of the second kind
        H0_2 = @(x) besselh(0, 2, x) ;
        K0 = @(x) besselk(0, x) ;

        % simplified coefficients
        Beta = (D / (rho * h))^(1/4) ; % sqrt(w)/k
        Q = F0 / (rho * h) ;

        steps = [1e-3 1e-3] ;
        Lx = round(10e-2 / steps(1)) ; Ly = round(10e-2 / steps(2)) ;
        xgrid = (0:Lx-1) * steps(1) ; ygrid = (0:Ly-1) * steps(2) ;
        [Ygrid, Xgrid] = meshgrid(ygrid, xgrid) ;

        % Green's function of an infinite plate (Time harmonic)
        Gr_fcn = @(f, X0, Y0) -Q / (4*pi * Beta^2 * w_fcn(f)) * (pi*1j/2 * H0_2(R(Xgrid-X0, Ygrid-Y0) * sqrt(w_fcn(f)) / Beta) +...
                                                                  K0(R(Xgrid-X0, Ygrid-Y0) * sqrt(w_fcn(f)) / Beta)       ) ; 
        
        vect_y0 = linspace(-100*Ly, -0.01, 50);
        S = zeros([size(Xgrid), length(vect_y0)]) ;
        for i_y0 = 1:length(vect_y0)
            % Construction of the signal
            f = 5e3 ; X0 = 0*steps(1) ; Y0 = vect_y0(i_y0)*steps(2) ;
            S(:,:,i_y0) = Gr_fcn(f, X0, Y0) ;
        end
        w = repmat(w_fcn(f), 1, size(S, 3)) ;

        theta = linspace(0,2*pi,10).' ;
        k_synth = sqrt(w_fcn(f))/Beta .* [cos(theta) sin(theta)].*steps ;
        k_synth = repmat(k_synth, 1, 1, length(vect_y0)) ;

        decayThreshold = inf ;

end

%% Mechanical parameters estimation

clear D_est

method = {'HRWA', 'HRWA', 'FAT', 'FAT', 'CFAT'} ;
model = {'Thin-Isotropic-Homogeneous', 'Thin-Anisotropic-Homogeneous', 'Thin-Isotropic-Homogeneous', 'Thin-Anisotropic-Homogeneous', 'Thin-Isotropic-Homogeneous'} ;

for i_method = 1:length(method)
    switch method{i_method}
        case 'HRWA'
            [k, decay, c, lbda, S_reco] = HRWA(steps, [1 1], w, S, decayThreshold) ;
    
            if Ly == 1
                S(:,2,:) = [] ; S_reco(:,2,:) = [] ;
            elseif Lx == 1
                S(2,:,:) = [] ; S_reco(2,:,:) = [] ;
            end
        
            switch model{i_method}
                case 'Thin-Isotropic-Homogeneous'
                    A = @(theta) ones(size(theta)) ;
                case 'Thin-Anisotropic-Homogeneous'
                    A = @(theta) [cos(theta).^4, sin(theta).^4, 2*cos(theta).^2.*sin(theta).^2, 4*cos(theta).^3.*sin(theta), 4*cos(theta).*sin(theta).^3] ;
            end
            D_est{i_method} = arrayfun(@(kk) A(angle(real(k{kk}./steps)*[1; 1j])) ...
                                             \ (vecnorm(real(k{kk})./steps, 2, 2).^4 / rho/h/w(kk)^2).^(-1), 1:length(k), 'Uni', 0) ;
    
        case {'FAT', 'CFAT'}
    
            tol = 1e-1 ; maxIter = 100 ;
            for i_w = 1:length(w)
                [Coefs{i_w}, k{i_w}] = C_FAT(steps, S(:,:,i_w), model{i_method}, method{i_method}, tol, maxIter) ;
                
                % normalization of wavevectors with respect to the spatial steps
                k{i_w} = k{i_w} .* steps ;

                % Coefficients of the differential equation (structural parameters)
                D_est{i_method}{i_w} = rho*h*w(i_w).^2 .* Coefs{i_w} ; clear Coefs
            end
    end
    D_est{i_method} = cat(2, D_est{i_method}{:}) ;
end


figure('Name', 'Structural parameters')
yline(D, '--')
hold on
plot(vect_y0*steps(2), real(D_est{1}(1,:)), 'g-')
plot(vect_y0*steps(2), real(D_est{2}(2,:)), 'r-')
plot(vect_y0*steps(2), real(D_est{3}(1,:)), 'g:')
plot(vect_y0*steps(2), real(D_est{4}(2,:)), 'r:')
plot(vect_y0*steps(2), real(D_est{5}(1,:)), 'g-.')
xlabel('Distance (m)') ; ylabel('D_{22} (N.m)') ;
set(gca, 'FontSize', 20)
set([gca().Children], 'LineWidth', 2)
legend(["Theoretical" string(method) + " " + string(model)],  'Location', 'SouthEast', 'FontSize', 14)
axis 'tight'

%% 2D FFT of signal

if Ly==1 || Lx==1
    S_FFT = fft(S) ;
else
    S_FFT = fft2(S) ;
end
S_FFT = fftshift(fftshift(S_FFT), 3) ;

% return;
%% Plot of the results

fig = findobj('Name', 'Input signal vs reconstruction') ;
if isempty(fig)
    figure('Name', 'Input signal vs reconstruction') ;
else
    figure(fig); clf
end

ax1 = nexttile ;
surf(ax1, Xgrid, Ygrid, db(real(S(:,:,1))), 'ZDataSource', 'db(real(S(:,:,i_w)))')
shading(ax1, 'interp'); view(ax1, [0 90]);
% imagesc(ax1, real(S))
axis(ax1, 'equal', 'tight')
subtitle(ax1, 'Input signal')
xlabel(ax1, 'x (m)') ; ylabel(ax1, 'y (m)')

ax2 = nexttile ; 
surf(ax2, Xgrid, Ygrid, db(real(S_reco(:,:,1))), 'ZDataSource', 'db(real(S_reco(:,:,i_w)))')
shading(ax2, 'interp'); view(ax2, [0 90]);
% imagesc(ax2, real(S_reco))
axis(ax2, 'equal', 'tight')
subtitle(ax2, 'Reconstructed signal')
xlabel(ax2, 'x (m)') ; ylabel(ax2, 'y (m)')

linkprop([ax1 ax2], {'View' 'Clim'}) ;

fig = findobj('Name', 'Expected wavenumbers vs estimated') ;
if isempty(fig)
    figure('Name', 'Expected wavenumbers vs estimated') ;
else
    figure(fig); clf
end

ax3 = nexttile ; hold(ax3, 'on') ;
ax4 = nexttile ; hold(ax4, 'on') ;

ph1 = pcolor(ax3, linspace(-pi, pi-rem(size(S,1)+1,2)*pi/size(S,1), size(S,1)), linspace(-pi, pi-rem(size(S,2)+1,2)*pi/size(S,2), size(S,2)), db(abs(S_FFT(:,:, 1))).') ; shading(ax3, 'interp')

try
    % plot3(ax3, real(k_synth(:,1,1)), real(k_synth(:,2,1)), w(1)*ones(size(k_synth,1)), 'b+', 'XDataSource', 'real(k_synth(:,1,i_w))',...
    %                                                                                          'YDataSource', 'real(k_synth(:,2,i_w))',...
    %                                                                                          'ZDataSource', 'w(i_w)*ones(size(k_synth,1))')
    
    % plot3(ax4, imag(k_synth(:,1,1)), imag(k_synth(:,2,1)), w(1)*ones(size(k_synth,1)), 'b+', 'XDataSource', 'imag(k_synth(:,1,i_w))',...
    %                                                                                          'YDataSource', 'imag(k_synth(:,2,i_w))',...
    %                                                                                          'ZDataSource', 'w(i_w)*ones(size(k_synth,1))') 

    plot(ax3, real(k_synth(:,1,1)), real(k_synth(:,2,1)), 'b+', 'XDataSource', 'real(k_synth(:,1,i_w))',...
                                                                'YDataSource', 'real(k_synth(:,2,i_w))')

    plot(ax4, imag(k_synth(:,1,1)), imag(k_synth(:,2,1)), 'b+', 'XDataSource', 'imag(k_synth(:,1,i_w))',...
                                                                'YDataSource', 'imag(k_synth(:,2,i_w))')
end

% plot3(ax3, real(k{1}(:,1)), real(k{1}(:,2)), w(1)*ones(size(k{1},1)), 'ro', 'XDataSource', 'real(k{i_w}(:,1))',...
%                                                                             'YDataSource', 'real(k{i_w}(:,2))',...
%                                                                             'ZDataSource', 'w(i_w)*ones(size(k{i_w},1))')

% plot3(ax4, imag(k{1}(:,1)), imag(k{1}(:,2)), w(1)*ones(size(k{1},1)), 'ro', 'XDataSource', 'imag(k{i_w}(:,1))',...
%                                                                             'YDataSource', 'imag(k{i_w}(:,2))',...
%                                                                             'ZDataSource', 'w(i_w)*ones(size(k{i_w},1))')
% zlim(ax3, [w(1) w(end)]) ;
% zlim(ax4, [w(1) w(end)]) ;

plot(ax3, real(k{1}(:,1)), real(k{1}(:,2)), 'ro', 'XDataSource', 'real(k{i_w}(:,1))',...
                                                  'YDataSource', 'real(k{i_w}(:,2))')

plot(ax4, imag(k{1}(:,1)), imag(k{1}(:,2)), 'ro', 'XDataSource', 'imag(k{i_w}(:,1))',...
                                                  'YDataSource', 'imag(k{i_w}(:,2))')

xlabel(ax3, 'Re(k_x)'); ylabel(ax3, 'Re(k_y)'); view(ax3, [0 90]) ; 
xlabel(ax4, 'Im(k_x)'); ylabel(ax4, 'Im(k_y)'); view(ax4, [0 90]) ; 

linkprop([ax3 ax4], {'View'}) ;

axis(ax3, 'equal', 'square', 'tight')
axis(ax4, 'equal', 'square')

%% Animating the plots along the 3rd dimension of S

for i_w = 1:size(S, 3)
    
    refreshdata([ax1 ax2 ax3 ax4])
    ax1.Title.String = [num2str(w(i_w)/2/pi) 'Hz'] ;

    set(ph1, 'CData', db(abs(S_FFT(:,:, i_w))).')

    drawnow ; pause(0.01)
end

%% Reconstruction error

% Estimation error between theoretical and estimated wavenumber
Err_k(:,1) = abs(vecnorm(real(k_synth(1,:,1))./steps,2,2) - vecnorm(real(cat(1, k{:})./steps),2,2)) ./ vecnorm(real(k_synth(1,:,1))./steps,2,2);
Err_k(:,2) = abs(vecnorm(imag(k_synth(1,:,1))./steps,2,2) - vecnorm(imag(cat(1, k{:})./steps),2,2)) ./ vecnorm(imag(k_synth(1,:,1))./steps,2,2);

if ~isscalar(Err_k)
    fig = findobj('Name', 'Estimation error') ;
    if isempty(fig)
        figure('Name', 'Estimation error') ;
    else
        figure(fig); 
    end
    yyaxis left
    plot(Err_k(:,1)) % vect_y0*steps(2), 
    ylabel('Estimation error (m^{-1})')
    yyaxis right
    plot(Err_k(:,2))
    xlabel('Source position along Y axis (m)')
    legend('Real part', 'Imaginary part')
end