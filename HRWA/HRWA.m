function [k, decay, c, lbda, S_reco] = HRWA(steps, R, w, S, decayThreshold)
    %% High Resolution Wavenumber Analysis of the signal S
    %
    % Inputs :
    %       steps          : 2D grid steps
    %       R              : bounds of the estimation orders for the 2D ESPRIT analysis
    %       w              : vector of natural frequencies corresponding to the 3rd dimension of S
    %       S              : signal to be analyzed (size X_samples*Y_samples*F_samples)
    %       decayThreshold : maximum spatial decay of estimated poles of the signal S
    %
    % Outputs :
    %       k              : Estimated wavenumbers (normalized between -pi and pi)
    %       decay          : Spatial decay of the estimated poles (wavenumbers)
    %       c              : Phase velocity of the estimated poles
    %       lbda           : Wavelength of the estimated poles
    %       S_reco         : Reconstructed signal from estimated poles

%%
    % Initialisation of the reconstructed signal array
    S_reco = zeros(size(S)) ;

    % Initialisation of wavenumber, phase velocity, wavelength and spatial decay ratio vectors
    k = cell(size(w)) ; c = k ; lbda = k ; decay = k ;

    % for each frequency
    for i_w = 1:length(w)

        % Wavenumber vector
        [k{i_w}, decay{i_w}, S_reco(:,:, i_w)] = esprit2D(R, S(:,:,i_w), decayThreshold, 'ESTER') ;

        % Wavenumber vector normalized by the grid steps size
        k_norm = k{i_w} ./ steps ;

        % Phase velocity vector
        c{i_w} = (w(i_w) / vecnorm(k_norm, 2, 2).^2) * conj(k_norm) ;
    
        % Wavelength vector
        lbda{i_w} = (2*pi / vecnorm(real(k_norm), 2, 2).^2) * real(k_norm) ;
    end
end