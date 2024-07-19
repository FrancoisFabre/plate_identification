function [pn, decay, S_reco] = esprit2D(N_comp, S, decayThreshold, Criterion)
    %% inputs: 
    %       N_comp         : potential numbers of components in the signal
    %       S              : signal to be analyzed
    %       decayThreshold : maximum decay accepted for the estimated poles of the signal
    %       Criterion      : pole selection criterion, either 'ESTER' or 'Stabchart'
    % 
    % outputs:
    %       pn             : extracted continuous poles of the signal
    %       decay          : decay of the estimated poles
    %       S_reco         : reconstructed 2D signal
    
%% Hankel block Hankel matrix construction

    % Dimensions of the grid in both dimensions
    L1 = size(S,1) ;
    L2 = size(S,2) ;

    % total number of components looked for in the noisy signal
    K1 = floor(L1/2) ;
    K2 = floor(L2/2) ;

    % "row" dimensions of the hankel and block-hankel matrices
    N1 = L1 - K1 + 1 ;
    N2 = L2 - K2 + 1 ;

    % indices for the hankel and block-hankel matrices
    hankel_1 = hankel((1:N1).', N1:L1) ;
    hankel_2 = hankel((1:N2).', N2:L2) ;

    % preallocate the hankel matrix
    h = zeros(size(hankel_2,1), size(hankel_2,2), L1) ;

    % for each line of S
    for i_row = 1:L1
        
        % construct the hankel matrix
        h(:,:,i_row) =  reshape(S(i_row, hankel_2), size(hankel_2)) ;
    end

    % construct the block-hankel matrix (from pages of the hankel matrix)
    HbH = reshape(h(:,:,hankel_1), [size(hankel_2), size(hankel_1)]) ;
    HbH = permute(HbH, [1 3 2 4]) ;
    HbH = reshape(HbH, size(hankel_2) .* size(hankel_1)) ;

%% Auto-covariance construction and eigen-value decomposition

    % Auto-covariance matrix
    Css = 1/L1/L2 * (HbH * HbH') ;

    % Signal subspace construction
    Rmax = N_comp(end) ;
    [W, ~] = eigs(Css, Rmax) ; 

%% Poles extraction

    [F1, F2] = F1_F2_computation(W, N_comp) ;

    for i_order = 1:length(F1)
        K = 0.55 * F1{i_order} + (1-0.55) * F2{i_order} ;
    
        % Identify the transfer matrix between W and V (Vandermonde matrix made of discrete poles)
        [T, ~] = eig(K) ;
    
        % Estimation of discrete poles
        Z1{i_order} = diag(T \ F1{i_order} * T) ;
        Z2{i_order} = diag(T \ F2{i_order} * T) ;
    
        % Estimation of continuous poles (wavenumbers or natural frequencies)
        pn{i_order} = -1i * [log(Z1{i_order}), log(Z2{i_order})] ;
    end

    if strcmp(Criterion, 'ESTER')
        
        pn = pn{1} ;
        Z1 = Z1{1} ;
        Z2 = Z2{1} ;

    elseif strcmp(Criterion, 'Stabchart')
        pn_glob = cat(1, pn{:}) ;
        Z1_glob = cat(1, Z1{:}) ;
        Z2_glob = cat(1, Z2{:}) ;

        % fig = findobj('Name', 'Stabilization chart') ;
        % if isempty(fig)
        %     figure('Name', 'Stabilization chart')
        % else
        %     figure(fig) ; clf ;
        % end
        % ax1 = nexttile; hold on ; ax2 = nexttile; hold on ;
        % for ii = 1:length(F1)
        %     scatter3(ax1, abs(real(pn{ii}) * [1; 1j]), angle(real(pn{ii}) * [1; 1j]), ii, 'b', 'SizeData', 12)
        %     scatter3(ax2, abs(imag(pn{ii}) * [1; 1j]), angle(imag(pn{ii}) * [1; 1j]), ii, 'r', 'SizeData', 12)
        % 
        %     xlabel(ax1, 'Re(k) (m^{-1})') ; xlabel(ax2, 'Im(k) (m^{-1})') ; 
        %     ylabel(ax1, 'arg(Re(k)) (rad)') ; ylabel(ax2, 'arg(Im(k)) (rad)') ; 
        %     zlabel([ax1 ax2], 'order') ; view([ax1 ax2], [0 0]) ; axis([ax1 ax2], 'square')
        %     linkprop([ax1 ax2], 'view') ;
        % end

        % Tolerances to determine the stability of the poles between consecutive orders
        tol = [1e-2 1e-2 5e-1 5e-1] ; 

        % Number of repetition of a pole to be extracted from the stabilization chart
        Nrep = round(length(pn) / 3) ;

        pn = Stabchart_selection(pn, tol, Nrep, 'Norm-Angle') ;

        Z1 = exp(1i*pn(:, 1)) ; Z2 = exp(1i*pn(:, 2)) ;
    end

    % Spatial decay ratio vector
    decay = vecnorm(imag(pn), 2, 2) ./ vecnorm(real(pn), 2, 2) ;

    % Keeping poles with a decay ratio below the threshold
    Z1 = Z1(decay <= decayThreshold) ;
    Z2 = Z2(decay <= decayThreshold) ;
    pn = pn(decay <= decayThreshold, :) ;
    decay = decay(decay <= decayThreshold) ;
   
%% Reconstruction of the signal from estimated poles

    % Construction of the Vandermonde matrix from estimated discrete poles
    V1 = (Z1.^(0:L1-1)).' ;
    V2 = (Z2.^(0:L2-1)).' ;
    V = Katri_Rao(V1, V2) ;

    % Computation of poles amplitudes
    Amp = V \ reshape(S.', [], 1) ;
    
    % Reconstruction of the signal
    S_reco = reshape(V * Amp, size(S.')).' ;

%% findSignalOrder function definition

    function [F1, F2] = F1_F2_computation(W, N_comp)

        % If a single estimation order is prescribed, duplicate it
        if isscalar(N_comp) ; N_comp = [N_comp N_comp] ; end

        % Initialisation of F1, F2 and  error vectors
        F1 = cell(1, diff(N_comp)+1) ;
        F2 = F1 ;
        Err_1 = zeros(1, diff(N_comp)+1) ;
        Err_2 = Err_1 ;

        for r = 1:(diff(N_comp)+1)
            W_R = W(:, 1:(N_comp(1)-1 + r)) ;

            % Indices of rows of W_R used to access blocks corresponding to both dimensions
            idx = reshape(1:N1*N2, N2, N1) ;

            % Computes matrices for rotationnal invariance property along the 1st dimension
            W_up_1 = W_R(idx(:,1:end-1), :) ;
            W_down_1 = W_R(idx(:,2:end), :) ;

            % Computes matrices for rotationnal invariance property along the 2nd dimension
            W_up_2 = W_R(idx(1:end-1,:), :) ;
            W_down_2 = W_R(idx(2:end,:), :) ;

            % Computes F1 and F2
            F1{r} = W_up_1 \ W_down_1;
            F2{r} = W_up_2 \ W_down_2 ;

            % Computes the error of estimation along both dimensions
            Err_1(r) = norm(W_down_1 * F1{r} - W_up_1, 2) ;
            Err_2(r) = norm(W_down_2 * F2{r} - W_up_2, 2) ;
        end
        
        if strcmp(Criterion, 'ESTER')

            % Find the order which minimizes the 2D ESTER criterion
            ESTER_crit = (Err_1 .* Err_2).^(-1/2) ;
            ESTER_thresh = 0.5 * max(ESTER_crit) ;
            % [~, R] = max(ESTER_crit) ;
            R = find(ESTER_crit>= ESTER_thresh, 1, 'last') ;
    
            % Select the corresponding cells of F1 and F2
            F1 = F1(R) ;
            F2 = F2(R) ;
    
            % R = N_comp(R)
    
            % fig = findobj('Name', 'ESTER criterion') ;
            % if isempty(fig)
            %     figure('Name', 'ESTER criterion') ;
            % else
            %     figure(fig); clf
            % end
            % plot(ESTER_crit, '-o', 'MarkerIndices', R); hold on;
            % yline(ESTER_thresh, ':'); hold off;
            % drawnow;
            % pause(0.01)
    
        end
    end
end