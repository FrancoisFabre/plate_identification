function [Coefs, k] = C_FAT(Steps, S, Model, Solver, tolerance, maxIter)
    %% Corrected Force Analysis Technique of the signal S with associated plate model Model
    %
    % Inputs :
    %       Steps     : Grid size along both dimensions
    %       S         : Signal to be analyzed
    %       Model     : Type of plate model
    %       Solver    : 'FAT' or 'CFAT'
    %       Tolerance : threshold for the relative change in solution of the
    %                   Gauss-Legendre iterative solver (only for CFAT)
    %       maxIter   : maximum number of iterations for the Gauss-Legendre iterative
    %                   solver (only for CFAT)
    % Outputs :
    %       Coefs  : Structural parameters corresponding to coefficients of the equation of motion
    %       k      : Estimated wavevectors 
    
%% Definition of bilaplacian estimator 

    % Renaming the grid steps
    dx = Steps(1) ; dy = Steps(end) ; clear Steps

    % Discretized stiffness operator of the system
    switch Model
        case 'Thin-Isotropic-Homogeneous'

            dw = @(i_x, i_y) cat(3, dw_4x(i_x, i_y), dw_4y(i_x, i_y), 2*dw_2x2y(i_x, i_y)) ;

            % Initialisation of the discrete bilaplacian estimator
            der_bounds = [2 2];
            dW = zeros([size(S) - der_bounds, 3]) ;

            k_th = @(Coefs, theta) (([1 1 2].*Coefs) * squeeze(prod([cos(theta) sin(theta)].^[4 0; ...
                                                                                              0 4; ...
                                                                                              2 2], 2))).^(-1/4);
    
        case 'Thin-Anisotropic-Homogeneous'
            dw = @(i_x, i_y) cat(3, dw_4x(i_x, i_y), dw_4y(i_x, i_y), 2*dw_2x2y(i_x, i_y), 4*dw_3x1y(i_x, i_y), 4*dw_1x3y(i_x, i_y)) ;

            % Initialisation of the discrete bilaplacian estimator
            der_bounds = [4 4];
            dW = zeros([size(S) - der_bounds, 5]) ;  

            k_th = @(Coefs, theta) (([1 1 2 4 4].*Coefs) * squeeze(prod([cos(theta) sin(theta)].^[4 0; ...
                                                                                                  0 4; ...
                                                                                                  2 2; ...
                                                                                                  3 1; ...
                                                                                                  1 3], 2))).^(-1/4);
    end

    % Displacement field (truncated to account for unavailable derivatives) with
    % contracted spatial dimensions x and y
    W = reshape(S(1+der_bounds(1):end-der_bounds(2), 1+der_bounds(1):end-der_bounds(2)), [], 1) ;

    % constuction of the discrete bilaplacian estimator
    for i_1 = 1+der_bounds(1):size(S, 1)-der_bounds(2) % along the 1st dimension
        for i_2 = 1+der_bounds(1):size(S, 2)-der_bounds(2) % along the 2nd dimension
            dW(i_1, i_2,:) = dw(i_1, i_2) ;
        end
    end

    dW = dW(1+der_bounds(1):end, 1+der_bounds(1):end, :) ;
    dW = reshape(permute(dW, [3 1 2]), size(dW, 3), []).' ;

%% Estimation of structural parameters and wavenumbers

    switch Solver
        case 'FAT'

            switch Model
                case 'Thin-Isotropic-Homogeneous'
                    Coefs = sum(dW, 2) \ W ;
                case 'Thin-Anisotropic-Homogeneous'
                    Coefs = dW \ W ;
            end
            
        case 'CFAT'
         
            % Gauss-Legendre iterative solving of the nonlinear LS problem
            % Coefs = NLsolver ;
            Coefs = NLsolver2 ;            
    end

    % Wanumber computation
    theta = permute(linspace(0, 2*pi, 50 + 1), [1 3 2]) ; theta = theta(1,1,1:end-1) ;
    k = k_th(Coefs.', theta).' .* [cos(theta(:)) sin(theta(:))] ;

%% Utilitary functions

    function A = dw_4x(i_x, i_y)
        % 4th derivative along x

        switch Model
            case 'Thin-Isotropic-Homogeneous'
                A = [1 -4 6 -4 1] * S(i_x + (-2:2), i_y) / dx^4 ;
            case 'Thin-Anisotropic-Homogeneous'
                A = [1 -4 6 -4 1] * S(i_x + 2 * (-2:2), i_y) / (2*dx)^4 ;
        end
        
    end

    function A = dw_4y(i_x, i_y)
        % 4th derivative along y

        switch Model
            case 'Thin-Isotropic-Homogeneous'
                A = S(i_x, i_y + (-2:2)) * [1 -4 6 -4 1].' / dy^4 ;
            case 'Thin-Anisotropic-Homogeneous'
                A = S(i_x, i_y + 2 * (-2:2)) * [1 -4 6 -4 1].' / (2*dy)^4 ;
        end

    end

    function A = dw_2x2y(i_x, i_y)
        % 2nd derivative along x and y

        coef = [ 1  -2  1;
                -2   4 -2;
                 1  -2  1] ;

        switch Model
            case 'Thin-Isotropic-Homogeneous'
                A = sum(coef .* S(i_x + (-1:1), i_y + (-1:1)), "all") / dx^2 / dy^2 ;
            case 'Thin-Anisotropic-Homogeneous'
                A = sum(coef .* S(i_x + 2 * (-1:1), i_y + 2 * (-1:1)), "all") / (2*dx)^2 / (2*dy)^2 ;
        end

    end

    function A = dw_3x1y(i_x, i_y)
        % 3rd derivative along x and 1st derivative along y

        coef = [ 1   0  -1;
                 3   0  -3;
                -3   0   3;
                 1  -2   1] ;

        switch Model
            case 'Thin-Isotropic-Homogeneous'
                A = sum(coef .* S(i_x + (-1.5:1:1.5), i_y + (-0.5:0.5:0.5)), "all") / dx^3 / dy ; % USELESS
            case 'Thin-Anisotropic-Homogeneous'
                A = sum(coef .* S(i_x + (-3:2:3), i_y + (-1:1)), "all") / (2*dx)^3 / (2*dy) ;
        end

    end

    function A = dw_1x3y(i_x, i_y)
        % 1st derivative along x and 3rd derivative along y

        coef = [ 1   0  -1;
                 3   0  -3;
                -3   0   3;
                 1  -2   1].' ;

        switch Model
            case 'Thin-Isotropic-Homogeneous'
                A = sum(coef .* S(i_x + (-0.5:0.5:0.5), i_y + (-1.5:1:1.5)), "all") / dx / dy^3 ; % USELESS
            case 'Thin-Anisotropic-Homogeneous'
                A = sum(coef .* S(i_x + (-1:1), i_y + (-3:2:3)), "all") / (2*dx) / (2*dy)^3 ;
        end
    end




    function A = dw_4x2(i_x, i_y)
        % 4th derivative along x

        A = tensorprod([1 -4 6 -4 1], S(i_x + (-2:2), i_y), 2, 1) / dx^4 ;
    end

    function A = dw_4y2(i_x, i_y)
        % 4th derivative along y

        A = tensorprod([1 -4 6 -4 1], S(i_x, i_y + (-2:2)), 2, 2 ) / dy^4 ;
    end

    function A = dw_2x2y2(i_x, i_y)
        % 2nd derivative along x and y

        coef = [ 1  -2  1;
                -2   4 -2;
                 1  -2  1] ;
        A = tensorprod(coef, S(i_x + (-1:1), i_y + (-1:1)), [1 2], [1 2]) / dx^2 /dy^2 ;
    end




    function out = NLsolver
        % Solves a nonlinear LS problem A = argmin ||A*B(A) - w||^2
        
        B = @(A) dW * [1;
                       1;
                       Xi(Zeta(A))/ 2] ;

        res = @(A) A * B(A) - W ;
        dres_dA = @(A) B(A) + A * dZeta_dA(A) * dXi_dX(Zeta(A)) * dW(:,3) / 2 ;

        % Initialisation of the nonlinear LS problem (from the FAT solution)
        A(1) = sum(dW, 2) \ W ;

        % A_glob = A(1) ;

        % Gauss-Legendre iterative solver
        Iter = 0 ;
        rel_Err = inf ;
        while all(rel_Err >= tolerance) && Iter<maxIter
            % Increment the iteration count
            Iter = Iter + 1 ;

            % Compute A_n+1
            A(2) = A(1) - dres_dA(A(1)) \ res(A(1)) ;

            % Relative error between A_n+1 and A_n (ckecking for convergence)
            rel_Err = abs(diff(A)) / abs(A(1)) ;

            % A_n becomes A_n+1
            A(1) = A(2) ;

            % A_glob(end+1) = A(2) ;
        end

        % plot(abs(A_glob));

        out = dx^2 * dy^2 * Zeta(A(2))^(-4) ;

        function Y = Xi(X)
            % custom function allowing the rewritting of the discrete operator (CFAT)
            
            Y = (1 - cos(X)).^2 ./ (1 - cos(X/sqrt(2))).^2 - 2 ;
        end
    
        function Y = dXi_dX(X)
            % 1st derivative of Xi w.r.t M
    
            Y = 2 * ((1 - cos(X)) .* (1 - cos(X / sqrt(2))) .* sin(X)   -   1/sqrt(2) .* sin(X / sqrt(2)) .* (1 - cos(X)).^2) ...
                    / (1 - cos(X / sqrt(2))).^3 ;
        end
    
        function Y = Zeta(A)
            % custom function allowing the rewritting of the discrete operator (CFAT)
    
            Y = acos(1 - dx*dy / 2/sqrt(A)) ;
        end
    
        function Y = dZeta_dA(A)
            % 1st derivative of Zeta w.r.t A
    
            Y = - dx*dy * A^(-3/2) / sqrt(1 - (1 - dx*dy/ 2/sqrt(A))^2) / 4 ;
        end
    end

    function out = NLsolver2
        % Solves a nonlinear LS problem A = argmin ||A*B(A) - w||^2
        
        % Initialisation of the nonlinear LS problem (from the FAT solution)
        Ratio_Init = dW \ W ;
        
        switch Model
            case 'Thin-Isotropic-Homogeneous'
                A(:,1) = sum(dW,2) \ W ;
            case 'Thin-Anisotropic-Homogeneous'
                A(:,1) = Ratio_Init ;
        end    

        % A_glob = A(:,1) ;

        % Gauss-Legendre iterative solver
        Iter = 0 ;
        rel_Err = inf ;
        while all(rel_Err >= tolerance) && Iter<maxIter
            % Increment the iteration count
            Iter = Iter + 1 ;

            % Update correction coefficients
            Ntheta = 1000 ;
            theta_tmp = linspace(0,pi, Ntheta).' ;
            X = 1/dx * sin(k_th(A(:,1).' , permute(theta_tmp, [2 3 1])).' .* cos(theta_tmp) * dx) ; 
            Y = 1/dy * sin(k_th(A(:,1).' , permute(theta_tmp, [2 3 1])).' .* sin(theta_tmp) * dy) ; 
            % X = pi * k_th(A(:,1).' , permute(theta_tmp, [2 3 1])).' .* cos(theta_tmp) .* sinc(k_th(A(:,1).' , permute(theta_tmp, [2 3 1])).' .* cos(theta_tmp) * dx / pi) ; 
            % Y = pi * k_th(A(:,1).' , permute(theta_tmp, [2 3 1])).' .* sin(theta_tmp) .* sinc(k_th(A(:,1).' , permute(theta_tmp, [2 3 1])).' .* sin(theta_tmp) * dy / pi) ; 

            switch Model
                case 'Thin-Isotropic-Homogeneous'
                    mu = ((X.^[4 0 2] .* Y.^[0 4 2]) .* ([1 1 2].*A(:,1).')) \ ones(size(theta_tmp)) ;

                    % Compute A_n+1
                    A(:,2) = mu \ Ratio_Init ;
                    % A(:,2) = (dW * mu) \ W ;

                case 'Thin-Anisotropic-Homogeneous'
                    mu = ((X.^[4 0 2 3 1] .* Y.^[0 4 2 1 3]) .* ([1 1 2 4 4].*A(:,1).')) \ ones(size(theta_tmp)) ;
                    
                    % Compute A_n+1
                    A(:,2) = Ratio_Init ./ mu ;
                    % A(:,2) = (dW .* mu.') \ W ;
            end           

            % Relative error between A_n+1 and A_n (ckecking for convergence)
            rel_Err = abs(diff(A,1,2)) ./ abs(A(:,1)) ;

            % A_n becomes A_n+1
            A(:,1) = A(:,2) ;

            % A_glob(:,end+1) = A(:,2) ;
        end

        % figure(1)
        % plot(abs(A_glob));
        out = A(:,2) ;
    end
end