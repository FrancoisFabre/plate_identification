function pn = Stabchart_selection(pn)
    %
    % Input :
    %       pn :
    %
    % Output
    %       pn : 
    

%% Determination of stable poles between consecutive orders

    % tolerances to determine the stability of the poles between consecutive orders
    tol_abs_Re_pn = 1e-2 ; 
    tol_arg_Re_pn = 1e-2 ;
    tol_abs_Im_pn = 1e-1 ;
    tol_arg_Im_pn = 1e-1 ;

    abs_Re_pn = cellfun(@(x) abs(real(x) * [1; 1j]), pn, 'UniformOutput', false) ;
    arg_Re_pn = cellfun(@(x) angle(real(x) * [1; 1j]), pn, 'UniformOutput', false) ;
    abs_Im_pn = cellfun(@(x) abs(imag(x) * [1; 1j]), pn, 'UniformOutput', false) ;
    arg_Im_pn = cellfun(@(x) angle(imag(x) * [1; 1j]), pn, 'UniformOutput', false) ;

    % Initialisation of the stable poles celle array
    pn_st = cell(size(pn)) ;
    for i_order = 2:length(pn)
        % for each order, compare the poles with the ones of the previous order
    
        [diff_abs_Re, idiff_abs_Re] = min(abs(abs_Re_pn{i_order} - abs_Re_pn{i_order-1}.') ./ abs_Re_pn{i_order}, [], 2) ;
        diff_arg_Re = abs(arg_Re_pn{i_order} - arg_Re_pn{i_order-1}(idiff_abs_Re)) ./ arg_Re_pn{i_order} ;
    
        % [diff_abs_Im, idiff_abs_Im] = min(abs(abs_Im_pn{i_order} - abs_Im_pn{i_order-1}.') ./ abs_Im_pn{i_order}, [], 2) ;
        diff_abs_Im = abs(abs_Im_pn{i_order} - abs_Im_pn{i_order-1}(idiff_abs_Re)) ./ abs_Im_pn{i_order} ;
        diff_arg_Im = abs(arg_Im_pn{i_order} - arg_Im_pn{i_order-1}(idiff_abs_Re)) ./ arg_Im_pn{i_order} ;
    
        % Determine stable poles with respect to the previous order (according to tolerances)
        StabCond =   (diff_abs_Re <= tol_abs_Re_pn)...
                   & (diff_arg_Re <= tol_arg_Re_pn)...
                   & (diff_abs_Im <= tol_abs_Im_pn)...
                   & (diff_arg_Im <= tol_arg_Im_pn) ;

        % Stable poles at the current order
        pn_st{i_order} = pn{i_order}(StabCond,:) ;

    end

%% Automatic selection of stable poles across all order
    
    nbOrd = length(pn_st) ;
    nbPol = cellfun(@(x) size(x,1) , pn_st) ;
    
    % from cell to matrix
    pn_st_m = zeros(nbOrd, max(nbPol), 2) ;
    for i=1:nbOrd
        pn_st_m(i, 1:nbPol(i),:) = pn_st{i} ;
    end
    
    io = nbOrd ;
    nbit = round(nbOrd/3) ;

    % Stable poles in frequency and damping at selected order
    pn_st_io = pn_st{io} ;

    % Norm of the real part of stable poles
    abs_Re_pn_m = abs(tensorprod(real(pn_st_m), [1; 1j], 3, 1)) ;
    abs_Re_pn_io = abs(real(pn_st_io) * [1; 1j]).' ;
        
    % Proximity condition between the selected order and all smaller orders
    proximityCond = (abs(reshape(abs_Re_pn_m(1:io, :), [], 1) - abs_Re_pn_io) ./ abs_Re_pn_io) < tol_abs_Re_pn ;

    % Number of repetition of the stable poles verifiying the proximity condition
    nbRep = sum(proximityCond) ;
    
    % Get the final list of stale poles
    pn = pn_st{io}(nbRep >= nbit, :) ;

end