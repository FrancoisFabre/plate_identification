function [pn, delta_pn] = Stabchart_selection(pn, tol, nbit, criteria)
    %%
    % Input :
    %       pn       : 
    %       tol      :
    %       nbit     : 
    %       criteria : 
    %
    % Output
    %       pn           : 
    %       delta_pn : 
    
%% Poles clusters construction

    switch criteria
        case 'Norm-Angle'
            poles = cellfun(@(x) [abs(real(x) * [1; 1j]) rem(angle(real(x) * [1; 1j]) + 2*pi, 2*pi)...
                                  abs(imag(x) * [1; 1j]) rem(angle(imag(x) * [1; 1j]) + 2*pi, 2*pi)], pn, 'UniformOutput', false) ;
        case 'Real-Imag'
            poles = cellfun(@(x) [real(x) imag(x)], pn, 'UniformOutput', false) ;
    end

    % Initialisation of the cell array containing clusters of poles
    clusters = mat2cell(pn{1}, ones(size(pn{1}, 1), 1), 2) ;
    for i_order = 2:length(pn)
        clear diff_poles
        
        % for each order, compare the poles with the existing clusters
    
        % compute the "center" of each cluster
        switch criteria
            case 'Norm-Angle'
                centroids = cellfun(@(x) mean([abs(real(x) * [1; 1j]) rem(angle(real(x) * [1; 1j]) + 2*pi, 2*pi)...
                                               abs(imag(x) * [1; 1j]) rem(angle(imag(x) * [1; 1j]) + 2*pi, 2*pi)], 1), clusters, 'UniformOutput', false) ;
            case 'Real-Imag'
                centroids = cellfun(@(x) mean([real(x) imag(x)], 1), clusters, 'UniformOutput', false) ;
        end
        centroids = cat(1, centroids{:}) ;

        % compare all poles of the current order with the clusters' centroids
        diff_poles(:,:,1) = abs(abs(poles{i_order}(:,1) - centroids(:,1).') ./ centroids(:,1).') ;
        diff_poles(:,:,3) = abs(abs(poles{i_order}(:,3) - centroids(:,3).') ./ centroids(:,3).') ;
            
        tmp = abs(poles{i_order}(:,2) - centroids(:,2).') ;
        if tmp <= pi & strcmp(criteria, 'Real-Imag')
            diff_poles(:,:,2) = abs(tmp ./ centroids(:,2).') ;
        else
            diff_poles(:,:,2) = abs(((-1).^(tmp>pi).*tmp + (tmp>pi).*2*pi) / 2/pi) ;
        end

        tmp = abs(poles{i_order}(:,4) - centroids(:,4).') ;
        if tmp <= pi & strcmp(criteria, 'Real-Imag')
            diff_poles(:,:,4) = abs(tmp ./ centroids(:,4).') ;
        else
            diff_poles(:,:,4) = abs(((-1).^(tmp>pi).*tmp + (tmp>pi).*2*pi) ./ 2/pi) ;
        end
    
        % Determine stable poles with respect to the existing clusters (according to tolerances)
        StabCond =   (diff_poles(:,:,1) <= tol(1))...
                   & (diff_poles(:,:,2) <= tol(2))...
                   & (diff_poles(:,:,3) <= tol(3))...
                   & (diff_poles(:,:,4) <= tol(4)) ;

        % Check if multiple poles verifies the stability criteria for the same cluster
        %    or if a pole verifies the stability criteria for multiple clusters
        Conflict_poles_idx = find(sum(StabCond,1)>1) ;
        Conflict_clusters_idx = find(sum(StabCond,2)>1).' ;

        if ~isempty(Conflict_poles_idx)
            disp('Error: Multiple poles for one cluster.')

            for iii = Conflict_poles_idx

                tmp = find(StabCond(:, iii)) ;
                [~, i_min] = min(sum(diff_poles(tmp, iii, :).^2, 3)) ;

                % Attribute the pole to the closest clusters
                StabCond(:, iii) = 0 ; 
                StabCond(tmp(i_min), iii) = 1 ; 
            end
        end

        if ~isempty(Conflict_clusters_idx)
            disp('Error: One pole for multiple clusters.')
            for iii = Conflict_clusters_idx

                tmp = find(StabCond(iii, :)) ;
                [~, i_min] = min(sum(diff_poles(iii, tmp, :).^2, 3)) ;

                % Attribute the pole to the closest clusters
                StabCond(iii, :) = 0 ; 
                StabCond(iii, tmp(i_min)) = 1 ; 
            end
        end

        % Update clusters (stable poles are added to existing clusters, 
        %                  others make a new cluster)
        [i_r, i_c] = find(StabCond) ;
        for iii = 1:length(i_r)
            clusters{i_c(iii)} = [clusters{i_c(iii)};
                                  pn{i_order}(i_r(iii), :)] ;
        end
        [i_r, ~] = find(~sum(StabCond, 2)) ;
        clusters(end + (1:length(i_r))) = mat2cell(pn{i_order}(i_r,:), ones(length(i_r), 1), 2) ;
    end

%% Clusters filtering and best poles extraction

    % Filtering out clusters with less than Nrep repetitions
    Nrep_cond = cellfun(@(x) size(x,1), clusters) ;
    if any(Nrep_cond>=nbit)
        clusters = clusters(Nrep_cond>=nbit) ;
    else
        clusters = clusters(Nrep_cond==max(Nrep_cond)) ;
    end

    % Extraction of the centroid
    centroids = cellfun(@(x) mean(x, 1), clusters, 'UniformOutput', false) ;
    pn = cat(1, centroids{:}) ;

    % Compute centroids absolute and relative uncertainty (95 % confidence)
    switch criteria
        case 'Norm-Angle'
            delta_pn = cellfun(@(x) 2*std([abs(real(x) * [1; 1j]) rem(angle(real(x) * [1; 1j]) + 2*pi, 2*pi)...
                                           abs(imag(x) * [1; 1j]) rem(angle(imag(x) * [1; 1j]) + 2*pi, 2*pi)], 0, 1) / sqrt(size(x, 1)), clusters, 'UniformOutput', false) ;
        case 'Real-Imag'
            delta_pn = cellfun(@(x) 2*std([real(x) imag(x)], 0, 1) / sqrt(size(x, 1)), clusters, 'UniformOutput', false) ;
    end
    % Absolute uncertainty
    delta_pn = cat(1, delta_pn{:}) ;

end