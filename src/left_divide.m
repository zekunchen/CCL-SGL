function [V] = left_divide(E,I,tol,~,V)
%[V] = LEFT_DIVIDE(E,I,tol,pp,V);
% 
% Implements left division for symmetric positive definite system solves
% such as the sparse forward solve and dense solve for a GN descent
% direction. LEFT_DIVIDE is optimised for symmetric matrices and overcomes
% small inefficiencies of matlab's mldivide. For non-symmetric solves 
% please use mldivide.
%
% Also uses conjugate gradients (for large problems).
%
% E   = The full rank system matrix
% I   = The currents matrix (RHS)
% tol = The tolerance in the forward solution, e.g. 1e-5
%
% pp,V are old options from previous solver. tilde used in arguments list
% to ignore pp and keep matlab's code analyzer happy

% (c) N. Polydorides 2003 % Copying permitted under terms of GNU GPL
% $Id: left_divide.m 6007 2019-06-29 13:32:48Z aadler $

if ischar(E) && strcmp(E,'UNIT_TEST'); do_unit_test; return; end

if ~exist('tol','var'); tol = 1e-8; end

[n_nodes,n_stims] = size(I);


    % V= E\I;
    % This takes MUCH longer when you have  more vectors in I,
    %  even if they are repeated. There must be some way to simplify
    %  this to speed it up. Matlab's sparse operators really should
    %  do this for you.
    
    % TODO: 
    % 1. change from QR implementation to basis implementation
    % 2. implement selection for required nodal values
    % 3. cache basis solve
    % 4. possibly change to itterative for successive solves on the same
    %    mesh
    if issparse(E)
        
% This should speed up, and help issue with octave on QR
        inotzeros = logical(any(I,2));
      if exist('OCTAVE_VERSION') == 5 % v 4.4 has problems with sparse qr
        [Qi,R] = qr(full(I(inotzeros,:)),0);
      else
        [Qi,R] = qr(I(inotzeros,:),0);
      end
        rnotzeros = logical(any(R,2));
        R= R(rnotzeros,:);
        Q = sparse(size(I,1), size(R,1));
        Q(inotzeros,:) = Qi(:,rnotzeros);
%        [Q,R] = qr(I,0);
%        rnotzeros = any(R~=0,2);
%        Q= Q(:,rnotzeros);
%        R= R(rnotzeros,:);
        V= (E \ Q)*R;
        
    else
        if isreal(E)
            try
                % for dense solve of tikhonov regularised least squares
                % matrix E is symmetric if it is of the form
                % (J.'*W*J + hp^2*R.'R) and is real
                opts.SYM=true;
                opts.POSDEF=true;
                
                V= linsolve(E,I,opts);
            catch Mexcp
                
                % error handling 
                if(strcmp(Mexcp.identifier,'MATLAB:posdef'))
                    
                    warning('EIDORS:leftDivideSymmetry',...
                        ['left_divide is optimised for symmetric ',...
                        'positive definite matrices.']);
                    
                else 
                    warning(['Error with linsolve in left_divide, trying backslash.\n',...
                        'Error identifier: ',Mexcp.identifier]);
                end
                
                % continue solve with backslash
                V=E\I;
            end
        else
            % cholesky only works for real valued system matrices
            V=E\I;
        end
    end
    
    


