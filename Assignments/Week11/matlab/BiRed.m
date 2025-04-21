function [ A_out, t_out, r_out ] = BiRed( A, t, r )
  [ ATL, ATR, ...
    ABL, ABR ] = FLA_Part_2x2( A, ...
                               0, 0, 'FLA_TL' );
  [ tT, ...
    tB ] = FLA_Part_2x1( t, ...
                         0, 'FLA_TOP' );
  [ rT, ...
    rB ] = FLA_Part_2x1( r, ...
                         0, 'FLA_TOP' );

  while ( size( ATL, 1 ) < size( A, 1 ) -1 ) % Changed loop condition
    
    [ A00,  a01,     A02,  ...
      a10t, alpha11, a12t, ...
      A20,  a21,     A22 ] = FLA_Repart_2x2_to_3x3( ATL, ATR, ...
                                                    ABL, ABR, ...
                                                    1, 1, 'FLA_BR' );
    [ t0, ...
      tau1, ...
      t2 ] = FLA_Repart_2x1_to_3x1( tT, ...
                                    tB, ...
                                    1, 'FLA_BOTTOM' );
    [ r0, ...
      rho1, ...
      r2 ] = FLA_Repart_2x1_to_3x1( rT, ...
                                    rB, ...
                                    1, 'FLA_BOTTOM' );

    %------------------------------------------------------------%
    
    % Householder from the left
    x = [ alpha11; a21 ];
    [ u, tau1 ] = Housev1( x );
    alpha11 = u(1);
    a21 = u(2:end); %implicitly zero after householder? householder stored over zeros
    
    % Apply H = I - (1/tau1) * u * u' to [a12t; A22]
    %u_full = [1; a21];
    W = (1/tau1) * (u' * [a12t; A22]);
    temp = [a12t; A22] - u * W;
    a12t = temp(1, :);
    A22 = temp(2:end, :);
    
    % Householder from the right (only if more than 1 column remains)
    % Step 1: Form the Householder vector for the right transformation
    
    x = a12t';  % Turn row into column vector
    [ u12, rho1 ] = Housev1( x );  % Compute Householder reflector
    
    % Step 2: Store Householder vector back in a12t
    a12t = u12';  % Store transposed (row vector)
    
    % Step 3: Build full Householder vector with first element = 1
    v = a12t';       % column form of stored a12t
    v(1) = 1;        % First element explicitly set to 1
    
    % Step 4: Apply H = I - rho1 * v * v' from the right to A22
    %         A22 := A22 * H = A22 - (A22 * v) * ((1/rho1) * v)'
    w = A22 * v;     % Compute A22 * v
    A22 = A22 - w * ((1/rho1) * v)';  % Apply the Householder from the right

    %if size( a12t, 2 ) > 0
      % x = [alpha11, a12t ];
      % [ v, rho1 ] = Housev1( x' ); %shouldn't this be just a12t'?
      % a12t = v(2:end)'; 
      % a12t(2:end) = 0; %implicitly zeros except first element?
      % alpha11 = v(1);
      % 
      % % Apply G = I - rho1 * v * v' to [a21, A22] from the right
      % v_full = [1; a12t']
      % W = rho1 * ([a21, A22] * v_full);
      % temp = [a21, A22] - W * v_full';
      % a21 = temp(:, 1);
      % A22 = temp(:, 2:end);
    % else
    %   rho1 = 0;
    % end

    %------------------------------------------------------------%
    [ ATL, ATR, ...
      ABL, ABR ] = FLA_Cont_with_3x3_to_2x2( A00,  a01,     A02,  ...
                                             a10t, alpha11, a12t, ...
                                             A20,  a21,     A22, ...
                                             'FLA_TL' );
    [ tT, ...
      tB ] = FLA_Cont_with_3x1_to_2x1( t0, ...
                                       tau1, ...
                                       t2, ...
                                       'FLA_TOP' );
    [ rT, ...
      rB ] = FLA_Cont_with_3x1_to_2x1( r0, ...
                                       rho1, ...
                                       r2, ...
                                       'FLA_TOP' );

  end

  A_out = [ ATL, ATR
            ABL, ABR ];
  t_out = [ tT
            tB ];
  r_out = [ rT
            rB ];

  % Extract the diagonal and first superdiagonal
  %d = diag(A_out);
  %sd = diag(A_out(1:end-1, 2:end));
  %A_out = [d, [sd; 0]]; % Store diagonal and first superdiagonal in A_out

return