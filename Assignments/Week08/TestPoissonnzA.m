function TestPoissonnzA()
    %N = 8; % Example grid size
    %[nzA, ir, ic] = Create_Poisson_nzA(N);
    [nzA,ir,ic] = Create_mdim1_Poisson_problem_nzA(2, 3)
    % Convert CRS format back into a sparse matrix
    % n = (N+2) * (N*2);
    % A_crs = sparse([], [], [], n, n, length(nzA));
    % 
    % for row = 1:n
    %     for idx = ir(row):(ir(row+1)-1)
    %         col = ic(idx);
    %         A_crs(row, col) = nzA(idx);
    %     end
    % end

    % Generate reference Poisson matrix using MATLAB's built-in spdiags
    % e = ones(n, 1);
    % main_diag = 4 * e;
    % off_diag = -1 * e;
    % 
    % A_ref = spdiags([off_diag, off_diag, main_diag, off_diag, off_diag], ...
    %                 [-N, -1, 0, 1, N], n, n);
    % 
    % Compare the matrices
    %disp(full(A_crs))
    %disp(full(A_ref));
    %assert(isequal(A_crs, A_ref), 'Test failed: Matrices do not match.');
    
    %disp('Test passed: The generated matrix matches the expected Poisson matrix.');
    
end
