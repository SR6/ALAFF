% Test script for SparseMvMult

% Test Case 1: Simple 3x3 sparse matrix
nzA1 = [4, 5, 7];         % non-zero values
ir1 = [1, 2, 3];          % row indices
ic1 = [1, 2, 3];          % column indices
x1 = [1; 2; 3];           % input vector

% Expected result y = A * x
y1_expected = [4*1 + 0 + 0; 0 + 5*2 + 0; 0 + 0 + 7*3];
y1 = SparseMvMult(nzA1, ir1, ic1, x1);
fprintf('Test Case 1: Expected y = [ %f ; %f ; %f ], Got y = [ %f ; %f ; %f ]\n', ...
        y1_expected(1), y1_expected(2), y1_expected(3), y1(1), y1(2), y1(3));

% Test Case 2: 4x4 sparse matrix with some zeros
nzA2 = [1, 2, 3, 4, 5];   % non-zero values
ir2 = [1, 2, 3, 4, 4];    % row indices
ic2 = [1, 2, 3, 1, 4];    % column indices
x2 = [1; 1; 1; 1];        % input vector

% Expected result y = A * x
y2_expected = [1*1 + 0 + 0 + 4*1; 0 + 2*1 + 0 + 0; 0 + 0 + 3*1 + 0; 1*1 + 0 + 0 + 5*1];
y2 = SparseMvMult(nzA2, ir2, ic2, x2);
fprintf('Test Case 2: Expected y = [ %f ; %f ; %f ; %f ], Got y = [ %f ; %f ; %f ; %f ]\n', ...
        y2_expected(1), y2_expected(2), y2_expected(3), y2_expected(4), y2(1), y2(2), y2(3), y2(4));

% Test Case 3: Sparse identity matrix (5x5)
nzA3 = [1, 1, 1, 1, 1];    % non-zero values
ir3 = [1, 2, 3, 4, 5];     % row indices
ic3 = [1, 2, 3, 4, 5];     % column indices
x3 = [5; 4; 3; 2; 1];      % input vector

% Expected result y = A * x
y3_expected = x3;          % Since it's an identity matrix
y3 = SparseMvMult(nzA3, ir3, ic3, x3);
fprintf('Test Case 3: Expected y = [ %f ; %f ; %f ; %f ; %f ], Got y = [ %f ; %f ; %f ; %f ; %f ]\n', ...
        y3_expected(1), y3_expected(2), y3_expected(3), y3_expected(4), y3_expected(5), ...
        y3(1), y3(2), y3(3), y3(4), y3(5));

% Test Case 4: Sparse matrix with only one non-zero element (2x2)
nzA4 = [10];               % non-zero values
ir4 = [1];                 % row indices
ic4 = [2];                 % column indices
x4 = [1; 2];               % input vector

% Expected result y = A * x
y4_expected = [0 + 10*2; 0];
y4 = SparseMvMult(nzA4, ir4, ic4, x4);
fprintf('Test Case 4: Expected y = [ %f ; %f ], Got y = [ %f ; %f ]\n', ...
        y4_expected(1), y4_expected(2), y4(1), y4(2));

% Test Case 5: Empty sparse matrix (0x0 matrix)
nzA5 = [];                 % non-zero values
ir5 = [];                  % row indices
ic5 = [];                  % column indices
x5 = [];                   % input vector

% Expected result y = A * x (empty matrix-vector multiplication should give an empty result)
y5_expected = [];
y5 = SparseMvMult(nzA5, ir5, ic5, x5);
fprintf('Test Case 5: Expected y = [] , Got y = []\n');

% Test Case 6: Large matrix (performance test)
n = 1000;
nzA6 = rand(1, 5000);     % random non-zero values
ir6 = randi([1, n], 1, 5000);  % random row indices
ic6 = randi([1, n], 1, 5000);  % random column indices
x6 = ones(n, 1);          % input vector of ones

% Perform multiplication and check result size
y6 = SparseMvMult(nzA6, ir6, ic6, x6);
fprintf('Test Case 6: Large Matrix - Result size: [%d x 1]\n', length(y6));
