B = diag(randn(10,1)) + diag(randn(9,1),1);  % Random bidiagonal matrix
[Sc, Uc, Vc] = SVD_BiDiag_ImpShift(B);
[Um, Sm, Vm] = svd(B);
S = svd(B);

fprintf("Singular values comparison (custom vs builtin):\n");
disp([Sc, S]);

B;
B_reconstructed = Uc * diag(Sc) * Vc';

fprintf("Frobenius norm of the reconstruction error: %.2e\n", norm(Uc*diag(Sc)*Vc' - B, 'fro'));
fprintf("Angle between first left singular vectors (U): %.2f degrees\n", rad2deg(acos(abs(Uc(:,1)' * Um(:,1)))));
fprintf("Angle between first right singular vectors (V): %.2f degrees\n", rad2deg(acos(abs(Vc(:,1)' * Vm(:,1)))));
fprintf("Is U orthogonal? %.2e\n", norm(Uc'*Uc - eye(size(Uc)), 'fro'));
fprintf("Is V orthogonal? %.2e\n", norm(Vc'*Vc - eye(size(Vc)), 'fro'));