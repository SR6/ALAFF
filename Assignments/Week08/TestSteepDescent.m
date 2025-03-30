function TestSteepDescent()
A = [[2, -1,0],
    [-1,4, -2],
    [0, -2, 6]]
x0 = [[0],
    [0],
    [0]];
b = [[0],
    [1],
    [14]];
[x, niters] = Method_of_Steepest_Descent(A,b,x0)
[x, niters] = CG(A,b,x0)