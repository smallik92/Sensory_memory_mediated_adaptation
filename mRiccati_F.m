function dXdt = mRiccati_F(t, X, A, B, R, Q)

% This function solves the matrix Riccati equation

X = reshape(X, size(A)); %Convert from "n^2"-by-1 to "n"-by-"n"
dXdt =( A'*X + X*A - X*B*inv(R)*B'*X + Q); %Determine derivative
dXdt = dXdt(:); %Convert from "n"-by-"n" to "n^2"-by-1    

end