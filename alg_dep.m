alg_dep := function(r,n)
    v := [1 | k in [0..n]];
    for k in [2..n] do
        v[k] := ZZ! (v[k]*r);
    end for;
    p := Modulus(r);
    Qp := Parent(r);
    
    M = ZeroMatrix(n+1,n);
    for i:=1 to n-1 do
        M[1,i] := -v[i+1];
	M[i+1,i] := 1;
    end for;
    n = Precision(r);
    pn = p^n;
    pnId = pn * IdentityMatrix(ZZ,n);
    M := HorizontalJoin(M,pnId);
    M, T := HermiteForm(M);
    M, T, r := LLL(HorizontalJoin(X,pnId);
    return M[1];
end for;