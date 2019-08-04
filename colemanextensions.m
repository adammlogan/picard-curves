// The following functions are modifications from the coleman.m
// package. Main added funcitonality: suppose C is a curve such that
// a divisor D in C(Q) is non-torsion, but C does not have any 
// points Pi,Pj in C(Q) such that Pi-Pj is non-torsion. Suppose
// L is an extension of Q such that the points in the support of D
// are in C(L). Suppose p is a prime which splits completely in L.
// We compute p-adic vanishing differential using D by integrating
// from a base point b in C(Q) to each point in the support of D 
// and summing up. 

// Most changes from the coleman.m package are to accommodate 
// divisors, created by set_bad_div. The output of set_bad_div
// is a list of x-coords instead of a single x-coord. The 
// associated functions (div_local_data, div_local_coord) 
// are the functions local_data, local_coord
// except that they return lists, one item for each
// x-coord in the divisor set by set_bad_div. 


// Notes to Sachi: set_bad_div is kept, we also use set_bad_point. 
// reason: set_bad_point is needed for p-adic information,
    // we use it to set points in the support of a divisor
    // and integrate to each of these points when we want to
    // integrate from P to a divisor D. Seems OK to me to 
    // keep both functions around - parts of the code are purely
    // p-adic, in which case there isn't harm in distinguishing
    // between point coming from Q point or point coming from
    // a rational divisor whose support is defined over L when
    // p splits completely in L, since then L embeds into Qp.
    

first_split_prime:=function(f)
  // returns first rational prime which splits completely
  // in splitting field of f
  // 
  
  L := SplittingField(f);
  OL := MaximalOrder(f);
  p := 5;
  while not IsTotallySplit(p,OL) do
      p := NextPrime(p);
  end while;
  return p;
end function;
  




mat_W0_basechange:=function(Q,L)

  // Compute the matrix W0 using MaximalOrderFinite over 
  // the number field L

  
  Qt:=RationalFunctionField(L);
  K:=ext<Qt|Q>;
  b0:=Basis(MaximalOrderFinite(K));
  d:=Degree(Q);
  mat:=ZeroMatrix(Qt,d,d);
  for i:=1 to d do
    for j:=1 to d do
      mat[i,j]:=Eltseq(K!b0[i])[j];
    end for;
  end for;
  return mat;

end function;



div_is_bad:=function(P,data)
  // check whether the point P is bad
  // changes from is_bad: performs evaluation of r at
  // all x-coordinates of points in support of P
  x0:=P`x; r:=data`r;
  if P`inf then // infinite point
    return true;
  end if;    
  
  for i:=1 to #x0 do
    if Valuation(Evaluate(r,x0[i])) gt 0 then // finite bad point
      return true;
    end if; 
  end for;
  return false;
  
end function;

div_is_very_bad:=function(P,data)

  // check whether the point P is very bad
  // changes from is_bad: - performs evaluation of r (if finite) or
  // checks valuation (if infinite) at all x-coordinates of points in support of P
  x0:=P`x; r:=data`r; N:=data`N;
  for i:=1 to #x0 do
      if P`inf then // infinite point
          if Valuation(x0[i]) ge N then // infinite very bad point
              return true;
          end if;
      else // finite point
          if Valuation(Evaluate(r,x0[i])) ge N then // finite very bad point
              return true;
          end if;
      end if;
  end for;

  return false;

end function;


set_bad_div:=function(x,b,inf,data,L,prime_over_p,divisor)
  // changes from set_bad_point: also stores the number field L 
  // such that P is defined over L and boolean divisor. In addition,
  // x is now a list of x-coordinates of points in support of P,
  // b is now a list of b-values for points in support of P
  // TO DO: replace Trace with proper coercion from ext of Qp to Qp

  Q:=data`Q; p:=data`p; N:=data`N; 
  Qp:=pAdicField(p,N); d:=Degree(Q);

  format:=recformat<x,b,inf,xt,bt,index,divisor,L>;
  P:=rec<format|>;
  support := [];
  P`inf:=inf;
  P`divisor := divisor;
  P`L := L;
  if divisor then
    Lp, embed := Completion(L,prime_over_p:Precision:=N);
    P`b := [[Qp!Trace(embed(b[i][j])) : j in [1..d]] : i in [1..#x]];
    for i in [1..#x] do
      support := Append(support,Qp!Trace(embed(x[i])));
    end for;
    P`x := support;
  else
    P`x:=[Qp!x];
    P`b:=[[Qp!b[i]:i in [1..d]]];
  end if;
  

  return P; 

end function;


vanishing_differentials_with_divs:=function(points,data:e:=1);

  // Compute the regular one forms of which the 
  // integrals vanish between all points in points.
  // NEW: If D is a divisor in points, then computes
  // integrals by integrating from base point to each 
  // point in the suport of D and summing.

  Q:=data`Q; p:=data`p;
  
  g:=genus(Q,p);

  IP1Pi:=[];
  NIP1Pi:=[];
  P1 := set_bad_point(points[1]`x[1],points[1]`b[1],points[1]`inf,data);
  for i:=1 to #points-1 do
    if not points[i+1]`divisor then
      xi := (points[i+1]`x)[1];
      bi := (points[i+1]`b)[1];
      P := set_bad_point(xi,bi,points[i+1]`inf,data);
      Ci,Ni:=coleman_integrals_on_basis(P1,P,data:e:=e);
      IP1Pi[i]:=Ci;
      NIP1Pi[i]:=Ni;
    else
      x1 := (points[i+1]`x)[1];
      b := (points[i+1]`b)[1];
      P := set_bad_point(x1,b,false,data);
      Ci, Ni := coleman_integrals_on_basis(P1,P,data:e:=e);
      for j in [2..#points[i+1]`x] do
        xj := (points[i+1]`x)[j];
        b := (points[i+1]`b)[j];
        P := set_bad_point(xj,b,false,data);
        Cij, Nij := coleman_integrals_on_basis(P1,P,data:e:=e);
        Ci := Ci + Cij;
        Ni := Minimum(Ni,Nij);
      end for;
      IP1Pi[i] := Ci;
      NIP1Pi[i] := Ni;
    end if;
  end for;

  Nint:=Minimum(NIP1Pi);

  K:=pAdicField(p,Nint);
  M:=ZeroMatrix(K,g,#points-1);
  for i:=1 to g do
    for j:=1 to #points-1 do
      M[i,j]:=K!reduce_mod_pN_Q(Rationals()!IP1Pi[j][i],p,Nint);
    end for;
  end for;

  v:=basis_kernel(M);

  return v,IP1Pi,NIP1Pi;
  
end function;
Q_divs:=function(data,bound,f)

  // Returns a list (not guaranteed to be complete) of Q-rational points 
  // upto height bound on the curve given by data, along with a 
  // divisor on the curve whose support is generated by the roots 
  // of f in a splitting field L of f. 

  Q:=data`Q; p:=data`p; N:=data`N; r:=data`r; W0:=data`W0; Winf:=data`Winf;
  d:=Degree(Q);

  
  pointlist:=[];

  A2:=AffineSpace(RationalField(),2);
  Qxy:=PolynomialRing(RationalField(),2);
  
  QA2:=Qxy!0;
  C:=Coefficients(Q);
  for i:=1 to #C do
    D:=Coefficients(C[i]);
    for j:=1 to #D do
      QA2:=QA2+D[j]*(Qxy.1)^(j-1)*(Qxy.2)^(i-1);
    end for;
  end for;
  
  X:=Scheme(A2,QA2);
  pts:=PointSearch(X,bound);
  xvalues:=[];
  for i:=1 to #pts do
    if not pts[i][1] in xvalues then
      xvalues:=Append(xvalues,pts[i][1]);
    end if;
  end for;

  Qx:=RationalFunctionField(RationalField());
  Qxy:=PolynomialRing(Qx);
  Qfun:=Qxy!0;
  C:=Coefficients(Q);
  for i:=1 to #C do
    D:=Coefficients(C[i]);
    for j:=1 to #D do
      Qfun:=Qfun+D[j]*(Qx.1)^(j-1)*(Qxy.1)^(i-1);
    end for;
  end for;
  FF:=FunctionField(Qfun);

  b0fun:=[];
  for i:=1 to d do
    bi:=FF!0;
    for j:=1 to d do
      bi:=bi+W0[i,j]*FF.1^(j-1);
    end for;
    b0fun[i]:=bi;
  end for;

  binffun:=[];
  for i:=1 to d do
    bi:=FF!0;
    for j:=1 to d do
      bi:=bi+Winf[i,j]*FF.1^(j-1);
    end for;
    binffun[i]:=bi;
  end for;
  
  divisor := false;
  L := RationalField();
  prime_over_p := p;
  for i:=1 to #xvalues do
    places:=Decomposition(FF,Zeros(Qx.1-xvalues[i])[1]);
    if Valuation(xvalues[i],p) ge 0 then
      for j:=1 to #places do
        if Degree(places[j]) eq 1 then
          x:=xvalues[i];
          b:=[];
          for k:=1 to d do
            b[k]:=Evaluate(b0fun[k],places[j]);
          end for;
          P:=set_bad_div(x,b,false,data,L,prime_over_p,divisor);
          pointlist:=Append(pointlist,P);
        end if;
      end for;
    else
      for j:=1 to #places do
        if Degree(places[j]) eq 1 then
          x:=1/xvalues[i];
          b:=[];
          for k:=1 to d do
            b[k]:=Evaluate(binffun[k],places[j]);
          end for;
          P:=set_bad_div(x,b,true,data,L,prime_over_p,divisor);
          pointlist:=Append(pointlist,P);
        end if;
      end for;
    end if;
  end for; 

  places:=InfinitePlaces(FF);
  for i:=1 to #places do
    if Degree(places[i]) eq 1 then
      x:=0;
      b:=[];
      for j:=1 to d do
        b[j]:=Evaluate(binffun[j],places[i]);
      end for;
      P:=set_bad_div(x,b,true,data,L,prime_over_p,divisor);
      pointlist:=Append(pointlist,P);
    end if;
  end for; 

     
  L,roots := SplittingField(f);
  OL := MaximalOrder(L);
  p_splitting := Factorization(p*OL);
  prime_over_p := p_splitting[1][1];
  Qx := RationalFunctionField(L);
  Qxy:=PolynomialRing(Qx);
  Qfun:=Qxy!0;
  C:=Coefficients(Q);
  for i:=1 to #C do
    D:=Coefficients(C[i]);
    for j:=1 to #D do
      Qfun:=Qfun+D[j]*(Qx.1)^(j-1)*(Qxy.1)^(i-1);
    end for;
  end for;
  FF:=FunctionField(Qfun);
  W0 := mat_W0_basechange(Q,L);
  b0fun:=[];
  for i:=1 to d do
    bi:=FF!0;
    for j:=1 to d do
      bi:=bi+FF!W0[i,j]*FF.1^(j-1);
    end for;
    b0fun[i]:=bi;
  end for;


  support := [];
  b_list := [];
  divisor := true;
  for root in roots do
    places:=Decomposition(FF,Zeros(Qx.1-root)[1]);
    if Valuation(root,prime_over_p) ge 0 then
      for j:=1 to #places do
        if Degree(places[j]) eq 1 then
          support := Append(support,root);
          b:=[];
          for k:=1 to d do
            b[k]:=Evaluate(b0fun[k],places[j]);
          end for;
          b_list := Append(b_list,b);
        end if;
      end for;
    end if;
  end for;
  P:=set_bad_div(support,b_list,false,data,L,prime_over_p,divisor);
  pointlist:=Append(pointlist,P);
  return pointlist;

end function;

div_local_data:=function(P,data)

  // For a point P, returns the ramification index of the map x on the residue disk at P
  // If P corresponds to a divisor, returns the ramification
  // indices for the map x on the points in the support of the 
  // map x.
  
  Q:=data`Q; p:=data`p; W0:=data`W0; Winf:=data`Winf; x0:=P`x; b:=P`b; d:=Degree(Q);
  L := P`L;
  ep := [];
  index := [];
  if not div_is_bad(P,data) then
    for i in [1..#P`x] do
      eP:=Append(ep,1);
      index := Append(index,0);
    end for;
    return eP,index;
  else
    Fp:=FiniteField(p); Fpx:=RationalFunctionField(Fp); Fpxy:=PolynomialRing(Fpx);
    f:=Fpxy!0;
    for i:=0 to d do
      for j:=0 to Degree(Coefficient(Q,i)) do
        f:=f+(Fp!Coefficient(Coefficient(Q,i),j))*Fpxy.1^i*Fpx.1^j;
      end for;
    end for;  
    FFp:=FunctionField(f); // function field of curve mod p
    bmodp := [];
    eP :=[];
    for i in [1..#x0] do
      if P`inf then
        places:=InfinitePlaces(FFp); // infinite places of function field of curve mod p
        W:=Winf;
      else
        Px0:=Zeros(Fpx.1-Fp!x0[i])[1]; 
        places:=Decomposition(FFp,Px0); // places of function field of curve mod p lying over x0 mod p
                                      
        if not P`divisor then
          W:=W0;
        else
          W := mat_W0_basechange(Q,L);
        end if;
      end if;

      bmodpi:=[]; // elements of b mod p, where b is either b^0 or b^inf
      for j:=1 to d do
        f:=FFp!0;
        for k:=1 to d do
          f:=f+(Fpx!W[j,k])*FFp.1^(k-1);
        end for;
        bmodpi[j]:=f;
      end for;
      bmodp := Append(bmodp,bmodpi);
      done:=false;

      for j:=1 to #places do
        same:=true;
        for k:=1 to d do
          if Evaluate(bmodp[i][k],places[j]) ne Fp!b[i][k] then
            same:=false;
          end if;
        end for;    
        if same then
          place:=places[j];
          done:=true;
        end if;
      end for;

      if not done then
        error "Point does not lie on curve";
      end if;
      eP:=RamificationIndex(place);

      if eP eq 1 then
        index[i]:=0;
      else
        done:=false;
        m:=1;
        while not done do
          ord:=Valuation(bmodp[i][m]-Evaluate(bmodp[i][m],place),place);
          if ord eq 1 then
            index[i]:=m;
            done:=true;
          end if;
          m:=m+1;
        end while;
      end if;

    end for;
  return eP,index,place,bmodp;
  end if;

end function;


div_local_coord:=function(P,prec,data);

  // Computes powerseries expansions of x and
  // the b^0_i or b^infty_i (depending on whether
  // P is infinite or not) in terms of the local
  // coordinate computed by local_data.
  // If P corresponds to a divisor, returns a 
  // list of powerseries expensions of x and b^0_i
  // or b^infty_i for each point in the support of P.
      
  if assigned P`xt then 
      xt := P`xt;
      for i:=1 to #xt do
          if Precision(Parent(P`xt[i])) ge prec then
              xt:=P`xt;
              bt:=P`bt;
              index:=P`index;
              return xt,bt,index;
          end if;
      end for;              
  end if;

  if div_is_bad(P,data) and not div_is_very_bad(P,data) then
    error "Cannot compute local parameter at a bad point which is not very bad";
  end if;

  x0:=P`x; Q:=data`Q; p:=data`p; N:=data`N; Winf:=data`Winf;
  d:=Degree(Q); b:=P`b;
  if P`divisor then
    L := P`L;
    W0 := data`W0; // is this ok? 
  else
    W0:=data`W0;
  end if;
  L := P`L;
  K:=Parent(x0[1]); Kt<t>:=PowerSeriesRing(K,prec); Kty:=PolynomialRing(Kt);
  Qt:=RationalFunctionField(L); Qty:=PolynomialRing(Qt);
  Fp:=FiniteField(p);

  f:=Qty!0;
  for i:=0 to d do
    for j:=0 to Degree(Coefficient(Q,i)) do
      f:=f+Coefficient(Coefficient(Q,i),j)*Qty.1^i*Qt.1^j;
    end for;
  end for;  
  FF:=FunctionField(f); // function field of curve

  if not div_is_bad(P,data) then // finite good point
    xt := [];
    bt := [];
    index:=[];
    for i :=1 to #x0 do
        
      xt := Append(xt,t+x0[i]);

      W0invx0:=Transpose(Evaluate(W0^(-1),x0[i]));
      ypowers:=Vector(b[i])*W0invx0;
      y0:=ypowers[2];

      C:=Coefficients(Q);
      D:=[];
      for j:=1 to #C do
        D[j]:=Evaluate(C[j],xt[i]); 
      end for;
      fy:=Kty!D;
      derfy:=Derivative(fy);

      yt:=hensel_lift(fy,Kt!y0);

      ypowerst:=[];
      ypowerst[1]:=FieldOfFractions(Kt)!1;
      ypowerst[2]:=yt;
      for i:=3 to d do
        ypowerst[i]:=ypowerst[i-1]*yt;
      end for;
      bti:=Eltseq(Vector(ypowerst)*Transpose(Evaluate(W0,xt[i])));
      btnew:=[];
      for j:=1 to d do
        btnew[j]:=Kt!bti[j];
      end for;
      bt:=Append(bt,btnew);

      index[i]:=0;
    end for;
  elif P`inf then // infinite point

    eP,index,place,bmodp:=div_local_data(P,data);
    FFp:=Parent(bmodp[1][1]);
    Fpx:=BaseRing(FFp);

    bfun:=[];
    for i:=1 to d do
      bi:=FF!0;
      for j:=1 to d do
        bi:=bi+Winf[i,j]*FF.1^(j-1);
      end for;
      bfun[i]:=bi;
    end for;
    
    if eP eq 1 then // P is an infinite point that is not ramified
      
      xt:=[t+x0[i] : i in [1..#x0]];
      bt:=[];
      for k := 1 to #x0 do
        btk := [];
        for i:=1 to d do

          if assigned data`minpolys and data`minpolys[2][1,i+1] ne 0 then
            poly:=data`minpolys[2][1,i+1]; 
          else 
            poly:=minpoly(FF!(1/Qt.1),bfun[i]);
          end if;

          C:=Coefficients(poly);
          D:=[];
          for j:=1 to #C do
            D[j]:=Evaluate(C[j],xt[k]); 
          end for;
          fy:=Kty!D;
          derfy:=Derivative(fy);

          modpprec:=mod_p_prec(fy);

          if assigned P`bt and Precision(Parent(P`bt[k][i])) ge modpprec then
            approxroot:=P`bt[k][i];
          else
            tmodp:=1/Fpx.1-Fp!x0;
            expamodp:=mod_p_expansion(bmodp[k][i],place,tmodp,modpprec);
            approxroot:=approx_root(fy,b[k][i],modpprec,expamodp);
          end if;

          btki:=hensel_lift(fy,approxroot);
          btk[i]:=btki;

        end for;
        bt := Append(bt,btk);
      end for;
    else // P is an infinite point that is ramified
      bt := [];
      xt := [];
      for k := 1 to #x0 do
        if assigned data`minpolys and data`minpolys[2][index[k]+1,1] ne 0 then
          poly:=data`minpolys[2][index[k]+1,1];
        else
          poly:=minpoly(bfun[index[k]],FF!1/(Qt.1));
        end if;

        C:=Coefficients(poly);

        D:=[];
        for j:=1 to #C do
          D[j]:=Evaluate(C[j],t+b[k][index[k]]); 
        end for;
        fy:=Kty!D;
        derfy:=Derivative(fy);

        modpprec:=mod_p_prec(fy);

        if assigned P`xt and Precision(Parent(P`xt[k])) ge modpprec then
          approxroot:=P`xt[k];
        else
          tmodp:=bmodp[k][index[k]]-Fp!b[k][index[k]];
          expamodp:=mod_p_expansion(FFp!1/Fpx.1,place,tmodp,modpprec);
          approxroot:=approx_root(fy,x0[k],modpprec,expamodp);
        end if;

        xtk:=hensel_lift(fy,approxroot);

        btk:=[];
        for i:=1 to d do 

          if i eq index[k] then
            btk[i]:=t+b[k][i];
          else

            if assigned data`minpolys and data`minpolys[2][index[k]+1,i+1] ne 0 then
              poly:=data`minpolys[2][index[k]+1,i+1];
            else
              poly:=minpoly(bfun[index[k]],bfun[i]);
            end if;

            C:=Coefficients(poly);
            D:=[];
            for j:=1 to #C do
              D[j]:=Evaluate(C[j],t+b[k][index[k]]); 
            end for;

            fy:=Kty!D;
            derfy:=Derivative(fy);

            modpprec:=mod_p_prec(fy);

            if assigned P`bt and Precision(Parent(P`bt[k][i])) ge modpprec then
              approxroot:=P`bt[k][i];
            else
              tmodp:=bmodp[k][index[k]]-Fp!b[k][index[k]];
              expamodp:=mod_p_expansion(bmodp[k][i],place,tmodp,modpprec);
              approxroot:=approx_root(fy,b[k][i],modpprec,expamodp);
            end if;

            btki:=hensel_lift(fy,approxroot);
            btk[i]:=btki;

          end if;

        end for;
        xt := Append(xt,xtk);
        bt := Append(bt,btk);
      end for;
          
    end if;

  else // finite bad point

    eP,index,place,bmodp:=div_local_data(P,data); 
    FFp:=Parent(bmodp[1][1]);
    Fpx:=BaseRing(FFp);

    bfun:=[];
    for i:=1 to d do
      bi:=FF!0;
      for j:=1 to d do
        bi:=bi+W0[i,j]*FF.1^(j-1);
      end for;
      bfun[i]:=bi;
    end for;

    if eP eq 1 then // P is a finite point that is not ramified
      xt := [t + x0[k] : k in [1..#x0]];
      bt := [];
      for k := 1 to #x0 do
        bt:=[];
        for i:=1 to d do

          if assigned data`minpolys and data`minpolys[1][1,i+1] ne 0 then
            poly:=data`minpolys[1][1,i+1];
          else
            poly:=minpoly(FF!Qt.1,bfun[i]);
          end if;

          C:=Coefficients(poly);
          D:=[];
          for j:=1 to #C do
            D[j]:=Evaluate(C[j],xt[k]); 
          end for;
          fy:=Kty!D;
          derfy:=Derivative(fy);

          modpprec:=mod_p_prec(fy);

          if assigned P`bt and Precision(Parent(P`bt[k][i])) ge modpprec then
            approxroot:=P`bt[k][i];
          else
            tmodp:=Fpx.1-Fp!x0;
            expamodp:=mod_p_expansion(bmodp[k][i],place,tmodp,modpprec);
            approxroot:=approx_root(fy,b[k][i],modpprec,expamodp);
          end if;

          btki:=hensel_lift(fy,approxroot);
          bt[k][i]:=btki;

        end for;
      end for;  
    else // P is a finite point that ramifies
      xt := [];
      bt := [];
      for k:=1 to #x0 do
        if assigned data`minpolys and data`minpolys[1][index[k]+1,1] ne 0 then
          poly:=data`minpolys[1][index[k]+1,1];
        else
          poly:=minpoly(bfun[index[k]],FF!Qt.1);
        end if;

        C:=Coefficients(poly);
        D:=[];
        for j:=1 to #C do
          D[j]:=Evaluate(C[j],t+b[k][index[k]]); 
        end for;
        fy:=Kty!D;
        derfy:=Derivative(fy);

        modpprec:=mod_p_prec(fy);

        if assigned P`xt and Precision(Parent(P`xt[k])) ge modpprec then
          approxroot:=P`xt[k];
        else
          tmodp:=bmodp[k][index[k]]-Fp!b[k][index[k]];
          expamodp:=mod_p_expansion(FFp!Fpx.1,place,tmodp,modpprec);
          approxroot:=approx_root(fy,x0,modpprec,expamodp);
        end if;

        xtk:=hensel_lift(fy,approxroot);  
        xt := Append(xt,xtk);
        btk:=[];
        for i:=1 to d do 

          if i eq index[k] then
            btk[i]:=t+b[k][i];
          else

            if assigned data`minpolys and data`minpolys[1][index[k]+1,i+1] ne 0 then
              poly:=data`minpolys[1][index[k]+1,i+1];
            else
              poly:=minpoly(bfun[index[k]],bfun[i]);
            end if;

            C:=Coefficients(poly);
            D:=[];
            for j:=1 to #C do
              D[j]:=Evaluate(C[j],t+b[k][index[k]]);
            end for;

            fy:=Kty!D;

            derfy:=Derivative(fy);

            modpprec:=mod_p_prec(fy);

            if assigned P`bt and Precision(Parent(P`bt[k][i])) ge modpprec then
              approxroot:=P`bt[k][i];
            else
              tmodp:=bmodp[k][index[k]]-Fp!b[k][index[k]];
              expamodp:=mod_p_expansion(bmodp[k][i],place,tmodp,modpprec);
              approxroot:=approx_root(fy,b[k][i],modpprec,expamodp);
            end if;

            btki:=hensel_lift(fy,approxroot);
            btk[i]:=btki;
            bt := Append(bt,btk);

          end if;

        end for;
      end for;
    end if;

  end if;

  return xt,bt,index;

end function;



effective_chabauty_with_Qdiv:=function(data:Qpoints:=[],bound:=0,e:=1);

  // Carries out effective Chabauty for the curve given by data.
  // First does a point search up to height bound. Then uses the
  // points found to determine the vanishing differentials. 
  // NEW: uses vanishing_differentials_with_divs to compute
  // vanishing differentials for 1-forms which vanish on .                    
  // Finally
  // goes over all residue disks mapping to points on the reduction
  // mod p and finds all common zeros of the vanishing differentials.
 
  if #Qpoints eq 0 then
    if bound eq 0 then
      error "have to specify either Qpoints or a bound for search";
    end if;
    Qpoints:=Q_points(data,bound);
  end if; 
   
  for i:=1 to #Qpoints do
    _,index:=div_local_data(Qpoints[i],data);
    for j:=1 to #Qpoints[i]`x do
      data:=update_minpolys(data,Qpoints[i]`inf,index[j]);
      if div_is_bad(Qpoints[i],data) then
        if div_is_very_bad(Qpoints[i],data) then
          xt,bt,index:=div_local_coord(Qpoints[i],tadicprec(data,e),data);
          Qpoints[i]`xt:=xt;
          Qpoints[i]`bt:=bt;
          Qpoints[i]`index:=index; 
        end if;
      else
        xt,bt,index:=div_local_coord(Qpoints[i],tadicprec(data,1),data);
        Qpoints[i]`xt:=xt;
        Qpoints[i]`bt:=bt;
        Qpoints[i]`index:=index; 
      end if;
    end for;
  end for;

  v,IP1Pi,NIP1Pi:=vanishing_differentials_with_divs(Qpoints,data:e:=e);

  Qppoints,data:=Qp_points(data:points:=Q_points(data,1000));
  for i:=1 to #Qppoints do
    if is_bad(Qppoints[i],data) then
      xt,bt,index:=local_coord(Qppoints[i],tadicprec(data,e),data);
    else
      xt,bt,index:=local_coord(Qppoints[i],tadicprec(data,1),data);
    end if;
    Qppoints[i]`xt:=xt;
    Qppoints[i]`bt:=bt;
    Qppoints[i]`index:=index;
  end for;

  pointlist:=[];
  for i:=1 to #Qppoints do
    k:=0;
    for j:=1 to #Qpoints do
      for m := 1 to #Qpoints[j]`x do
        if (Qppoints[i]`x eq Qpoints[j]`x[m]) and (Qppoints[i]`b eq Qpoints[j]`b[m]) and (Qppoints[i]`inf eq Qpoints[j]`inf) then
          k:=j;
        end if;
      end for;
    end for;
    x1 := (Qpoints[1]`x)[1];
    b := (Qpoints[1]`b)[1];
    inf := Qpoints[1]`inf;
    P1 := set_bad_point(x1,b,Qpoints[1]`inf,data);

    if k lt 2 then
      pts:=zeros_on_disk(P1,Qppoints[i],v,data:e:=e);
    else
      pts:=zeros_on_disk(P1,Qppoints[i],v,data:e:=e,integral:=[*IP1Pi[k-1],NIP1Pi[k-1]*]);
    end if;
    for j:=1 to #pts do
      pointlist:=Append(pointlist,pts[j]);
    end for;
  end for;

  return pointlist, v;

end function;


/** Computes integrals from infinity to each point in Qpoints and ...
    Lpoints 
**/
    

compute_integrals := function(data,rat_points,L,Lpoints,e)
  Q:=data`Q; p:=data`p; N:=data`N; r:=data`r; W0:=data`W0; Winf:=data`Winf;
  d:=Degree(Q);

  Qp := pAdicField(p,N);
  
  xvalues:=[];
  for i:=1 to #rat_points do
    if not rat_points[i][1] in xvalues then
      xvalues:=Append(xvalues,rat_points[i][1]);
    end if;
  end for;

  Qx:=RationalFunctionField(RationalField());
  Qxy:=PolynomialRing(Qx);
  Qfun:=Qxy!0;
  C:=Coefficients(Q);
  for i:=1 to #C do
    D:=Coefficients(C[i]);
    for j:=1 to #D do
      Qfun:=Qfun+D[j]*(Qx.1)^(j-1)*(Qxy.1)^(i-1);
    end for;
  end for;
  FF:=FunctionField(Qfun);

  b0fun:=[];
  for i:=1 to d do
    bi:=FF!0;
    for j:=1 to d do
      bi:=bi+W0[i,j]*FF.1^(j-1);
    end for;
    b0fun[i]:=bi;
  end for;

  binffun:=[];
  for i:=1 to d do
    bi:=FF!0;
    for j:=1 to d do
      bi:=bi+Winf[i,j]*FF.1^(j-1);
    end for;
    binffun[i]:=bi;
  end for;
  
  pointlist := [];
  /* make point at infinity, assuming there's only 1 */
  x:=0;
  b:=[];
  places:=InfinitePlaces(FF);
  for j:=1 to d do
      b[j]:=Evaluate(binffun[j],places[1]);
  end for;
  infty:=set_bad_point(x,b,true,data);
  pointlist := Append(pointlist,infty);
  

  for i:=1 to #xvalues do
    places:=Decomposition(FF,Zeros(Qx.1-xvalues[i])[1]);
    if Valuation(xvalues[i],p) ge 0 then
      for j:=1 to #places do
        if Degree(places[j]) eq 1 then
          x:=xvalues[i];
          b:=[];
          for k:=1 to d do
            b[k]:=Evaluate(b0fun[k],places[j]);
          end for;
          P:=set_bad_point(x,b,false,data);
          pointlist:=Append(pointlist,P);
        end if;
      end for;
    else
      for j:=1 to #places do
        if Degree(places[j]) eq 1 then
          x:=1/xvalues[i];
          b:=[];
          for k:=1 to d do
            b[k]:=Evaluate(binffun[k],places[j]);
          end for;
          P:=set_bad_point(x,b,true,data);
          pointlist:=Append(pointlist,P);
        end if;
      end for;
    end if;
  end for; 




  OL := MaximalOrder(L);
  p_splitting := Factorization(p*OL);
  prime_over_p := p_splitting[1][1];
  Qx := RationalFunctionField(L);
  Qxy:=PolynomialRing(Qx);
  Qfun:=Qxy!0;
  C:=Coefficients(Q);
  for i:=1 to #C do
    D:=Coefficients(C[i]);
    for j:=1 to #D do
      Qfun:=Qfun+D[j]*(Qx.1)^(j-1)*(Qxy.1)^(i-1);
    end for;
  end for;
  FF:=FunctionField(Qfun);
  W0 := mat_W0_basechange(Q,L);
  b0fun:=[];
  for i:=1 to d do
    bi:=FF!0;
    for j:=1 to d do
      bi:=bi+FF!W0[i,j]*FF.1^(j-1);
    end for;
    b0fun[i]:=bi;
  end for;


  Lp, embed := Completion(L,prime_over_p:Precision:=N);
  for P in Lpoints do
    x := P[1];
    places:=Decomposition(FF,Zeros(Qx.1-x)[1]);
    if Valuation(x,prime_over_p) ge 0 then
      for j:=1 to #places do
        if Degree(places[j]) eq 1 then
          b:=[];
          for k:=1 to d do
            b[k]:=Evaluate(b0fun[k],places[j]);
          end for;
          x := Qp!Trace(embed(x));
          b := [Qp!Trace(embed(b[j])) : j in [1..d]];
          Append(~pointlist,set_bad_point(x,b,false,data));
        end if;
      end for;
    end if;
  end for;
      
  for i:=1 to #pointlist do
    _,index:=local_data(pointlist[i],data);
    data:=update_minpolys(data,pointlist[i]`inf,index);
    if is_bad(pointlist[i],data) then
        if is_very_bad(pointlist[i],data) then
            xt,bt,index:=local_coord(pointlist[i],tadicprec(data,e),data);
            pointlist[i]`xt:=xt;
            pointlist[i]`bt:=bt;
            pointlist[i]`index:=index; 
        end if;
    else
        xt,bt,index:=local_coord(pointlist[i],tadicprec(data,1),data);
        pointlist[i]`xt:=xt;
        pointlist[i]`bt:=bt;
        pointlist[i]`index:=index; 
    end if;
  end for;
  
      
  IP1Pi := [];
  NIP1Pi := [];
  for i:=1 to #pointlist-1 do
      Ci,Ni:=coleman_integrals_on_basis(pointlist[1],pointlist[i+1],data:e:=e);
      IP1Pi[i]:=Ci;
      NIP1Pi[i]:=Ni;
  end for;
  return IP1Pi;
  
end function;