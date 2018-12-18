load "coleman.m";
cc_parameters := AssociativeArray();
cc_parameters["height"] := 1000;
cc_parameters["precision"] := 20;
cc_parameters["precision_inc"] := 5;
cc_parameters["e"] := 50;
cc_parameters["e_inc"] := 10;

/*
Sorts output of effective_chabauty() into known rational 
points and possibly extra points. 
*/
function compare(f, p, cc_parameters)
    height := cc_parameters["height"];
    prec := cc_parameters["precision"];
    e := cc_parameters["e"];
    data := coleman_data(y^3 - f, p, prec);
    Qpoints := Q_points(data, height);
    point_coords := [Qpoints[i]`x : i in [1..#Qpoints]];
    L,v := effective_chabauty(data:Qpoints:=Qpoints, e:=e);
    candidates := [L[i]`x : i in [1..#L]];    
    for xP in point_coords do
        for a in candidates do
            if Integers()!xP eq Integers()!a then
                Remove(~candidates, Index(candidates,a));
            end if;
        end for;
    end for;
    return point_coords, candidates;
end function;

/*
Runs compare() on increasing precision and e until 
effective_chabauty() doesn't throw an assert error.
*/

function compare_errors(f, p, cc_parameters)
    precision := cc_parameters["precision"];
    e := cc_parameters["e"];
    precision_inc := cc_parameters["precision_inc"];
    e_inc := cc_parameters["e_inc"];
    try 
        matches, extras := compare(f, p, cc_parameters);
    catch err
        cc_parameters["precision"] := precision + precision_inc;
        cc_parameters["e"] := e + e_inc;
        matches, extras := compare_errors(f, p, cc_parameters);
    end try;
    return matches, extras;
end function;





/*
Runs Chabauty and finds the extra padic points.

-f: polynomial of degree 4
-prime_list: list of good primes to use for effective_chabauty
-height: bound for PointSearch
-precision: initial precision (I use 20)
-e: (I use 50)
-precision_increment: I use between 2-5
-e_increment: I use between 5-15
-fout: directory for output file
*/

procedure extra_points(curve, prime_list, cc_parameters, fout)
    //new_precision := precision;
    //new_e := e;
    disc := Integers()!curve[1];
    f := curve[2];
    Write(fout, "[*");
    for p in prime_list do
        if not (disc mod p eq 0) then 
            p_matches, p_extras := compare_errors(f, p, cc_parameters);
            for a in p_extras do
                Write(fout, Sprint([Integers()!a,p,Precision(a)])*",");
            end for;
        end if;    
    end for;
    Write(fout, "*]");
end procedure;            

/*
Returns true if the point is torsion.
    
-f: degree 4 polynomial
-a: x-coord output by effective_chabauty() not corresponding 
    to rational point 
-points_height: bound for PointSearch
-relns_height: bound for torsion

Let P = (a,a^(1/3),1) be a point on C over K(a^(1/3)). 
We check if N(P-infty) is torsion for N<reln_height.
*/
function torsion_test(f,a,points_height,relns_height)
    P2 := ProjectiveSpace(Rationals(),2);
    C := Curve(P2,Numerator(Evaluate(f,P2.1/P2.3)*P2.3^4-P2.3*P2.2^3));
    points := PointSearch(C,points_height);
    a := Rationals()!a;
    K := NumberField(x^3 - Evaluate(f,a));
    C := BaseChange(C,K);
    places := [Place(C![Q[i] : i in [1..3]]) : Q in points];
    infty := Place(C![0,1,0]);
    a := K!a;
    P := Place(C![a,K.1,1]);
    for N in [1..relns_height] do 
        if IsPrincipal(N*P - N*infty) then
            return true;
        end if;
    end for;
    return false;
end function;

/* 
Returns true if a is padic zero of f.
*/
function ws_test(f,a)
    Qp := Parent(a);
    return Evaluate(f,a) eq Qp!0;
end function;


/*
Runs tests on points in extras_file.

f is the degree 4 polynomial for the picard curve.
extras_file is the output of extras_
points_height: bound for PointSearch
reln_height: bound for torsion search
fout: directory for output

This will test whether the x-coordinates, which are 
outputs of effective_chabauty(), are either:
 - a weierstrass point
 - a torsion point
*/
procedure test_extras(curve,extras_file,points_height,reln_height,fout)
    f := curve[2];
    extras := eval(Read(extras_file));
    for data in extras do
        Qp := pAdicField(data[2],data[3]);
        a := Qp!data[1];
        if ws_test(f,a) then
            Write(fout,Sprint(data)*", Weierstrass point");
        elif torsion_test(f,a,points_height,reln_height) then
            Write(fout,Sprint(data)*", torsion");
        else
            Write(fout,Sprint(data)*", unexplained");
        end if;
    end for;
end procedure;
        
