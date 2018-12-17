load "coleman.m";

/*
Sorts output of effective_chabauty() into known rational 
points and possibly extra points. 
*/
function compare(f,height,p,precision,e)
    data := coleman_data(y^3 - f, p, precision);
    Qpoints := Q_points(data, height);
    L,v := effective_chabauty(data:Qpoints:=Qpoints, e:=e);
    matches := [**];
    extras := [**];
    for i in [1..#L] do
        matched := false;
        for j in [1..#Qpoints] do
            if L[i]`inf then
                if Qpoints[j]`inf then
                    Append(~matches,"infty");
                    matched := true;
                    break;
                end if;
            else
                prec1 := Precision(L[i]`x);
                prec2 := Precision(Qpoints[j]`x);
                if prec1 eq 0 or prec2 eq 0 then
                    if prec1 eq 0 and prec2 eq 0 then
                        if not Qpoints[j]`inf then
                            Append(~matches,Integers()!L[i]`x);
                            matched := true;
                            break;
                        end if;
                    end if;
                else    
                    minimum_precision := Min(prec1,prec2);
                    x_Qpoints := ChangePrecision(Qpoints[j]`x, minimum_precision);
                    x_chabauty := ChangePrecision(L[i]`x,minimum_precision);
                    if x_Qpoints eq x_chabauty then
                        Append(~matches,Integers()!(L[i]`x));
                        matched := true;
                        break;
                    end if;
                end if;
            end if;    
        end for;
        if not matched then
            Append(~extras,L[i]`x);
        end if;
    end for;
    return matches, extras;
end function;

/*
Runs compare() on increasing precision and e until 
effective_chabauty() doesn't throw an assert error.
*/

function compare_errors(f,
                        height,
                        p,
                        precision,
                        e,
                        precision_increment,
                        e_increment)
    new_precision := precision;
    new_e := e;
    try 
        extras, matches := compare(f, height, p, precision, e);
    catch err
        new_precision := precision + precision_increment;
        new_e := e + e_increment;
        extras, matches := compare_errors(f, height, p, new_precision, 
                                          new_e, precision_increment, 
                                          e_increment);
    end try;
    return extras, matches, new_precision, new_e;
end function;



procedure compare_primes(f, prime_list, height, precision, e, precision_increment, e_increment)
    new_precision := precision;
    new_e := e;
    extras := AssociativeArray();
    for p in prime_list do
        matches, p_extras, new_precision, new_e := compare_errors(f, height, p, new_precision, new_e, precision_increment, e_increment);
        extras[p] := p_extras;
    end for;
    printf "Rational points are %o\n", matches;
    for p in Keys(extras) do
        printf "%o: extras are %o\n", p, extras[p];
    end for;
end procedure;

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

procedure extra_points(f,
                        prime_list,
                        height,
                        precision,
                        e,
                        precision_increment,
                        e_increment,
                        fout)
    //new_precision := precision;
    //new_e := e;
    Write(fout, "[*");
    for p in prime_list do
        new_precision := precision;
        new_e := e;
        matches, p_extras, new_precision, new_e := 
        compare_errors(f,
                       height,
                       p,
                       new_precision,
                       new_e,
                       precision_increment,
                       e_increment);
        for a in p_extras do
            Write(fout, Sprint([Integers()!a,p,Precision(a)])*",");
        end for;
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
procedure test_extras(f,extras_file,points_height,reln_height,fout)
    extras := eval(Read(extras_file));
    for data in extras do
        Qp := pAdicField(data[2],data[3]);
        a := Qp!data[1];
        if ws_test(f,a) then
            Write(fout,Sprint(data)*", Weierstrass point");
        elif torsion_test(f,a,points_height,reln_height) then
            Write(fout,Sprint(data)*", torsion");
        end if;
    end for;
end procedure;
        
    