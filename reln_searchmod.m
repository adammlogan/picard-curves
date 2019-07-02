



/*
This is from the example list of extra points. We find a relation
between the point R not defined over the rationals and a
point of infinite order on our curve. 

_<x>:=PolynomialRing(Rationals());
F := x^4+x^3-3*x^2-2*x;
K<a> := NumberField(x^3+3*x^2+1);
R := [a,a^2+2*a,1];
reln_search(F,R,K,100,10);

_<x>:=PolynomialRing(Rationals());
F:= x^4 + 2*x^3 + 5*x^2 + 4*x + 1;
K<a> := NumberField(x^3 - 1/16);
R := [-1/2,a,1];
reln_search(F,R,K,100,5);
*/

/* Seaches for relations over a box of size search_bnd 
    in J(C) between a point R defined over 
    the number field L, and points in C(Q) up to 
    height. C is the Picard curve given by y^3-F=0. 
    
    
*/
/*

function reln_search(F,R,L,height,search_bnd)
    P2 := ProjectiveSpace(Rationals(),2);
    C := Curve(P2,Numerator(Evaluate(F,P2.1/P2.3)*P2.3^4-P2.3*P2.2^3));
    I := [-search_bnd..search_bnd];
    points := PointSearch(C,height);
    fin_pts := Remove(points,Index(points,C![0,1,0]));
    C := BaseChange(C,L);
    R := Place(C!R);
    fin_pts := [C![P[i] : i in [1..3]] : P in fin_pts];
    infty := Place(C![0,1,0]);
    fin_pts := [Place(P) : P in fin_pts];
    non_tor := [];
    for P in fin_pts do                 
        if IsPrincipal(3*(P-infty)) then
            continue;
        else
            Append(~non_tor,P);
        end if;
    end for;
    comps := [I : i in [1..#non_tor]];
    Append(~comps,[1..search_bnd]);
    box := CartesianProduct(comps);
    relns := [];
    for a in box do
        vector := [a[i] : i in [1..#a]];
        if Gcd(vector) eq 1 then
            D := vector[#vector]*R;
            weight := vector[#vector];
            for i in [1..#vector-1] do
                D := D + vector[i]*fin_pts[i];
                weight := weight + vector[i];
            end for;
            D := D - weight*infty;
            if IsPrincipal(D) then 
                Append(~non_tor,R);
                return vector, non_tor;
                break;
            end if;
        end if;
    end for;
    return relns;
end function;
*/


function reln_search(F,R,L,height,search_bnd)
    P2 := ProjectiveSpace(Rationals(),2);
    C := Curve(P2,Numerator(Evaluate(F,P2.1/P2.3)*P2.3^4-P2.3*P2.2^3));
    points := PointSearch(C,height);
    fin_pts := Remove(points,Index(points,C![0,1,0]));
    C := BaseChange(C,L);
    R := Place(C!R);
    fin_pts := [C![P[i] : i in [1..3]] : P in fin_pts];
    infty := Place(C![0,1,0]);
    fin_pts := [Place(P) : P in fin_pts];
    non_tor := [];
    for P in fin_pts do                 
        if IsPrincipal(3*(P-infty)) then
            continue;
        else
            Append(~non_tor,P);
        end if;
    end for;
    relns := [];
    for P in non_tor do
        for i in [-search_bnd..search_bnd] do
            for j in [1..search_bnd] do
                D := i*P + j*R - (i+j)*infty;
                if IsPrincipal(D) then
                    return [i,j],P,R;
                end if;
            end for;
        end for;
    end for;
    return relns;
end function;

/*
inputs: F is deg 4 poly in x defining picard curve y^3-F=0
        L is a number field Q(a)
        R is a point [R1,R2,R3] defined over L
        MWgen is a tuple [P1,P2,P3] st (P1:P2:P3)-inf is non-tor
        search_bnd a bound on how high you want to search for a
            relation M*R + N*MWgen - (M+N)*infty is zero in J(C) 
            where |M|,|N| < search_bnd
minimal working example:
_<x>:=PolynomialRing(Rationals());
F:= x^4 + x^3 - 3*x^2 - 2*x;
K<a> := NumberField(x^3 +3*x^2 +1);
R := [a,a^2+2*a,1];
MWgen := [-1,1,1];
reln_search(F,R,K,100,10)
            
*/

function division_test(F,R,MWgen,L,search_bnd)
    P2 := ProjectiveSpace(Rationals(),2);
    C := Curve(P2,Numerator(Evaluate(F,P2.1/P2.3)*P2.3^4-P2.3*P2.2^3));
    C := BaseChange(C,L);
    MWgen := Place(C!MWgen);
    R := Place(C!R);
    infty := Place(C![0,1,0]);
    relns := [];
    for i in [-search_bnd..search_bnd] do
        for j in [1..search_bnd] do
            D := i*R + j*MWgen - (i+j)*infty;
            if IsPrincipal(D) then
                return [i,j],R,MWgen;
            end if;
        end for;
    end for;
    return relns;
end function;