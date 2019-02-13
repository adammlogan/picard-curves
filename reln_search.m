
/* 
    fin_places: a list of the non-infinite rational places
    infty: the place above infinity = [0,1,0]
    Q: point on our curve not defined over rationals
    m: bound on the box size you want to search over
    
    This creates the box [-m,m]*...*[-m,m]*[1,m] and 
    checks whether  a_i*P_i +i*Q - w*infty is principal,
    where a_i in [-m,m],
    i in [1,m], 
    and w = i + sum(a_i).
    
    It returns a list of all tuples of relations. 
    
*/
function reln_search(fin_places,infty,Q,m)
    I := [-m..m];
    comps := [I : i in [1..#fin_places]];
    Append(~comps,[1..m]);
    box := CartesianProduct(comps);
    relns := [];
    for a in box do
        vector := [a[i] : i in [1..#a]];
        if Gcd(vector) eq 1 then
            D := vector[#vector]*Q;
            weight := vector[#vector];
            for i in [1..#vector-1] do
                D := D + vector[i]*fin_places[i];
                weight := weight + vector[i];
            end for;
            D := D - weight*infty;
            if IsPrincipal(D) then 
                Append(~relns,vector);
            end if;
        end if;
    end for;
    return relns;
end function;


/*
This is from the example list of extra points. We find a relation
between the point Q not defined over the rationals and a
point of infinite order on our curve. One remark: in this
example, I removed the nontorsion points from my list of 
rational points to speed up the search just for this example. 
*/

f := x^4+x^3-3*x^2-2*x;
K<a> := NumberField(x^3+3*x^2+1);
P2 := ProjectiveSpace(Rationals(),2);
C := Curve(P2,Numerator(Evaluate(f,P2.1/P2.3)*P2.3^4-P2.3*P2.2^3));
points := PointSearch(C,100);
C := BaseChange(C,K);
Q := Place(C![a,a^2+2*a,1]);
fin_pts := points[2..#points];
fin_pts := [C![P[i] : i in [1..3]] : P in fin_pts];
infty := Place(C![points[1][i]: i in [1..3]]);
fin_pts := [Place(P) : P in fin_pts];
non_tor := [fin_pts[2],fin_pts[4]];
reln_search(non_tor,infty,Q,3);