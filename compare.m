load "coleman.m";
cc_parameters := AssociativeArray();
cc_parameters["height"] := 1000;
cc_parameters["precision"] := 20;
cc_parameters["precision_inc"] := 5;
cc_parameters["e"] := 50;
cc_parameters["e_inc"] := 10;
cc_parameters["max_prec"] := 40; 

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
    return candidates;
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
        extras := compare(f, p, cc_parameters);
        success := true;
    catch err
        new_prec := precision + precision_inc;
        new_e := e + e_inc;
        if new_prec le cc_parameters["max_prec"] then
            cc_parameters["precision"] := new_prec;
            cc_parameters["e"] := new_e;
            success, extras := compare_errors(f, p, cc_parameters);
        else
            success := false;
            extras := [p,new_prec,new_e];
        end if;
    end try;
    return success, extras;
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

function extra_points(curve, prime_list, cc_parameters)
    //new_precision := precision;
    //new_e := e;
    disc := Integers()!curve[1];
    f := curve[2];
    extras := [**];
    for p in prime_list do
        if not (disc mod p eq 0) then 
            success, results := compare_errors(f, p, cc_parameters);
            if success then 
                p_extras := results;
                for a in p_extras do
                    Append(~extras,[Integers()!a,p,Precision(a)]);
                end for;
            else
                Append(~extras,results);
            end if;
        end if;    
    end for;
    return extras;
end function;            

/*
Returns true if the point is torsion.
    
-f: degree 4 polynomial
-a: x-coord output by effective_chabauty() not corresponding 
    to rational point 
-points_height: bound for PointSearch

Let Pt = (a, a^(1/3)) in Qp. We check if the Coleman integrals
integral from infinty to Pt on basis differentials are zero

NEED TO FIX: WHAT PRECISION IS IT AFTER INTEGRATING

*/
function torsion_test(f,extras_data)
    p := extras_data[2];
    N := extras_data[3];
    data:=coleman_data(y^3-f,p,N);
    Qp := pAdicField(p,N);
    infty := set_bad_point(0,[1,0,0],true,data);
    Pt :=set_point(Qp!extras_data[1],Root(Evaluate(f,Qp!extras_data[1]),3),data);
    I := coleman_integrals_on_basis(infty,Pt,data:e:=100);
    return Valuation(I[1]) ge Max(N-5,3) and Valuation(I[2]) ge Max(N-5,3) and Valuation(I[3]) ge Max(N-5,3);
end function;

/* 
Returns true if a is padic zero of f.
*/
function ws_test(f,extras_data)
    Qp := pAdicField(extras_data[2],extras_data[3]);
    a := Qp!extras_data[1];
    return Valuation(Evaluate(f,a)) gt 0;
end function;


/*
Returns N if the point is torsion of order N.
    
-f: degree 4 polynomial
-a: x-coord output by effective_chabauty() not corresponding 
    to rational point but with x-coord rational
-points_height: bound for PointSearch
-relns_height: bound for torsion

Let P = (a,a^(1/3),1) be a point on C over K(a^(1/3)). 
We check if N(P-infty) is torsion for N<reln_height.
*/
function rational_torsion_test(f,a,relns_height)
    P2 := ProjectiveSpace(Rationals(),2);
    C := Curve(P2,Numerator(Evaluate(f,P2.1/P2.3)*P2.3^4-P2.3*P2.2^3));
    a := Rationals()!a;
    K := NumberField(x^3 - Evaluate(f,a));
    C := BaseChange(C,K);
    infty := Place(C![0,1,0]);
    a := K!a;
    P := Place(C![a,K.1,1]);
    for N in [1..relns_height] do 
        if IsPrincipal(N*P - N*infty) then
            return N;
        end if;
    end for;
    return -1;
end function;


/*
Runs tests to see if the data has torsion points or has failed, of so prints to file.

Should be run on files that have already gone through cc_file_io
*/

procedure parse_data(fin, fout)
    data := eval(Read(fin));
    for curve in data do
        c := curve[1];
        b := curve[3];
        pts := curve[4];
        if b eq false then
            Failed := true;     
        else
            Failed := false;    
        end if;
        for c in curve[4] do
            if #c ge 3 then
                t:= #c[3];
                if t ne 0 then
                    Torsion := true;
                    break c;
                else Torsion := false;
                end if;
            else Torsion := false;
            end if;
        end for;
        for c in curve[4] do
            if #c ge 4 then
                u:= #c[4];
                if u ne 0 then
                    Unexplained := true;
                    break c;
                else Unexplained := false;
                end if;
            else Unexplained := false;
            end if;
        end for;
    if Failed and Unexplained and Torsion then
        Write(fout,Sprint(curve[1])*"Failed, Torsion");
    elif Failed and Unexplained then
        Write(fout,Sprint(curve[1])*"Failed");
    elif Torsion then
        Write(fout,Sprint(curve[1])*"Torsion");
    end if;
    end for;
end procedure;

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
        if ws_test(f,data) then
            Write(fout,Sprint(data)*", Weierstrass point");
        elif torsion_test(f,data,points_height) then
            Write(fout,Sprint(data)*", torsion");
        else
            Write(fout,Sprint(data)*", unexplained");
        end if;
    end for;
end procedure;
       
function sort_cc_data(p_list,cc_extras_output)
    data := cc_extras_output;
    curve := data[1];
    f := curve[2];
    disc := Integers()!curve[1];
    P2 := ProjectiveSpace(Rationals(),2);
    C := Curve(P2,Numerator(Evaluate(f,P2.1/P2.3)*P2.3^4-P2.3*P2.2^3));
    points := PointSearch(C,1000);
    points := [Coordinates(P) : P in points];
    extras := data[2];
    test_results := [**];
    points_found := false;
    for p in p_list do
        p_results := [*p*];
        p_extras_ws := [];
        p_extras_tor := [];
        p_unexplained := [];
        p_failure := [];
        if disc mod p ne 0 then
            for cc_result in extras do
                if cc_result[2] eq p then
                    if ws_test(f,cc_result) then
                        Append(~p_extras_ws,cc_result);
                    elif torsion_test(f,cc_result) then
                        Append(~p_extras_tor,cc_result);
                    else
                        Append(~p_unexplained,cc_result);
                    end if;
                elif cc_result eq [p,45,100] then
                    Append(~p_failure,cc_result);
                end if;
            end for;
            if #p_unexplained eq 0 and #p_failure eq 0 then
                points_found := true;
            end if;
            Append(~p_results,p_extras_ws);
            Append(~p_results,p_extras_tor);
            Append(~p_results,p_unexplained);
            Append(~p_results,p_failure);
        end if;
        Append(~test_results,p_results);
    end for;
    sorted_data := [**];
    Append(~sorted_data,curve);
    Append(~sorted_data,points);
    Append(~sorted_data,points_found);
    Append(~sorted_data,test_results);
    return sorted_data;
end function;
                
procedure cc_file_io(extras_file,fout);
    cc_data := eval(Read(extras_file));
    for data in cc_data do
        sorted_data := sort_cc_data([5,7,11,13],data);
        Write(fout,Sprint(sorted_data)*",");
    end for;
end procedure;



        
procedure batch_extras(ind_start,ind_end,p_list,curves_file,fout)
    curves := eval(Read(curves_file));
    for i in [ind_start..ind_end] do
        curve := curves[i];
        results := [*curve*];
        extras := extra_points(curve,p_list,cc_parameters);
        Append(~results,extras);
        Write(fout,Sprint(results)*",");
    end for;
end procedure;  