load "coleman.m";
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

function compare_errors(f, height, p, precision, e, precision_increment, e_increment)
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
    p_matches := [**];
    p_extras := [**];
    for p in prime_list do
        matches, extras, new_precision, new_e := compare_errors(f, height, p, new_precision, new_e, precision_increment, e_increment);
        Append(~p_extras,[*p,extras*]);
    end for;
    printf "Rational points are %o\n", matches;
    for pair in p_extras do
        printf "%o: extras are %o\n", pair[1], pair[2];
    end for;
end procedure;
            
            