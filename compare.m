function compare(f,height,p,precision,e)
    load "coleman.m";
    data := coleman_data(y^3 - f, p, precision);
    Qpoints := Q_points(data, height);
    L,v := effective_chabauty(data:Qpoints:=Qpoints, e:=e);
    matches := [];
    extras := [];
    for i in [1..#L)] do
        for i in [1..Len(Qpoints)] do
            if L[i]`inf then
                if Qpoints[i]`inf then
                    Append(~matches,(0:1:0));
                else
                    continue;
                end if;
            else
                if L[i]`x == Qpoints[i]`x then
                    Append(~matches,(L[i]`x,Qpoints[i]`x));
                    continue;
                else
                    Append(~extras,L[i]`x);
                    continue;
                end if;
           end if;
        end for;
    end for;
    return matches, extras;
end function;
                