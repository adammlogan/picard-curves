/* Runs RankBounds on a file fin of monic curves. Picks out lower bound number of generators. Prints the result to a file fout.  */

procedure rankboundsbatch(fin,fout)
	R<x> := PolynomialRing(RationalField());
	data := eval(Read(fin));
	for curve in data do
		d:= curve[1];
		f := curve[2][1];
		try 
			l, u, gens := RankBounds(f,3:ReturnGenerators);
			fact := Factorization(f);
			factorlist = [**];
			for factor in fact do
				Append(~factorlist,factor[1]);
			end for;
			for x in fact do
				if x in gens then
					Exclude(~S,x);
				end if;
			end for;
			sgens:=Set(gens);
			gens:=[s: s in sgens];
			glist := gens[1..l];
			Write(fout,Sprint(d)*":"*"["*Sprint(f)*"]"*Sprint(l)*","*Sprint(u)*"," * Sprint(glist));
		catch e
			Write(fout,Sprint(d)*":"*"["*Sprint(f)*"]"*e`Object);
		end try;
	end for;
end procedure;
