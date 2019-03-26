/* Runs RankBounds on a file fin of monic curves. Picks out lower bound number of generators. Prints the result to a file fout.  */

procedure rkboundswdiv(fin,fout)
	R<x> := PolynomialRing(RationalField());
	data := eval(Read(fin));
	for curve in data do
		d:= curve[1];
		f := curve[2];
		try 
			l, u, gens := RankBounds(f,3:ReturnGenerators);
			fact := Factorization(f);
			print fact;
			factorlist := [];
			for factor in fact do
				Append(~factorlist,factor[1]);
			end for;
			seqgens :=[a : a in gens];
			sgens:=Set(seqgens);
			gens:=[s: s in sgens];
			for elt in factorlist do
				if elt in gens then
					Exclude(~gens,elt);
				end if;
			end for;
			glist := gens[1..l];
			Write(fout,Sprint(d)*":"*"["*Sprint(f)*"]"*Sprint(l)*","*Sprint(u)*"," * Sprint(glist));
		catch e
			Write(fout,Sprint(d)*":"*"["*Sprint(f)*"]"*e`Object);
		end try;
	end for;
end procedure;
