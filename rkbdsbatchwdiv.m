/* Runs RankBounds on a file fin of monic curves of rank 1. Picks out one generator. Prints the result to a file fout.  */

/*can easily be modified to give all nontorsion gens for higher rank curves, didn't print rank because we know it's rank 1*/

procedure rkboundswdiv(fin,fout)
	R<x> := PolynomialRing(RationalField());
	data := eval(Read(fin));
	Write(fout,"[*");
	for curve in data do
		d:= curve[1];
		f := curve[2];
		try 
			l, u, gens := RankBounds(f,3:ReturnGenerators);
			fact := Factorization(f);
			factorlist := [**];
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
			glist := gens[1];
			Write(fout,"["*Sprint(d)*","*Sprint(f)*"," * Sprint(glist)* "]"*",");
		catch e
			Write(fout,"["*Sprint(d)*","*Sprint(f)*"]"*e`Object);
		end try;
	end for;
	Write(fout,"*]");
end procedure;
