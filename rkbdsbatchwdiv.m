/* Runs RankBounds on a file fin of monic curves of rank 1. Picks out one generator. Prints the result to a file fout.  */

/*can easily be modified to give all nontorsion gens for higher rank curves, didn't print rank because we know it's rank 1*/

procedure rkboundswdiv(fin,fout1,fout2)
	R<x> := PolynomialRing(RationalField());
	data := eval(Read(fin));
	Write(fout1,"[");
	Write(fout2,"[");

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
			Write(fout1,"["*Sprint(d)*","*Sprint(f)*"," * Sprint(l) *Sprint(u)* Sprint(gens)* "]"*",");
			if l eq u and l eq 1 then
				Write(fout2,"["*Sprint(d)*","*Sprint(f)*"," * Sprint(gens[1])* "]"*",");
			end if;
		catch e
			Write(fout1,"["*Sprint(d)*","*Sprint(f)*"]"*e`Object);
		end try;
	end for;
	Write(fout1,"]");
	Write(fout2,"]");
end procedure;
