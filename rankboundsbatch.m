/* Runs RankBounds on a file fin of curves. Prints the result to a file fout.  */

procedure rankboundsbatch(fin,fout)
	R<x> := PolynomialRing(RationalField());
	data := eval(Read(fin));
	for curve in data do
		d:= curve[1];
		f := curve[2];
		try 
			l, u := RankBounds(f,3);
			Write(fout,Sprint(d)*":"*"["*Sprint(f)*"]"*Sprint(l)*","*Sprint(u));
		catch e
			Write(fout,Sprint(d)*":"*"["*Sprint(f)*"]"*e`Object);
		end try;
	end for;
end procedure;
