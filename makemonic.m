/*Makes monic models for picard curves */

procedure makemonic(fin,fout)
	R<x> := PolynomialRing(RationalField());
	data := eval(Read(fin));
	Write(fout, "[");
	for curve in data do
		d:= curve[1];
		f := curve[2];
		g := MonicModel(f,3);
		Write(fout,"["*Sprint(d)*","*Sprint(g)*"],");
	end for;
	Write(fout, "]");
end procedure;