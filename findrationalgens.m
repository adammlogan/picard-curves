procedure findrationalgens(fin,fout)
	R<x> := PolynomialRing(RationalField());
	data := eval(Read(fin));
	Write(fout,"[");
	for curve in data do
		d:=curve[1];
		f:=curve[2];
		gens := curve[3];
		for g in gens do
			if Degree(g) eq 1 then
				Write(fout,"["*Sprint(d)*","*Sprint(f)*"," * Sprint(g)* "]"*",");
				break g;
			end if;
		end for;
	end for;
	Write(fout,"]");
end procedure;

