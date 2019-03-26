#takes in a file of curves on which magma has run RankBounds and finds the rank n curves, prints to a file
#you should manually remove the last comma by hand before running through magma

def findRankOne(fin,fout):
	out = open(fout,"w")
	curves = open(fin,"r")
	out.write('[')
	for c in curves:
		c=c.strip()
		try:
			if int(c[-1:])==2 and int(c[-3:-2]) ==2:
				c=c.replace(':',',')
				c=c.replace('[','')
				c=c.replace(']','')
				out.write("["+c[:-3]+"]"+","+"\n")
		except:
			continue
	out.write(']')
findRankOne("runtime_errorsoutput.txt","rank2fromrte.txt")