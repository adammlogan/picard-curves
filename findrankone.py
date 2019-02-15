#takes in a file of curves on which magma has run RankBounds and finds the rank one curves, prints to a file
#you should manually remove the last comma by hand before running through magma

def findRankOne(fin,fout):
	out = open(fout,"w")
	curves = open(fin,"r")
	out.write('[')
	for c in curves:
		c=c.strip()
		try:
			if int(c[-1:])==1 and int(c[-3:-2]) ==1:
				c=c.replace(':',',')
				c=c.replace('[','')
				c=c.replace(']','')
				out.write("["+c[:-3]+"]"+","+"\n")
		except:
			continue
	out.write(']')
findRankOne("weird_rankboundsoutput.txt","rank1fromwrb.txt")