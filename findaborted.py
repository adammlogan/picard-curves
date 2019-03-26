#takes in a file of curves on which magma has run RankBounds and finds the aborted curves, prints to files separating aborted from reasonable


def findaborted(fin,fout1,fout2):
	out1 = open(fout1,"w")
	out2 = open(fout2,"w")
	curves = open(fin,"r")
	out.write('[')
	for c in curves:
		c=c.strip()
		try:
			if str(c[-28:])== "RankBounds aborted after 120":
				c=c.replace(':',',')
				c=c.replace('[','')
				c=c.replace(']','')
				out1.write("["+c+"]"+","+"\n")
			elif c[-5:]=="Error":
				out2.write("["+c[:-29]+"]"+","+"\n")
			else:
				out2.write("["+c[:-3]+"]"+","+"\n")
		except:
			continue
	out.write(']')
findaborted("RankBounds_data","abortedcurves.txt","allpicardcurves.txt")