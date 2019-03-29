#takes in a file of curves on which magma has run RankBounds and finds the aborted curves, prints to files separating aborted and no valid from reasonable


def findaborted(fin,fout1,fout2):
	out1 = open(fout1,"w")
	out2 = open(fout2,"w")
	curves = open(fin,"r")
	for c in curves:
		c=c.strip()
		try:
			if str(c[-28:])== "RankBounds aborted after 120":
				c=c.replace(':',',')
				c=c.replace('[','')
				c=c.replace(']','')
				out1.write("["+c+"]"+","+"\n")
			elif c[-8:]=="elements":
				c=c.replace(':',',')
				c=c.replace('[','')
				c=c.replace(']','')
				out1.write("["+c+"]"+","+"\n")
			elif c[-6:]=="image.":
				c=c.replace(':',',')
				c=c.replace('[','')
				c=c.replace(']','')
				out1.write("["+c+"]"+","+"\n")			
			else:
				c=c.replace(':',',')
				c=c.replace('[','')
				c=c.replace(']','')
				out2.write("["+c[:-4]+"]"+","+"\n")
		except:
			continue
findaborted("master_rankboundsdata.txt","badcurves.txt","goodpicardcurves.txt")