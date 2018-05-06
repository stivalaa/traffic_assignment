zin=open("512.graph.rmat.txt.0based")
zout=open("512.graph.rmat.txt","w")
for line in zin:
    sline = line.split()
    zout.write("%d\t%d\t%s\n" % ( int(sline[0])+1, int (sline[1])+1, sline[2]))

