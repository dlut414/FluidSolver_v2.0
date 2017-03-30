name = "26"
with open(name+".vtk", mode='w') as fw:
    fw.write("# vtk DataFile Version 3.0\n")
    fw.write("point data\n")
    fw.write("ASCII\n")
    fw.write("DATASET POLYDATA\n")
    with open(name+".out", mode='r') as fr:
        fr.readline()
        fr.readline()
        a=fr.readline()
        np=a.split()[0]
        fw.write("POINTS " + np + " float" + "\n")
        for line in fr:
            fw.write(line.split()[1]+" "+line.split()[2]+" 0\n")
        fw.write("VERTICES" + " " + np + " " + str(2*int(np)) + "\n")
        loop = 0
        while loop < int(np):
            fw.write("1 " + str(loop) + "\n")
            loop=loop+1
        fw.write("POINT_DATA " + np + "\n")
        fw.write("VECTORS U float\n")
##        fw.write("LOOKUP_TABLE default\n")
        fr.seek(0)
        fr.readline()
        fr.readline()
        fr.readline()
        for line in fr:
            fw.write(line.split()[3]+" "+line.split()[4]+" 0\n")

