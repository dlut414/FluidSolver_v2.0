name = "26"
with open(name, mode='w') as fw:
    with open(name+".out", mode='r') as fr:
        time=float(fr.readline().split()[0])
        dp=float(fr.readline().split()[0])
        np=int(fr.readline().split()[0])
        fw.write(str(time) + "\n")
        for line in fr:
            t=line.split()[0]
            x=line.split()[1]
            y=line.split()[2]
            u=line.split()[3]
            v=line.split()[4]
            fw.write(x+" "+y+" "+t+" "+u+" "+v+" 0 0 0 1.0 1.0 "+str(float(dp*3.1/2.0))+"\n")
