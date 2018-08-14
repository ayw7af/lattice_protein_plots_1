from src.AAdict import aa as aa, imat as imat
def getCoord (conformation): 
    x = []
    y = []
    z = []
    for i in range(len(conformation)+1):
        x.append(0)
        y.append(0)
        z.append(0)
    for i in range(len(conformation)):
        if conformation[i]=="A":
            for j in range(i,len(conformation)):
                y[j+1] = y[j+1] - 1
        if conformation[i]=="B":
            for j in range(i,len(conformation)):
                x[j+1] = x[j+1] + 1
        if conformation[i]=="C":
            for j in range(i,len(conformation)):
                y[j+1] = y[j+1] + 1
        if conformation[i]=="D":
            for j in range(i,len(conformation)):
                x[j+1] = x[j+1] - 1
        if conformation[i]=="E":
            for j in range(i,len(conformation)):
                z[j+1] = z[j+1] + 1
        if conformation[i]=="F":
            for j in range(i,len(conformation)):
                z[j+1] = z[j+1] - 1
    return x,y,z

def energy(conformation, seq, native_list):
    
    x,y,z = getCoord(conformation)
   
    energy = 0
    contact = 0
    native_contact = 0

    for (i, j) in native_list:
        if abs(x[i]-x[j])+abs(y[i]-y[j])+abs(z[i]-z[j])==1:
            native_contact = native_contact + 1

    # Self-avoiding checks.
    for i in range(len(conformation)-1):
        for j in range(i+2, len(conformation)+1):
            if (abs(x[i]-x[j])**2 + abs(y[i]-y[j])**2 + abs(z[i]-z[j])**2)==0:
                return float('nan'), 0, 0

    for i in range(len(conformation)-2):
        for j in range(i+3, len(conformation)+1):
            if abs(x[i]-x[j])+abs(y[i]-y[j])+abs(z[i]-z[j])==1:
                energy = energy + imat[aa.index(seq[i])][aa.index(seq[j])]
                contact = contact + 1

    #energy and number of native contacts for 1 conformation
    return energy, contact, native_contact

