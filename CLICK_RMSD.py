import math

f = open("str1_1ok8-str2_1svb.1.pdb")
g = open("str2_1svb-str1_1ok8.1.pdb")
h = open("str1_1ok8-str2_1svb.pdb.1.clique")
i = open("ResRMSD.txt", "xt")

PDB1ResNum = []
PDB2ResNum = []
PDB1ResName = []
PDB2ResName = []
PDB1AtomList = []
PDB2AtomList = []
ResNumList = []
x1 = []
x2 = []
y1 = []
y2 = []
z1 = []
z2 = []
coordinate_distance = []
coordinate_distanceSUM = 0
RMSD = 0

for x in range (0,7):
    h.readline()
for line in h:
    ResList = line.split()
    PDB1ResNum.append(int (ResList[1]))
    PDB2ResNum.append(int (ResList[5]))

for res in PDB1ResNum:
    for line in f:
        PDB1AtomList = line.split()
        if res == int(PDB1AtomList[5]) and PDB1AtomList[2] == 'CA':
            PDB1ResName.append(PDB1AtomList[3])
            x1.append(float(PDB1AtomList [6]))
            y1.append(float(PDB1AtomList [7]))
            z1.append(float(PDB1AtomList [8]))
            f.seek(0,0)
            break

for res in PDB2ResNum:
    for line in g:
        PDB2AtomList = line.split()       
        if res == int(PDB2AtomList[5]) and PDB2AtomList[2] == 'CA':
            PDB2ResName.append(PDB2AtomList[3])
            x2.append(float(PDB2AtomList [6]))
            y2.append(float(PDB2AtomList [7]))
            z2.append(float(PDB2AtomList [8]))
            g.seek(0,0)
            break
        

#collect per residue RMSD
for n in range (0, len(x1)):
    coordinate_distance.append( math.sqrt (((x1[n]-x2[n])**2) + ((y1[n]-y2[n])**2) + ((z1[n]-z2[n])**2)))

#calculate total RMSD
coordinate_sub = []
for n in range (0, len(x1)):
    coordinate_sub.append(((x1[n]-x2[n])**2) + ((y1[n]-y2[n])**2) + ((z1[n]-z2[n])**2))
    coordinate_distanceSUM = coordinate_distanceSUM + coordinate_sub[n]

RMSD = math.sqrt (coordinate_distanceSUM/len(x1))


#print table
i.write('''Below is a table summarizing the matched atoms and their distance. The columns are as follows:
Bracket, PDB1ResNum[n], PDB1ResName[n], x1[n], y1[n], z1[n], PDB2ResNum[n], PDB2ResName[n], x2[n], y2[n], z2[n], coordinate_distance[n], Bracket
''')
  
for n in range (0, len(x1)):
    table1 = ('',PDB1ResNum[n], PDB1ResName[n], x1[n], y1[n], z1[n], PDB2ResNum[n], PDB2ResName[n], x2[n], y2[n], z2[n], coordinate_distance[n],'')
    stringed1= str(table1)
    print (stringed1)
    i.write(stringed1 + '\n')
print (RMSD)
q = str(RMSD)

i.write('The RMSD is ' + q + '\n \n \n')

#protlen = max(ResNumList)

full_coordinate_distance = coordinate_distance

for line in f:
    PDB1AtomList = line.split()
    if PDB1AtomList[2] == 'CA':
        ResNumList.append(PDB1AtomList[5])
    
for x in ResNumList:
    if int(x) not in PDB1ResNum:
        full_coordinate_distance.insert(int(x)-1,15)
        
i.write('''Below is a table with coordinate distances for matched residues and a
default distance of 15 for unmatched residues. The columns are:
Bracket, PDB1ResNum, full_coordinate_distance, Bracket \n''')
for n in range (0, len(ResNumList)):
    table2 = ('', ResNumList[n], full_coordinate_distance[n],'')
    stringed2= str(table2)
    print (stringed2)
    i.write(stringed2 + '\n')

            

i.close 

