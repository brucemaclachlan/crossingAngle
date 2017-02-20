import numpy

### File loader ###

def getFile():
    import sys
    fileFromTerm = ''
    for arg in sys.argv[1:]:
        fileFromTerm = arg
    if len(fileFromTerm) == 0:
        print 'No file inputed'
        return None
    else:
        print '\n \nOpening file: ' + fileFromTerm
        return fileFromTerm

def readFile(file):
    if file == None:
        print 'No file to read'
    else:
        print 'Reading file: ' + str(file)
        return open(file, "r")

def terminalImport():
    IMPORTED = getFile()
    FILE = readFile(IMPORTED)
    return FILE

def fileToList(inFile):
    allLines=[]
    for line in inFile:
        allLines.append(line)
    return allLines

### PDB Atom Parser ###

def atomParser(allLines):
    atomList=[]
    elseList=[]
    for line in allLines:
        if "ATOM" in line[0:5]:
           #print lines
           atomList.append(line)
        else:
           elseList.append(line)

    print str(len(atomList)) + " atom lines detected"
    return atomList, elseList
    
def atomRowParser(atomRow):
    atomRowParsed=[]
    atomRowParsed.append(atomRow[0:6])           #0 record name 
    atomRowParsed.append(int(atomRow[6:11]))     #1 atom serial numeber
    atomRowParsed.append(atomRow[12:16])         #2 atom name
    atomRowParsed.append(atomRow[16])            #3 alternate location indicator
    atomRowParsed.append(atomRow[17:20])         #4 residue name
    atomRowParsed.append(atomRow[21])            #5 chain identifier
    atomRowParsed.append(int(atomRow[22:26]))    #6 residue sequence number
    atomRowParsed.append(atomRow[26])            #7 code for insertion of residues
    atomRowParsed.append(float(atomRow[30:38]))  #8 coordinate x
    atomRowParsed.append(float(atomRow[38:46]))  #9 coordinate y
    atomRowParsed.append(float(atomRow[46:54]))  #10 coordinate z
    atomRowParsed.append(float(atomRow[54:60]))  #11 occupancy
    atomRowParsed.append(float(atomRow[60:66]))  #12 temperature b factor
    atomRowParsed.append(atomRow[76:78])         #13 element symbol
    atomRowParsed.append(atomRow[78:-1])         #14 charge on atom
    return atomRowParsed
            
def atomMatrix(atomList):
    atomMatrix=[]
    print '\nCreating atom matrix...\n'
    for row in atomList:
        atomMatrix.append(atomRowParser(row))
    return atomMatrix
    
### Robust Cysteine pair finder ###

def SSBondParser(inFile):
    print "Finding disulphide bridges..."
    inFile
    SSBondList=[]
    for line in inFile:
        if "SSBOND" in line[0:6]:
           #print line
           parsedLine=[]
           parsedLine.append(line[15])            #Chain 1
           parsedLine.append(int(line[18:21]))    #Residue 1
           parsedLine.append(line[29])            #Chain 2
           parsedLine.append(int(line[32:35]))         #Residue 1           
           #print parsedLine
           SSBondList.append(parsedLine)
           
    print str(len(SSBondList)) + " SSBOND lines detected!\n"
    
    return SSBondList


def SSBondvalidator(SSBondList, chain):
    SSBondsinchain=[]
    for row in SSBondList:
        if row[15] == chain:
            SSBondsinchain.append(row)
    print SSBondsinchain
    if len(SSBondsinchain) == 2:
        return True
    else:
        return False
        
def TCRpairFinder(SSBondList):
    
    '''
    Takes in all the SSBonds and outputs two lists.. 
    
    the alpha chain cys and the beta chain cys
    '''
    print "Finding TCRa and TCRb cysteine pairs...\n "
    
    alphaCys=[]
    betaCys=[]
    for row in SSBondList:
        if row[0] == 'D':
            if 15 <= row[1] <= 35:
                if row[2] == 'D':
                    if 85 <= row[3] <= 100:
                        alphaCys=row
                                        
    for row in SSBondList:

        if row[0] == 'E':
            if 15 <= row[1] <= 35:
                if row[2] == 'E':
                    if 85 <= row[3] <= 100:
                        betaCys = row
            
    print alphaCys
    print betaCys
    print "\n"
    return alphaCys, betaCys
    

def TCRpairAtomFinder(TCRpairLocation, atomMatrix):
    print "Finding SG atoms..." 
    coordinate1=[]
    coordinate2=[]       
    chain=''
    cysResidue1=0
    cysResidue2=0
    chain=TCRpairLocation[0]
    cysResidue1=TCRpairLocation[1]
    cysResidue2=TCRpairLocation[3]
    
    for atoms in atomMatrix:
        if atoms[5] == chain:
            if atoms[6] == cysResidue1:
                if "SG" in atoms[2]:
                    if atoms[3] == '' or 'A':
                        coordinate1=atoms
                
    for atoms in atomMatrix:
        if atoms[5] == chain:
            if atoms[6] == cysResidue2:
                if "SG" in atoms[2]:
                    if atoms[3] == '' or 'A':
                        coordinate2=atoms
                
    return  coordinate1, coordinate2
    
    
#Cys pair Method 2#    
    
    
    
def cysPairFinder(atomMatrix, chain):
    cysteines = []
    cyspairs = []
    cysPairsSG = []
    for row in atomMatrix:   
        if row[5] == chain and "CYS" in row[4]:
            cysteines.append(row)
    for row in cysteines:
        if 15 <= row[6] <= 35 or 85 <= row[6] <= 105:
            if row[3] == '' or 'A':
                cyspairs.append(row)
    for x in cyspairs:
        print x             
    print '\nNumber of cys atoms found in chain ' + str(chain) + ': ' + str(len(cyspairs))
    for row in cyspairs:
        if 'SG' in row[2]:
            cysPairsSG.append(row)
        
    print '\nNumber of cys SG atoms found in chain ' + str(chain) + ': ' + str(len(cysPairsSG)) + "\n"
    for x in cysPairsSG:
        print x   
    print "\n"
    return cysPairsSG
    
def SGcoords(row):
    x=''
    y=''
    z=''
    SGcoords=[]

    x=row[8]
    y=row[9]
    z=row[10]
    SGcoords=[x, y, z]
    return SGcoords
    
            
### MHC coordinates finder ###

def MHCaxis(atomMatrix):
    Calphas = []
    helixCalphas = []
    for row in atomMatrix:   
        if row[5] == 'A' and "CA" in row[2]:
            Calphas.append(row)
    for row in Calphas:
        if 46 <= row[6] <= 78:
            if row[3] == '' or 'A':
                helixCalphas.append(row)
    Calphas=[]
    for row in atomMatrix:   
        if row[5] == 'B' and "CA" in row[2]:
            Calphas.append(row)
    for row in Calphas:
        if 54 <= row[6] <= 64:
            if row[3] == '' or 'A':
                helixCalphas.append(row) 
        if 67 <= row[6] <= 91:
            if row[3] == '' or 'A':
                helixCalphas.append(row)                              
        
    for x in helixCalphas:
        print x             
    print '\nNumber of Calpha atoms in MHC axis: ' + str(len(helixCalphas)) + "\n"
    return helixCalphas
                       
def MHCaxisCoords(helixCalphas):
    MHCx=[]
    MHCy=[]
    MHCz=[]
    MHCcoords=[]
    for row in helixCalphas:
        MHCx.append(row[8])
        MHCy.append(row[9])
        MHCz.append(row[10])
    MHCcoords=[MHCx, MHCy, MHCz]
    return MHCcoords         
    
def findMHCaxisCoords(atomMatrix):
    print 'Establishing MHC axis...'
    return MHCaxisCoords(MHCaxis(atomMatrix))
    
### Matrix M ###

def generateMatrixM(MHCcoords):
    matrixMlist=[]
    MHCxbar= sum(MHCcoords[0]) / float(len(MHCcoords[0]))
    MHCybar=sum(MHCcoords[1]) / float(len(MHCcoords[1]))
    MHCzbar=sum(MHCcoords[2]) / float(len(MHCcoords[2]))
    
    MHCxminusxbar=[]
    for values in MHCcoords[0]:
        MHCxminusxbar.append(values - MHCxbar)
        
    MHCyminusybar=[]
    for values in MHCcoords[1]:
        MHCyminusybar.append(values - MHCybar)
        
    MHCzminuszbar=[]
    for values in MHCcoords[2]:
        MHCzminuszbar.append(values - MHCzbar)
        
    matrixMlist = [MHCxminusxbar, MHCyminusybar, MHCzminuszbar]  
    return numpy.array(matrixMlist, dtype=float)

### Matrix A ###

def generateMatrixA(matrixM):
    print "\nGenerating matrix A...\n"
    matrixMt = matrixM.transpose()
    matrixMtt = matrixMt.transpose()
    print "Matrix A: \n"
    print numpy.dot(matrixMtt,matrixMt)
    return numpy.dot(matrixMtt,matrixMt)
    
### MHC vector ###

def generateEigens(matrixA):
    
    eigenvecval = numpy.linalg.eigh(matrixA)
    eigenvalues = eigenvecval[0]
    eigenvectors = eigenvecval[1]
    print "\nLargest Eigenvalue:"
    print eigenvalues[-1] 
    print "\nLargest Eigenvector"
    print eigenvectors[:,2]*-1
    print "\n"
    return (eigenvalues[-1], eigenvectors[:,2]*-1)
    
    
    
### Crossing angle ###

def generateCrossingAngle(l1A, l1B, l2A, l2B):
    p=   (l1B[0]-l1A[0]) * (l2B[0]-l2A[0])     +     (l1B[1]-l1A[1]) * (l2B[1]-l2A[1])     +     (l1B[2]-l1A[2]) * (l2B[2]-l2A[2]) 
    print "p = " + str(p)    
    q=((l1B[0]-l1A[0])**2) + ((l1B[1]-l1A[1])**2) + ((l1B[2]-l1A[2])**2)
    print "q = " + str(q)    
    r=((l2B[0]-l2A[0])**2) + ((l2B[1]-l2A[1])**2) + ((l2B[2]-l2A[2])**2)
    print "r = "+ str(r)    
    s=p/((q**0.5)*(r**0.5))
    print "s = "+ str(s)
    crossingAngle=(180/numpy.pi)*numpy.arccos(s)
    return crossingAngle

        
### Initialiser ###

print('     ~  Running crossingAngle.py v0.6 ALPHA  ~')
inFile = terminalImport()
inFileName = getFile()
if type(inFileName) != str:
    raise IOError('No file was loaded. Please view usage and provide a valid .pdb to be processed')
    

inFileNameProc = inFileName.split('.')
allLines=fileToList(inFile)
atomList, elseList=atomParser(allLines)
atomMatrix = atomMatrix(atomList)


#outFile = open(str(inFileNameProc[0]) + '_noh2o.pdb' , 'w')
#logFile = open(str(inFileNameProc[0]) + '_log.txt' , 'w')
#output=''


### Body ###


SSLines = SSBondParser(elseList)
TCRaCys,TCRbCys = TCRpairFinder(SSLines)
TCRaCys1,TCRaCys2 = TCRpairAtomFinder(TCRaCys, atomMatrix)
print TCRaCys1
print TCRaCys2
TCRbCys1,TCRbCys2 = TCRpairAtomFinder(TCRbCys, atomMatrix)
print TCRbCys1
print TCRbCys2

cysAlpha1=SGcoords(TCRaCys1)
cysAlpha2=SGcoords(TCRaCys2)
cysBeta1=SGcoords(TCRbCys1)
cysBeta1=SGcoords(TCRbCys2)








MHCcoords=findMHCaxisCoords(atomMatrix)
matrixM=generateMatrixM(MHCcoords)
matrixA=generateMatrixA(matrixM)
Evalue,Evector=generateEigens(matrixA)







### Cheat inputs ###

### HA 1.7  cheat input ###

cysAlpha1=[]
cysAlpha1=[-16.920,  1.00, 36.05]
cysAlpha2=[]
cysAlpha2=[-16.069,  1.00, 46.72]

#### Cys pairs ###


centroidAlpha=[]
centroidAlpha.append((cysAlpha1[0]+cysAlpha2[0])/2)
centroidAlpha.append((cysAlpha1[1]+cysAlpha2[1])/2)
centroidAlpha.append((cysAlpha1[2]+cysAlpha2[2])/2)
#print centroidAlpha


### HA 1.7  cheat input ###

cysBeta1=[]
cysBeta1=[-5.460,  1.00, 41.66]
#print cysBeta1
cysBeta2=[]
cysBeta2=[-6.060,  1.00, 43.54]
#print cysBeta2


centroidBeta=[]
centroidBeta.append((cysBeta1[0]+cysBeta2[0])/2)
centroidBeta.append((cysBeta1[1]+cysBeta2[1])/2)
centroidBeta.append((cysBeta1[2]+cysBeta2[2])/2)
#print centroidBeta
 
    
### Crossing angle calculator ###

l1A=[]
l1B=[]
l2A=[]
l2B=[]

l1A=centroidAlpha
l1B=centroidBeta
l2A=[0,0,0]
for n in Evector:
    l2B.append(n)


print l1A
print l1B
print l2A
print l2B


p=   (l1B[0]-l1A[0]) * (l2B[0]-l2A[0])     +     (l1B[1]-l1A[1]) * (l2B[1]-l2A[1])     +     (l1B[2]-l1A[2]) * (l2B[2]-l2A[2]) 

print "p = "
print p

q=((l1B[0]-l1A[0])**2) + ((l1B[1]-l1A[1])**2) + ((l1B[2]-l1A[2])**2)

r=((l2B[0]-l2A[0])**2) + ((l2B[1]-l2A[1])**2) + ((l2B[2]-l2A[2])**2)

s=p/((q**0.5)*(r**0.5))
print "s = "
print s

crossingAngle=(180/numpy.pi)*numpy.arccos(s)

print "Crossing angles = "
print crossingAngle
print 180-crossingAngle
