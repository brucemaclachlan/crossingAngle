import numpy
### File loader ###

def getFile():
    import sys
    fileFromTerm = ''
    args=sys.argv
    fileFromTerm = args[1]
    if len(fileFromTerm) == 0:
        print 'No file inputed'
        return None
    else:
        print '\n \nOpening file: ' + fileFromTerm
        return fileFromTerm
        
def MHCclass():
    import sys
    classFromTerm = ''
    args=sys.argv
    classFromTerm=args[2]
    if len(classFromTerm) == 0:
        print 'No MHC class information inputed. Please see usage'
        return None
    else:
        print '\n \nInput flagged as MHC Class ' + classFromTerm
        return classFromTerm
        
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
    print "Done!"
    return atomMatrix
    
### Robust Cysteine pair finder ###

def SSBondParser(inFile):
    print "Finding disulphide bridges..."
    inFile
    SSBondList=[]
    for line in inFile:
        if "SSBOND" in line[0:6]:
           parsedLine=[]
           parsedLine.append(line[15])            #Chain 1
           parsedLine.append(int(line[18:21]))    #Residue 1
           parsedLine.append(line[29])            #Chain 2
           parsedLine.append(int(line[32:35]))    #Residue 1           
           SSBondList.append(parsedLine)
    print str(len(SSBondList)) + " SSBOND lines detected!\n"
    if len(SSBondList)>0:
        return SSBondList
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
                    if 85 <= row[3] <= 110:
                        alphaCys=row
        if row[0] == 'D':
            if 15 <= row[3] <= 35:
                if row[2] == 'D':
                    if 85 <= row[1] <= 110:
                        alphaCys=row                                        
    for row in SSBondList:
        if row[0] == 'E':
            if 15 <= row[1] <= 35:
                if row[2] == 'E':
                    if 85 <= row[3] <= 110:
                        betaCys = row
        if row[0] == 'E':
            if 15 <= row[3] <= 35:
                if row[2] == 'E':
                    if 85 <= row[1] <= 110:
                        betaCys = row          
    print alphaCys
    print betaCys
    print "\n"
    if len(alphaCys) > 0 and len(betaCys) > 0:
        return alphaCys, betaCys
    else:
        return False, False
    
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
    
def cysRobust(inFile, atomMatrix):
    SSLines=SSBondParser(inFile)
    if SSLines==False:        
        print "No SSBOND lines were detected in the PDB!\n"
        return False, False, False, False                  
    a,b=TCRpairFinder(SSLines)
    if a == False or b == False:
        print "Could not find TCR cysteine pairs in the SSBOND list!\n"
        return False, False, False, False
    A1atoms, A2atoms=TCRpairAtomFinder(a, atomMatrix)
    B1atoms, B2atoms=TCRpairAtomFinder(b, atomMatrix)
    if A1atoms == False or A2atoms == False or B1atoms == False or B2atoms == False:
        print "Could not find TCR cysteine atoms from the cysteine pairs listed in SSBOND list!\n"
        return False, False, False, False 
    return A1atoms, A2atoms, B1atoms, B2atoms    
    
# Weak Cys pair finder #    
    
def cysPairFinder(atomMatrix, chain):
    print "Finding CYS residues between residues 15 to 35 and between 85 to 105..."
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
    print '\nNumber of cys atoms found in chain ' + str(chain) + ' around the S-S bridge region: ' + str(len(cyspairs))
    for row in cyspairs:
        if 'SG' in row[2]:
            cysPairsSG.append(row)    
    print '\nNumber of cys SG atoms found in chain ' + str(chain) + ' around the S-S bridge region: ' + str(len(cysPairsSG)) + "\n"
    for x in cysPairsSG:
        print x   
    print "\n"
    return cysPairsSG
    
def cysWeak(atomMatrix):
    a = cysPairFinder(atomMatrix, 'D')
    a1,a2 = a[0], a[1]
    b = cysPairFinder(atomMatrix, 'E')
    b1,b2 = b[0], b[1]
    return a1, a2, b1, b2
            
### MHC coordinates finder ###

def MHCaxisI(atomMatrix):

    print "\nMHC axis atoms are:\n"
    Calphas = []
    helixCalphas = []
    for row in atomMatrix:   
        if row[5] == 'A' and "CA" in row[2]:
            Calphas.append(row)
    for row in Calphas:
        if 50 <= row[6] <= 86 or 140 <= row[6] <= 176:
            if row[3] == '' or 'A':
                helixCalphas.append(row)
    for x in helixCalphas:
        print x             
    print '\nNumber of Calpha atoms in MHC axis: ' + str(len(helixCalphas)) + "\n"
    return helixCalphas

def MHCaxisII(atomMatrix):
    print "\nMHC axis atoms are:\n"
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
    print "Extracting MHC axis coordinates..."
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
    
def findMHCaxisCoords(atomMatrix,MHCclass):
    
    if MHCclass == "I":
        print 'Establishing MHC I axis...'
        return MHCaxisCoords(MHCaxisI(atomMatrix))
    if MHCclass == "II":
        print 'Establishing MHC II axis...'
        return MHCaxisCoords(MHCaxisII(atomMatrix))        
            
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
    p=   (l1B[0]-l1A[0])*(l2B[0]-l2A[0])+(l1B[1]-l1A[1]) * (l2B[1]-l2A[1])+(l1B[2]-l1A[2])*(l2B[2]-l2A[2]) 
    print "p = " + str(p)    
    q=((l1B[0]-l1A[0])**2) + ((l1B[1]-l1A[1])**2) + ((l1B[2]-l1A[2])**2)
    print "q = " + str(q)    
    r=((l2B[0]-l2A[0])**2) + ((l2B[1]-l2A[1])**2) + ((l2B[2]-l2A[2])**2)
    print "r = "+ str(r)    
    s=p/((q**0.5)*(r**0.5))
    print "s = "+ str(s)
    crossingAngle=(180/numpy.pi)*numpy.arccos(s)
    return crossingAngle
  
########### BODY    ###############
  
### Initialiser ###

print('     ~  Running crossingAngle_panClass.py v0.9 ALPHA  ~')
inFile = terminalImport()
inFileName = getFile()
if type(inFileName) != str:
    raise IOError('No file was loaded. Please view usage and provide a valid .pdb to be processed')
inFileNameProc = inFileName.split('.')
allLines=fileToList(inFile)
atomList, elseList=atomParser(allLines)
atomMatrix = atomMatrix(atomList)


#outFile = open(str(inFileNameProc[0]) + '_crossingAngle.txt' , 'w')
#logFile = open(str(inFileNameProc[0]) + '_crossingAngle_log.txt' , 'w')
#output=''


### Cysteine axis ###

cysWARNING=False
TCRaCys1,TCRaCys2,TCRbCys1,TCRbCys2=cysRobust(elseList, atomMatrix)
if TCRaCys1==False or TCRaCys2==False or TCRbCys1==False or TCRbCys2 == False:
    cysWARNING=True
    TCRaCys1,TCRaCys2,TCRbCys1,TCRbCys2 = cysWeak(atomMatrix)
print "\nTCRa and TCRb cysteine pair SG atoms are:\n"
print TCRaCys1
print TCRaCys2
print TCRbCys1
print TCRbCys2
cysAlpha1=SGcoords(TCRaCys1)
cysAlpha2=SGcoords(TCRaCys2)
cysBeta1=SGcoords(TCRbCys1)
cysBeta2=SGcoords(TCRbCys2)
print "\nTCRa and TCRb cysteine pair SG coordinates are [x,y,z]:\n"
print cysAlpha1
print cysAlpha2
print cysBeta1
print cysBeta1
centroidAlpha=[]
centroidAlpha.append((cysAlpha1[0]+cysAlpha2[0])/2)
centroidAlpha.append((cysAlpha1[1]+cysAlpha2[1])/2)
centroidAlpha.append((cysAlpha1[2]+cysAlpha2[2])/2)
centroidBeta=[]
centroidBeta.append((cysBeta1[0]+cysBeta2[0])/2)
centroidBeta.append((cysBeta1[1]+cysBeta2[1])/2)
centroidBeta.append((cysBeta1[2]+cysBeta2[2])/2)
print "\nTCRa cysteine centroid is [x,y,z]:"
print centroidAlpha
print "\nTCRa cysteine centroid is [x,y,z]:"
print centroidBeta
print "\n\n"

# MHC axis #
MHCclass=MHCclass()
MHCcoords=findMHCaxisCoords(atomMatrix,MHCclass)
matrixM=generateMatrixM(MHCcoords)
matrixA=generateMatrixA(matrixM)
Evalue,Evector=generateEigens(matrixA)
    
# Crossing angle calculator #

l1A=[]
l1B=[]
l2A=[]
l2B=[]
l1A=centroidAlpha
l1B=centroidBeta
l2A=[0,0,0]
for n in Evector:
    l2B.append(n)
print "Line 1A:"
print l1A
print "\nLine 1B:"
print l1B
print "\nLine 2A:"
print l2A
print "\nLine 2B:"
print l2B
print "\nCalculating crossing angle...\n"
p=   (l1B[0]-l1A[0]) * (l2B[0]-l2A[0])     +     (l1B[1]-l1A[1]) * (l2B[1]-l2A[1])     +     (l1B[2]-l1A[2]) * (l2B[2]-l2A[2]) 
q=((l1B[0]-l1A[0])**2) + ((l1B[1]-l1A[1])**2) + ((l1B[2]-l1A[2])**2)
r=((l2B[0]-l2A[0])**2) + ((l2B[1]-l2A[1])**2) + ((l2B[2]-l2A[2])**2)
s=p/((q**0.5)*(r**0.5))
crossingAngle=(180/numpy.pi)*numpy.arccos(s)

### ANSWER ###
print "#########    ANSWER    ############"
print "Crossing angles = "
print crossingAngle
print 180-crossingAngle

print "\nIf the TCR binds using conventional polarity, then the smallest number is your crossing angle. If reverse polarity as described in Beringer et al. Nature Immunology 2015, the larger value is the crossing angle. \
Later versions will be written to account for TCR docking polarity.\n"

if cysWARNING==True:
    print "WARNING!!!!!.. this value was calculated by attempting to find the disulphide bridge in each TCR chain and may cause errors if there are CYS residues within the CDR3 region.\n \
    Please check the output to confirm that the 'TCRa and TCRb cysteine pair SG atoms' are indeed forming a disulphide bridge in your structure."
    
if cysWARNING==False:
    print "Calculated using the robust Cys detection method"

print('\n     ~  End of crossingAngle_panClass.py v0.9 ALPHA  ~')
