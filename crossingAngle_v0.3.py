import numpy

# File loader

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

def atomParser(inFile):
    inFileYo=inFile
    atomList=[]
    for line in inFileYo:
       if "ATOM" in line:
           #print lines
           atomList.append(line)
    print str(len(atomList)) + " atom lines detected"
    return atomList
    
def atomRowParser(atomRow):
    atomRowYo=atomRow
    atomRowParsed=[]
    for x in atomRowYo.split('\t'):
        atomRowParsed.append(x)
    return atomRowParsed
        
   

### Initialiser ###

print('     ~  Running crossingAngle.py v0.3 ALPHA  ~')
inFile = terminalImport()
inFileName = getFile()
if type(inFileName) != str:
    raise IOError('No file was loaded. Please view usage and provide a valid .pdb to be processed')
    

inFileNameProc = inFileName.split('.')
atomList=atomParser(inFile)

#outFile = open(str(inFileNameProc[0]) + '_noh2o.pdb' , 'w')
#logFile = open(str(inFileNameProc[0]) + '_log.txt' , 'w')
#output=''


    



### Cheat inputs ###

cysAlpha1=[]
cysAlpha1=[6.735,15.867, 72.983]
cysAlpha2=[]
cysAlpha2=[7.266, 13.815, 73.171]


MHCx=[]
MHCy=[]
MHCz=[]


MHCx=[-6.316, -4.79, -2.975, -6.254, -7.187, -4.398, -5.769, -5.275, -2.099, -0.391, -1.07, -0.115, 3.151, 3.996, 3.112, 5.33, 8.318, 8.017, 7.79, 10.769, 12.542, 11.811, 13.451, 16.701, 16.79, 16.428, 19.372, 21.335, 20.981, 22.632, 25.27, 26.283, 26.273, 28.88, 30.863, 31.244, 31.571, 30.847, 33.681, 32.784, 29.139, 30.454, 32.453, 29.394, 27.139, 29.659, 30.134, 26.29, 25.819, 23.355, 24.224, 23.097, 19.65, 19.241, 20.193, 17.674, 15.028, 15.915, 16.272, 13.228, 10.548, 10.951, 10.46, 7.681, 5.669, 6.469, 5.024, 2.003, 1.372, 1.72, -0.924, -3.323, -2.847, -2.668]
MHCy=[-6.125, -2.693, -3.491, -4.636, -0.98, -0.519, -0.209, -2.17, -0.313, -1.068, -4.813, -5.034, -3.259, -5.436, -8.573, -7.568, -6.718, -9.944, -12.132, -10.324, -11.061, -14.799, -14.837, -13.289, -15.591, -18.513, -17.218, -16.651, -20.262, -21.379, -18.706, -19.517, -23.252, -22.476, -20.119, -22.914, -25.768, -12.907, -10.662, -11.42, -10.397, -7.273, -6.263, -6.68, -4.673, -1.888, -1.547, -1.511, 1.128, -1.038, 0.58, 4.028, 2.632, 0.825, 3.847, 6.176, 3.564, 2.821, 6.515, 7.658, 5.282, 2.574, 4.95, 6.875, 3.703, 2.415, 5.71, 5.264, 1.676, 2.741, 5.387, 3.281, 0.295, 2.299]
MHCz=[102.925, 102.632, 99.41, 97.821, 97.286, 94.715, 91.178, 87.968, 86.959, 90.219, 89.962, 86.316, 87.161, 90.088, 88.166, 85.311, 87.52, 89.521, 86.375, 84.846, 88.109, 87.817, 84.364, 85.559, 88.536, 86.045, 83.883, 87.047, 88.217, 84.918, 85.009, 88.614, 87.965, 85.276, 87.447, 90.039, 87.484, 92.84, 91.656, 88.071, 88.821, 90.612, 87.594, 85.351, 87.743, 88.179, 84.433, 84.153, 86.866, 88.684, 92.014, 90.963, 89.885, 93.103, 95.176, 93.434, 94.18, 97.832, 98.835, 96.81, 95.498, 98.128, 101.028, 99.305, 98.766, 102.234, 103.665, 101.442, 102.66, 106.369, 105.885, 103.839, 106.107, 109.338]


#### Cys pairs ###


centroidAlpha=[]
centroidAlpha.append((cysAlpha1[0]+cysAlpha2[0])/2)
centroidAlpha.append((cysAlpha1[1]+cysAlpha2[1])/2)
centroidAlpha.append((cysAlpha1[2]+cysAlpha2[2])/2)
#print centroidAlpha

cysBeta1=[]
cysBeta1=[13.909, -8.319, 65.094]
#print cysBeta1
cysBeta2=[]
cysBeta2=[12.958, -6.916, 66.234]
#print cysBeta2


centroidBeta=[]
centroidBeta.append((cysBeta1[0]+cysBeta2[0])/2)
centroidBeta.append((cysBeta1[1]+cysBeta2[1])/2)
centroidBeta.append((cysBeta1[2]+cysBeta2[2])/2)
#print centroidBeta

### MHC axis ###


#print "Averages calculated by sum over len:"

MHCxbar= sum(MHCx) / float(len(MHCx))
MHCybar=sum(MHCy) / float(len(MHCy))
MHCzbar=sum(MHCz) / float(len(MHCz))

#print MHCxbar
#print MHCybar
#print MHCzbar


MHCxminusxbar=[]
for values in MHCx:
    MHCxminusxbar.append(values - MHCxbar)
    
MHCyminusybar=[]
for values in MHCy:
    MHCyminusybar.append(values - MHCybar)
    
MHCzminuszbar=[]
for values in MHCz:
    MHCzminuszbar.append(values - MHCzbar)
    
#print "\n"    
#print MHCxminusxbar
    


#print "\n"    
#print MHCyminusybar



#print "\n"    
#print MHCzminuszbar


### Matrix multiplication ###

nestedlistmatrixM=[MHCxminusxbar, MHCyminusybar, MHCzminuszbar]
matrixM=numpy.array(nestedlistmatrixM, dtype=float)


#matrixM=numpy.array((MHCxminusxbar, MHCyminusybar, MHCzminuszbar), dtype=float)
#print "\n Numpy array:"
print matrixM
print type(matrixM)

matrixMt = matrixM.transpose()
print matrixMt
matrixMtt = matrixMt.transpose()
#print matrixMtt

matrixA=numpy.dot(matrixMtt,matrixMt)
#print type(matrixA)


print "\n Dot product: \n"

print matrixA
print type(matrixA)
print "\n"

matrixAt = matrixA.transpose()
eigenvecval = numpy.linalg.eigh(matrixA)

print "Eigenvectors 1: "
print eigenvecval[1]

eigenvecval2 = numpy.linalg.eigh(-matrixA)
print "Eigenvectors 2: "
print eigenvecval2[1]

eigenvecval3 = numpy.linalg.eigh(matrixAt)
print "Eigenvectors 3: "
print eigenvecval3[1]

eigenvecval4 = numpy.linalg.eigh(-matrixAt)
print "Eigenvectors 4: "
print eigenvecval4[1]

eigenvecval5 = numpy.linalg.eigh(-matrixAt)
print "Eigenvectors 5: "
print eigenvecval4[1]*-1

print "Eigenvalues: \n"
eigenvalues = eigenvecval[1]
print eigenvalues



print "Eigenvectors: \n"
eigenvectors = eigenvecval[1]
print eigenvectors


print "\n Largest Eigenvalue \n"
print eigenvalues[-1]
print "\n Largest Eigenvector \n"
print eigenvectors[:,2]
print "\n"
 
    
### Crossing angle calculator ###

l1A=[]
l1B=[]
l2A=[]
l2B=[]


l1A=centroidAlpha
l1B=centroidBeta
l2A=[0,0,0]
for n in eigenvectors[:,2]:
    l2B.append(n*-1)


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
