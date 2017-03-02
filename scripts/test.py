import sys,os

def parse(ss):
    bps = []
    stack = []
    for i,c in enumerate(ss):
        if c=="(":
            stack.append(i)
        elif c==")":
            j=stack.pop()
            bps.append((j,i))
    return bps

BPs = set([('A','U'),('U','A'),('C','G'),('G','C'),('G','U'),('U','G')])

def generateSS(n):
    if n == 0:
        return [""]
    else:
        res = ["."+s for s in generateSS(n-1)]
        for i in range(0,n-1):
            for s1 in generateSS(i):
                for s2 in generateSS(n-i-2):
                    res .append("("+s1+")"+s2)
        return res

def generate(n=0):
    if n==0:
        return [""]
    res = []
    rec = generate(n-1)
    res += ['A'+s for s in rec]
    res += ['C'+s for s in rec]
    res += ['G'+s for s in rec]
    res += ['U'+s for s in rec]
    return res

def compatible(seq,bps):
    for i,j in bps:
        if (seq[i],seq[j]) not in BPs:
            return False
    return True

if __name__ == "__main__":
    l = int(sys.argv[1])
    allstructs = generateSS(l)
    ns = len(allstructs)
    print allstructs
    numS = 0
    for i1 in range(len(allstructs)):        
        for i2 in range(i1+1,len(allstructs)):
            for i3 in range(i2+1,len(allstructs)):
                s1,s2,s3 = allstructs[i1],allstructs[i2],allstructs[i3]
                structs = [s1,s2,s3]
                os.chdir("."+os.sep+"bin")
                cmd = "."+os.sep+"RNARedPrint %s"%(" ".join(["\"%s\""%(s) for s in structs]))
                #print cmd
                sys.stdout.flush()
                os.system(cmd)
                os.chdir("..")
                bps = []
                for s in structs:
                    bps += parse(s) 
                acc = 0
                tot = 0
                
                for seq in generate(l):
                    if compatible(seq,bps):
                        acc+=1
                    tot+=1
                print "%s/%s\n"%(acc,tot)
                sys.stdout.flush()
                sys.stderr.write("Structs: %s/%s\n"%(numS+1,ns*(ns-1)*(ns-2)/6.))
                numS+=1
