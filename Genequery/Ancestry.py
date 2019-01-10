import myvariant
mv = myvariant.MyVariantInfo(url='http://myvariant.info/v1')

#to open file
file = open("Ancestry_small.txt","r")
fileb = open("AncestryDNA_404.txt","w")
#get the chrosome, position and genotype and put them into the things suitable for mv commend
#print the information in the file
count = 0
lineNumber = 0
nonrs = 0
for lines in file:
    lineNumber = lineNumber +1
    #print("line",count)
    if(lines[0]=='r' and lines[1]=='s'):
        #file.readline()
        # for every element get them into a list
        count = count +1
        wordslist = []
        words = lines.split("\t")
        for item in words:
            wordslist.append(item)
        #print(wordslist)
        rsid = wordslist[0]
        #chromosome = wordslist[1]
        #position = wordslist[2]
        #genotype = wordslist[3]
        #originalBase = genotype[0]
        #postBase = genotype[1]
        queryinfo = ("dbsnp.rsid:"+ rsid)
        #print(queryinfo)
        result = mv.query(queryinfo,fields='all')
        print(result)
        if(count >10):
        	break
    elif(count>0):
       	nonrs = nonrs + 1
        fileb.write(lines) 
print("line number is",lineNumber)
print("rsid number is",count)
print("non-rs line is:",nonrs)
# to close the file
file.close()
fileb.close()
