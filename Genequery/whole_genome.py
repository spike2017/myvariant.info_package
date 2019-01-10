import myvariant
mv = myvariant.MyVariantInfo(url='http://localhost:8000/v1')
from pyliftover import LiftOver
lo = LiftOver('hg38','hg19')

#to open file
file = open("var-GS000037347-ASM.tsv","r")
fileb = open("whole_genome_404.txt","w")
#get the chrosome, position and genotype and put them into the things suitable for mv commend
#print the information in the file
count = 0
lineNumber = 0
nonrs = 0
'''
ccc = "chr7:g.140453134T>C"
testt = mv.getvariant(ccc)
#print(testt)
aa = []
aa = lo.convert_coordinate('chr12',1889144)
bb = str(aa[0])
ss = bb.split(",")
print(ss)
print("this is first line")
#test2 = mv.getvariant(ddd)
print("\nthis is first line result")
#print(test2)
'''
print ("\ntest end, real begin")
for lines in file:
    lineNumber = lineNumber +1
    #print(lines(0))
    #print("line",count)
    if(lines[0]!='#' and lines[0] != '>'):
        #file.readline()
        # for every element get them into a list
        
        wordslist = []
        words = lines.split("\t")
        for item in words:
            wordslist.append(item)
        #print(wordslist)
	#length = len(wordslist)
        if(len(wordslist)>=5):
            #rsid = wordslist[0]
            chromosome = wordslist[3]
            position = wordslist[5]
            vartype = wordslist[6]
            originalBase = wordslist[7]
            postBase = wordslist[8]
            if (vartype != 'snp'):
                continue
            count = count + 1
            chro = chromosome
            print("mark"+chro)
            print(position)
            convert = []
            position = str(position)
            #convert = lo.convert_coordinate(chro,position)
            #print(convert)
            #resultt = str(position)
            #ss = resultt.split(",")
            #cc = ss[0]
            #dd = []
            #dd = cc.split("'")
            #chromosome = dd[1]
            print(chromosome)
            #position = ss[1]
            #position = position.strip()
            print(position)
        else:
            continue	
        #queryinfo = ('\"'+chromosome+':g.'+position+originalBase+'>'+postBase+'\"')
        queryinfo = (chromosome+':g.'+position+originalBase+'>'+postBase)
        print(queryinfo)
        result = mv.getvariant(queryinfo)
        if(result != ""):
            print('get result {}'.format(result))
        if(count >10):
            break
    elif(count>0):
       	nonrs = nonrs + 1
        #fileb.write(lines) 
print("line number is",lineNumber)
print("identiable line number is",count)
print("non-identiable line number is:",nonrs)
# to close the file
file.close()
fileb.close()
