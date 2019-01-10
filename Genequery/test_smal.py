from pyliftover import LiftOver
import myvariant
import math
import re
import csv
import datetime


class query:
    # start
    mv = myvariant.MyVariantInfo(url='http://myvariant.info/v1')
    '''
    parameters = input('Enter parameters :')
    filetype = input('file type :')
    hgversion = #input('NCBI snp version (eg. :38): ')
    lo = hgVersionJudge(hgversion)
    '''
    def __init__(self,fileName, outputfile_name, type, outputtype, version = 19, fieldKeyWords = 'None', lineBegin = -1, lineEnd = math.inf):
        self.parameters = fileName
        self.outputfile_name = outputfile_name
        self.hgversion = version
        self.filetype = type
        self.outputtype = outputtype
        self.lo = self.hgVersionJudge(version)
        self.fieldKeyWords = fieldKeyWords
        self.lineBegin = lineBegin
        self.lineEnd = lineEnd

    def infosearch(self,line):
        p = re.compile('\bthis\b')
        print(p.search('no class at all'))
        print(re.search(line))

    def hgVersionJudge(self, nowVersion):
        if (int(nowVersion) != 19):
            strs = 'hg' + str(nowVersion)
            lo = LiftOver(strs, 'hg19')
            return lo
        else:
            return 0

    def whole_genomeProcessor(self, wordslist, hgVersionNow, lo):
        if (len(wordslist) < 8):
            print(len(wordslist))
            return 'no'
        chromosome = wordslist[3]
        position = wordslist[5]
        vartype = wordslist[6]
        originalBase = wordslist[7]
        postBase = wordslist[8]
        if (vartype != 'snp'):
            return 'no'
        chro = chromosome
        print("mark" + chro)
        print(position)
        print(chromosome)
        queryinfos = (chromosome + ':g.' + position + originalBase + '>' + postBase)
        return queryinfos

    # def hgVersion_ChrPosConvert(self,lo, hgVersionNow,chro,position):
    #    return
    def AncestryAndmeProcessor(self,wordslist, hgVersionNow, lo):
        rsid = wordslist[0]
        # chromosome = wordslist[1]
        # position = wordslist[2]
        # genotype = wordslist[3]
        # originalBase = genotype[0]
        # postBase = genotype[1]
        queryinfom = (rsid)
        return queryinfom

    def vcfFileProcessor(self, wordslist, hgVersionNow, lo):
        # rsid = wordslist[0]
        chromosome = wordslist[0]
        position = wordslist[1]
        # genotype = wordslist[3]
        originalBase = wordslist[3]
        postBase = wordslist[4]
        chro = "chr" + chromosome
        print("mark" + chro)
        print(position)
        convert = []
        position = int(position)
        if (hgVersionNow == 19):
            convert = lo.convert_coordinate(chro, position)
            print(convert)
            resultt = str(convert[0])
            ss = resultt.split(",")
            cc = ss[0]
            dd = []
            dd = cc.split("'")
            chromosome = dd[1]
            print(chromosome)
            position = ss[1]
            position = position.strip()
            print(ss)
        position = str(position)
        queryinfo = 'chr' + chromosome + ':g.' + position + originalBase + '>' + postBase
        return queryinfo

    def expansion(self, dict1, dict0, key1):
        #print(type(dict1))
        if (isinstance(dict1, dict)):
            for key2 in dict1.keys():
                #print('keys: ' + key2)
                if (key1 != ''):
                    key3 = str(key1) + '.' + str(key2)
                else:
                    key3 = str(key2)
                self.expansion(dict1.get(key2), dict0, key3)
        else:
            dict0[key1] = str(dict1)

    def queries(self, queryinfo,wordslist):
        if (self.fieldKeyWords == 'None'):
        # print('130: ' + queryinfo)
            result_collection = self.mv.query(queryinfo)
        else:
        # print('133: ' + queryinfo + ',fields=' + self.fieldKeyWords)
            result = self.mv.query(queryinfo, fields=self.fieldKeyWords)
        if (self.outputtype == 'vcf'):
            output = wordslist[1] + '\t' + wordslist[2] + '\t' + wordslist[0] + '\t' + wordslist[3][0] + '\t' + \
                     wordslist[3][1] + '\t.\t.\t'
        elif (self.outputtype == 'csv'):
            print(result)
            result = result['hits'][0]



    def genequery(self):
        # get the chrosome, position and genotype and put them into the things suitable for mv commend
        # print the information in the file
        # to open file
        result = ''
        file = open(self.parameters, "r")
        fileb = open(self.outputfile_name, "w", encoding='utf-8')
        count = 0
        lineNumber = 0
        nonrs = 0
        title = {}
        outputs = []
        multiple_query = []
        list_of_wordlist = []
        # print("\ntest end, real begin")
        if (self.outputtype != 'csv') :
            fileb.write('##fileformat=VCFv4.1 \n##fileDate=' + str(datetime.datetime.now()) + '\n##version=hg19 \n##CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO')
        for lines in file:
            lineNumber = lineNumber + 1
            if (lines[0] != '#' and lines[0] != '>' and lineNumber >= self.lineBegin and lineNumber <= self.lineEnd):
                count = count + 1
                wordslist = []
                words = lines.split("\t")
                queryinfo = ''
                output = ''
                for item in words:
                    wordslist.append(item)
                if (len(wordslist) >= 1):
                    if (self.filetype == 'vcf'):
                        queryinfo = self.vcfFileProcessor(wordslist, self.hgversion, self.lo)
                        if (self.fieldKeyWords == 'None') :
                            print('140 : queryInfo' + self.filetype)
                            print('141: ' + queryinfo)
                            result = self.mv.getvariant(queryinfo)
                        else:
                            result = self.mv.getvariant(queryinfo ,fields =self.fieldKeyWords)
                    elif (self.filetype == '23andme' or self.filetype == 'ancestry'):
                        queryinfo = self.AncestryAndmeProcessor(wordslist, self.hgversion, self.lo)
                        #print(queryinfo)
                        multiple_query.append(queryinfo)
                        list_of_wordlist.append(wordslist)
                    elif (self.filetype == 'whole_genome'):
                        queryinfo = self.whole_genomeProcessor(wordslist, self.hgversion, self.lo)
                        if (self.fieldKeyWords == 'None') :
                            result = self.mv.getvariant(queryinfo)
                        else:
                            result = self.mv.getvariant(queryinfo ,fields =self.fieldKeyWords)
                    else:
                        print('possible type : vcf, 23andme, ancestry, whole_genome')
                        break
                else:
                    continue
                if (queryinfo == 'no'):
                    continue
                #print(queryinfo)
                if (result != "" and self.outputtype =='vcf'):
                    lineNumber = lineNumber
                    #print(str(lines) + '\t' + str(result))
                    #fileb.write(str(output)+ '\t' + str(result) + '\n' )
                elif (result != "" and self.outputtype == 'csv') :
                    test_output = {}
                    self.expansion(result, test_output, '')
                    print(type(test_output))
                    print(type(outputs))
                    outputs.append(test_output)
                    for keys in test_output.keys() :
                        if (test_output.get(keys, 'a') != 'a'):
                            title[keys] = 1
                    #w = csv.DictWriter(fileb, test_output.keys())
                    #print(test_output)
                    #w.writeheader()
                    #w.writerow(test_output)
                if (count > 2000):
                    break
            elif (count >= 0):
                nonrs = nonrs + 1
                if (self.outputtype != 'csv' and (lines[0] == '#' or lines[0] == '>')):
                    #fileb.write(lines)
                    1 + 1
        if (self.outputtype == 'csv'):
            print('197')
            all_keys = title.keys()
            dict_writer = csv.DictWriter(fileb, title)
            dict_writer.writeheader()
            dict_writer.writerows(outputs)
        print("line number is: ", lineNumber)
        print("identiable line number is: ", count)
        print("non-identiable line number is: ", nonrs)
        # to close the file
        file.close()
        fileb.close()
        return multiple_query

import myvariant
import time
mv = myvariant.MyVariantInfo(url='http://myvariant.info/v1')
list2 = []
demo = query('23andme_small.txt','test2.csv','23andme','vcf', 19,'None',29,999)
list2 = demo.genequery()
print('len of list2 is : ',len(list2))
#print(list2)
start_time = time.time()
result2 = mv.querymany(list2, 'dbsnp.rsid',returnall=True)
print("test 1 --- %s seconds ---" % (time.time() - start_time))
start_time = time.time()
#for count in range(0,999) :
    #result2 = mv.query('rs58991260')
print("test 2 --- %s seconds ---" % (time.time() - start_time))

#print(result2)




#'emv.clinvar_rcv, emv.egl_classification, emv.egl_protein, emv.egl_variant, emv.gene, emv.hgvs'
