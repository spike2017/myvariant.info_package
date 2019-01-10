#from guppy import hpy
#import time
#fileb = open('test.txt','w',encoding = 'utf-8')
#fileb.write("%s"%time.time())

import myvariant
import time
import os
import psutil
from pyliftover import LiftOver
import math
import re
import csv
import datetime
import pandas as pd

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

    def AncestryAndmeProcessor_vcf_title(self,wordslist, hgVersionNow, lo):
        # CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
        rsid = wordslist[0]
        chromosome = wordslist[1]
        position = wordslist[2]
        genotype = wordslist[3]
        originalBase = genotype[0]
        postBase = genotype[1]
        return (chromosome + '\t' + position + '\t' + rsid + '\t' + originalBase + '\t' + postBase )

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
        elif (isinstance(dict1,list)):
            for item in dict1 :
                self.expansion(item,dict0,key1)
        else :
            dict0[key1] = str(dict1)


    def genequery(self):
        # get the chrosome, position and genotype and put them into the things suitable for mv commend
        # print the information in the file
        # to open file
        #print('line133 success')
        file = open(self.parameters, "r")
        fileb = open(self.outputfile_name, "a", encoding='utf-8')
        count = 0
        lineNumber = 0
        nonrs = 0
        title = {}
        outputs = []
        queryinfo_list = []
        vcfinfo_list = []
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
                result = ''
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
                        #print(lineNumber)
                        queryinfo_list.append(queryinfo)
                        if (self.outputtype == 'vcf'):
                            titleinfo = self.AncestryAndmeProcessor_vcf_title(wordslist, self.hgversion, self.lo)
                            vcfinfo_list.append(titleinfo)
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
                #got a query list but no result
                #if (count > 20):
                    #break
            elif (count >= 0):
                nonrs = nonrs + 1
                #if (self.outputtype != 'csv' and (lines[0] == '#' or lines[0] == '>')):
                    #fileb.write(lines)
        #for loop end here

        #return queryinfo_list

        if (self.outputtype =='vcf'):
            count = 0
            dul_count = 0
            test_out = ''
            results = self.mv.querymany(queryinfo_list, 'dbsnp.rsid', returnall=True)
            for result in results['out']:
                # print('210',result)
                # result = result['hits'][0]
                test_output = {}
                if (True):
                    self.expansion(result, test_output, '')
                    strall = ''
                    for key in test_output.keys():
                        strall = strall + str(key) + ' : ' + test_output.get(key) + ' '
                    test_out = vcfinfo_list[count] +strall + '\n'
                    fileb.write(test_out)
                    for keys in test_output.keys():
                        if (test_output.get(keys, 'a') != 'a'):
                            title[keys] = 1
                            # detect whether there is dup in dup[] list
                    if (dul_count == 0):
                        for item in results['dup']:
                            # print('223',queryinfo_list[count])
                            # print(item)
                            if queryinfo_list[count] in item:
                                dul_count = item[1] - 1
                    else:
                        dul_count = dul_count - 1
                    if (dul_count == 0):
                        count = count + 1
            #count = 0
            #dul_count = 0
            #results = self.mv.querymany(queryinfo_list, 'dbsnp.rsid', returnall=True)
            #print(str(lines) + '\t' + str(result))
            #fileb.write(str(output)+ '\t' + str(result) + '\n'

        if (self.outputtype == 'csv'):
            count = 0
            dul_count = 0
            start_time = time.time()
            results = self.mv.querymany(queryinfo_list, 'dbsnp.rsid', returnall=True)
            return results
            print("--- time " + "%s sconds ---"%(time.time() - start_time))
            fileb.write("%s \n"%(time.time() - start_time))
            for result in results['out']:
                #print('210',result)
                #result = result['hits'][0]
                if (result != ""):
                    test_output = {'info': queryinfo_list[count]}
                    self.expansion(result, test_output, '')
                    #print(type(test_output))
                    #print(type(outputs))
                    outputs.append(test_output)
                    for keys in test_output.keys():
                        if (test_output.get(keys, 'a') != 'a'):
                            title[keys] = 1
                # detect whether there is dup in dup[] list
                    if (dul_count == 0) :
                        for item in results['dup'] :
                            #print('223',queryinfo_list[count])
                            #print(item)
                            if queryinfo_list[count]  in item :
                                dul_count = item[1] - 1
                    else :
                        dul_count = dul_count - 1
                    if (dul_count == 0):
                        count = count + 1
        if (self.outputtype == 'csv'):
            #print('197')
            all_keys = title.keys()
            #print('all keys : ' , all_keys)
            #dict_writer = csv.DictWriter(fileb, title, restval='Nan',)
            #dict_writer.writeheader()
            #dict_writer.writerows(outputs)
        #print("line number is: ", lineNumber)
        #print("identiable line number is: ", count)
        #print("non-identiable line number is: ", nonrs)
        # to close the file
        file.close()
        fileb.close()







def get_process_memory():
    pid = os.getpid()
    print(pid)
    ps = psutil.Process(pid)
    memoryUse = ps.memory_info()
    print('memoryUse.rss : ',memoryUse.rss)
    return memoryUse.rss

def test_func(num1,num2):

    mv = myvariant.MyVariantInfo(url='http://myvariant.info/v1')
    demo = query('23andme_large.txt','time_record.txt','23andme','csv', 19,'None',num1,num2)
    result = demo.genequery()
    return result

def read_in_chunks(file_object, chunk_size=1024):
    """Lazy function (generator) to read a file piece by piece.
    Default chunk size: 1k."""
    while True:
        data = file_object.read(chunk_size)
        if not data:
            break
        yield data


f = open('23andme_large.txt')
for piece in read_in_chunks(f):
    print(len(piece))


fileb = open('out.txt', "a", encoding='utf-8')
num = 2500
while num < 1000000:
    mem_before = get_process_memory()
    start = time.time()
    result = test_func(num,num+2500)
    elapsed_time = time.time() - start
    mem_after = get_process_memory()
    print(num)
    print("{}: memory before: {:,}, after: {:,}, consumed: {:,}; exec time: {}".format(
                'test_func',
                mem_before, mem_after, mem_after - mem_before,
                elapsed_time))
    fileb.write("{}: memory before: {:,}, after: {:,}, consumed: {:,}; exec time: {}".format(
                'test_func',
                mem_before, mem_after, mem_after - mem_before,
                elapsed_time))
    fileb.write("\n")
    num = num + 2500
