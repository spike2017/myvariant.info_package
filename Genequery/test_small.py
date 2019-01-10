#parallele
import myvariant
import time
from multiprocessing import Pool



if __name__ == '__main__':
    list_test = []
    list_all = []
    for i in range(1,1001):
        list_test.append('chr1:g.866422C>T')
    print('len of list',len(list_test))

    mv = myvariant.MyVariantInfo(url='http://myvariant.info/v1')

    for i in range(1,11):
        list_all.append(list_test)

    start = time.time()
    for item in list_all:
        result = mv.getvariants(item)
    #print(len(result))
    print('%s'%(time.time() - start))
#65.52 sec

    print('test1 end')
    p = Pool(2)
    start = time.time()
    result = p.map(mv.getvariants,list_all)
    print('%s'%(time.time() - start))
#20.96

'''
import time
from multiprocessing import Pool

def f(x):
    return x*x

list= [1,2,3,4,5,6,7,8,9,10]

def change(x):
    time.sleep(1)
    return x


    start = time.time()
    for item in list:
        print(change(item))
    print('%s'%(time.time() - start))

    print('process down')

    p = Pool()
    start = time.time()
    result = p.map(change,list)
    print('%s'%(time.time() - start))
'''
