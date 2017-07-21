#coding:utf-8
from decode import *
from numpy import *
import xlrd
import random
import time
'''
   * If you want the source code and the notions, please go to the
     'https://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders#BCH_error_correction'
     for more details.
   **The project is also available on the
     'https://github.com/tomerfiliba/reedsolomon'
'''

prim = 0x11d
n = 52
k = 32
group_num =2
init_tables(prim)
colValues = [[]]*100
RepetCodeSize = 3

def fileRead(fileName,index):
    RightData = xlrd.open_workbook(fileName)
    table1 = RightData.sheets()[0]
    nrows1 = table1.nrows  # 行数
    ncols1 = table1.ncols  # 列数
    for i in range(0,ncols1-1):
        colValues[i] = table1.col_values(i+1)[index + 1 : index + 257]
        while '' in colValues[i]:
            colValues[i].remove('')
    return colValues


def readNameAndIndex():
    Data = xlrd.open_workbook('index.xlsx')
    table = Data.sheets()[0]
    nrows1 = table.nrows  # 行数
    colValues = {}
    name = []
    for i in range(nrows1):
        colValues[str(table.row_values(i)[1])]=int(table.row_values(i)[0])
        name.append(str(table.row_values(i)[1]))
    return colValues,name


def dataSet(fileName,index):
    return fileRead(fileName,index)


def extrcRightData(rawValue):
    correct_data = []
    for i in range(256):
        count = 0
        for j in range(0, 100):
            if int(rawValue[j][i]) == 1:
                count += 1
        if count >= 50:
            correct_data.append(1)
        else:
            correct_data.append(0)
    return correct_data


def bin2dec(data):
    dataset = []
    for i in range(len(data)/4):
        temp = 0
        for j in range(4):
            temp += data[i * 4 + j]*(2**j)
        dataset.append(temp)
    return dataset



def dec2bin(data):
    dataset = []
    for i in range(len(data)):
        for j in range(4):
            dataset.append((data[i] >> j) & 1)
    return dataset


def list2matrix(list,size = 16):
    matrix = []
    for i in range(len(list)/size):
        matrix.append(list[size * i : size * (i + 1)])
    matrix = array(matrix)
    return matrix


def randindex():
    return random.randint(0,99)


def reproData(dataset):
    '''重新生成PUF数据'''
    PUFOutput = [[]]*RepetCodeSize
    for reNum in range(RepetCodeSize):
        PUFOutput[reNum] = [int(i) for i in dataset[randindex()]]
    return PUFOutput


def repetitionCode(rawValue):
    codeWord = []
    numRead = RepetCodeSize
    for i in range(256):
        count1 = 0
        count0 = 0
        for j in range(numRead):
            if rawValue[j][i] == 1:
                count1 += 1
            else:
                count0 += 1
        if count1 >= count0:
            codeWord.append(1)
        else:
            codeWord.append(0)
    return codeWord


def encode(message):
    messageSet = [int(i) for i in message]
    messageMatrix = [[]] * group_num
    checkbit = [[]] * group_num
    for j in range(group_num):
        messageMatrix[j] = messageSet[j * 32: (j + 1) * 32]
        mesecc = rs_encode_msg(messageMatrix[j], n - k)
        checkbit[j] = mesecc[k:]
    return checkbit


def decode(receivebit,checkbit):
    receBitAfterRepet = bin2dec(receivebit)
    receiveBitMatrix = [[]] * group_num
    codedWord=[]
    errlocation=[]
    for j in range(group_num):
        receiveBitMatrix[j] = receBitAfterRepet[j * k: (j + 1) * k]
        receiveBitMatrix[j].extend(checkbit[j])
        corrected_message,  err_pos = rs_correct_msg(receiveBitMatrix[j], n - k)
        err_pos = [i + 1 + 32 * j  for i in err_pos]
        codedWord += corrected_message
        errlocation += err_pos
    errlocation.sort()
    return codedWord,errlocation




if __name__ == '__main__':
    print ('Ready to test!')
    flag = 1

    while flag == 1:
        rightStructName = 1
        while rightStructName == 1:
            Structure_name = raw_input('=> 请输入结构名 (structure):\n(8T_V1~8T_V4)')
            rightStructName = 0
            rightDieName = 1
            if Structure_name in ['8T_V1', '8T_V2', '8T_V3', '8T_V4']:
                print (
                '* 晶圆序号可选择有：(1)150328  (2)150934  (3)151057  (4)151220  (5)151343  (6)151506  (7)151629  (8)151752 ' \
                '(9)151915  (10)152038  (11)152201  (12)152324  (13)152446  (14)152609  \n'
                '(15)152732  (16)152855  (17)153018  (18)153141  (19)153304  (20)153427  (21)153550  (22)153713  (23)153836  (24)153959 ' \
                '(25)154122  (26)154245  (27)154408  (28)154531  \n'
                '(29)154654  (30)154817  (31)154940  (32)155103  (33)155226  (34)155349  (35)155512  (36)155635  (37)155758  (38)155921  '
                '(39)160044  (40)160207  (41)160330')
                while rightDieName == 1:
                    Die_name = raw_input('=> 请输入晶圆序号 (die): \n(logfile161226_xxxxxx)')
                    rightDieName = 0
                    if 'logfile161226_' + Die_name in readNameAndIndex()[1]:
                        print ('请稍等，正在读取数据......')
                        start = time.time()
                        rawValue = fileRead(Structure_name + '.xlsx',
                                            readNameAndIndex()[0]['logfile161226_' + Die_name])
                        rightDataBin = extrcRightData(rawValue)
                        rightDataDec = bin2dec(rightDataBin)
                        checkBit = encode(rightDataDec)
                        end = time.time()
                        print ('读取成功！')
                        print '耗时：%ss' % (end - start)

                        mode = raw_input('=> 请输入检测模式：\n1.依次检测(100组全部进行纠错检测，不使用重复码) \n2.随机检测(模拟PUF采样模式，并使用重复码)')
                        if mode == '1':
                            for j in range(100):
                                testData = [int(i) for i in rawValue[j]]
                                print '\n', '第%s组数据为(十进制)：' % (j + 1), '\n', list2matrix(bin2dec(testData))
                                decodeWordDec, err_pos = decode(testData, checkBit)
                                print '第%s组数据经过纠正后为(十进制)：' % (j + 1), '\n', list2matrix(decodeWordDec)
                                if decodeWordDec == rightDataDec:
                                    print '第%s组数据纠正通过' % (j + 1), '错误发生位置为：', err_pos, '错误个数：', len(err_pos)
                        # #下面代码是用来检验代码是否真的纠了错
                        # #十进制
                        # for j in range(10):
                        #     c = 0
                        #     l = []
                        #     testData = [int(i) for i in rawValue[j]]
                        #     testDataDec = bin2dec(testData)
                        #     for k in range(64):
                        #         if testDataDec[k] != rightDataDec[k]:
                        #             c += 1
                        #             l.append(k+1)
                        #     print j,c, l
                        # #二进制
                        # for j in range(10):
                        #     c = 0
                        #     l = []
                        #     testData = [int(i) for i in rawValue[j]]
                        #     for k in range(256):
                        #         if testData[k] != rightDataBin[k]:
                        #             c += 1
                        #             l.append(k+1)
                        #     print j,c, l
                        elif mode == '2':
                            PUFdata = reproData(rawValue)
                            dataAfterRepet_Bin = repetitionCode(PUFdata)
                            dataAfterRepet_Dec = bin2dec(dataAfterRepet_Bin)
                            print 'PUF数据经过重复码纠正后为(十进制)：', '\n', list2matrix(dataAfterRepet_Dec)
                            decodeWordDec, err_pos = decode(dataAfterRepet_Bin, checkBit)
                            print 'PUF数据经过纠错码纠正后为(十进制)：', '\n', list2matrix(decodeWordDec)
                            if decodeWordDec == rightDataDec:
                                print '纠错成功！！', '错误发生位置为：', err_pos, '错误个数：', len(err_pos)
                    else:
                        print ('不存在的晶圆!!')
                        rightDieName = 1
            else:
                print ('不存在的结构！！')
                rightStructName = 1

        flag = input('是否继续测试？1：继续 2：结束')

    print '测试结束！'


#-------------------------------------------------------------
# def repeCodeTest1(cyc):
#     '''重复码测试：测试重复发送次数与输出的关系,输入参数为重复发送的次数，输出为每个cell在这个码下的错误率'''
#     data = dataSet()
#     g = []
#     for f in range(100):
#         count = 0
#         for h in range(100):
#             count1 = 0
#             count2 = 0
#             for a in range(cyc):
#                 b = randIndex()
#                 if data[b][f] == 1:
#                     count1 += 1
#                 else:
#                     count2 += 1
#             if count1 > count2:
#                 g.append(1)
#                 count += 1
#             elif count1 < count2:
#                 g.append(0)
#             else:
#                 pass
#         if count >= 50:
#             count = 100 - count
#         print '%d: %d%s'%(f,count,'%')
# repeCodeTest(3)

# def repeCodeTest12():
#     '''重复码测试：测试在重复码纠正下错误发生位置和错误数量与未经过重复码的随机测试的对比'''
#     a = extrcRightData(dataSet())
#     b = repetitionCode(dataSet())
#     c = dataSet()[randIndex()]
#     total1 = 0
#     total2 = 0
#     err_loc1 = []
#     err_loc2 = []
#     for i in range(256):
#         if a[i] != b[i]:
#             total1 += 1
#             err_loc1.append(i)
#     print 'With repetiton code:','\n','  error location:',err_loc1,'\n','  number of errors:', total1
#     for i in range(256):
#         if a[i] != c[i]:
#             total2 += 1
#             err_loc2.append(i)
#     print 'Without repetiton code:','\n','  error location:',err_loc2,'\n','  number of errors:', total2
# repeCodeTest12()
