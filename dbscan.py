from TriNode import *
#from rtree import index

#using MultilMatch.GeoStruct;
#using MultilMatch.BaseGeoAlgorithm;
class GCPoint(TriNode):
    def __init__(self,trinode):
        TriNode.__init__(self)
        self.ID = trinode.ID
        self.X = trinode.X
        self.Y = trinode.Y
        self.FeatureType =trinode.FeatureType
        self._isKey = False
        self._isClassed = False

    def IsKey(self):
        return self._isKey

    def SetKey(self, _isKey):
        self._isKey = _isKey
        self._isClassed = True

    def IsClassed(self):
        return self._isClassed

    def SetClassed(self, _isClassed):
        self._isClassed = _isClassed

class Utility():
    def __init__(self):
        self._isok = False

    @staticmethod
    def Euclidist(p,q):
        return math.sqrt((p.X - q.X)*(p.X - q.X) + (p.Y -q.Y)*(p.Y - q.Y))

    @staticmethod
    def isKeyPoint(lst, pp, ee, minp):
        count = 0
        tmpLst = []
        i = 0
        while i < len(lst) - 1:
            p = lst[i]
            q = lst[i] #lst[i + 1]
            if Utility.Euclidist(pp, q) <= ee:
                count += 1
                # if not tmpLst.contains(q):
                #     tmpLst.append(q)
                if q not in tmpLst:
                    tmpLst.append(q)
            i += 1
        if count >= minp:
            pp.SetKey(True)  #如果以pp为中心，半径为r，邻域内其他点的数量大于等于minpts，则该点就是核心点
            print(count)
            return tmpLst
        return None

    @staticmethod
    def MergeList(a, b):
        merge = False
        if a == None or b == None:
            return False
        index = 0
        while index < len(b):
            p = b[index]
            if p.IsKey() and p in a:
                merge = True
                break
            index += 1
        if merge:
            index = 0
            while index < len(b):
                if not b[index] in a:
                    a.append(b[index])
                index += 1
        return merge

class DBSCANCLUSTERING():
    def __init__(self):
        self._isclusterfinished = False

    @staticmethod
    def cluster(pointlist ,ee_value,minp_value):
        '''根据dbscan，进行聚类'''
        ee = ee_value #650
        minp = minp_value #24
        resultList =[]
        if len(pointlist)<3:
            print('输入要素小于三个，无法聚类，退出')
        else:
            for i in range(0,len(pointlist)):
                tmpLst = []
                point = pointlist[i]
                if (point._isClassed == True):
                    continue
                temLst = Utility.isKeyPoint(pointlist,point,ee,minp)
                if tmpLst is not None:
                    resultList.append(temLst)
            length = len(resultList)
            for i in range(0,length):
                for j in range(0,length):
                    if i !=j:
                        if Utility.MergeList(resultList[i],resultList[j]):
                            resultList[j].clear()

        return resultList
