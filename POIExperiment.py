from MapObject import *
from TriNode import *
from Util import *
from dbscan import *
import math
import os
import numpy as np
import ogr
#from rtree import index
import networkx as nx
#import osmnx as onx
from sklearn.cluster import DBSCAN



from anytree import Node,RenderTree,AsciiStyle
from LineAlgorithm import *


line_algorithm = LineAlgorithm()

def ReadPolygon(filename):
	driver = ogr.GetDriverByName("ESRI Shapefile")
	#filename = #'F:\\POI数据\\第一批POI数据\\成果数据\\01-北京\\deleteonepoint\\鱼or食.shp'
	dataSource = driver.Open(filename,0)
	if(dataSource is not None):
		layer = dataSource.GetLayer(0)
		n = layer.GetFeatureCount()
		extent = layer.GetExtent()

		Points = [TriNode]
		for fea in layer:
			geom = fea.GetGeometryRef()#GetGeometryRef  #GetGeometry
			if(geom is not None):
				geo_polygon = geom.GetGeometryRef(0)
				no_of_fea = geo_polygon.GetPointCount()
				for vertex in range(0,no_of_fea):
					x,y,z = geo_polygon.GetPoint(vertex)
					#print(x,y,z)
					p = TriNode(x,y,vertex)
					Points.append(p)

		Points.pop(0)
		#p = index.Property()
		idx = index.Index(property = Points)


		#ProxiNode.Create_WriteProxiNodes2Shp('E:/PycharmProjects/rastertest/test2','testExport',Points)
		print('success export')
		return Points


def ReadPolyline(fullfilename):
	driver = ogr.GetDriverByName("ESRI Shapefile")
	filename = fullfilename #E:\travelWork\TrafficTrajectory\RoadRecom\data2\OSMBJP.shp  E:\travelWork\TrafficTrajectory\RoadRecom\data2\data\OSMhandgeneral.shp
	dataSource = driver.Open(filename,0)
	if(dataSource is not None):
		layer = dataSource.GetLayer(0)
		n = layer.GetFeatureCount()
		extent = layer.GetExtent()
		Roads = []
		count = 0
		for fea in layer:
			coordinates = []
			geom = fea.GetGeometryRef()  # GetGeometryRef  #GetGeometry
			if (geom is not None):
				geo_polyline = geom #.GetGeometryRef(1)
				no_of_fea = geo_polyline.GetPointCount()
				for vertex in range(0, no_of_fea):
					x, y, z = geo_polyline.GetPoint(vertex)
					# print(x,y,z)
					p = TriNode(x, y, vertex)
					coordinates.append(p)
				polyline_i = PolylineObject(count,coordinates)
				Roads.append(polyline_i)
				count += 1
		return Roads


def ReadPoi(fullfilename,field_num,whetherToUTM = False):
	#print(gdal.__version__)
	driver = ogr.GetDriverByName("ESRI Shapefile")
	filename = fullfilename # 'E:\\PycharmProjects\\PoiExperiment2\\鱼or食.shp'#Export_Output.shp'
	dataSource = driver.Open(filename,0)
	if(dataSource is not None):
		layer = dataSource.GetLayer(0)
		n = layer.GetFeatureCount()
		extent = layer.GetExtent()

		Points = []
		count = 0
		for fea in layer:
			geom = fea.geometry()#GetGeometryRef()#GetGeometryRef  #GetGeometry
			value = fea.GetField(field_num)
			if(geom is not None):
				geo_point = geom#.GetGeometryRef()#这一步取0就是多边形
				x,y,z = geo_point.GetPoint()
				if whetherToUTM is False:
					#print(x,y,z)
					p = TriNode(x,y,count,value)
					Points.append(p)
					count += 1
					fea.Destroy()
				else:
					utmx ,utmy = laglon_to_UTMxy(y,x,117)
					p = TriNode(utmx,utmy,count,value)
					Points.append(p)
					count += 1
					fea.Destroy()

		Points.pop(0)
		#p = index.Property() #ProxiNode.Create_WriteProxiNodes2Shp('E:/PycharmProjects/rastertest/test2','testExport',Points)
		print('success export to POIS')
		return extent,Points

def ReadPoints(fullfilename,whetherToUTM = False):
	#print(gdal.__version__)
	driver = ogr.GetDriverByName("ESRI Shapefile")
	filename = fullfilename # 'E:\\PycharmProjects\\PoiExperiment2\\鱼or食.shp'#Export_Output.shp'
	dataSource = driver.Open(filename,0)
	if(dataSource is not None):
		layer = dataSource.GetLayer(0)
		n = layer.GetFeatureCount()
		extent = layer.GetExtent()

		Points = []
		count = 0
		for fea in layer:
			geom = fea.geometry()#GetGeometryRef()#GetGeometryRef  #GetGeometry
			value = fea.GetField(0)
			if(geom is not None):
				geo_point = geom#.GetGeometryRef()#这一步取0就是多边形
				x,y,z = geo_point.GetPoint()
				if whetherToUTM is False:
					#print(x,y,z)
					p = TriNode(x,y,count)
					Points.append(p)
					count += 1
					fea.Destroy()
				else:
					utmx ,utmy = laglon_to_UTMxy(y,x,117)
					p = TriNode(utmx,utmy,count)
					Points.append(p)
					count += 1
					fea.Destroy()

		#Points.pop(0)
		#p = index.Property() #ProxiNode.Create_WriteProxiNodes2Shp('E:/PycharmProjects/rastertest/test2','testExport',Points)
		print('success export to POIS')
		return extent,Points

def get_fieldname_index(target_name,full_file_name):
	daShapefile = full_file_name#r"E:\PycharmProjects\PoiExperiment2\鱼or食.shp"

	dataSource = ogr.Open(daShapefile)
	daLayer = dataSource.GetLayer(0)
	layerDefinition = daLayer.GetLayerDefn()

	#print('Name  -  Type  Width  Precision')
	for i in range(layerDefinition.GetFieldCount()):
		fieldName = layerDefinition.GetFieldDefn(i).GetName()
		if (fieldName == target_name):
			return i
		# fieldTypeCode = layerDefinition.GetFieldDefn(i).GetType()
		# fieldType = layerDefinition.GetFieldDefn(i).GetFieldTypeName(fieldTypeCode)
		# fieldWidth = layerDefinition.GetFieldDefn(i).GetWidth()
		# GetPrecision = layerDefinition.GetFieldDefn(i).GetPrecision()
		#
		# print(fieldName + " - " + fieldType + " " + str(fieldWidth) + " " + str(GetPrecision))

def ConvertRoadNetworkToG(junctions,Roads):
	nodes = {}
	paths = {}
	for juc in junctions:
		key = juc.ID
		nodes[key] = juc

	for path in Roads:
		key = path.ID
		paths[key] = path
	return nodes, paths

def laglon_to_UTMxy(lat,lon,center = 117):
	'lat 纬度，lon 经度，center 中央经度'
	b = center
	c = lat
	d = lon
	e = int(c)+(int(c * 100)-int(c)*100)/60 + (c*10000-int(c * 100)*100)/3600
	f = int(d)+(int(d *100) -int(d)*100)/60+(d * 10000-int(d *100)*100)/3600
	g = f - b
	h = g / 57.2957795130823
	i = math.tan(math.radians(e))
	j = math.cos(math.radians(e))
	k = 0.006738525415 * j * j
	l = i * i
	m = 1 + k
	n = 6399698.9018/math.sqrt(m)
	o = h * h * j * j
	p = i * j
	q = p * p
	r =(32005.78006+q*(133.92133+q*0.7031))
	x = 6367558.49686*e/57.29577951308-p*j*r+((((l-58)*l+61)*o/30+(4*k+5)*m-l)*o/12+1)*n*i*o/2
	y = ((((l-18)*l-(58*l-14)*k+5)*o/20+m-l)*o/6+1)*n*(h*j)
	return x,y

def partitionPOIwithtype(POIS):
	E1 = []
	E2 = []
	E3,E4,E5,E6 = [],[],[],[]
	for i in range(POIS):
		p = POIS[i]
		if p.type == 1:
			E1.append(p)
		elif p.type ==2:
			E2.append(p)
		elif p.type ==3:
			E3.append(p)
		elif p.type == 4:
			E4.append(p)
		elif p.type == 5:
			E5.append(p)
		elif p.type == 6:
			E6.append(p)
		else:
			continue
	return E1,E2,E3,E4,E5,E6


def prepare_road_edge_G(roads,intersectsND,name):
	'''roads 是shapefile格式的道路网数据，没有拓扑节点
	intersects是建立了拓扑节点的交叉路口点
	这样返回的edges有重合的路段，是子路段的更长覆盖长度的路段情况'''
	intersects = intersectsND[1]
	edges = []
	G = nx.MultiDiGraph(name = name)
	nodes = {}
	paths = {}
	key = 0
	imax = len(intersects)
	kmax = len(roads)
	for i in range(0,imax):
		j = i+1
		while j < len(intersects):
			p1 = intersects[i] #[i]
			p2 = intersects[j] #[j]
			segment = TriEdge.TriEdge(p1,p2)
			for k in range(0,kmax):
				#比较segment和roads每条线，是否有一条近似平行，如果有就是roads上的一个片段路段
				segment_whether_in_polyline = line_algorithm.PolylineInRectangle(segment,line_algorithm.PolylineBuffer(roads[k],1))#[k],1))
				if segment_whether_in_polyline is True:

					edges.append(segment)
					nodes[key] = p1
					paths[key] = segment
					key += 1
			j += 1
	for node,data in nodes.items():
		G.add_node(node,**data)
	for path in paths.values():
		G.add_path(G,path,False)
	return G,edges

def delete_longer_edges(edges):
	'''删除包含了小闭合环的大闭合环，其实arcgis建立了拓扑关系后
	不需要这一步删除操作'''
	return



#得到shapefile文件属性表中的属性
#indexId = get_fieldname_index('type',r"E:\PycharmProjects\PoiExperiment2\鱼or食.shp")

# E:\travelWork\TrafficTrajectory\RoadRecom\data2\data\OSMhandgeneral.shp
#
#
# E:\travelWork\TrafficTrajectory\RoadRecom\data2\data\OSMhandgeneral_ND_Junctions.shp
# Geometry Type:	Point

#distest = Util.Euclidist(Util,TriNode(1,1),TriNode(5,7))

#建立路网的拓扑结构图(0,20)和9号路
#intersect_nodes = ReadPoints(r"C:\Study\travelWork\TrafficTrajectory\RoadRecom\data2\data\OSMhandgeneral_ND_Junctions.shp",False)
#roads = ReadPolyline(r"C:\Study\travelWork\TrafficTrajectory\RoadRecom\data2\data\OSMhandgeneral.shp")
#这里图的拓扑组织没写好，有bug
#G,edges = prepare_road_edge_G(roads,intersect_nodes,'Road_G')

#Roads = ReadPolyline('E:\\travelWork\\TrafficTrajectory\\RoadRecom\\data2\\data\\OSMhandgeneral.shp')
#x 是经度，y是纬度

#Extent, POIS = ReadPoi(r"E:\travelWork\TrafficTrajectory\RoadRecom\data2\BJALL.shp",True)#E:\travelWork\TrafficTrajectory\RoadRecom\data2\BJALL.shp

# check_company = nx.read_shp(r"C:\Study\travelWork\TrafficTrajectory\RoadRecom\data\check_company.shp")
# A = nx.adjacency_matrix(check_company)
# numpy_K = nx.convert_matrix.to_numpy_matrix(A)
#nx.draw(check_company)#这个时间很慢

Extent, POIS = ReadPoi(r"D:\POI城市功能区\data\dbscantest.shp",False)
# extentx_low_utm,extenty_low_utm = laglon_to_UTMxy(Extent[2],Extent[0])
# extentx_up_utm,extenty_up_utm = laglon_to_UTMxy(Extent[3],Extent[1]) #我写的这个坐标转换没有用，还是得用arcgis的坐标系转换和定义投影操作
del_x_extent = Extent[2] - Extent[0] #extentx_up_utm - extentx_low_utm
del_y_extent = Extent[3] - Extent[1] #extenty_up_utm - extenty_low_utm

radius = math.fabs(del_x_extent / 8)
n_r_y = del_y_extent / radius

#尝试一下，ski kit learn的聚类
# x_list = []
# for i in range(0,len(POIS)):
# 	x_list.append([POIS[i].X,POIS[i].Y,POIS[i].ID])
# X = np.array(x_list)
# clustering2 = DBSCAN(eps = radius,min_samples= 100).fit(X)

#自己写的dbscan聚类
A = []
for i in range(0,len(POIS)):
	A.append(GCPoint(POIS[i]))
print('all poi trans to GCPoints')
#聚类
clustering = DBSCANCLUSTERING.cluster(A,radius,3)

#划分poi类别和签到类别后，建立树
#e1,e2,e3,e4,e5,e6 = partitionPOIwithtype(POIS)




G_extent = nx.MultiDiGraph(name = 'extent')
G_extent.add_nodes_from(Extent)
#G_extent_toUTM = onx.project_graph(G_extent)


POIS_with_Index = index.Index(property = POIS)

road_path_name = 'E:\travelWork\TrafficTrajectory\RoadRecom\BJR.shp'
G = nx.read_shp(road_path_name)

#G = nx.MultiDiGraph(name='POIeat', crs=settings.default_crs)

#G.add_nodes_from(POIS)
nodes = {}
paths = {}

print('success poi read')
