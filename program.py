# Library imports
import os
from networkx import astar_path, DiGraph, NetworkXNoPath
from osmnx import graph_from_xml
from rtree import index
from pyproj import Geod
from sys import exit
from shapely.geometry import LineString, Point
from helpers import Simplifier
from geopandas import GeoSeries
from matplotlib.pyplot import subplots
import matplotlib.pyplot as plt
from fiona import crs

# gets the current working directory
path = os.getcwd()
print(f'Current working directory :: {path}')

# set the working directory to the same file that the current script is running from
os.chdir(os.path.dirname(os.path.abspath(__file__)))


# function to simplify a list of coordinates
def simplify(coordinates, number=None, ratio=None, threshold=None):
    return Simplifier(coordinates).simplify(number=number, ratio=ratio, threshold=threshold).tolist()


# function for ellipsoidal distance
# Calculate the 'as the crow flies' distance between two locations along an ellipsoid using
# the Inverse Vincenty method, via the PyProj library.
def ellipsoidal_distance(a, b):
    # set the ellipsoid to be used - EXPLAIN WHY THIS ONE
    # this one is a pretty safe bet for global stuff
    g = Geod(ellps='WGS84')

    # extract nodes take the id of the nodes to be measured (a and b) use the grpah.nodes (data = True) function of
    # networkx to retrieve a list of nodes in the graph data=True means that the node list contains the data stored
    # inside each node, including coords, rather than just a list of ids extract the specific nodes in which we are
    # interested using [a] and [b] to retrieve the node with the corresponding id from the nodes list

    start = graph.nodes(data=True)[a]
    end = graph.nodes(data=True)[b]

    # compute forward and back azimuths, plus distance
    azf, azb, distance = g.inv(start['x'], start['y'], end['x'], end['y'])

    # return only the distance we are only interested in the distance
    return distance


# LOAD DATA IN XML (.osm) FORMAT AND CREATE GRAPH (from the xml dataset)
# downloaded as .osm file for greater manchester (GM)
# renamed .osm with .xml before loading data
# convert the graph into a DiGraph (to permit duplicate edges - EXPLAIN WHY)
# store the result in the variable 'graph'
# use relative file path to load the data into the graph
graph = DiGraph(graph_from_xml("mcr_osm.xml"))

# print to check that the .xml file loaded correctly (remove this later) - doesnt work
# # print(graph)
plt.plot(graph)
plt.show()


# CALCULATE A SPATIAL INDEX FROM THE DIGRAPH build a spatial index - needed to take locations and find the nearest
# node in the graph (one to define the start point and one to define the end point of the route)
idx = index.Index()
for _id, data in list(graph.nodes(data=True)):
    # convert the id to a number using int(id)
    idx.insert(int(_id), (data['x'], data['y'], data['x'], data['y']))

# CALCULATE NEAREST NODE TO THE START AND END POINT Just assign two random points for purposes of development
# calculate the 'from' and 'to' node as the nearest to the specified coordinates rtree spatial indexes think about
# the world in rectangles - so need to put bottom-left and top-right of rectangle in (but are the same coordinates as
# it is a point) can then use the rtree.nearest() function, which uses the spatial index to quickly locate the
# closest features (nodes in this case) in this case, as bottom-left and top-right share the same coord, ony want the
# nearest one coordinate to access the result, must convert it a list using the built-in python method list,
# and as it is only one item, it can be extracted directly with [0]
# define the start and end node coordinates (to be added back into the simplified list later to ensure simplified
# route is complete )
startPoint = Point(-2.23397, 53.48167)
endPoint = Point(-2.25337, 53.47556)
fromNode = list(idx.nearest((-2.23397, 53.48167, -2.23397, 53.48167), 1))[0]
toNode = list(idx.nearest((-2.25337, 53.47556, -2.25337, 53.47556), 1))[0]

# print the 'from' and 'to' nodes to the console
print(graph.nodes()[fromNode])
print(graph.nodes()[toNode])

# CALCULATE THE ROUTE using networkx A* algorithm to calculate the shortest path between the two nodes A* algorithm
# us one of the most popular and widely used algorithms function for A* routing in networkx: nx.astar_path(graph,
# source, target, heuristic) where graph (G), source node (fromNode), target node (toNode) and an heuristic function
# 'Heuristic' common in optimisation algorithms (like finding the shortest path) - simply, a heuristic function is
# one that ranks alternative approaches according to some sort of measure (e.g. distance) TRY-EXCEPT STATEMENT TO
# PREVENT USER SEEING 'normal' ERRORS avoid the user seeing a common NetworkXNoPath that occurs when networkx cannot
# calculate a route (usually because its not possible) using try-except statement
# use try statement to catch exceptions
try:
    # calculate the shortest path across the network and extract from graph astar function in the try block as this
    # has the potential to raise the NetworkXNoPath Exception) heuristic function used - enables measurement of the
    # geographical distance between each node when deciding which way to go
    shortestPath = astar_path(graph, source=fromNode, target=toNode, heuristic=ellipsoidal_distance)

# catch exception for no path available
# if exception is raised, block of code will stop what iit is doing and run some different code (a nice error message)
except NetworkXNoPath:
    print("Sorry, there is no path between those locations in the provided network")
    exit()

print(f'Shortest path in dataset: {shortestPath}')

# CONVERT LIST OF NODES TO A LineString GEOMETRY to draw route onto map, need to get the coordinates out of the graph
# object to make a LineString to save to a shape file
# loop through each node in the shortest path and store in an empty list
# (where each element refers to the id of the node )
line = []
for n in shortestPath:
    # get the relevant node from the graph with lat lng data, getting the node object using its id
    node = graph.nodes(data=True)[n]

    # line.append(tuple([-2.25337, 53.47556]))

    # load the lat lng data coordinate pair from the node to the list (into the lineString )
    # append adds the each set of coordinates to the list in order
    line.append([node['x'], node['y']])

# WHY DOES THE VW ALGORITHM NOT WORK IF I STORE IT AS A LINESTRING, BUT DOES IF I I USE THE LIST (AND IS THIS LIST
# KEOT IN ORDER) store as a MultiLineString after the loop is finished
original_route_lineString = LineString(line)

print(f'LineString of shortest path linestring (coords): {original_route_lineString}')
print(f'LineString of shortest path (coords): ({len(line)}')  # 31

vw_coords = simplify(line)
print(f'No. nodes in vw simplified linestring: ({len(vw_coords)} ')  # 27

# returns the linestring, simplified, keeping 90% of nodes (by percentage of points to keep)
vw_coords_2 = simplify(line, ratio=0.5)
print(f'No. nodes in vw simplified linestring ratio 0.5: ({len(vw_coords_2)} ')  # 15

# returns the linestring, simplified, using area threshold
vw_coords_3 = simplify(line, threshold=0.01)
print(f'No. nodes in vw simplified area threshold: ({len(vw_coords_3)} ')  # 15

# simplify the list of coordinates, with the ratio set to 0.26
vw_coords_2 = simplify(line, ratio=0.26)
print(f'No. nodes in vw simplified linestring ratio 0.5: ({len(vw_coords_2)} ')  # 15

# convert the list of simplified vw coordinates to a linestring to be plotted
vw_route_lineString = LineString(vw_coords_2)

# PLOT THE MAP - (CHECK IT WORKS)
# convert LineString to a Geoseries so that there is access to the .plot() function from GeoPandas pass to GeoSeries
# constructor and set the coord ref system  - use WGS84 , as this is the CRS of the underlying OpenStreetMap data can
# achieve this using the crs.from_epsg() function from the fiona library, which allows to retrive the CRS details
# simply by passing its EPSG code (4326 for WGS84) then project it to the Universal Traverse Mercator projection zone
# 34N, which is a projection appropriate to use in this location (check this is right because this was for Kalingrad)


# convert linestring to GeoSeries and project to UTM zone 34
utm34 = "+proj=utm +zone=34 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
path_1 = GeoSeries(original_route_lineString, crs=crs.from_epsg(4326)).to_crs(utm34)

# convert linestring to GeoSeries and project to UTM zone 34
utm34 = "+proj=utm +zone=34 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
path_2 = GeoSeries(vw_route_lineString, crs=crs.from_epsg(4326)).to_crs(utm34)

# create map axis object, remove axes, set title
fig, ax = subplots(1, 1, figsize=(15, 8))
ax.axis('off')
ax.set(title=" Map of route")

# set bounds
buffer = 100
ax.set_xlim([path_2.geometry.iloc[0].bounds[0] - buffer, path_2.geometry.iloc[0].bounds[2] + buffer])
ax.set_ylim([path_2.geometry.iloc[0].bounds[1] - buffer, path_2.geometry.iloc[0].bounds[3] + buffer])

# add the path
path_1.plot(
    ax=ax,
    color='blue',
    linewidth=2,
)

# add the path
path_2.plot(
    ax=ax,
    color='red',
    linewidth=2,
)

#plt.show()
