# HistSFC

HistSFC is a novel framework for querying nD point clouds considering skewed data distributions. It converts a multidimensional point — using Morton code — to a one-dimensional key. Then, HistSFC organizes and indexes the keys using a B+-tree. To perform a window query or an aribitrary geoemtrical query (e.g., the polytope query), it utilizes an nD-histogram to recursively decompose the data extent to approaximate the query geometry by one-dimensional ranges. Then, it searches the ranges in the B+-tree organized table to retrieve the data. Afterwards, a second filter conductes the decoding and point-wise filtering to derive the accurate result.

![mortonrange-1](https://user-images.githubusercontent.com/35140221/125162461-3cfdcd00-e188-11eb-9491-6aa7742a12c0.png)

Figure 1. Executing a window query on a uniformly distributed 2D point set

![queryflow (4)-1](https://user-images.githubusercontent.com/35140221/125162724-b813b300-e189-11eb-92a1-1f8b1cae4a78.png)

Figure 2. The loading and querying procedure of the Morton index-organized table approach

![histtree-1](https://user-images.githubusercontent.com/35140221/125191321-12704a80-e242-11eb-84a4-8d9742a4f212.png)

Figure 3. A 2D HistogramTree example, where the threshold is 100; left: point counting in Morton quadrants, middle: pointer structure of HistogramTree, with each node storing a Morton key and number of points, right: structure of a HistTree node

Publications:

An efficient nD-point data structure for querying flood risk,
https://www.researchgate.net/publication/352893491_AN_EFFICIENT_ND-POINT_DATA_STRUCTURE_FOR_QUERYING_FLOOD_RISK 

HistSFC: Optimization for nD Massive Spatial Points Querying, https://www.researchgate.net/publication/342500095_HISTSFC_OPTIMIZATION_FOR_ND_MASSIVE_SPATIAL_POINTS_QUERYING 

An optimized SFC approach for nD window querying on point clouds, https://www.researchgate.net/publication/342500474_AN_OPTIMIZED_SFC_APPROACH_FOR_ND_WINDOW_QUERYING_ON_POINT_CLOUDS
