# HistSFC

HistSFC is a novel framework for querying nD point clouds considering skewed distribution. It converts a multidimensional point using Morton code into a one-dimensional key. Then, HistSFC organizes and indexes the keys using a B+-tree. To perform a window query or an aribitrary geoemtrical query (e.g., the polytope query), it utilizes an nD-histogram to recursively decomposed the data extent to approaximate the query geometry by one-dimensional ranges. Then, it searches the ranges in the B+-tree organized table to retrive data. Afterwards, a second filter conductes the decoding and point-wise filtering to derive the accurate result.

![mortonrange-1](https://user-images.githubusercontent.com/35140221/125162461-3cfdcd00-e188-11eb-9491-6aa7742a12c0.png)

![queryflow (4)](https://user-images.githubusercontent.com/35140221/125162545-ce6d3f00-e188-11eb-88fe-b2c34ee18312.png)

Publications:
An efficient nD-point data structure for querying flood risk
https://www.researchgate.net/publication/352893491_AN_EFFICIENT_ND-POINT_DATA_STRUCTURE_FOR_QUERYING_FLOOD_RISK 

HistSFC: Optimization for nD Massive Spatial Points Querying, https://www.researchgate.net/publication/342500095_HISTSFC_OPTIMIZATION_FOR_ND_MASSIVE_SPATIAL_POINTS_QUERYING 

An optimized SFC approach for nD window querying on point clouds, https://www.researchgate.net/publication/342500474_AN_OPTIMIZED_SFC_APPROACH_FOR_ND_WINDOW_QUERYING_ON_POINT_CLOUDS
