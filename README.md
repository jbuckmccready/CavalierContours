# Summary
C++ header only library for offsetting open and closed 2D curves, including self intersecting curves. Supports polylines defined by straight line and constant radius arc segments, other curves may be offset by first approximating them as a polyline. For interactive UI and development go to the development project [here](https://github.com/jbuckmccready/CavalierContoursDev).

<img src="https://raw.githubusercontent.com/jbuckmccready/CavalierContoursDoc/master/images/pretty_examples/example1.png" width="400"/> <img src="https://raw.githubusercontent.com/jbuckmccready/CavalierContoursDoc/master/images/pretty_examples/example2.png" width="400"/>

<img src="https://raw.githubusercontent.com/jbuckmccready/CavalierContoursDoc/master/images/pretty_examples/example6.png" width="400"/> <img src="https://raw.githubusercontent.com/jbuckmccready/CavalierContoursDoc/master/images/pretty_examples/example7.png" width="400"/>

# Table of Contents
- [Summary](#summary)
- [Table of Contents](#table-of-contents)
- [Quick Code Example](#quick-code-example)
- [Polyline Structure](#polyline-structure)
- [Algorithm and Stepwise Example](#algorithm-and-stepwise-example)
  - [Original input polyline, *pline* in blue, vertexes in red](#original-input-polyline-pline-in-blue-vertexes-in-red)
  - [Raw offset segments generated in purple (Step 1)](#raw-offset-segments-generated-in-purple-step-1)
  - [Raw offset polyline created from raw offset segments, *pline1* (in green) (Step 2)](#raw-offset-polyline-created-from-raw-offset-segments-pline1-in-green-step-2)
  - [Raw offset polyline self intersects (dark cyan) (Step 4)](#raw-offset-polyline-self-intersects-dark-cyan-step-4)
  - [Valid open polyline slices created from self intersects (in green, red, and blue) (Step 5 & 6)](#valid-open-polyline-slices-created-from-self-intersects-in-green-red-and-blue-step-5--6)
  - [Open polyline slices stitched together (in red and blue) (Step 7)](#open-polyline-slices-stitched-together-in-red-and-blue-step-7)
- [Interactively Exploring the Algorithm](#interactively-exploring-the-algorithm)
- [Performance](#performance)
- [Implementation Notes and Variations](#implementation-notes-and-variations)
  - [Float Comparing and Thresholding](#float-comparing-and-thresholding)
  - [Joining Raw Offset Segments](#joining-raw-offset-segments)
  - [Stitching Open Polylines](#stitching-open-polylines)
- [Development](#development)
- [API Stability](#api-stability)
- [Project Background](#project-background)
- [Algorithm Complexity and 2D Spatial Indexing](#algorithm-complexity-and-2d-spatial-indexing)
  - [Packed Hilbert R-Tree](#packed-hilbert-r-tree)
- [References](#references)


# Quick Code Example
```c++
#include "cavc/polyline.h"

// input polyline
cavc::Polyline<double> input;
// add vertexes as (x, y, bulge)
input.addVertex(0, 25, 1);
input.addVertex(0, 0, 0);
input.addVertex(2, 0, 1);
input.addVertex(10, 0, -0.5);
input.addVertex(8, 9, 0.374794619217547);
input.addVertex(21, 0, 0);
input.addVertex(23, 0, 1);
input.addVertex(32, 0, -0.5);
input.addVertex(28, 0, 0.5);
input.addVertex(39, 21, 0);
input.addVertex(28, 12, 0);
input.isClosed() = true;

// compute the resulting offset polylines, offset = 3
std::vector<cavc::Polyline<double>> results = cavc::parallelOffset(input, 3);
```
# Polyline Structure
Polylines are defined by a sequence of vertexes and a bool indicating whether the polyline is closed or open. Each vertex has a 2D position (x and y) as well as a bulge value. Bulge is used to define arcs, where `bulge = tan(theta/4)`. `theta` is the arc sweep angle from the starting vertex position to the next vertex position. If the polyline is closed then the last vertex connects to the first vertex, otherwise it does not (and the last vertex bulge value is unused). See [[2]](#references) for more details regarding bulge calculations.

# Algorithm and Stepwise Example
1. Generate raw offset segments from the input polyline, *pline*.
2. Create the raw offset polyline, *pline1*, by trimming/joining raw offset segments acquired in step 1.
3. If the input polyline, *pline*, has self intersections or is an open polyline then repeat steps 1 and 2 with the offset negated (e.g. if the offset was 0.5 then create raw offset polyline with offset of -0.5), this is known as *pline2*.
4. Find all self-intersects of *pline1*. If step 3 was performed then also find all intersects between *pline1* and *pline2*. If *pline* is an open polyline then also find intersects between *pline1* and circles at the start and end vertex points of *pline* with radius equal to the offset.
5. Create a set of open polylines by slicing *pline1* at all of the intersect points found in step 4. 
6. Discard all open polyline slices whose minimum distance to *pline* is less than the offset.
7. Stitch together the remaining open polyline slices found in step 6, closing the final stitched results if *pline* is closed.

The algorithm is mostly based on Liu et al. [[1]](#references) with some differences since the algorithm they describe for GCPP (general closest point pair) clipping fails for certain inputs with large offsets (or I atleast I am unable to make their algorithm work).

The key clarifications/differences are:
- When raw offset segments are extended to form a raw offset polyline they are always joined by an arc to form a rounded constant distance from the input polyline.
- Dual offset clipping is only applied if input polyline is open or has self intersects, it is not required for a closed polyline with no self intersects.
- If the polyline is open then a circle is formed at each end point with radius equal to the offset, the intersects between those circles and the raw offset polyline are included when forming slices.
- GCPP (general closest point pair) clipping is never performed and instead slices are formed from intersects, then they are discarded if too close to the original polyline, and finally stitched back together.
- No special handling is done for adjacent segments that overlap (it is not required given the slice and stitch method).
- Collapsing arc segments (arcs whose radius is less than the offset value) are converted into a line and specially marked for joining purposes.

Here is example code and visualizations of the algorithm operating on a closed polyline with no self intersects as input. 
## Original input polyline, *pline* in blue, vertexes in red
```c++
#include "cavc/polyline.h"
using namespace cavc;
// Create input polyline
Polyline<double> input;
// Add vertexes as (x, y, bulge)
input.addVertex(0, 25, 1);
input.addVertex(0, 0, 0);
input.addVertex(2, 0, 1);
input.addVertex(10, 0, -0.5);
input.addVertex(8, 9, 0.374794619217547);
input.addVertex(21, 0, 0);
input.addVertex(23, 0, 1);
input.addVertex(32, 0, -0.5);
input.addVertex(28, 0, 0.5);
input.addVertex(39, 21, 0);
input.addVertex(28, 12, 0);
input.isClosed() = true;
```
![Input Polyline](https://raw.githubusercontent.com/jbuckmccready/CavalierContoursDoc/master/images/algorithm_steps/input_polyline.png)

## Raw offset segments generated in purple (Step 1)
```c++
// Step 1. Create untrimmed offset segments
double offset = 3;
std::vector<PlineOffsetSegment<double>> rawOffsets = createUntrimmedOffsetSegments(input, offset);
```
![Raw Offset Segments](https://raw.githubusercontent.com/jbuckmccready/CavalierContoursDoc/master/images/algorithm_steps/raw_offset_segments.png)

## Raw offset polyline created from raw offset segments, *pline1* (in green) (Step 2)

```c++
// Step 2. Create raw offset polyline
// NOTE The code is factored such that createRawOffsetPline handles creating the untrimmed 
// offsets done in step 1
Polyline<double> pline1 =  createRawOffsetPline(input, offset);
```
![Raw Offset Polyline](https://raw.githubusercontent.com/jbuckmccready/CavalierContoursDoc/master/images/algorithm_steps/raw_offset_polyline.png)

## Raw offset polyline self intersects (dark cyan) (Step 4)

```c++
// Step 4. Find all self intersects of pline1 (Step 3 skipped for this example case)
// spatial index used for finding self intersects
StaticSpatialIndex<double> rawOffsetPlineSpatialIndex = createApproxSpatialIndex(pline1);
std::vector<PlineIntersect<double>> selfIntersects;
// selfIntersects passed by reference to be filled by allSelfIntersects
allSelfIntersects(pline1, selfIntersects, rawOffsetPlineSpatialIndex);
```
![Raw Offset Polyline Intersects](https://raw.githubusercontent.com/jbuckmccready/CavalierContoursDoc/master/images/algorithm_steps/raw_offset_polyline_intersects.png)

## Valid open polyline slices created from self intersects (in green, red, and blue) (Step 5 & 6)

```c++
// Step 5 and 6. Create valid open polyline slices
// NOTE The code is factored such that sliceAtIntersects handles finding the self intersects in step 4
// and discarding invalid slices in step 6
std::vector<OpenPolylineSlice<double>> slices = sliceAtIntersects(input, pline1, offset);
```
![Valid Slices](https://raw.githubusercontent.com/jbuckmccready/CavalierContoursDoc/master/images/algorithm_steps/valid_slices.png)

## Open polyline slices stitched together (in red and blue) (Step 7)

```c++
// Step 7. Stitch valid slices together
std::size_t origMaxIndex = pline1.size() - 1;
std::vector<Polyline<double>> finalResult = stitchSlicesTogether(slices, input.isClosed(), origMaxIndex)
```
![Final Output Polylines](https://raw.githubusercontent.com/jbuckmccready/CavalierContoursDoc/master/images/algorithm_steps/output.png)


# Interactively Exploring the Algorithm
An interactive UI app (implemented using Qt and QML) is available to visualize and explore in real time the offset algorithm. The app was used to generate all the images in this markdown. Link to project [here](https://github.com/jbuckmccready/CavalierContoursDev)

# Performance
The implementation is not entirely geared around performance but some profiling has been done to determine where to reserve memory up front and avoid trig functions where possible. I suspect there is quite a bit of performance to be gained in using a custom memory allocator, and there may be ways to modify the algorithm to prune invalid segments faster.

Additionally the structures and control flow are factored for readability/maintainability and may be limiting how much is compiled to simd instructions.

# Implementation Notes and Variations
## Float Comparing and Thresholding
When comparing two float values for equality, `a` and `b`, this implementation uses a simple fuzzy comparing method of `abs(a - b) < epsilon`, where `epsilon` is a very small number, e.g. 1e-8 or 1e-5 depending on the context. These numbers were picked through anecdotal use case trial/error where input values are typically between 0.1 and 1000.0. If finer comparisons are required, or if the input values get quite large or small then numerical stability issues may arise and the means of fuzzy comparing will need to be adjusted. See [mathutils.h](include/cavc/mathutils.h) for c++ implementation.

## Joining Raw Offset Segments
When joining raw offset segments together the current implementation always uses an arc to connect adjacent segments if they do not intersect, the arc maintains a constant distance from the original input polyline. Alternatively these joins could be done by extending the original line/arc segments until they intersect or by some other means entirely. Note that any other type of join will result in the offset polyline not being at a constant offset distance from the original input.

## Stitching Open Polylines
When stitching the open polylines together the current implementation only looks to stitch end points to start points. This assumes that the slices will always stitch together this way, this is a valid assumption if all slices are from the same raw offset polyline, but in the case that you want to stitch slices from different offsets together the implementation must look at both start and end points since directionality may change.

When stitching open polylines together the current implementation attempts to stitch the slices into the longest polylines possible (there are multiple possibilities when slices become coincident or tangent at end points with one another). Alternatively one could implement it in a way to stich slices into the shortest polylines possible. This may be useful if in the case of coincident stretches the result should be marked or discarded.

# Development
Pull requests, feature requests/ideas, issues, and bug reports are welcome. Please attempt to follow the code style and apply clang-format (using the .clang-format file) before making a pull request.

# API Stability
There is not an official release yet - all functions and structures are subject to change. The purpose of this repository for now is to serve as an implementation reference that is easy to understand and possibly transcribe to other programming languages. The code can be used as is but there is no guarantee that future development will maintain the same functions and structures. Ideas/pull requests for what a stable API interface should look like are welcome.

# Project Background
There are many papers on offsetting curves for CAD/CAM uses [[1][9][10][16]](#references), but there is often no reference implementation. I have found that most algorithms described are dense and difficult to reproduce. Issues such as numeric stability, how to handle coincident segments, etc. are often not mentioned and algorithmic description detail is inconsistent (e.g. it may be very clear how to perform some of the steps but other steps are quickly glossed over and it becomes unclear how to go about implementing). All of these issues would not be much of a problem if an open source reference implementation was supplied but for all the papers I have read not a single implementation was given.

In addition to papers being difficult to use as a pragmatic tool, most papers focus on offsetting straight segment polylines or polygons (sometimes referred to as point sequence curves) [[9][10][16]](#references). And there are a few notable open source libraries that work only on straight segment polylines as well [[11][12][13]](#references). Unfortunately if only straight segments are supported then all curves, even simple constant radius arcs, must be approximated using straight segments, which is very inefficient due to the memory footprint required. Additionally if arcs must be reconstructed from points then that adds to the code complexity. Constant radius arcs are very common in CAD/CAM applications (tool compensation, tool offsetting for cleanout, part sizing, etc.), and arcs may be used to approximate other curves more memory efficiently than straight line segments, e.g. for Bezier curves [[14]](#references).

There are a few papers on offsetting curves with arc segments [[1][15]](#references), some involve a Voronoi diagram approach [[15]](#references) and others are more similar to the approach this library takes by processing self-intersects and applying a clipping algorithm [[1]](#references). These papers do not provide an opensource reference implementation, and many of them are difficult to understand, and even more difficult to reproduce.

The goal of this project is to provide a simple, direct, and pragmatic algorithm and reference implementation for polyline offsetting (supporting both lines and arcs, open and closed) for use in CAD/CAM applications. It includes the minimum computational geometry building blocks required (vectors, spatial indexing, etc.) to accomplish this goal. The code is written with the intent that it can be easily read, modified, and transcribed to other programming languages. The [interactive UI app](https://github.com/jbuckmccready/CavalierContoursDev) can be used to help understand the algorithm visually and make the parts of the code easier to digest.

# Algorithm Complexity and 2D Spatial Indexing
The algorithm requires finding self-intersects, testing the distance between points and the original polyline, testing for intersects between new segments and the original polyline, and stitching open polylines together end to end. The naïve approach to such steps typically results in O(n<sup>2</sup>) algorithmic complexity where each segment must be tested against every other segment and each point must be tested against every segment.

The approach used by this library is to use a packed Hilbert R-Tree [[3]](#references). See [staticspatialindex.h](include/cavc/staticspatialindex.h) for c++ implementation. See [[4]](#references) and [[5]](#references) for more implementation references. This results in an algorithm complexity of O(n log n) for typical inputs, with worse case O(n<sup>2</sup>) for pathological inputs.

## Packed Hilbert R-Tree
Here is an image of a closed polyline approximating a circle using 100 line segments (blue lines and red vertexes) with spatial index bounding boxes made visible (magenta, orange, and light green boxes). The root of the R-Tree is the light green box, its children are the orange boxes, and its grand children are the magenta boxes.

![Simple Spatial Index](https://raw.githubusercontent.com/jbuckmccready/CavalierContoursDoc/master/images/spatial_index/simple_example.png)

The packed Hilbert R-Tree is very fast to build and query, but requires rebuilding the tree anytime the data changes – however for computational geometry algorithms such as those performed in this library we start with no indexed data and our data does not change for the duration of the algorithm, making the tradeoff ideal for the use case.

There are alternative approaches, e.g. by combining an Interval Tree [[6]](#references) and binary heap [[7]](#references) a sweep line algorithm [[8]](#references) approach can be taken to avoid O(n<sup>2</sup>) complexity.

Note that there are pathological input cases that will still result in O(n<sup>2</sup>) behavior even with a spatial index or sweep line algorithm. E.g. a polyline for which all segments have axis-aligned bounding boxes that overlap. But for common inputs seen in CAD/CAM applications this is not the case, and a run time complexity of O(n log n) will result.

# References
[1] Liu, X.-Z., Yong, J.-H., Zheng, G.-Q., & Sun, J.-G. (2007). An offset algorithm for polyline curves. Computers in Industry, 58(3), 240–254. doi:10.1016/j.compind.2006.06.002

[2] Bulge conversions: http://www.lee-mac.com/bulgeconversion.html

[3] https://en.wikipedia.org/wiki/Hilbert_R-tree#Packed_Hilbert_R-trees

[4] JavaScript spatial index implementation: https://github.com/mourner/flatbush/

[5] Fast 2D to 1D Hilbert curve mapping used for spatial index: https://github.com/rawrunprotected/hilbert_curves

[6] https://en.wikipedia.org/wiki/Interval_tree

[7] https://en.wikipedia.org/wiki/Binary_heap

[8] https://en.wikipedia.org/wiki/Sweep_line_algorithm

[9] Lin, Z., Fu, J., He, Y., & Gan, W. (2013). A robust 2D point-sequence curve offset algorithm with multiple islands for contour-parallel tool path. Computer-Aided Design, 45(3), 657–670. doi:10.1016/j.cad.2012.09.002

[10] Kim, D.-S. (1998). Polygon offsetting using a Voronoi diagram and two stacks. Computer-Aided Design, 30(14), 1069–1076. doi:10.1016/s0010-4485(98)00063-3

[11] Clipper library: http://www.angusj.com/delphi/clipper.php

[12] CGAL library for offsetting polylines: https://doc.cgal.org/latest/Straight_skeleton_2/index.html

[13] Boost geometry: https://www.boost.org/doc/libs/1_53_0/libs/geometry/doc/html/index.html

[14] Bezier biarc approximating: https://github.com/domoszlai/bezier2biarc

[15] Held, M., & Huber, S. (2009). Topology-oriented incremental computation of Voronoi diagrams of circular arcs and straight-line segments. Computer-Aided Design, 41(5), 327–338. doi:10.1016/j.cad.2008.08.004 

[16] Kim, H.-C., Lee, S.-G., & Yang, M.-Y. (2005). A new offset algorithm for closed 2D lines with Islands. The International Journal of Advanced Manufacturing Technology, 29(11-12), 1169–1177. doi:10.1007/s00170-005-0013-1
