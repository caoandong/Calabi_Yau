# Machine Learning on Sasaki-Einstein Manifolds Raw Data

The content of this repository can be separated into two parts: triangulation (data generation) and machine learning. All of the data obtained from both parts are located in the `/output` directory. Some directories in the `/output` folder are outdated, such as `output/vol` and `output/failed/triang`. Most of the files in the root directory, such as `hilbert.sage` and `train.py`, are outdated as well.
The triangulation algorithm contains everything from generating the polytopes to calculating the hilbert series to finding the minimum volume.

## Triangulation

The triangulation dataset contains:
* Cutting corners from a 3D cubic toric diagram
    * Algorithm: `triangulate/cube.sage`
    * Data:
        * Volume and the coordinates of the vertices: `output/train/cube/cube_nxn.txt`. Note: the `n` indicates the number of vertices on the cube's side. For those geometries that failed to show a volume, the volume is -1. (I should probably extract these geoemtries out and save them to the `output/failed` directory)
        * Volume and topological quantities (number of corner pts, edge pts, face pts, and body pts): `output/train/cube/count_nxn.txt`
        * Hilbert series: `output/series/cube`. Note: some of the Hilbert series files in the `output/series` folder, such as `series_cube_n.txt`, are outdated. For many of the cube's data, I forgot to save the Hilbert series; I can regenerate these Hilbert Series if necessary.
        * Vertices and topologial quantities (the inputs to the triangulation algorithm): `output/polygon/cube` and `output/topology/cube`
* Triangular and square prisms
    * Algorithm:
        * Triangular base prism: `triangulate/Triangulation.sage`
        * Sqaure base prism: triangulate/square_prism.sage`
    * Data:
        * Volume and the coordinates of the vertices: `output/train/cylinder/tri_nxn.txt` or `output/train/cylinder/sq_nxn.txt` or `output/train/cylinder/lift_1_to_50.txt`. Note: for those geometries that failed to show a volume, the volume is -1. The files `tri_1_to_50_2.txt` and `sq_1_to_50_2.txt` contain data for all triangular/square prisms from total height = 1 to total height = 50.
        * Volume and topological quantities (number of corner pts, edge pts, face pts, and body pts): `output/train/cube/count_nxn.txt`
        * Hilbert series: `output/series/cylinder`.
        * Vertices (the inputs to the triangulation algorithm): `output/polygon/cylinder`

## Machine Learning

The machine learning algorithm has two parts: fully-connected linear regression or PointNet
* Fully connected linear regression: machine_leanring/fully_connected
`** Note: the fully connected linear regressions contain: 1 layer, 2 layers, and 3 layers models. The cost functions is the mean average of the squared distance between the prediction and the actual value. To preprocess the input data, I multiply them to gether to produce terms like x^2, y^2, z^2, xy, xz, yz, x^2y, x^2z, xy^2, xz^2, etc., in order to approximate the Taylor expansion of an equation.
* PointNet: machine_leanring/pointnet

A sample output of the machine learning algorithms are located at `output/train/cylinder/loss`. I am too embarassed to save all the result, and I am trying to obtain better ones.
