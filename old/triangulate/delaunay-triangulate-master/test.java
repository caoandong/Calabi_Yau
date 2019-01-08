te = require("delaunay-triangulate")

var points = [
  [0, 1],
  [1, 0],
  [1, 1],
  [0, 0],
  [0.5, 0.5]
]

var triangles = triangulate(points)

console.log(triangles)
