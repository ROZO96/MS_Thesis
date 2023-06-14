rbox 5000 c D2 > points.txt
cat points.txt | qdelaunay Fv Qt TO triangles.txt
pythonw order_triangles.py