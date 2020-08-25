import scipy.spatial,numpy,time

# Create some dummy data
y_array = numpy.random.random(10000).reshape(100,100)
x_array = numpy.random.random(10000).reshape(100,100)

points = numpy.random.random(10000).reshape(2,5000)
# Shoe-horn existing data for entry into KDTree routines
combined_x_y_arrays = numpy.dstack([y_array.ravel(),x_array.ravel()])[0]
points_list = list(points.transpose())


def do_kdtree(combined_x_y_arrays,points):
    mytree = scipy.spatial.cKDTree(combined_x_y_arrays)
    dist, indexes = mytree.query(points)
    return indexes

start = time.time()
results2 = do_kdtree(combined_x_y_arrays,points_list)
end = time.time()
print results2
print 'Completed in: ',end-start
