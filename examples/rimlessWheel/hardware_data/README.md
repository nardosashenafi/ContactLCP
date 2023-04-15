To convert rosbag to files
rosbag play mybag.bag
rostopic echo /foo/position/x > output_x.txt
rostopic echo /foo/position/y > output_y.txt
rostopic echo /foo/position/z > output_z.txt