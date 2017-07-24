from imageprep import imageprep

shift_file = 'night1_roughshifts.txt'

imageprep('AWI0000fgw', 'H', shift_file, 108, 5)
imageprep('AWI0000enp', 'H', shift_file, 115, 5)
imageprep('AWI0000enp', 'Y', shift_file, 122, 5)
imageprep('AWI00008ye', 'Y', shift_file, 129, 5)
imageprep('AWI00008ye', 'H', shift_file, 135, 5)
imageprep('AWI00062l8', 'H', shift_file, 143, 5)
imageprep('AWI00062l8', 'Y', shift_file, 149, 5)
imageprep('AWI000601n', 'Y', shift_file, 156, 5)
imageprep('AWI000601n', 'H', shift_file, 162, 5)
imageprep('AWI0005wd8', 'H', shift_file, 169, 5)
imageprep('AWI0005wd8', 'Y', shift_file, 175, 5)
imageprep('AWI0005yog', 'Y', shift_file, 182, 5)
imageprep('AWI0005yog', 'H', shift_file, 188, 5)
imageprep('AWI0005yp9', 'H', shift_file, 196, 5)
imageprep('AWI0005yp9', 'Y', shift_file, 202, 5)
imageprep('AWI0005ypo', 'Y', shift_file, 209, 5)
imageprep('AWI0005ypo', 'H', shift_file, 215, 5)
imageprep('AWI0005ypv', 'H', shift_file, 221, 5)
imageprep('AWI0005ypv', 'Y', shift_file, 227, 5)
imageprep('AWI0005wg8', 'Y', shift_file, 234, 5)
imageprep('AWI0005wg8', 'H', shift_file, 240, 5)
imageprep('AWI0005whl', 'H', shift_file, 246, 5)
imageprep('AWI0005whl', 'Y', shift_file, 253, 5)
imageprep('AWI0005wjy', 'Y', shift_file, 260, 5)
imageprep('AWI0005wjy', 'H', shift_file, 266, 5)

print "Done"
