from string import find
def _read_triangulation(fd):
    """Read the generated triangulation, NOT the outline."""

    delimiter = " "

    # loading the point info
    line = fd.readline()
    fragments = line.split()
    if fragments == []:
        NumOfVertices = 0
        NumOfVertAttributes = 0
    else:
        NumOfVertices = fragments[0]
        NumOfVertAttributes = fragments[1]
    points = []
    pointattributes = []
    for index in range(int(NumOfVertices)):
        fragments = fd.readline().split()
        fragments.pop(0)        # pop off the index

        # pop the x & y off so we're left with a list of attributes
        vert = [float(fragments.pop(0)), float(fragments.pop(0))]
        points.append(vert)
        apointattributes  = []
        for fragment in fragments:
            apointattributes.append(float(fragment))
        if apointattributes != []:
            pointattributes.append(apointattributes)

    # loading the point title info
    line = fd.readline()
    vertTitle = []
    for index in range(int(NumOfVertAttributes)):
        fragments = fd.readline().strip()
        vertTitle.append(fragments)

    # loading the triangle info
    line = fd.readline()
    fragments = line.split()
    NumOfTriangles = fragments[0]
    triangles = []
    triangleattributes = []
    triangleneighbors = []
    for index in range(int(NumOfTriangles)):
        line = fd.readline()
        line.strip()            # so we can get the region string
        fragments = line.split()
        fragments.pop(0)        # pop off the index

        tri = [int(fragments[0]), int(fragments[1]), int(fragments[2])]
        triangles.append(tri)
        neighbors = [int(fragments[3]), int(fragments[4]), int(fragments[5])]
        triangleneighbors.append(neighbors)
        stringtag = " ".join(fragments[6:]) # after pop the first
        triangleattributes.append(stringtag)

    # loading the segment info
    line = fd.readline()
    fragments = line.split()
    NumOfSegments = fragments[0]
    segments = []
    segmenttags = []
    for index in range(int(NumOfSegments)):
        line = fd.readline()
        line.strip()            # to get the segment string
        fragments = line.split()
        fragments.pop(0)        #pop off the index
        seg = [int(fragments[0]), int(fragments[1])]
        segments.append(seg)
        stringtag = " ".join(fragments[2:])
        segmenttags.append(stringtag)

    meshDict = {}
    meshDict['vertices'] = points
    if pointattributes == []:
        meshDict['vertex_attributes'] = None
    else:
        meshDict['vertex_attributes'] = pointattributes
    meshDict['triangles'] = triangles
    meshDict['triangle_tags'] = triangleattributes
    meshDict['triangle_neighbors'] = triangleneighbors
    meshDict['segments'] = segments
    meshDict['segment_tags'] = segmenttags
    meshDict['vertex_attribute_titles'] = vertTitle

    return meshDict

def _read_triangulation_old(fd):
    """Read the generated triangulation, NOT the outline."""

    delimiter = " "

    # loading the point info
    line = fd.readline()
    fragments = line.split()
    if fragments == []:
        NumOfVertices = 0
        NumOfVertAttributes = 0
    else:
        NumOfVertices = fragments[0]
        NumOfVertAttributes = fragments[1]
    points = []
    pointattributes = []
    for index in range(int(NumOfVertices)):
        fragments = fd.readline().split()
        fragments.pop(0)        # pop off the index

        # pop the x & y off so we're left with a list of attributes
        vert = [float(fragments.pop(0)), float(fragments.pop(0))]
        points.append(vert)
        apointattributes  = []
        for fragment in fragments:
            apointattributes.append(float(fragment))
        if apointattributes != []:
            pointattributes.append(apointattributes)

    # loading the point title info
    line = fd.readline()
    vertTitle = []
    for index in range(int(NumOfVertAttributes)):
        fragments = fd.readline().strip()
        vertTitle.append(fragments)

    # loading the triangle info
    line = fd.readline()
    fragments = line.split()
    NumOfTriangles = fragments[0]
    triangles = []
    triangleattributes = []
    triangleneighbors = []
    for index in range(int(NumOfTriangles)):
        line = fd.readline()
        line.strip()            # so we can get the region string
        fragments = line.split()
        fragments.pop(0)        # pop off the index

        tri = [int(fragments[0]), int(fragments[1]), int(fragments[2])]
        triangles.append(tri)
        neighbors = [int(fragments[3]), int(fragments[4]), int(fragments[5])]
        triangleneighbors.append(neighbors)
        for x in range(7):  # remove index [<vertex #>] [<neigbouring tri #>]
            line = line[find(line, delimiter):]     # remove index
            line = line.lstrip()
        stringtag = line.strip()
        triangleattributes.append(stringtag)

    # loading the segment info
    line = fd.readline()
    fragments = line.split()
    NumOfSegments = fragments[0]
    segments = []
    segmenttags = []
    for index in range(int(NumOfSegments)):
        line = fd.readline()
        line.strip()            # to get the segment string
        fragments = line.split()
        fragments.pop(0)        #pop off the index
        seg = [int(fragments[0]), int(fragments[1])]
        segments.append(seg)
        line = line[find(line, delimiter):] # remove index
        line = line.lstrip()
        line = line[find(line, delimiter):] # remove x
        line = line.lstrip()
        line = line[find(line, delimiter):] # remove y
        stringtag = line.strip()
        segmenttags.append(stringtag)

    meshDict = {}
    meshDict['vertices'] = points
    if pointattributes == []:
        meshDict['vertex_attributes'] = None
    else:
        meshDict['vertex_attributes'] = pointattributes
    meshDict['triangles'] = triangles
    meshDict['triangle_tags'] = triangleattributes
    meshDict['triangle_neighbors'] = triangleneighbors
    meshDict['segments'] = segments
    meshDict['segment_tags'] = segmenttags
    meshDict['vertex_attribute_titles'] = vertTitle

    return meshDict


if __name__ == "__main__":
   fd=open("merimbula_43200.tsh")

   import profile
   profile.run("_read_triangulation(fd)", "newread.txt")
   fd.seek(0)

   import pstats
   p = pstats.Stats("newread.txt")
   p.sort_stats("cumulative").print_stats()

   profile.run("_read_triangulation_old(fd)", "oldread.txt")
   q = pstats.Stats("oldread.txt")
   q.sort_stats("cumulative").print_stats()

