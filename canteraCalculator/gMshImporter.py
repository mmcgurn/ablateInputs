import math

import gmsh

gmsh.initialize()

gmsh.open("/Users/mcgurn/chrestScratch/ablateInputs/geom/slabMotor.3D.msh")

# Load in the element
elementIds = gmsh.model.mesh.getElementsByCoordinates(0.044636811516239397, 0.0014524456734120389,-0.007080303862734838, 3, False)

# Print the element information
for elementId in elementIds:
    print("Element: ", elementId)
    # Get the element
    element = gmsh.model.mesh.getElement(elementId)
    # Print the nodes
    for n in element[1]:
        print("\t", str(gmsh.model.mesh.getNode(n)[0][0]), ", ", str(gmsh.model.mesh.getNode(n)[0][1]), ", ", str(gmsh.model.mesh.getNode(n)[0][2]))


element = gmsh.model.mesh.getElement(2419)

# get the element nodes
nodes = element[1]
print(nodes)

# define the outward facing faces
faces = [
    [1, 2, 6, 5],
    [0, 4, 7, 3],
    [2, 3, 7, 6],
    [1, 5, 4, 0],
    [5, 6, 7, 4],
    [0, 3, 2, 1]]

# Compute the area on each face
areas = [[]]
areaMags = []
centers = [[]]
areasFromNorms = [[]]
areasFromOrgs = [[]]
areasFromOrgsMags = []
norms = [[]]
areaSum = [0, 0, 0]
areaSumNorm = [0, 0, 0]
centroid = [0, 0, 0]

for n in nodes:
    node = gmsh.model.mesh.getNode(n)
    for d in range(3):
        centroid[d] += node[0][d] / (len(nodes))


def cross(a, b):
    return [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1]-a[1]*b[0]]

def norm(a):
    mag = math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])
    return [a[0] / mag, a[1] / mag, a[2] / mag]

for face in faces:
    print("Face: ", face)
    a = gmsh.model.mesh.getNode(nodes[face[0]])
    b = gmsh.model.mesh.getNode(nodes[face[1]])
    c = gmsh.model.mesh.getNode(nodes[face[2]])
    ab = [b[0][0] - a[0][0], b[0][1] - a[0][1], b[0][2] - a[0][2]]
    ac = [c[0][0] - a[0][0], c[0][1] - a[0][1], c[0][2] - a[0][2]]
    area0 = cross(ab, ac)

    # assume the norm is from first element
    norm0 = norm(area0)
    norms.append(norm0)

    a = gmsh.model.mesh.getNode(nodes[face[0]])
    b = gmsh.model.mesh.getNode(nodes[face[2]])
    c = gmsh.model.mesh.getNode(nodes[face[3]])
    ab = [b[0][0] - a[0][0], b[0][1] - a[0][1], b[0][2] - a[0][2]]
    ac = [c[0][0] - a[0][0], c[0][1] - a[0][1], c[0][2] - a[0][2]]
    area1 = cross(ab, ac)
    area = ([0.5*(area1[0]+area0[0]), 0.5*(area1[1]+area0[1]), 0.5*(area1[2]+area0[2])])

    areaSum[0] += area[0]
    areaSum[1] += area[1]
    areaSum[2] += area[2]
    areas.append(area)
    areaMag = math.sqrt(area[0]*area[0] + area[1]*area[1]+ area[2]*area[2])
    areaMags.append(areaMag)

    # Compute using the norm
    areaNorm = [norm0[0]*areaMag, norm0[1]*areaMag, norm0[2]*areaMag]
    areasFromNorms.append(areaNorm)
    areaSumNorm[0] += areaNorm[0]
    areaSumNorm[1] += areaNorm[1]
    areaSumNorm[2] += areaNorm[2]

    center = [0, 0, 0]
    for n in face:
        node = gmsh.model.mesh.getNode(nodes[n])
        for d in range(3):
            center[d] += node[0][d]/(len(face))

    centers.append(center)

    # Compute the area around summing each part
    areaFromOrg = [0, 0, 0]
    segments = [[0, 1], [1, 2], [2, 3], [3, 0]]
    for segment in segments:
        oA = gmsh.model.mesh.getNode(nodes[face[segment[0]]])[0]
        oB = gmsh.model.mesh.getNode(nodes[face[segment[1]]])[0]
        areaSegment = cross(oA, oB)
        areaFromOrg[0] += areaSegment[0]
        areaFromOrg[1] += areaSegment[1]
        areaFromOrg[2] += areaSegment[2]

    areasFromOrgs.append(areaFromOrg)
    areasFromOrgsMags.append(math.sqrt(areaFromOrg[0]*areaFromOrg[0] + areaFromOrg[1]*areaFromOrg[1]+ areaFromOrg[2]*areaFromOrg[2]))


print("areas: ", areas)
print("centers: ", centers)
print("areaSum: ", areaSum)
print("areaMags: ", areaMags)
print("areasFromOrgs: ", areasFromOrgs)
print("areasFromOrgsMags: ", areasFromOrgsMags)

print("areaSumNorm: ", areaSumNorm)
print("centroid: " , centroid)
