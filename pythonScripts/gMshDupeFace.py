import math

import gmsh

gmsh.initialize()

gmsh.open("/Users/mcgurn/chrestScratch/ablateInputs/geom/slabMotor.3D.hex.msh")

# Load in the element
[allElementTypes, allElementTags, allNodeTags] = gmsh.model.mesh.getElements(3)

elementType = allElementTypes[0]
elementTags = allElementTags[0]
nodeTags = allNodeTags[0]

# Get the type info
typeInfo = gmsh.model.mesh.getElementProperties(elementType)

faceHistory = {}

# define the outward facing faces
faces = [
    [1, 2, 6, 5],
    [0, 4, 7, 3],
    [2, 3, 7, 6],
    [1, 5, 4, 0],
    [5, 6, 7, 4],
    [0, 3, 2, 1]]

for eleTag in elementTags:
    ele = gmsh.model.mesh.getElement(eleTag)
    nodes = ele[1]

    for face in faces:
        faceNodes = [];
        for n in face:
            faceNodes.append(nodes[n])

        faceNodes.sort()
        nodeKey = str(faceNodes)

        if nodeKey in faceHistory:
            faceHistory[nodeKey] += 1
        else:
            faceHistory[nodeKey] = 1

for nodeKey in faceHistory:
    if faceHistory[nodeKey] > 2:
        print("gMshDupeFace Duplicate: ", nodeKey)

