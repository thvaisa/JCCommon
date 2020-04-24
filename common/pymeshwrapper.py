import pymesh
import numpy as np

#######################################
##Interacts with the PyMesh
class PyMeshObject:

    def __init__(self,fname):
        self.mesh = pymesh.meshio.load_mesh(fname)
        self.mesh.add_attribute("face_normal")
        self.mesh.add_attribute("face_area")
        self.mesh.add_attribute("face_centroid")

    def get_face_information(self,TYPE):
        faces = self.mesh.get_face_attribute("face_normal").astype(TYPE)
        areas = self.mesh.get_face_attribute("face_area").astype(TYPE)
        centroids = self.mesh.get_face_attribute("face_centroid").astype(TYPE)
        return faces,areas,centroids
       
    def get_vertices(self,TYPE):
        return self.mesh.vertices.astype(TYPE)

    def get_faces(self,TYPE):
        return self.mesh.faces.astype(TYPE)

    def n_faces(self):
        return self.mesh.get_face_attribute("face_normal").shape[0]

#######################################
