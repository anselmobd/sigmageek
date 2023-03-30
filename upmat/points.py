#!/usr/bin/env python3.8

import numbers
import numpy as np
import sympy as sp
from itertools import combinations
from math import (
    asin,
    cos,
    dist,
    pi,
    sin,
    sqrt,
)
from pprint import (
    pformat,
    pprint,
)


class Vector3D:

    def __init__(self, x=0, y=0, z=0):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

        self.round = 5

    def __str__(self) -> str:
        return pformat((
            round(self.x, self.round),
            round(self.y, self.round),
            round(self.z, self.round),
        ))

    __repr__ = __str__

    def __mul__(self, mult):
        if isinstance(mult, numbers.Number):
            return Vector3D(
                self.x * mult,
                self.y * mult,
                self.z * mult,
            )
        else:
            raise TypeError

    __rmul__ = __mul__

    def __add__(self, vector):
        if isinstance(vector, Vector3D):
            return Vector3D(
                self.x + vector.x,
                self.y + vector.y,
                self.z + vector.z,
            )
        else:
            raise TypeError

    def astuple(self) -> str:
        return (
            self.x,
            self.y,
            self.z,
        )

    def copy(self) -> str:
        return Vector3D(
            self.x,
            self.y,
            self.z,
        )

    def distance(self, p2):
        try:
            return dist(
                self.astuple(),
                p2.astuple(),
            )
        except Exception:
            raise TypeError
    
    def degree2rad(self, degree):
        """Degrees to radians"""
        return degree * pi / 180

    def unitvec_rads(self, yaw, pitch):
        """Unit vector generated from yaw and pitch angles"""
        try:
            return Vector3D(
                cos(yaw) * cos(pitch),
                sin(yaw) * cos(pitch),
                sin(pitch),
            )
        except Exception:
            raise TypeError

    def unitvec_degs(self, yaw, pitch):
        """unitvec_rads working with degrees"""
        return self.unitvec_rads(
            self.degree2rad(yaw),
            self.degree2rad(pitch),
        )

    def vec_rads_length(self, yaw, pitch, length=1):
        """Vector generated from unitvec_rads times length"""
        return self.unitvec_rads(yaw, pitch) * length
        
    def vec_degs_length(self, yaw, pitch, length=1):
        """Vector generated from unitvec_degs times length"""
        return self.unitvec_degs(yaw, pitch) * length


def deg2rad(degree):
    """Degrees to radians"""
    return degree * pi / 180


def rad2deg(radian):
    """Radians to degrees"""
    return radian * 180 / pi


def triangle_bcB(b, c, B):
    Br = deg2rad(B)
    Cr = asin(c * sin(B) / b)
    Ar = pi - Br - Cr
    a = b * sin(Ar) / sin(Br)
    C = rad2deg(Cr)
    A = rad2deg(Ar)
    return (a, b, c, A, B, C)


def triangle_center(va, vb, vc):
    v = []
    for a, b, c in zip(va.astuple(), vb.astuple(), vc.astuple()):
        v.append((a + b + c) / 3)
    return Vector3D(*v)

class RegularOctahedron():
    """Regular polyhedron formed by 12 edges, 6 vertices, and 8 faces
    """

    def __init__(self, x, y, z, edge_length):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.edge_length = float(edge_length)

        self.vertices = self.calc_vertices()
        self.faces = self.calc_faces()
        self.centers = self.calc_faces_centers()

    def calc_vertices(self):
        half = self.edge_length / 2
        v = {}
        v[1] = Vector3D(half, half, 0)
        v[2] = Vector3D(-half, half, 0)
        v[3] = Vector3D(-half, -half, 0)
        v[4] = Vector3D(half, -half, 0)

        origin_v1 = v[1].distance(Vector3D())
        v1_top = self.edge_length
        v1_top_angle = 90
        triang = triangle_bcB(v1_top, origin_v1, v1_top_angle)
        top_origin = triang[2]
        v[5] = Vector3D(0, 0, top_origin)

        v[6] = v[5] * -1
        return v

    def calc_faces(self):
        return {
            1: [self.vertices[1], self.vertices[2], self.vertices[5]],
            2: [self.vertices[2], self.vertices[3], self.vertices[5]],
            3: [self.vertices[3], self.vertices[4], self.vertices[5]],
            4: [self.vertices[4], self.vertices[1], self.vertices[5]],
            5: [self.vertices[1], self.vertices[2], self.vertices[6]],
            6: [self.vertices[2], self.vertices[3], self.vertices[6]],
            7: [self.vertices[3], self.vertices[4], self.vertices[6]],
            8: [self.vertices[4], self.vertices[1], self.vertices[6]],
        }

    def calc_faces_centers(self):
        centers = []
        for f_id in self.faces:
            centers.append(
                triangle_center(
                    *self.faces[f_id]))
        return centers


def get_plane(v1, v2, v3):
    A = np.array(v1.astuple())
    B = np.array(v2.astuple())
    C = np.array(v3.astuple())
    u = B - A
    v = C - A
    n = np.cross(u, v)
    (a, b, c) = n
    if not n.any():
        raise ValueError
    d = np.dot(n, A)
    result = (a, b, c, d)
    return result


def in_plane(plane, point):
    calc = 0
    for pl, po in zip(plane[:3], point.astuple()):
        calc += pl * po
    return round(calc, 3) == round(plane[3], 3)

def sorted_edges_lengths(t):
    lengths = []
    lengths.append(round(t[0].distance(t[1]), 3))
    lengths.append(round(t[0].distance(t[2]), 3))
    lengths.append(round(t[0].distance(t[3]), 3))
    lengths.append(round(t[1].distance(t[2]), 3))
    lengths.append(round(t[1].distance(t[3]), 3))
    lengths.append(round(t[2].distance(t[3]), 3))
    return tuple(sorted(lengths))

def count_tetrahedrons(points):
    tetra = set()
    count = 0
    for comb in combinations(points, 4):
        try:
            plane = get_plane(comb[0], comb[1], comb[2])
        except ValueError:
            continue
        if not in_plane(plane, comb[3]):
            count += 1
            key = sorted_edges_lengths(comb)
            tetra.add(key)
    return count, len(tetra)


if __name__ == '__main__':
    ro = RegularOctahedron(0, 0, 0, 1)

    points = [
        *ro.vertices.values(),
        *ro.centers,
    ]

    count, diff = count_tetrahedrons(points)
    print(count)
