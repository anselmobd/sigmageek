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
    # print('A B C')
    # pprint(A)
    # pprint(B)
    # pprint(C)
    u = B - A
    v = C - A
    # print('u v')
    # pprint(u)
    # pprint(v)
    n = np.cross(u, v)
    # print('multipliers')
    # pprint(n)
    (a, b, c) = n
    if not n.any():
        raise ValueError
    d = np.dot(n, A)
    # print('d')
    # pprint(d)
    result = (a, b, c, d)
    return result


def in_plane(plane, point):
    # print('inplane')
    # pprint(plane)
    # print(point)
    calc = 0
    # print(calc)
    for pl, po in zip(plane[:3], point.astuple()):
        calc += pl * po
        # print(calc)
    return round(calc, 3) == round(plane[3], 3)


def count_tetrahedrons(points):
    count = 0
    for comb in combinations(points, 4):
        try:
            plane = get_plane(comb[0], comb[1], comb[2])
        except ValueError:
            continue
        if not in_plane(plane, comb[3]):
            count += 1
    return count



if __name__ == '__main__':
    s1 = Vector3D(0, 0, 0)
    s2 = Vector3D(1, 1, 1)
    s3 = Vector3D(1, 2, 3)

    # print(s2)
    # s3 = s2.copy()
    # s3.x = 11
    # print(s3)
    # print(s2)

    # print('distance')
    # print(s1.distance(s2))

    # print('unitvec_rads')
    # print(0, 0, s1.unitvec_rads(0, 0))
    # print(30, 0, s1.unitvec_rads(pi/2/3, 0))
    # print(90, 0, s1.unitvec_rads(pi/2, 0))
    # print(90, 90, s1.unitvec_rads(pi/2, pi/2))
    # print('unitvec_degs')
    # print(90, 90, s1.unitvec_degs(90, 90))
    # print(45, 45, s1.unitvec_degs(45, 45))

    # print('vec_rads_length')
    # print(45, 45, 10, s1.vec_degs_length(45, 45, 10))

    ro = RegularOctahedron(0, 0, 0, 1)
    # print('vertices')
    # pprint(ro.vertices)
    # pprint(ro.faces)
    # pprint(ro.centers)

    points = [
        *ro.vertices.values(),
        *ro.centers,
    ]
    pprint(points)

    print(count_tetrahedrons(points))

    # try:
    #     plane = get_plane(s1, s2, s3)
    #     pprint(plane)
    # except ValueError:
    #     print('Not a plane')

    # print('face')
    # s1 = Vector3D(0.5, 0.5, 0.0)
    # s2 = Vector3D(-0.5, 0.5, 0.0)
    # s3 = Vector3D(0.0, 0.0, 0.70711)
    # s4 = Vector3D(0.0, 0.33333, 0.2357)
    # try:
    #     plane = get_plane(s1, s2, s3)
    #     pprint(plane)
    # except ValueError:
    #     print('Not a plane')
    # print(in_plane(plane, s4))
    # print(s4)

    # print('manual')
    # s1 = Vector3D(1, 0, 1)
    # s2 = Vector3D(0, 1, 1)
    # s3 = Vector3D(0, 0, 1)
    # s4 = Vector3D(1, 1, 0)
    # s5 = Vector3D(1, 1, 1)
    # try:
    #     plane = get_plane(s1, s2, s3)
    #     pprint(plane)
    # except ValueError:
    #     print('Not a plane')

    # print(in_plane(plane, s1))
    # print(in_plane(plane, s2))
    # print(in_plane(plane, s3))
    # print(in_plane(plane, s4))
    # print(in_plane(plane, s5))
