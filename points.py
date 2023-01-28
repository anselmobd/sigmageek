#!/usr/bin/env python3.8

import numbers
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

if __name__ == '__main__':
    s1 = Vector3D(0, 0, 0)
    s2 = Vector3D(1, 1, 1)

    print(s2)
    s3 = s2.copy()
    s3.x = 11
    print(s3)
    print(s2)

    print('distance')
    print(s1.distance(s2))

    print('unitvec_rads')
    print(0, 0, s1.unitvec_rads(0, 0))
    print(30, 0, s1.unitvec_rads(pi/2/3, 0))
    print(90, 0, s1.unitvec_rads(pi/2, 0))
    print(90, 90, s1.unitvec_rads(pi/2, pi/2))
    print('unitvec_degs')
    print(90, 90, s1.unitvec_degs(90, 90))
    print(45, 45, s1.unitvec_degs(45, 45))

    print('vec_rads_length')
    print(45, 45, 10, s1.vec_degs_length(45, 45, 10))

    ro = RegularOctahedron(0, 0, 0, 1)
    print('vertices')
    pprint(ro.vertices)
    pprint(ro.faces)
