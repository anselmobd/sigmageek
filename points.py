#!/usr/bin/env python3.8

import numbers
from math import (
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

    def tuple(self) -> str:
        return (
            self.x,
            self.y,
            self.z,
        )

    def distance(self, p2):
        try:
            return dist(
                self.tuple(),
                p2.tuple(),
            )
        except Exception:
            raise TypeError
    
    def unitary_angles(self, yaw, pitch):
        try:
            return Vector3D(
                cos(yaw) * cos(pitch),
                sin(yaw) * cos(pitch),
                sin(pitch),
            )
        except Exception:
            raise TypeError

    def angles_distance(self, yaw, pitch, distance=1):
        return self + self.unitary_angles(yaw, pitch) * distance
        

class RegularOctahedron():
    """Regular polyhedron formed by 12 edges, 6 vertices, and 8 faces
    """

    def __init__(self, x, y, z, edge_length):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.edge_length = float(edge_length)

        self.calc_vertices()

    def calc_vertices(self):
        half = self.edge_length / 2
        v1 = (half, half, 0)
        v2 = (-half, half, 0)
        v3 = (-half, -half, 0)
        v4 = (half, -half, 0)
        v5 = Vector3D(*v1).angles_distance(0, pi/2).tuple()
        v6 = (Vector3D(*v5) * -1).tuple()
        self.vertices= [
            v1,
            v2,
            v3,
            v4,
            v5,
            v6,
        ]


if __name__ == '__main__':
    s1 = Vector3D(0, 0, 0)
    s2 = Vector3D(1, 1, 1)

    print('distance')
    print(s1.distance(s2))

    print('unitary_angles')
    print(0, 0, s1.unitary_angles(0, 0))
    print(30, 0, s1.unitary_angles(pi/2/3, 0))
    print(90, 0, s1.unitary_angles(pi/2, 0))
    print(90, 90, s1.unitary_angles(pi/2, pi/2))

    print('angles_distance')
    print(0, 0, 10, s1.angles_distance(0, 0, 10))

    ro = RegularOctahedron(0, 0, 0, 1)
    print('vertices')
    pprint(ro.vertices)
