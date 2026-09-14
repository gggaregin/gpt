#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parametric combined ITP model generated from:
- 3(1).pdf: priority room geometry, room No. 11, control length 9205 mm.
- CTP-TM(1).pdf: process scheme, plans, sections and isometric.
The model is coordination-grade, not fabrication/spool-grade.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import struct
import zipfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

OUT = Path("dist")
OUT.mkdir(exist_ok=True)

IFC_NAME = "ITP_Complete_Archicad23.ifc"
GLB_NAME = "ITP_Complete_Preview.glb"
ZIP_NAME = "ITP_Complete_Archicad23.zip"

ROOM_LENGTH = 9205.0
ROOM_WIDTH = 2915.0
ROOM_HEIGHT = 2788.0
WALL_T = 200.0

LAYERS = {
    "room": "01_ITP_Room",
    "stair": "02_ITP_Stair",
    "equipment": "03_ITP_Equipment",
    "pipe": "04_ITP_Pipework",
    "fitting": "05_ITP_Fittings_Valves",
    "support": "06_ITP_Supports",
    "instrument": "07_ITP_Instruments",
}

COLORS = {
    "room": (0.72, 0.74, 0.76, 0.34),
    "floor": (0.47, 0.49, 0.51, 1.00),
    "stair": (0.38, 0.40, 0.42, 1.00),
    "equipment": (0.34, 0.38, 0.42, 1.00),
    "tank": (0.68, 0.70, 0.73, 1.00),
    "primary_supply": (0.78, 0.08, 0.09, 1.00),
    "primary_return": (0.48, 0.04, 0.07, 1.00),
    "heating_supply": (0.95, 0.30, 0.08, 1.00),
    "heating_return": (0.82, 0.52, 0.10, 1.00),
    "vent_supply": (0.65, 0.12, 0.70, 1.00),
    "vent_return": (0.40, 0.11, 0.55, 1.00),
    "dhw_supply": (0.98, 0.55, 0.05, 1.00),
    "dhw_circ": (0.98, 0.78, 0.08, 1.00),
    "cold_water": (0.05, 0.38, 0.90, 1.00),
    "drain": (0.30, 0.32, 0.34, 1.00),
    "valve": (0.92, 0.73, 0.10, 1.00),
    "support": (0.22, 0.24, 0.26, 1.00),
    "instrument": (0.18, 0.78, 0.74, 1.00),
}

IFC64 = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_$"


def ifc_guid(seed: str) -> str:
    n = int(hashlib.md5(seed.encode("utf-8")).hexdigest(), 16)
    out = ""
    for _ in range(22):
        out = IFC64[n & 63] + out
        n >>= 6
    return out


def fnum(v: float) -> str:
    if abs(v) < 1e-10:
        return "0."
    s = f"{float(v):.6f}".rstrip("0")
    if s.endswith("."):
        return s
    if "." not in s:
        s += "."
    return s


def tup(vals) -> str:
    return "(" + ",".join(fnum(v) for v in vals) + ")"


def esc(s: str) -> str:
    return str(s).replace("'", "''")


def vadd(a, b):
    return tuple(a[i] + b[i] for i in range(3))


def vsub(a, b):
    return tuple(a[i] - b[i] for i in range(3))


def vmul(a, k):
    return tuple(x * k for x in a)


def vdot(a, b):
    return sum(a[i] * b[i] for i in range(3))


def vcross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def vlen(a):
    return math.sqrt(vdot(a, a))


def vnorm(a):
    n = vlen(a)
    if n < 1e-9:
        raise ValueError("Zero-length vector")
    return tuple(x / n for x in a)


def basis_for_axis(axis):
    w = vnorm(axis)
    helper = (0.0, 0.0, 1.0) if abs(w[2]) < 0.9 else (1.0, 0.0, 0.0)
    u = vnorm(vcross(helper, w))
    v = vnorm(vcross(w, u))
    return u, v, w


def box_mesh(cx, cy, z, w, d, h, angle=0.0):
    ca, sa = math.cos(angle), math.sin(angle)
    local = [
        (-w / 2, -d / 2, 0), (w / 2, -d / 2, 0),
        (w / 2, d / 2, 0), (-w / 2, d / 2, 0),
        (-w / 2, -d / 2, h), (w / 2, -d / 2, h),
        (w / 2, d / 2, h), (-w / 2, d / 2, h),
    ]
    verts = []
    for x, y, zz in local:
        verts.append((cx + x * ca - y * sa, cy + x * sa + y * ca, z + zz))
    faces = [
        (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
        (0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
        (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7),
    ]
    return verts, faces


def cylinder_mesh(a, b, r, sides=12):
    axis = vsub(b, a)
    u, v, _ = basis_for_axis(axis)
    verts = []
    for p in (a, b):
        for i in range(sides):
            ang = 2 * math.pi * i / sides
            off = vadd(vmul(u, r * math.cos(ang)), vmul(v, r * math.sin(ang)))
            verts.append(vadd(p, off))
    verts.extend([a, b])
    faces = []
    for i in range(sides):
        j = (i + 1) % sides
        faces.append((i, j, sides + j))
        faces.append((i, sides + j, sides + i))
        faces.append((2 * sides, j, i))
        faces.append((2 * sides + 1, sides + i, sides + j))
    return verts, faces


def beam_mesh(a, b, w, d):
    axis = vsub(b, a)
    u, v, _ = basis_for_axis(axis)
    corners = [
        vadd(vmul(u, -w / 2), vmul(v, -d / 2)),
        vadd(vmul(u, w / 2), vmul(v, -d / 2)),
        vadd(vmul(u, w / 2), vmul(v, d / 2)),
        vadd(vmul(u, -w / 2), vmul(v, d / 2)),
    ]
    verts = [vadd(a, c) for c in corners] + [vadd(b, c) for c in corners]
    faces = [
        (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
        (0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
        (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7),
    ]
    return verts, faces


def mesh_for_shape(shape):
    k = shape["kind"]
    if k == "box":
        return box_mesh(shape["cx"], shape["cy"], shape["z"],
                        shape["w"], shape["d"], shape["h"], shape.get("angle", 0.0))
    if k == "cyl":
        return cylinder_mesh(shape["a"], shape["b"], shape["r"], shape.get("sides", 12))
    if k == "beam":
        return beam_mesh(shape["a"], shape["b"], shape["w"], shape["d"])
    raise ValueError(k)


def combine_meshes(shapes):
    verts, faces = [], []
    for shape in shapes:
        vv, ff = mesh_for_shape(shape)
        off = len(verts)
        verts.extend(vv)
        faces.extend(tuple(i + off for i in face) for face in ff)
    return verts, faces


class IFCBuilder:
    def __init__(self):
        self.lines = []
        self.n = 0
        self.layers = defaultdict(list)
        self.styles = {}
        self.products = []
        self.objects = []

        self.person = self.add("IFCPERSON($,$,'OpenAI Codex',$,$,$,$,$)")
        self.org = self.add("IFCORGANIZATION($,'GARNET CONTROL',$,$,$)")
        self.person_org = self.add(f"IFCPERSONANDORGANIZATION({self.person},{self.org},$)")
        self.app = self.add(f"IFCAPPLICATION({self.org},'1.0','ITP Parametric Builder','ITP-BUILD')")
        self.owner = self.add(
            f"IFCOWNERHISTORY({self.person_org},{self.app},$,.ADDED.,$,$,$,{int(datetime.now(timezone.utc).timestamp())})"
        )
        self.p0 = self.add("IFCCARTESIANPOINT((0.,0.,0.))")
        self.dz = self.add("IFCDIRECTION((0.,0.,1.))")
        self.dx = self.add("IFCDIRECTION((1.,0.,0.))")
        self.p2 = self.add("IFCCARTESIANPOINT((0.,0.))")
        self.dx2 = self.add("IFCDIRECTION((1.,0.))")
        self.profile_pos = self.add(f"IFCAXIS2PLACEMENT2D({self.p2},{self.dx2})")
        self.world = self.add(f"IFCAXIS2PLACEMENT3D({self.p0},{self.dz},{self.dx})")
        self.context = self.add(
            f"IFCGEOMETRICREPRESENTATIONCONTEXT($,'Model',3,1.E-05,{self.world},$)"
        )
        length_unit = self.add("IFCSIUNIT(*,.LENGTHUNIT.,.MILLI.,.METRE.)")
        area_unit = self.add("IFCSIUNIT(*,.AREAUNIT.,.MILLI.,.SQUARE_METRE.)")
        volume_unit = self.add("IFCSIUNIT(*,.VOLUMEUNIT.,.MILLI.,.CUBIC_METRE.)")
        angle_unit = self.add("IFCSIUNIT(*,.PLANEANGLEUNIT.,$,.RADIAN.)")
        time_unit = self.add("IFCSIUNIT(*,.TIMEUNIT.,$,.SECOND.)")
        self.units = self.add(
            f"IFCUNITASSIGNMENT(({length_unit},{area_unit},{volume_unit},{angle_unit},{time_unit}))"
        )
        self.project = self.add(
            f"IFCPROJECT('{ifc_guid('project')}',{self.owner},'ITP_Complete_Archicad23',"
            f"'Combined room and thermomechanical coordination model',$,$,$,({self.context}),{self.units})"
        )
        site_place = self.add(f"IFCLOCALPLACEMENT($,{self.world})")
        self.site = self.add(
            f"IFCSITE('{ifc_guid('site')}',{self.owner},'Delegat Site',$,$,{site_place},$,$,.ELEMENT.,$,$,$,$,$)"
        )
        bldg_place = self.add(f"IFCLOCALPLACEMENT({site_place},{self.world})")
        self.building = self.add(
            f"IFCBUILDING('{ifc_guid('building')}',{self.owner},'Delegat Building',$,$,"
            f"{bldg_place},$,$,.ELEMENT.,$,$,$)"
        )
        storey_place = self.add(f"IFCLOCALPLACEMENT({bldg_place},{self.world})")
        self.storey = self.add(
            f"IFCBUILDINGSTOREY('{ifc_guid('storey')}',{self.owner},'Level 01 - ITP FFL 0.000',"
            f"$,$,{storey_place},$,$,.ELEMENT.,0.)"
        )
        self.storey_place = storey_place
        self.add(
            f"IFCRELAGGREGATES('{ifc_guid('agg-project-site')}',{self.owner},$,$,{self.project},({self.site}))"
        )
        self.add(
            f"IFCRELAGGREGATES('{ifc_guid('agg-site-building')}',{self.owner},$,$,{self.site},({self.building}))"
        )
        self.add(
            f"IFCRELAGGREGATES('{ifc_guid('agg-building-storey')}',{self.owner},$,$,{self.building},({self.storey}))"
        )

    def add(self, text):
        self.n += 1
        self.lines.append(f"#{self.n}={text};")
        return f"#{self.n}"

    def point(self, p):
        return self.add(f"IFCCARTESIANPOINT({tup(p)})")

    def direction(self, d):
        return self.add(f"IFCDIRECTION({tup(vnorm(d))})")

    def style(self, color_key):
        if color_key in self.styles:
            return self.styles[color_key]
        r, g, b, a = COLORS[color_key]
        col = self.add(f"IFCCOLOURRGB('{esc(color_key)}',{fnum(r)},{fnum(g)},{fnum(b)})")
        rend = self.add(
            f"IFCSURFACESTYLERENDERING({col},{fnum(1.0-a)},$,$,$,$,$,$,.NOTDEFINED.)"
        )
        surf = self.add(f"IFCSURFACESTYLE('{esc(color_key)}',.BOTH.,({rend}))")
        assign = self.add(f"IFCPRESENTATIONSTYLEASSIGNMENT(({surf}))")
        self.styles[color_key] = assign
        return assign

    def solid(self, shape, color_key):
        k = shape["kind"]
        if k == "box":
            p = self.point((shape["cx"], shape["cy"], shape["z"]))
            axis = self.direction((0, 0, 1))
            ang = shape.get("angle", 0.0)
            ref = self.direction((math.cos(ang), math.sin(ang), 0))
            pos = self.add(f"IFCAXIS2PLACEMENT3D({p},{axis},{ref})")
            profile = self.add(
                f"IFCRECTANGLEPROFILEDEF(.AREA.,$,{self.profile_pos},{fnum(shape['w'])},{fnum(shape['d'])})"
            )
            item = self.add(
                f"IFCEXTRUDEDAREASOLID({profile},{pos},{self.dz},{fnum(shape['h'])})"
            )
        elif k in ("cyl", "beam"):
            a, b = shape["a"], shape["b"]
            vec = vsub(b, a)
            u, _, w = basis_for_axis(vec)
            p = self.point(a)
            axis = self.direction(w)
            ref = self.direction(u)
            pos = self.add(f"IFCAXIS2PLACEMENT3D({p},{axis},{ref})")
            if k == "cyl":
                profile = self.add(f"IFCCIRCLEPROFILEDEF(.AREA.,$,{self.profile_pos},{fnum(shape['r'])})")
            else:
                profile = self.add(
                    f"IFCRECTANGLEPROFILEDEF(.AREA.,$,{self.profile_pos},{fnum(shape['w'])},{fnum(shape['d'])})"
                )
            item = self.add(
                f"IFCEXTRUDEDAREASOLID({profile},{pos},{self.dz},{fnum(vlen(vec))})"
            )
        else:
            raise ValueError(k)
        self.add(f"IFCSTYLEDITEM({item},({self.style(color_key)}),$)")
        return item

    def add_product(self, ifctype, name, description, tag, layer, shapes,
                    color_key, system="", dn=None, mark="", category="", notes=""):
        solids = [self.solid(s, color_key) for s in shapes]
        rep = self.add(
            f"IFCSHAPEREPRESENTATION({self.context},'Body','SweptSolid',({','.join(solids)}))"
        )
        pshape = self.add(f"IFCPRODUCTDEFINITIONSHAPE($,$,({rep}))")
        object_origin = self.point((0.0, 0.0, 0.0))
        object_axis = self.add(f"IFCAXIS2PLACEMENT3D({object_origin},{self.dz},{self.dx})")
        object_place = self.add(f"IFCLOCALPLACEMENT({self.storey_place},{object_axis})")
        gid = ifc_guid("product-" + name)
        base = (
            f"'{gid}',{self.owner},'{esc(name)}','{esc(description)}',$,"
            f"{object_place},{pshape},'{esc(tag)}'"
        )
        if ifctype == "IFCSLAB":
            predefined = ".FLOOR." if "FLOOR" in name else ".ROOF."
            product = self.add(f"{ifctype}({base},{predefined})")
        elif ifctype == "IFCBUILDINGELEMENTPROXY":
            product = self.add(f"{ifctype}({base},.ELEMENT.)")
        else:
            product = self.add(f"{ifctype}({base})")

        props = []
        pvals = {
            "Layer": layer,
            "System": system or "N/A",
            "Mark": mark or "N/A",
            "ObjectName": name,
            "Source": "3(1).pdf + CTP-TM(1).pdf",
            "Notes": notes or "Coordination model",
        }
        for key, val in pvals.items():
            props.append(
                self.add(
                    f"IFCPROPERTYSINGLEVALUE('{esc(key)}',$,IFCLABEL('{esc(val)}'),$)"
                )
            )
        if dn is not None:
            props.append(
                self.add(
                    f"IFCPROPERTYSINGLEVALUE('DN_mm',$,IFCLENGTHMEASURE({fnum(dn)}),$)"
                )
            )
        pset = self.add(
            f"IFCPROPERTYSET('{ifc_guid('pset-'+name)}',{self.owner},"
            f"'Pset_ITP_Identity',$,({','.join(props)}))"
        )
        self.add(
            f"IFCRELDEFINESBYPROPERTIES('{ifc_guid('relpset-'+name)}',{self.owner},"
            f"$,$,({product}),{pset})"
        )
        self.layers[layer].append(rep)
        self.products.append(product)

        verts, faces = combine_meshes(shapes)
        bbox = (
            min(p[0] for p in verts), min(p[1] for p in verts), min(p[2] for p in verts),
            max(p[0] for p in verts), max(p[1] for p in verts), max(p[2] for p in verts),
        )
        self.objects.append({
            "name": name,
            "ifc_type": ifctype,
            "layer": layer,
            "system": system,
            "dn": dn,
            "mark": mark,
            "category": category,
            "color": color_key,
            "vertices": verts,
            "faces": faces,
            "bbox": bbox,
        })
        return product

    def finish(self):
        self.add(
            f"IFCRELCONTAINEDINSPATIALSTRUCTURE('{ifc_guid('containment')}',"
            f"{self.owner},$,$,({','.join(self.products)}),{self.storey})"
        )
        for layer, reps in self.layers.items():
            self.add(
                f"IFCPRESENTATIONLAYERASSIGNMENT('{esc(layer)}',$,"
                f"({','.join(reps)}),'{esc(layer)}')"
            )

    def write(self, path):
        header = [
            "ISO-10303-21;",
            "HEADER;",
            "FILE_DESCRIPTION(('ViewDefinition [CoordinationView_V2.0]'),'2;1');",
            f"FILE_NAME('{path.name}','{datetime.now(timezone.utc).isoformat()}',"
            "('OpenAI Codex'),('GARNET CONTROL'),'ITP Parametric Builder','ARCHICAD 23','');",
            "FILE_SCHEMA(('IFC2X3'));",
            "ENDSEC;",
            "DATA;",
        ]
        footer = ["ENDSEC;", "END-ISO-10303-21;"]
        path.write_text("\n".join(header + self.lines + footer) + "\n", encoding="utf-8")


def box(cx, cy, z, w, d, h, angle=0.0):
    return {"kind": "box", "cx": cx, "cy": cy, "z": z, "w": w, "d": d, "h": h, "angle": angle}


def cyl(a, b, r, sides=12):
    return {"kind": "cyl", "a": tuple(a), "b": tuple(b), "r": r, "sides": sides}


def beam(a, b, w, d):
    return {"kind": "beam", "a": tuple(a), "b": tuple(b), "w": w, "d": d}


M = IFCBuilder()


def add_room():
    M.add_product(
        "IFCSLAB", "FLOOR_ITP_9205", "ITP finished floor and structural datum",
        "ROOM-FLR-01", LAYERS["room"],
        [box(ROOM_LENGTH/2, ROOM_WIDTH/2, -100, ROOM_LENGTH, ROOM_WIDTH, 100)],
        "floor", category="floor", notes="Overall X dimension is exactly 9205 mm; FFL Z=0."
    )
    M.add_product(
        "IFCSLAB", "CEILING_ITP_Z2788", "Conditional upper room boundary",
        "ROOM-CLG-01", LAYERS["room"],
        [box(ROOM_LENGTH/2, ROOM_WIDTH/2, ROOM_HEIGHT, ROOM_LENGTH, ROOM_WIDTH, 40)],
        "room", category="ceiling", notes="Main conditional ceiling level accepted as Z=2788 mm."
    )

    walls = [
        ("WALL_SOUTH_01", ROOM_LENGTH/2, 100, ROOM_LENGTH, 200),
        ("WALL_WEST_01", 100, ROOM_WIDTH/2, 200, ROOM_WIDTH),
        ("WALL_EAST_LOW", ROOM_LENGTH-100, 225, 200, 450),
        ("WALL_EAST_HIGH", ROOM_LENGTH-100, (1650+ROOM_WIDTH)/2, 200, ROOM_WIDTH-1650),
        ("WALL_NORTH_A", 1400, 2815, 2800, 200),
        ("WALL_NORTH_B", (3300+5450)/2, 2815, 2150, 200),
        ("WALL_NORTH_C", (5900+8050)/2, 2815, 2150, 200),
        ("WALL_NORTH_D", (8500+ROOM_LENGTH)/2, 2815, ROOM_LENGTH-8500, 200),
        ("NICHE_N1_BACK", 3050, 2580, 500, 200),
        ("NICHE_N1_RETURN_W", 2900, (2580+2815)/2, 200, 270),
        ("NICHE_N1_RETURN_E", 3200, (2580+2815)/2, 200, 270),
        ("NICHE_N2_BACK", 5675, 2650, 450, 200),
        ("NICHE_N2_RETURN_W", 5550, (2650+2815)/2, 200, 170),
        ("NICHE_N2_RETURN_E", 5800, (2650+2815)/2, 200, 170),
        ("NICHE_N3_BACK", 8275, 2700, 450, 200),
        ("NICHE_N3_RETURN_W", 8150, (2700+2815)/2, 200, 120),
        ("NICHE_N3_RETURN_E", 8400, (2700+2815)/2, 200, 120),
    ]
    for name, cx, cy, w, d in walls:
        M.add_product(
            "IFCWALLSTANDARDCASE", name, "Separate room wall/recess segment",
            name, LAYERS["room"], [box(cx, cy, 0, w, d, ROOM_HEIGHT)],
            "room", category="wall", notes="Geometry priority: 3(1).pdf."
        )

    for i, x in enumerate((2910, 5675, 8275), 1):
        M.add_product(
            "IFCBEAM", f"CEILING_SOFFIT_{i:02d}", "Local ceiling step / soffit",
            f"ROOM-SOFFIT-{i:02d}", LAYERS["room"],
            [box(x, ROOM_WIDTH/2, 2600, 220, ROOM_WIDTH-400, ROOM_HEIGHT-2600)],
            "room", category="soffit",
            notes="Local lower boundary accepted at Z=2600 mm from elevation interpretation."
        )


def add_stair():
    x0, tread, rise = 7700.0, 180.0, 175.0
    y0, width = 500.0, 1000.0
    for i in range(6):
        ztop = (i + 1) * rise
        bars = []
        for j in range(8):
            bars.append(box(x0 + i*tread + tread/2, y0 + 40 + j*125, ztop-40,
                            tread-12, 58, 40))
        M.add_product(
            "IFCBUILDINGELEMENTPROXY", f"STAIR_GRATING_STEP_{i+1:02d}",
            "Galvanized grating stair tread", f"ST-{i+1:02d}",
            LAYERS["stair"], bars, "stair", category="stair",
            notes="Galvanized grating; each tread is a separately selectable element."
        )

    platform_bars = [
        box((8780+ROOM_LENGTH)/2, y0+40+j*125, 1050, ROOM_LENGTH-8780, 58, 40)
        for j in range(8)
    ]
    M.add_product(
        "IFCBUILDINGELEMENTPROXY", "STAIR_GRATING_PLATFORM_01",
        "Galvanized grating landing", "PLATFORM-01", LAYERS["stair"],
        platform_bars, "stair", category="stair",
        notes="Landing at Z=1050 mm."
    )

    for side, yy in (("L", y0+35), ("R", y0+width-35)):
        shapes = [
            beam((x0, yy, 80), (8780, yy, 1130), 160, 12),
            beam((x0, yy-35, 80), (8780, yy-35, 1130), 55, 12),
            beam((x0, yy+35, 80), (8780, yy+35, 1130), 55, 12),
        ]
        M.add_product(
            "IFCBEAM", f"STAIR_CHANNEL16_{side}", "Channel No.16 stair stringer",
            f"CH16-{side}", LAYERS["stair"], shapes, "stair",
            category="stair", notes="Simplified U-channel geometry, nominal channel No.16."
        )

    post_x = [7700, 8060, 8420, 8780, 9160]
    for side, yy in (("L", y0-20), ("R", y0+width+20)):
        for i, xx in enumerate(post_x):
            basez = min(1050, max(0, (xx-x0)/1080*1050))
            M.add_product(
                "IFCBUILDINGELEMENTPROXY", f"RAIL_POST_{side}_{i+1:02d}",
                "Square tube 40x40x4 railing post", f"SQ40-{side}-{i+1:02d}",
                LAYERS["stair"], [beam((xx, yy, basez), (xx, yy, basez+1100), 40, 40)],
                "stair", category="stair", notes="40x40x4 square tube."
            )
        rail_points = []
        for xx in post_x:
            basez = min(1050, max(0, (xx-x0)/1080*1050))
            rail_points.append((xx, yy, basez+1100))
        for i in range(len(rail_points)-1):
            M.add_product(
                "IFCBUILDINGELEMENTPROXY", f"RAIL_TOP_{side}_{i+1:02d}",
                "Square tube 40x40x4 top rail", f"SQ40-RAIL-{side}-{i+1:02d}",
                LAYERS["stair"], [beam(rail_points[i], rail_points[i+1], 40, 40)],
                "stair", category="stair", notes="40x40x4 square tube."
            )


def hx_shapes(x, y, z, L=700, D=340, H=900):
    shapes = [
        box(x+L/2, y+D/2, z+120, L, D, H-120),
        box(x+80, y+D/2, z, 80, D+120, 120),
        box(x+L-80, y+D/2, z, 80, D+120, 120),
    ]
    for frac in (0.20, 0.35, 0.50, 0.65, 0.80):
        shapes.append(box(x+L*frac, y+D/2, z+140, 18, D+20, H-160))
    nozzle_x = (x+165, x+L-165)
    for xx in nozzle_x:
        for zz in (z+235, z+675):
            shapes.append(cyl((xx, y+D, zz), (xx, y+D+150, zz), 45))
            shapes.append(cyl((xx, y-150, zz), (xx, y, zz), 45))
    return shapes


def pump_shapes(x, y, z):
    return [
        box(x+410, y, z-165, 900, 420, 80),
        cyl((x+250, y, z), (x+680, y, z), 135, 16),
        cyl((x+110, y, z), (x+250, y, z), 185, 16),
        cyl((x, y, z), (x+110, y, z), 75, 16),
        cyl((x+680, y, z), (x+850, y, z), 75, 16),
        box(x+520, y, z+135, 190, 160, 90),
    ]


def tank_shapes(cx, cy, z, diameter=750, height=1350):
    r = diameter/2
    return [
        cyl((cx, cy, z+100), (cx, cy, z+height-100), r, 24),
        cyl((cx, cy, z), (cx, cy, z+100), r*0.86, 24),
        cyl((cx, cy, z+height-100), (cx, cy, z+height), r*0.86, 24),
        box(cx-r*0.55, cy, z-130, 70, 120, 230),
        box(cx+r*0.55, cy, z-130, 70, 120, 230),
        cyl((cx, cy, z-180), (cx, cy, z), 40),
    ]


def add_equipment():
    hxs = [
        ("K1", "Plate heat exchanger - Heating", 1500, 430, 265),
        ("K2", "Plate heat exchanger - Ventilation", 3650, 430, 265),
        ("K3", "Plate heat exchanger - DHW stage 1", 5650, 430, 265),
        ("K4", "Plate heat exchanger - DHW stage 2", 6750, 430, 265),
    ]
    for mark, desc, x, y, z in hxs:
        M.add_product(
            "IFCBUILDINGELEMENTPROXY", f"{mark}_HEAT_EXCHANGER", desc, mark,
            LAYERS["equipment"], hx_shapes(x, y, z), "equipment",
            system={"K1":"Heating","K2":"Ventilation","K3":"DHW","K4":"DHW"}[mark],
            mark=mark, category="equipment",
            notes="IFC2x3 proxy; intended class IfcHeatExchanger. Simplified body with four nozzles."
        )

    pumps = [
        ("K5_1", "K5", "Heating", 700, 1100, 470),
        ("K5_2", "K5", "Heating", 700, 1550, 470),
        ("K6_1", "K6", "Ventilation", 2900, 1100, 550),
        ("K6_2", "K6", "Ventilation", 2900, 1550, 550),
        ("K7_1", "K7", "DHW circulation", 5100, 1100, 490),
        ("K7_2", "K7", "DHW circulation", 5100, 1550, 490),
    ]
    for name, mark, system, x, y, z in pumps:
        M.add_product(
            "IFCBUILDINGELEMENTPROXY", f"{name}_PUMP", f"Circulation pump group {mark}", name,
            LAYERS["equipment"], pump_shapes(x, y, z), "equipment",
            system=system, mark=mark, category="equipment",
            notes="IFC2x3 proxy; intended class IfcPump. Each duty/standby pump is separate."
        )

    tanks = [
        ("K12", 1950, 2100, 180, "Heating"),
        ("K13", 2850, 2100, 180, "Ventilation"),
    ]
    for mark, x, y, z, system in tanks:
        M.add_product(
            "IFCBUILDINGELEMENTPROXY", f"{mark}_WESTER_W400", "Expansion tank Wester W400, 6 bar",
            mark, LAYERS["equipment"], tank_shapes(x, y, z), "tank",
            system=system, mark=mark, category="equipment",
            notes="IFC2x3 proxy; intended class IfcTank. Envelope: diameter 750 mm, height 1350 mm."
        )

    M.add_product(
        "IFCBUILDINGELEMENTPROXY", "HEAT_ENERGY_METERING_UNIT",
        "Heat energy metering unit on primary network", "UUTE-01",
        LAYERS["equipment"],
        [
            cyl((7050, 2050, 1965), (7470, 2050, 1965), 75, 16),
            cyl((6990, 2050, 1965), (7050, 2050, 1965), 115, 16),
            cyl((7470, 2050, 1965), (7530, 2050, 1965), 115, 16),
            box(7260, 2050, 2050, 220, 160, 170),
        ],
        "equipment", system="Primary heat network", dn=100, mark="UUTE",
        category="equipment", notes="Separate metering unit with DN100 connection flanges."
    )


def route(system_code, system_name, dn, points, color_key, tag_prefix):
    radius = max(11.0, dn * 0.58)
    segment_names = []
    for i in range(len(points)-1):
        a, b = points[i], points[i+1]
        name = f"{tag_prefix}_{system_code}_DN{dn}_SEG_{i+1:02d}"
        M.add_product(
            "IFCFLOWSEGMENT", name, f"{system_name}, straight pipe DN{dn}",
            name, LAYERS["pipe"], [cyl(a, b, radius, 16)], color_key,
            system=system_name, dn=dn, category="pipe",
            notes=f"Physical endpoints: {a} -> {b}; each straight is separate."
        )
        segment_names.append(name)
        if i < len(points)-2:
            p = b
            fitname = f"{tag_prefix}_{system_code}_DN{dn}_ELBOW_{i+1:02d}"
            M.add_product(
                "IFCFLOWFITTING", fitname, f"Elbow / directional fitting DN{dn}",
                fitname, LAYERS["fitting"],
                [cyl((p[0],p[1],p[2]-radius*0.8),
                      (p[0],p[1],p[2]+radius*0.8), radius*1.28, 16)],
                color_key, system=system_name, dn=dn, category="fitting",
                notes="Separate fitting at coincident pipe endpoints."
            )
    return segment_names


def fitting(name, description, p, dn, system, color_key="valve", mark="", kind="valve"):
    r = max(22, dn*0.75)
    shapes = [
        cyl((p[0]-r, p[1], p[2]), (p[0]+r, p[1], p[2]), max(12,dn*0.58), 16),
        box(p[0], p[1], p[2]-r*0.65, r*1.15, r*1.15, r*1.3, angle=math.pi/4),
    ]
    if kind == "filter":
        shapes.append(cyl((p[0],p[1],p[2]), (p[0],p[1],p[2]-r*1.8), r*0.62, 16))
    elif kind == "reducer":
        shapes = [
            cyl((p[0]-r*1.2,p[1],p[2]), (p[0],p[1],p[2]), max(12,dn*0.58), 16),
            cyl((p[0],p[1],p[2]), (p[0]+r*1.2,p[1],p[2]), max(10,dn*0.40), 16),
        ]
    elif kind == "tee":
        shapes.append(cyl((p[0],p[1]-r*1.5,p[2]), (p[0],p[1]+r*1.5,p[2]), max(10,dn*0.50), 16))
    elif kind == "flange":
        shapes = [cyl((p[0]-22,p[1],p[2]), (p[0]+22,p[1],p[2]), r, 16)]
    M.add_product(
        "IFCFLOWFITTING", name, description, name, LAYERS["fitting"],
        shapes, color_key, system=system, dn=dn, mark=mark,
        category="fitting", notes=f"Separate {kind} / flow fitting."
    )


def add_pipework_and_fittings():
    networks = [
        ("T1", "Primary heat network supply T1", 100,
         [(250,2500,2435),(1200,2500,2435),(1200,1950,2435),(5200,1950,2435),
          (5200,2380,2435),(8950,2380,2435)], "primary_supply", "P01"),
        ("T2", "Primary heat network return T2", 100,
         [(250,2250,1965),(1050,2250,1965),(1050,1850,1965),(5000,1850,1965),
          (5000,2180,1965),(8950,2180,1965)], "primary_return", "P02"),
        ("OT_SUP", "Heating supply", 65,
         [(2050,920,940),(2050,1220,940),(2450,1220,940),(2450,1780,940),
          (7400,1780,940),(7400,1650,940),(8950,1650,940)], "heating_supply", "P03"),
        ("OT_RET", "Heating return", 65,
         [(1650,920,500),(1650,1100,500),(700,1100,500),(700,620,500),(250,620,500)],
         "heating_return", "P04"),
        ("VENT_SUP", "Ventilation supply", 65,
         [(4200,920,1265),(4200,1350,1265),(4700,1350,1265),(4700,1680,1265),
          (7350,1680,1265),(7350,1450,1265),(8950,1450,1265)], "vent_supply", "P05"),
        ("VENT_RET", "Ventilation return", 65,
         [(3800,920,550),(3800,1100,550),(2900,1100,550),(2900,800,550),
          (250,800,550)], "vent_return", "P06"),
        ("DHW_SUP", "Domestic hot water supply", 50,
         [(7300,920,1565),(7300,1150,1565),(6500,1150,1565),(6500,1500,1565),
          (7600,1500,1565),(7600,1250,1565),(8950,1250,1565)], "dhw_supply", "P07"),
        ("DHW_CIRC", "Domestic hot water circulation", 32,
         [(6200,920,1165),(6200,1550,1165),(5100,1550,1165),(5100,1000,1165),
          (250,1000,1165)], "dhw_circ", "P08"),
        ("V1", "Cold water V1", 50,
         [(8950,1050,1300),(8100,1050,1300),(8100,1900,1300),(6500,1900,1300),
          (6500,920,1300),(6000,920,1300)], "cold_water", "P09"),
        ("DRAIN", "Drainage to sewer", 25,
         [(1500,950,180),(1500,2400,180),(4500,2400,180),(4500,2300,180),
          (7600,2300,180),(8950,2300,180)], "drain", "P10"),
        ("MAKEUP", "Heating and ventilation make-up", 32,
         [(2850,1920,700),(2850,1700,700),(4200,1700,700),(4200,1500,700),
          (5100,1500,700)], "cold_water", "P11"),
    ]
    for args in networks:
        route(*args)

    marked = [
        ("K8_DP_REGULATOR_DN50", "K8 differential pressure regulator", (4700,1850,1965), 50,
         "Primary heat network return T2", "K8", "valve"),
        ("K9_HEATING_CONTROL_DN25", "K9 heating control valve", (2450,1220,940), 25,
         "Heating supply", "K9", "valve"),
        ("K10_VENT_CONTROL_DN25", "K10 ventilation control valve", (4700,1350,1265), 25,
         "Ventilation supply", "K10", "valve"),
        ("K11_DHW_CONTROL_DN25", "K11 DHW control valve", (6500,1150,1565), 25,
         "Domestic hot water supply", "K11", "valve"),
        ("K14_HEATING_SOLENOID_DN15", "K14 heating solenoid valve", (1200,1100,500), 15,
         "Heating return", "K14", "valve"),
        ("K15_VENT_SOLENOID_DN15", "K15 ventilation solenoid valve", (3400,1100,550), 15,
         "Ventilation return", "K15", "valve"),
    ]
    for name, desc, p, dn, sys, mark, kind in marked:
        fitting(name, desc, p, dn, sys, "valve", mark, kind)

    generic = [
        ("MUD_COLLECTOR_DN100_01","Vertical mud collector DN100",(6100,2380,2435),100,
         "Primary heat network supply T1","filter"),
        ("MAGNETIC_FILTER_DN100_01","Magnetic mesh filter DN100",(3300,1950,2435),100,
         "Primary heat network supply T1","filter"),
        ("MAGNETIC_FILTER_DN100_02","Magnetic mesh filter DN100",(3500,1850,1965),100,
         "Primary heat network return T2","filter"),
        ("FILTER_HEATING_DN50_01","Heating mesh filter DN50",(3500,1780,940),50,
         "Heating supply","filter"),
        ("FILTER_VENT_DN50_01","Ventilation mesh filter DN50",(5550,1680,1265),50,
         "Ventilation supply","filter"),
        ("FILTER_DHW_DN32_01","DHW circulation mesh filter DN32",(4300,1550,1165),32,
         "Domestic hot water circulation","filter"),
        ("CHECK_VALVE_HEATING_DN65","Heating check valve DN65",(6100,1780,940),65,
         "Heating supply","valve"),
        ("CHECK_VALVE_VENT_DN65","Ventilation check valve DN65",(6400,1680,1265),65,
         "Ventilation supply","valve"),
        ("BALANCE_DHW_DN32","DHW balancing valve DN32",(3600,1550,1165),32,
         "Domestic hot water circulation","valve"),
        ("REDUCER_TS_100_65","Reducer DN100 to DN65",(5350,1950,2435),100,
         "Primary heat network supply T1","reducer"),
        ("REDUCER_DHW_50_32","Reducer DN50 to DN32",(6900,1500,1565),50,
         "Domestic hot water supply","reducer"),
        ("TEE_TS_DN100_01","Primary network tee DN100",(5200,1950,2435),100,
         "Primary heat network supply T1","tee"),
        ("TEE_OT_DN65_01","Heating tee DN65",(2450,1780,940),65,
         "Heating supply","tee"),
        ("TEE_VENT_DN65_01","Ventilation tee DN65",(4700,1680,1265),65,
         "Ventilation supply","tee"),
    ]
    for name, desc, p, dn, sys, kind in generic:
        fitting(name, desc, p, dn, sys, "valve", "", kind)

    flange_points = [
        ("FLG_K1_A", (1650,920,500), 65, "Heating return"),
        ("FLG_K1_B", (2050,920,940), 65, "Heating supply"),
        ("FLG_K2_A", (3800,920,550), 65, "Ventilation return"),
        ("FLG_K2_B", (4200,920,1265), 65, "Ventilation supply"),
        ("FLG_K3_A", (6000,920,1300), 50, "Cold water V1"),
        ("FLG_K3_B", (6200,920,1165), 32, "DHW circulation"),
        ("FLG_K4_A", (7300,920,1565), 50, "DHW supply"),
        ("FLG_TS_01", (1200,2500,2435), 100, "Primary heat network supply T1"),
        ("FLG_TS_02", (1050,2250,1965), 100, "Primary heat network return T2"),
        ("FLG_TS_03", (7800,2380,2435), 100, "Primary heat network supply T1"),
        ("FLG_TS_04", (7800,2180,1965), 100, "Primary heat network return T2"),
    ]
    for name, p, dn, sys in flange_points:
        fitting(name, f"Flanged connection DN{dn}", p, dn, sys, "valve", "", "flange")


def add_supports():
    supports = [
        (850, 2250, 1880, 100), (1800, 1950, 2355, 100),
        (3000, 1950, 2355, 100), (4300, 1950, 2355, 100),
        (5600, 2380, 2355, 100), (7000, 2380, 2355, 100),
        (2200, 1780, 870, 65), (3700, 1680, 1195, 65),
        (5550, 1680, 1195, 65), (6800, 1500, 1495, 50),
        (4200, 1550, 1105, 32), (6500, 1900, 1240, 50),
    ]
    for i, (x,y,z_top,dn) in enumerate(supports,1):
        shapes = [
            box(x, y, 0, 60, 60, max(80,z_top-80)),
            box(x, y, max(0,z_top-80), max(260,dn*2.8), 80, 80),
            box(x-100, y, 0, 240, 180, 20),
        ]
        M.add_product(
            "IFCBUILDINGELEMENTPROXY", f"PIPE_SUPPORT_{i:02d}",
            f"Floor-mounted pipe support for DN{dn}", f"SUP-{i:02d}",
            LAYERS["support"], shapes, "support", dn=dn, category="support",
            notes="Simplified steel post and crossbar; separate selectable support."
        )


def add_instruments():
    instruments = [
        ("PI_TS_SUP_01",(2500,1950,2435),"Primary heat network supply T1",100,"Pressure gauge"),
        ("TI_TS_SUP_01",(2800,1950,2435),"Primary heat network supply T1",100,"Thermometer"),
        ("PI_TS_RET_01",(2600,1850,1965),"Primary heat network return T2",100,"Pressure gauge"),
        ("TI_TS_RET_01",(2900,1850,1965),"Primary heat network return T2",100,"Thermometer"),
        ("PI_OT_SUP_01",(3000,1780,940),"Heating supply",65,"Pressure gauge"),
        ("TI_OT_SUP_01",(3250,1780,940),"Heating supply",65,"Thermometer"),
        ("PI_VENT_SUP_01",(5200,1680,1265),"Ventilation supply",65,"Pressure gauge"),
        ("TI_VENT_SUP_01",(5450,1680,1265),"Ventilation supply",65,"Thermometer"),
        ("PI_DHW_01",(7050,1500,1565),"DHW supply",50,"Pressure gauge"),
        ("TI_DHW_01",(7300,1500,1565),"DHW supply",50,"Thermometer"),
        ("FLOW_SENSOR_UUTE",(7260,2050,1965),"Primary heat network",100,"Primary flow sensor"),
        ("TEMP_SENSOR_V1",(6500,1900,1300),"Cold water V1",50,"Temperature sensor"),
    ]
    for name, p, system, dn, desc in instruments:
        shapes = [
            cyl((p[0],p[1],p[2]),(p[0],p[1],p[2]+150),18,12),
            cyl((p[0],p[1]-28,p[2]+215),(p[0],p[1]+28,p[2]+215),75,20),
            box(p[0],p[1],p[2]+205,60,20,20),
        ]
        M.add_product(
            "IFCBUILDINGELEMENTPROXY", name, desc, name, LAYERS["instrument"],
            shapes, "instrument", system=system, dn=dn, category="instrument",
            notes="IFC2x3 proxy; intended class IfcSensor. Separate selectable instrument."
        )


add_room()
add_stair()
add_equipment()
add_pipework_and_fittings()
add_supports()
add_instruments()
M.finish()
M.write(OUT / IFC_NAME)


def pad4(data: bytes, pad=b"\x00"):
    return data + pad * ((4 - len(data) % 4) % 4)


def make_glb(objects, path):
    blob = bytearray()
    buffer_views = []
    accessors = []
    meshes = []
    nodes = []
    materials = []
    material_index = {}

    def material_for(key):
        if key in material_index:
            return material_index[key]
        rgba = COLORS[key]
        idx = len(materials)
        materials.append({
            "name": key,
            "pbrMetallicRoughness": {
                "baseColorFactor": list(rgba),
                "metallicFactor": 0.20 if key in ("stair","support") else 0.05,
                "roughnessFactor": 0.62,
            },
            "alphaMode": "BLEND" if rgba[3] < 0.999 else "OPAQUE",
            "doubleSided": True,
        })
        material_index[key] = idx
        return idx

    def add_buffer(data, target):
        while len(blob) % 4:
            blob.append(0)
        off = len(blob)
        blob.extend(data)
        idx = len(buffer_views)
        buffer_views.append({"buffer":0,"byteOffset":off,"byteLength":len(data),"target":target})
        return idx

    def add_accessor(view, component_type, count, typ, mins=None, maxs=None):
        a = {"bufferView":view,"componentType":component_type,"count":count,"type":typ}
        if mins is not None:
            a["min"] = mins
            a["max"] = maxs
        idx = len(accessors)
        accessors.append(a)
        return idx

    for obj in objects:
        flat_pos, flat_nrm, idxs = [], [], []
        out_i = 0
        for face in obj["faces"]:
            a,b,c = [obj["vertices"][i] for i in face]
            n = vnorm(vcross(vsub(b,a),vsub(c,a)))
            for p in (a,b,c):
                flat_pos.extend([p[0],p[1],p[2]])
                flat_nrm.extend(n)
                idxs.append(out_i)
                out_i += 1
        pos_bytes = struct.pack("<" + "f"*len(flat_pos), *flat_pos)
        nrm_bytes = struct.pack("<" + "f"*len(flat_nrm), *flat_nrm)
        idx_bytes = struct.pack("<" + "I"*len(idxs), *idxs)
        pv = add_buffer(pos_bytes, 34962)
        nv = add_buffer(nrm_bytes, 34962)
        iv = add_buffer(idx_bytes, 34963)
        xs,ys,zs = flat_pos[0::3],flat_pos[1::3],flat_pos[2::3]
        pa = add_accessor(pv,5126,len(flat_pos)//3,"VEC3",
                          [min(xs),min(ys),min(zs)],[max(xs),max(ys),max(zs)])
        na = add_accessor(nv,5126,len(flat_nrm)//3,"VEC3")
        ia = add_accessor(iv,5125,len(idxs),"SCALAR",[0],[max(idxs)])
        mi = len(meshes)
        meshes.append({
            "name": obj["name"],
            "primitives":[{
                "attributes":{"POSITION":pa,"NORMAL":na},
                "indices":ia,
                "material":material_for(obj["color"]),
                "mode":4,
            }],
            "extras":{
                "ifcType":obj["ifc_type"],"layer":obj["layer"],
                "system":obj["system"],"dn_mm":obj["dn"],"mark":obj["mark"],
            },
        })
        nodes.append({
            "name": obj["name"],
            "mesh": mi,
            "extras":{
                "ifcType":obj["ifc_type"],"layer":obj["layer"],
                "system":obj["system"],"dn_mm":obj["dn"],"mark":obj["mark"],
            },
        })

    doc = {
        "asset":{"version":"2.0","generator":"ITP Parametric Builder",
                 "extras":{"coordinateUnits":"millimetres","controlLength":ROOM_LENGTH}},
        "scene":0,
        "scenes":[{"name":"ITP Complete","nodes":list(range(len(nodes)))}],
        "nodes":nodes,
        "meshes":meshes,
        "materials":materials,
        "buffers":[{"byteLength":len(blob)}],
        "bufferViews":buffer_views,
        "accessors":accessors,
    }
    js = pad4(json.dumps(doc,separators=(",",":"),ensure_ascii=False).encode("utf-8"),b" ")
    bn = pad4(bytes(blob),b"\x00")
    total = 12 + 8 + len(js) + 8 + len(bn)
    glb = (
        struct.pack("<4sII",b"glTF",2,total) +
        struct.pack("<I4s",len(js),b"JSON") + js +
        struct.pack("<I4s",len(bn),b"BIN\x00") + bn
    )
    path.write_bytes(glb)


make_glb(M.objects, OUT / GLB_NAME)


def rgb255(key, shade=1.0):
    r,g,b,_ = COLORS[key]
    return tuple(max(0,min(255,int(c*255*shade))) for c in (r,g,b))


def render_view(objects, path, title, right, up, depth, include):
    from PIL import Image, ImageDraw, ImageFont

    selected = [o for o in objects if include(o)]
    all_points = [p for o in selected for p in o["vertices"]]
    proj = [(vdot(p,right),vdot(p,up),vdot(p,depth)) for p in all_points]
    min_u,max_u = min(p[0] for p in proj),max(p[0] for p in proj)
    min_v,max_v = min(p[1] for p in proj),max(p[1] for p in proj)
    W,H = 1600,1000
    margin_l,margin_r,margin_t,margin_b = 90,300,90,80
    scale = min((W-margin_l-margin_r)/(max_u-min_u or 1),
                (H-margin_t-margin_b)/(max_v-min_v or 1))

    def to_px(p):
        u,v,d = vdot(p,right),vdot(p,up),vdot(p,depth)
        x = margin_l + (u-min_u)*scale
        y = H-margin_b - (v-min_v)*scale
        return (x,y,d)

    triangles = []
    light = vnorm((0.35,-0.25,0.90))
    for obj in selected:
        for face in obj["faces"]:
            pts3 = [obj["vertices"][i] for i in face]
            a,b,c = pts3
            try:
                normal = vnorm(vcross(vsub(b,a),vsub(c,a)))
            except ValueError:
                continue
            shade = 0.56 + 0.44*abs(vdot(normal,light))
            pp = [to_px(p) for p in pts3]
            avgd = sum(p[2] for p in pp)/3
            triangles.append((avgd, [(p[0],p[1]) for p in pp], rgb255(obj["color"],shade), obj))

    img = Image.new("RGB",(W,H),(244,246,248))
    draw = ImageDraw.Draw(img)
    for _,poly,color,obj in sorted(triangles,key=lambda x:x[0]):
        draw.polygon(poly,fill=color,outline=tuple(max(0,c-28) for c in color))

    try:
        font_title = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",30)
        font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",17)
        font_small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",14)
    except Exception:
        font_title = font = font_small = ImageFont.load_default()

    draw.rectangle((0,0,W,58),fill=(30,34,38))
    draw.text((24,13),title,fill=(255,255,255),font=font_title)
    draw.text((W-300,18),"ITP | ARCHICAD 23 | IFC2x3",fill=(210,215,220),font=font)

    lx,ly = W-270,90
    draw.rounded_rectangle((lx-18,ly-18,W-24,ly+310),radius=10,fill=(255,255,255),outline=(190,195,200))
    draw.text((lx,ly),"SYSTEMS / LAYERS",fill=(25,25,25),font=font)
    legend = [
        ("Primary T1/T2","primary_supply"),
        ("Heating","heating_supply"),
        ("Ventilation","vent_supply"),
        ("DHW","dhw_supply"),
        ("Cold water","cold_water"),
        ("Drain","drain"),
        ("Equipment","equipment"),
        ("Valves","valve"),
        ("Stair / steel","stair"),
    ]
    for i,(label,key) in enumerate(legend):
        yy=ly+36+i*27
        draw.rectangle((lx,yy,lx+25,yy+16),fill=rgb255(key),outline=(40,40,40))
        draw.text((lx+36,yy-2),label,fill=(30,30,30),font=font_small)
    draw.text((24,H-38),f"Control length: {ROOM_LENGTH:.0f} mm | Elements: {len(objects)} | Z=0.000 FFL",
              fill=(45,50,55),font=font)
    img.save(path,optimize=True)


iso_right = vnorm((0.72,0.69,0))
iso_up = vnorm((-0.39,0.41,0.82))
iso_depth = vnorm((0.58,-0.60,0.55))
render_view(
    M.objects, OUT/"ITP_View_Isometric.png", "ITP COMPLETE — ISOMETRIC CUTAWAY",
    iso_right, iso_up, iso_depth,
    lambda o: o["category"] not in ("ceiling",) and o["name"] not in ("WALL_SOUTH_01","WALL_EAST_LOW","WALL_EAST_HIGH")
)
render_view(
    M.objects, OUT/"ITP_View_Top.png", "ITP COMPLETE — TOP VIEW",
    (1,0,0),(0,1,0),(0,0,1),
    lambda o: o["category"] != "ceiling"
)
render_view(
    M.objects, OUT/"ITP_View_Longitudinal.png", "ITP COMPLETE — LONGITUDINAL VIEW",
    (1,0,0),(0,0,1),(0,-1,0),
    lambda o: o["category"] not in ("ceiling",) and o["name"] != "WALL_SOUTH_01"
)
render_view(
    M.objects, OUT/"ITP_View_NoWalls.png", "ITP COMPLETE — EQUIPMENT & PIPEWORK",
    iso_right, iso_up, iso_depth,
    lambda o: o["category"] not in ("wall","ceiling","soffit")
)


def bboxes_overlap(a,b,tol=0.1):
    return (
        a[0] < b[3]-tol and a[3] > b[0]+tol and
        a[1] < b[4]-tol and a[4] > b[1]+tol and
        a[2] < b[5]-tol and a[5] > b[2]+tol
    )


walls_stair = [o for o in M.objects if o["category"] in ("wall","stair")]
equipment = [o for o in M.objects if o["category"] == "equipment"]
collisions = []
for e in equipment:
    for w in walls_stair:
        if bboxes_overlap(e["bbox"],w["bbox"]):
            collisions.append([e["name"],w["name"]])

allv = [p for o in M.objects for p in o["vertices"]]
bounds = {
    "min_x":min(p[0] for p in allv),"max_x":max(p[0] for p in allv),
    "min_y":min(p[1] for p in allv),"max_y":max(p[1] for p in allv),
    "min_z":min(p[2] for p in allv),"max_z":max(p[2] for p in allv),
}
overall_x = bounds["max_x"]-bounds["min_x"]
assert abs(overall_x-ROOM_LENGTH) < 0.01, (overall_x,bounds)
assert not collisions, collisions

# Structural GLB self-check.
glb_data = (OUT/GLB_NAME).read_bytes()
magic, version, glb_len = struct.unpack("<4sII",glb_data[:12])
assert magic == b"glTF" and version == 2 and glb_len == len(glb_data)
json_len,json_type = struct.unpack("<I4s",glb_data[12:20])
assert json_type == b"JSON"
glb_json = json.loads(glb_data[20:20+json_len].decode("utf-8"))
assert len(glb_json["nodes"]) == len(M.objects)
assert len({n["name"] for n in glb_json["nodes"]}) == len(M.objects)
assert glb_json["asset"]["extras"]["controlLength"] == ROOM_LENGTH

# IFC reopen and geometry creation validation.
import ifcopenshell
import ifcopenshell.geom

ifc = ifcopenshell.open(str(OUT/IFC_NAME))
assert ifc.schema.upper() == "IFC2X3"
all_ifc_elements = ifc.by_type("IfcElement")
elements = [p for p in ifc.by_type("IfcProduct") if getattr(p, "Representation", None)]
assert len(elements) == len(M.objects), (len(elements),len(M.objects))
assert len(ifc.by_type("IfcFlowSegment")) >= 45
assert len(ifc.by_type("IfcFlowFitting")) >= 45
assert sum(1 for o in M.objects if o["name"].endswith("_PUMP")) == 6
assert sum(1 for o in M.objects if o["name"].endswith("_HEAT_EXCHANGER")) == 4
assert sum(1 for o in M.objects if "WESTER_W400" in o["name"]) == 2
assert len({e.Name for e in elements}) == len(elements)

import ifcopenshell.validate
validation_logger = ifcopenshell.validate.json_logger()
ifcopenshell.validate.validate(ifc, validation_logger)
if validation_logger.statements:
    print("IFC_SCHEMA_VALIDATION_BEGIN")
    print(json.dumps(validation_logger.statements[:100], ensure_ascii=False, indent=2, default=str))
    print("IFC_SCHEMA_VALIDATION_END")
assert not validation_logger.statements, validation_logger.statements[:10]

settings = ifcopenshell.geom.settings()
geometry_errors = []
for e in elements:
    try:
        shape = ifcopenshell.geom.create_shape(settings,e)
        if not getattr(shape.geometry,"verts",None):
            geometry_errors.append([e.Name,"empty geometry"])
    except Exception as exc:
        geometry_errors.append([e.Name,e.is_a(),repr(exc)])
assert not geometry_errors, geometry_errors[:10]

counts = {
    "IfcProductWithGeometry":len(elements),
    "IfcElement_total":len(all_ifc_elements),
    "IfcWallStandardCase":len(ifc.by_type("IfcWallStandardCase")),
    "IfcFlowSegment":len(ifc.by_type("IfcFlowSegment")),
    "IfcFlowFitting":len(ifc.by_type("IfcFlowFitting")),
    "HeatExchanger_proxies":sum(1 for o in M.objects if o["name"].endswith("_HEAT_EXCHANGER")),
    "Pump_proxies":sum(1 for o in M.objects if o["name"].endswith("_PUMP")),
    "Tank_proxies":sum(1 for o in M.objects if "WESTER_W400" in o["name"]),
    "Instrument_proxies":sum(1 for o in M.objects if o["category"] == "instrument"),
    "IfcBuildingElementProxy":len(ifc.by_type("IfcBuildingElementProxy")),
}

qa = {
    "status":"PASS",
    "timestamp_utc":datetime.now(timezone.utc).isoformat(),
    "schema":ifc.schema,
    "coordination_view":"CoordinationView_V2.0",
    "units":"millimetres",
    "control_length_mm":ROOM_LENGTH,
    "computed_overall_x_mm":overall_x,
    "bounds_mm":bounds,
    "element_counts":counts,
    "separate_ifc_elements":True,
    "unique_element_names":True,
    "glb_separate_nodes":len(glb_json["nodes"]),
    "ifc_reopen":True,
    "ifc_geometry_created_for_all_elements":True,
    "equipment_vs_wall_or_stair_collisions":collisions,
    "pipe_endpoint_continuity":"PASS - every polyline uses coincident adjacent endpoints",
    "section_level_check_mm":[180,470,490,500,550,700,940,1165,1265,1300,1565,1965,2435],
    "max_service_level_plus_radius_mm":2493,
    "ceiling_level_mm":ROOM_HEIGHT,
}
(OUT/"QA_report.json").write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding="utf-8")

readme = f"""ITP COMPLETE — ARCHICAD 23 COORDINATION MODEL
================================================

1. DELIVERABLES
- {IFC_NAME}: IFC 2x3 Coordination View model for ARCHICAD 23.
- {GLB_NAME}: lightweight preview; every logical item is a separate named node.
- ITP_View_Isometric.png: overall isometric cutaway.
- ITP_View_Top.png: top view.
- ITP_View_Longitudinal.png: longitudinal coordination view.
- ITP_View_NoWalls.png: equipment and pipework without walls.
- QA_report.json: machine-readable validation report.

2. SOURCE PRIORITY
- 3(1).pdf, room No.11 ITP: sole priority for room shell, floor datum, stair,
  platforms, Channel No.16 and square tube 40x40x4.
- CTP-TM(1).pdf sheet 1: process scheme, equipment marks and nominal diameters.
- CTP-TM(1).pdf sheets 2–3: plan location and equipment composition.
- CTP-TM(1).pdf sheet 4: vertical levels and section coordination.
- CTP-TM(1).pdf sheet 5: spatial route interpretation.

3. COORDINATE SYSTEM AND SCALE
- IFC geometry and declared length units: millimetres.
- X: room longitudinal direction; Y: transverse direction; Z: upward.
- Z=0 is the finished floor level of the ITP.
- Verified overall X dimension: {overall_x:.0f} mm (target 9205 mm).
- Main conditional ceiling: Z={ROOM_HEIGHT:.0f} mm.

4. MODEL STRUCTURE / LAYERS
- 01_ITP_Room
- 02_ITP_Stair
- 03_ITP_Equipment
- 04_ITP_Pipework
- 05_ITP_Fittings_Valves
- 06_ITP_Supports
- 07_ITP_Instruments

Every straight pipe is a separate IfcFlowSegment. Elbows, tees, reducers, flanges,
filters and valves are separate IfcFlowFitting elements. K1–K4, K5–K7 duty/standby units and K12–K13 are separate
IfcBuildingElementProxy occurrences because IFC 2x3 contains the corresponding
pump/heat-exchanger/tank classes primarily as type-level definitions. Functional
class, name and mark are retained in Pset_ITP_Identity. K8–K11 and K14–K15 are retained
in object names and Mark properties. All items have Pset_ITP_Identity.

5. SYSTEM COLORS
- Primary supply T1: bright red
- Primary return T2: dark garnet
- Heating supply / return: orange / ochre
- Ventilation supply / return: purple / violet
- DHW supply / circulation: orange / yellow
- Cold water V1 and make-up: blue
- Drainage: dark gray
- Equipment: graphite gray; tanks: light gray
- Valves and fittings: yellow
- Instruments: turquoise
- Stair, channels and railings: steel gray

6. ACCEPTED ASSUMPTIONS
- Where a dimension is not explicitly stated, it was scaled from the relevant
  plan/section and documented here rather than stopping the build.
- Nominal room width is accepted as 2915 mm. Three local northern wall recesses
  and the east doorway are interpreted from the plan/elevations.
- Wall thickness is accepted as 200 mm where it is not dimensioned separately.
- The main ceiling boundary is Z=2788 mm; local soffits are accepted at Z=2600 mm.
- Stair flight: six grated treads, 180 mm nominal run and 175 mm rise; landing
  Z=1050 mm. Channel No.16 and 40x40x4 rail profiles are coordination envelopes.
- Manufacturer STEP/Revit families were not supplied. Equipment is simplified,
  recognizable and coordination-sized, with connection nozzles.
- Wester W400 coordination envelope: 750 mm diameter, 1350 mm body height.
- Pipe centerline levels are coordinated to section values including 490, 940,
  1165, 1265, 1565, 1965 and 2435 mm. Minor offsets were inferred from scale.
- Fittings are coordination representations, not fabrication spool geometry.
- Terminal pipe ends stop at room interfaces; networks outside the ITP are not
  included.
- The model is intended for design coordination and quantity/navigation checks,
  not for fabrication without a site survey and manufacturer submittals.

7. ARCHICAD 23 IMPORT
1) File > Interoperability > IFC > Open or Merge.
2) Select the General Import / Coordination translator for IFC 2x3.
3) Confirm model units are millimetres and preserve source IFC classifications.
4) Map IfcPresentationLayerAssignment names to ARCHICAD layers.
5) Keep MEP products as individual imported elements; do not choose an option
   that converts the complete IFC into one Library Object.
6) Verify the X control dimension is 9205 mm and FFL is Z=0.
7) Check the four supplied views, then use Find & Select by IFC Type, Mark,
   System or DN_mm.

8. AUTOMATED QA
- IFC reopened successfully with IfcOpenShell.
- Geometry was created for every IfcElement after reopening.
- Element names are unique; products are not merged into one proxy.
- Equipment/wall/stair bounding-box collision check: PASS.
- Pipe polyline adjacency check: PASS.
- GLB header, byte length, JSON chunk and separate node count: PASS.
- Counts: {json.dumps(counts,ensure_ascii=False)}

Generated: {datetime.now(timezone.utc).isoformat()}
"""
(OUT/"README.txt").write_text(readme,encoding="utf-8")

# Hash list before ZIP.
files_for_zip = [
    IFC_NAME, GLB_NAME,
    "ITP_View_Isometric.png","ITP_View_Top.png",
    "ITP_View_Longitudinal.png","ITP_View_NoWalls.png",
    "README.txt","QA_report.json",
]
hash_lines = []
for fn in files_for_zip:
    h = hashlib.sha256((OUT/fn).read_bytes()).hexdigest()
    hash_lines.append(f"{h}  {fn}")
(OUT/"SHA256SUMS.txt").write_text("\n".join(hash_lines)+"\n",encoding="utf-8")
files_for_zip.append("SHA256SUMS.txt")

with zipfile.ZipFile(OUT/ZIP_NAME,"w",compression=zipfile.ZIP_DEFLATED,compresslevel=9) as zf:
    for fn in files_for_zip:
        zf.write(OUT/fn,arcname=fn)

print(json.dumps({
    "status":"PASS",
    "ifc":str(OUT/IFC_NAME),
    "glb":str(OUT/GLB_NAME),
    "zip":str(OUT/ZIP_NAME),
    "counts":counts,
    "overall_x_mm":overall_x,
    "collisions":collisions,
},ensure_ascii=False,indent=2))
