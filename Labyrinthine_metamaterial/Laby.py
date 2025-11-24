import math
import gmsh
import sys

# =============================================================================
# 1. PARAMÈTRES
# =============================================================================

def get_params_PMMA():
    a = 50.0
    cols = 1
    rows = 1
    array_width = cols * a
    array_height = rows * a

    return {
        "a": a,
        "half_a": 0.5 * a,
        "tb": 1.25,
        "R_start": 23.75,
        "num_void_layers": 9,
        "cols": cols,
        "rows": rows,
        "array_width": array_width,
        "array_height": array_height,
        "tfoot": 1.25,
        "doutY": 4.0,
        "Lplate": array_width + 4.0,
        "tplate": 6.0,
        "t": 8.0,
        "din1": 1.0,
        "din3": 9.0,
        "n_spec_1mm": 3,
        "n_spec_9mm": 2,
        "doutZ": 4.0,
        "eps": 0.15,
        "seq_card": [1, 0, 1, 0, 1, 0, 1, 0, 1],
        "seq_diag": [0, 1, 0, 1, 0, 1, 0, 1, 0],

        # Maillage régulier (taille cible)
        "mesh_h": a / 80.0,

        # Transfinite sur les bords périodiques x=0/a et y=0/a
        "n_tf_x": 3,       # nb de points sur [0, a] en x (bords y=0 et y=a)
        "n_tf_y": 3,       # nb de points sur [0, a] en y (bords x=0 et x=a)

        # Transfinite sur les poutres externes (Top_Beam / Bot_Beam)
        "n_tf_beam_x": 41,  # le long de la longueur de la poutre
        "n_tf_beam_y": 5,   # à travers l'épaisseur de la poutre
    }

def get_params_3D(specimen_label="b3D"):
    a = 50.0
    cols = 1
    rows = 1
    array_width = cols * a
    array_height = rows * a

    return {
        "a": a,
        "half_a": 0.5 * a,
        "tb": 1.25,
        "R_start": 23.75,
        "num_void_layers": 9,
        "cols": cols,
        "rows": rows,
        "array_width": array_width,
        "array_height": array_height,
        "tfoot": 0.0,
        "doutY": 15.5,
        "Lplate": array_width + 4.0,
        "tplate": 8.0,
        "t": 73.0,
        "doutZ": 16.7,
        "eps": 0.15,
        "seq_card": [1, 0, 1, 0, 1, 0, 1, 0, 1],
        "seq_diag": [0, 1, 0, 1, 0, 1, 0, 1, 0],

        # Maillage régulier
        "mesh_h": a / 40.0,

        # Transfinite bords périodiques
        "n_tf_x": 5,
        "n_tf_y": 5,

        # Transfinite poutres externes
        "n_tf_beam_x": 41,
        "n_tf_beam_y": 5,
    }

# =============================================================================
# 2. GÉOMÉTRIE 2D
# =============================================================================

def create_ring_surface(cx, cy, r_outer, thickness):
    r_inner = r_outer - thickness
    if r_inner < 0.001:
        c = gmsh.model.occ.addCircle(cx, cy, 0, r_outer)
        loop = gmsh.model.occ.addCurveLoop([c])
        return gmsh.model.occ.addPlaneSurface([loop])
    else:
        c1 = gmsh.model.occ.addCircle(cx, cy, 0, r_outer)
        c2 = gmsh.model.occ.addCircle(cx, cy, 0, r_inner)
        loop_out = gmsh.model.occ.addCurveLoop([c1])
        loop_in = gmsh.model.occ.addCurveLoop([c2])
        return gmsh.model.occ.addPlaneSurface([loop_out, loop_in])

def create_bridge_surface(cx, cy, r_outer, r_inner, width, angle_deg, eps):
    length = (r_outer - r_inner) + 2 * eps
    start_r = r_inner - eps
    br = gmsh.model.occ.addRectangle(start_r, -width / 2.0, 0, length, width)
    gmsh.model.occ.rotate([(2, br)], 0, 0, 0, 0, 0, 1, math.radians(angle_deg))
    gmsh.model.occ.translate([(2, br)], cx, cy, 0)
    return br

def create_anchor_surface(cx, cy, angle_deg, params):
    rad = math.radians(angle_deg)

    R = params["R_start"]
    tan_x = cx + R * math.cos(rad)
    tan_y = cy + R * math.sin(rad)

    anchor_angle = rad + math.pi / 2.0

    length = 40.0
    width = params["tb"]

    b = gmsh.model.occ.addRectangle(-length / 2.0, -width / 2.0, 0, length, width)
    gmsh.model.occ.rotate([(2, b)], 0, 0, 0, 0, 0, 1, anchor_angle)
    gmsh.model.occ.translate([(2, b)], tan_x, tan_y, 0)

    return b

def create_plate_surface(x, y, width, height):
    return gmsh.model.occ.addRectangle(x, y, 0, width, height)

# =============================================================================
# 3. COMPOSANTS 2D
# =============================================================================

def generate_unit_cell_parts(cx, cy, params):
    parts = []

    # Ancrages diagonaux
    for ang in [45, 135, 225, 315]:
        parts.append(create_anchor_surface(cx, cy, ang, params))

    # Anneaux + poutres
    current_R = params["R_start"]
    tb = params["tb"]
    for k in range(params["num_void_layers"]):
        parts.append(create_ring_surface(cx, cy, current_R, tb))
        current_R -= tb

        void_outer = current_R
        void_inner = current_R - tb

        if params["seq_card"][k]:
            for ang in [0, 90, 180, 270]:
                parts.append(
                    create_bridge_surface(
                        cx, cy, void_outer, void_inner, tb, ang, params["eps"]
                    )
                )
        if params["seq_diag"][k]:
            for ang in [45, 135, 225, 315]:
                parts.append(
                    create_bridge_surface(
                        cx, cy, void_outer, void_inner, tb, ang, params["eps"]
                    )
                )
        current_R -= tb

    if current_R > 0.001:
        center_r = current_R + params["eps"]
        parts.append(create_ring_surface(cx, cy, center_r, center_r))

    return parts

def generate_array_parts(params):
    all_parts = []
    print(f"Generating array: {params['cols']} x {params['rows']} unit cells...")
    for i in range(params["cols"]):
        for j in range(params["rows"]):
            cx = i * params["a"] + params["half_a"]
            cy = j * params["a"] + params["half_a"]
            all_parts.extend(generate_unit_cell_parts(cx, cy, params))
    return all_parts

# =============================================================================
# 4. ASSEMBLAGE 2D
# =============================================================================

def build_fused_specimen_2d(array_parts, params):
    print("Building Body 1 (specimen + internal feet) in 2D (Safe Method)...")
    parts_to_fuse = [(2, p) for p in array_parts]

    tfoot = params["tfoot"]
    eps = params["eps"]

    # Pieds internes
    if tfoot > 0.0:
        foot_bot = create_plate_surface(0.0, -tfoot, params["array_width"], tfoot + eps)
        parts_to_fuse.append((2, foot_bot))

        foot_top_y = params["array_height"] - eps
        foot_top = create_plate_surface(
            0.0, foot_top_y, params["array_width"], tfoot + eps
        )
        parts_to_fuse.append((2, foot_top))

        y_min = -tfoot
        total_h = params["array_height"] + 2.0 * tfoot
    else:
        y_min = 0.0
        total_h = params["array_height"]

    gmsh.model.occ.synchronize()

    # 1) Fragmentation
    fragmented_shapes, _ = gmsh.model.occ.fragment(parts_to_fuse, [])

    # 2) Fusion
    if len(fragmented_shapes) > 1:
        fused_body_list, _ = gmsh.model.occ.fuse(
            [fragmented_shapes[0]], fragmented_shapes[1:]
        )
    else:
        fused_body_list = fragmented_shapes

    # 3) Recoupe propre avec une boîte
    bbox = gmsh.model.occ.addRectangle(0.0, y_min, 0.0, params["array_width"], total_h)
    final_body, _ = gmsh.model.occ.intersect(fused_body_list, [(2, bbox)])
    gmsh.model.occ.remove([(2, bbox)])

    # 4) Nettoyage
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()

    return final_body

def build_external_beams_2d(params):
    Lplate = params["Lplate"]
    tplate = params["tplate"]
    tfoot = params["tfoot"]

    print(f"Building external beams: Lplate = {Lplate}, tplate = {tplate}")

    beam_x = (params["array_width"] - Lplate) / 2.0

    if tfoot > 0.0:
        top_y = params["array_height"] + tfoot
        bot_y = -tfoot - tplate
    else:
        top_y = params["array_height"]
        bot_y = -tplate

    top_beam = create_plate_surface(beam_x, top_y, Lplate, tplate)
    bot_beam = create_plate_surface(beam_x, bot_y, Lplate, tplate)

    return [(2, top_beam)], [(2, bot_beam)]

# =============================================================================
# 4bis. RÉGLAGE MESH : RÉGULIER + TRI / QUAD
# =============================================================================

def setup_regular_mesh(elem_type, params):
    """
    - Maillage régulier : même taille partout (mesh_h).
    - elem_type = 'tri' ou 'quad' pour choisir triangulaire / quadrilatères (2D).
    """
    gmsh.model.occ.synchronize()

    h = params.get("mesh_h", params["a"] / 40.0)
    print(f"Target uniform mesh size h = {h:.4f}")

    # Taille uniforme sur tous les points géométriques
    pts = gmsh.model.getEntities(0)
    if pts:
        gmsh.model.mesh.setSize(pts, h)

    # Pour rester cohérent avec ton code précédent
    gmsh.option.setNumber("Mesh.MeshSizeMin", h)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)

    mt = elem_type.lower()
    if mt.startswith("q"):  # 'quad', 'quadrilatere', 'quadri', etc.
        print(" -> Using quadrilateral-dominant 2D mesh (recombine).")
        gmsh.option.setNumber("Mesh.RecombineAll", 1)
        # Algo frontal/quads
        gmsh.option.setNumber("Mesh.Algorithm", 8)  # Frontal-Delaunay for quads
        # On recombine toutes les surfaces 2D
        for dim, tag in gmsh.model.getEntities(2):
            gmsh.model.mesh.setRecombine(2, tag)
        # gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)
        # gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)


    else:
        print(" -> Using triangular 2D mesh.")
        gmsh.option.setNumber("Mesh.RecombineAll", 0)
        gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal

def apply_transfinite_and_periodic_2d(params,
                                      specimen_surfaces,
                                      top_beam_surfaces,
                                      bot_beam_surfaces):
    """
    - Transfinite sur les bords du carré [0,a]x[0,a] (x=0, x=a, y=0, y=a).
    - Périodicité Left<->Right, Bot<->Top.
    - Transfinite structuré sur Top_Beam et Bot_Beam.
    """

    gmsh.model.occ.synchronize()

    a = params["array_width"]
    b = params["array_height"]
    tol = 1e-6

    # -------------------------------------------------------------------------
    # 1) Récupérer les lignes de bord x=0, x=a, y=0, y=a
    # -------------------------------------------------------------------------
    all_lines = gmsh.model.getEntities(1)
    left_lines = []
    right_lines = []
    bot_lines = []
    top_lines = []

    for dim, tag in all_lines:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)

        # Lignes verticales x = 0 ou x = a, limitées à 0 <= y <= b
        if abs(xmin - 0.0) < tol and abs(xmax - 0.0) < tol and ymin >= -tol and ymax <= b + tol:
            left_lines.append(tag)
        elif abs(xmin - a) < tol and abs(xmax - a) < tol and ymin >= -tol and ymax <= b + tol:
            right_lines.append(tag)

        # Lignes horizontales y = 0 ou y = b, limitées à 0 <= x <= a
        if abs(ymin - 0.0) < tol and abs(ymax - 0.0) < tol and xmin >= -tol and xmax <= a + tol:
            bot_lines.append(tag)
        elif abs(ymin - b) < tol and abs(ymax - b) < tol and xmin >= -tol and xmax <= a + tol:
            top_lines.append(tag)

    # -------------------------------------------------------------------------
    # 2) Transfinite sur ces lignes
    # -------------------------------------------------------------------------
    n_tf_x = params["n_tf_x"]   # pour les bords horizontaux (Bot/Top)
    n_tf_y = params["n_tf_y"]   # pour les bords verticaux (Left/Right)

    def sort_by_center(tags, vertical):
        """Trie les lignes par centre en y (vertical) ou en x (horizontal)."""
        def center(tag):
            xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(1, tag)
            if vertical:
                return 0.5 * (ymin + ymax)
            else:
                return 0.5 * (xmin + xmax)
        return sorted(tags, key=center)

    left_sorted = sort_by_center(left_lines, vertical=True)
    right_sorted = sort_by_center(right_lines, vertical=True)
    bot_sorted = sort_by_center(bot_lines, vertical=False)
    top_sorted = sort_by_center(top_lines, vertical=False)

    # Transfinite courbes
    for tag in left_sorted + right_sorted:
        gmsh.model.mesh.setTransfiniteCurve(tag, n_tf_y)
    for tag in bot_sorted + top_sorted:
        gmsh.model.mesh.setTransfiniteCurve(tag, n_tf_x)

    # Groupes physiques pour les bords
    if left_lines:
        gmsh.model.addPhysicalGroup(1, left_lines, name="Left")
    if right_lines:
        gmsh.model.addPhysicalGroup(1, right_lines, name="Right")
    if bot_lines:
        gmsh.model.addPhysicalGroup(1, bot_lines, name="Bot")
    if top_lines:
        gmsh.model.addPhysicalGroup(1, top_lines, name="Top")

    # -------------------------------------------------------------------------
    # 3) Périodicité Left <-> Right, Bot <-> Top (pour les noeuds)
    # -------------------------------------------------------------------------
    # Matrice 4x4 (row-major) : translation
    # master -> slave : x_slave = x_master + dx, y_slave = y_master + dy

    # Gauche (master) -> Droite (slave)
    if left_sorted and right_sorted and len(left_sorted) == len(right_sorted):
        dx = a
        T_lr = [
            1, 0, 0, dx,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1,
        ]
        for sl, ma in zip(right_sorted, left_sorted):
            gmsh.model.mesh.setPeriodic(1, [sl], [ma], T_lr)

    # Bas (master) -> Haut (slave)
    if bot_sorted and top_sorted and len(bot_sorted) == len(top_sorted):
        dy = b
        T_bt = [
            1, 0, 0, 0,
            0, 1, 0, dy,
            0, 0, 1, 0,
            0, 0, 0, 1,
        ]
        for sl, ma in zip(top_sorted, bot_sorted):
            gmsh.model.mesh.setPeriodic(1, [sl], [ma], T_bt)

    # -------------------------------------------------------------------------
    # 4) Transfinite structuré sur les poutres Top_Beam et Bot_Beam
    # -------------------------------------------------------------------------
    n_tf_beam_x = params["n_tf_beam_x"]
    n_tf_beam_y = params["n_tf_beam_y"]

    def set_beam_transfinite(surfaces):
        for s in surfaces:
            # surface rectangulaire -> transfinite surface + recombine après
            gmsh.model.mesh.setTransfiniteSurface(s)
            boundary = gmsh.model.getBoundary([(2, s)], oriented=False)
            for dim, ltag in boundary:
                if dim != 1:
                    continue
                xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(1, ltag)
                dx = abs(xmax - xmin)
                dy = abs(ymax - ymin)
                if dx >= dy:  # segment "long" (direction x)
                    gmsh.model.mesh.setTransfiniteCurve(ltag, n_tf_beam_x)
                else:         # segment "court" (direction y = épaisseur)
                    gmsh.model.mesh.setTransfiniteCurve(ltag, n_tf_beam_y)

    set_beam_transfinite(top_beam_surfaces)
    set_beam_transfinite(bot_beam_surfaces)

# =============================================================================
# 5. EMPILAGE 3D
# =============================================================================

def extrude_and_stack_all_PMMA(body1_tags, body2_tags, body3_tags, params):
    gmsh.model.occ.synchronize()
    all_volumes_specimen = []
    all_volumes_beams = []

    h_spec = params["t"]
    g_small = params["din1"]
    g_large = params["din3"]
    n1 = params["n_spec_1mm"]
    nB_add = params["n_spec_9mm"]
    n_blocks = 1 + max(0, nB_add)

    if n1 == 1:
        H_block = h_spec
    else:
        H_block = n1 * h_spec + (n1 - 1) * g_small

    z_bottoms = []
    for b in range(n_blocks):
        z_block_bottom = b * (H_block + g_large)
        for i in range(n1):
            z_spec_bottom = z_block_bottom + i * (h_spec + g_small)
            z_bottoms.append(z_spec_bottom)

    z_stack_bottom = z_bottoms[0]
    z_stack_top = z_bottoms[-1] + h_spec
    stack_height = z_stack_top - z_stack_bottom

    for zb in z_bottoms:
        copy_list = gmsh.model.occ.copy(body1_tags)
        gmsh.model.occ.translate(copy_list, 0, 0, zb)
        ext = gmsh.model.occ.extrude(copy_list, 0, 0, h_spec)
        vols = [tag for dim, tag in ext if dim == 3]
        if vols:
            all_volumes_specimen.extend(vols)
        gmsh.model.occ.remove(copy_list)

    doutZ = params["doutZ"]
    over_z = 0.5 * doutZ
    beam_base_z = z_stack_bottom - over_z
    beam_height_z = stack_height + 2.0 * over_z

    for beam_tags in [body2_tags, body3_tags]:
        if beam_tags:
            b_copy = gmsh.model.occ.copy(beam_tags)
            gmsh.model.occ.translate(b_copy, 0, 0, beam_base_z)
            ext = gmsh.model.occ.extrude(b_copy, 0, 0, beam_height_z)
            vols = [tag for dim, tag in ext if dim == 3]
            if vols:
                all_volumes_beams.extend(vols)
            gmsh.model.occ.remove(b_copy)

    gmsh.model.occ.remove(body1_tags)
    gmsh.model.occ.remove(body2_tags)
    gmsh.model.occ.remove(body3_tags)
    gmsh.model.occ.synchronize()
    return all_volumes_specimen, all_volumes_beams

def extrude_single_specimen_with_beams(body1_tags, body2_tags, body3_tags, params):
    gmsh.model.occ.synchronize()
    all_volumes_specimen = []
    all_volumes_beams = []
    h_spec = params["t"]
    doutZ = params["doutZ"]

    ext = gmsh.model.occ.extrude(body1_tags, 0, 0, h_spec)
    vols = [tag for dim, tag in ext if dim == 3]
    all_volumes_specimen.extend(vols)

    over_z = 0.5 * doutZ
    beam_base_z = -over_z
    beam_height = h_spec + doutZ

    for beam_tags in [body2_tags, body3_tags]:
        if beam_tags:
            b_copy = gmsh.model.occ.copy(beam_tags)
            gmsh.model.occ.translate(b_copy, 0, 0, beam_base_z)
            ext = gmsh.model.occ.extrude(b_copy, 0, 0, beam_height)
            vols = [tag for dim, tag in ext if dim == 3]
            all_volumes_beams.extend(vols)
            gmsh.model.occ.remove(b_copy)

    gmsh.model.occ.remove(body1_tags)
    gmsh.model.occ.remove(body2_tags)
    gmsh.model.occ.remove(body3_tags)
    gmsh.model.occ.synchronize()

    return all_volumes_specimen, all_volumes_beams

# =============================================================================
# 6. MAIN
# =============================================================================

def main(model_dim="3D", specimen_family="b3D", elem_type="tri"):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("Labyrinth_Specimen")

    try:
        fam = specimen_family.upper()
        if fam in ("PMMA", "SP", "BP", "SPG", "BPG"):
            print("Using PMMA stacked parameters...")
            params = get_params_PMMA()
            use_stack = True
        else:
            print(f"Using 3D-printed parameters for {specimen_family}...")
            params = get_params_3D(specimen_family)
            use_stack = False

        # 1) Géométrie 2D
        array_parts = generate_array_parts(params)
        body1_2d = build_fused_specimen_2d(array_parts, params)
        body2_2d, body3_2d = build_external_beams_2d(params)

        dim_flag = model_dim.upper()

        if dim_flag == "2D":
            print("Meshing 2D model...")

            specimen_surfaces = [tag for dim, tag in body1_2d]
            top_beam_surfaces = [tag for dim, tag in body2_2d]
            bot_beam_surfaces = [tag for dim, tag in body3_2d]

            gmsh.model.addPhysicalGroup(2, specimen_surfaces, name="Specimen")
            gmsh.model.addPhysicalGroup(2, top_beam_surfaces, name="Top_Beam")
            gmsh.model.addPhysicalGroup(2, bot_beam_surfaces, name="Bot_Beam")

            # >>> NOUVEAU : transfinite + périodicité + transfinite beams
            apply_transfinite_and_periodic_2d(
                params,
                specimen_surfaces,
                top_beam_surfaces,
                bot_beam_surfaces,
            )

            # Maillage régulier + TRI/QUAD
            setup_regular_mesh(elem_type, params)

            gmsh.model.mesh.generate(2)
            gmsh.write(f"output_2D_{specimen_family}_{elem_type}.inp")


        else:
            print("Meshing 3D model...")
            if use_stack:
                vol_specimens, vol_beams = extrude_and_stack_all_PMMA(
                    body1_2d, body2_2d, body3_2d, params
                )
            else:
                vol_specimens, vol_beams = extrude_single_specimen_with_beams(
                    body1_2d, body2_2d, body3_2d, params
                )

            if vol_specimens:
                gmsh.model.addPhysicalGroup(3, vol_specimens, name="Specimen")
            if vol_beams:
                gmsh.model.addPhysicalGroup(3, vol_beams, name="Beams")

            # Maillage régulier (le choix TRI/QUAD agit sur les surfaces 2D)
            setup_regular_mesh(elem_type, params)

            gmsh.model.mesh.generate(3)
            gmsh.write(f"output_3D_{specimen_family}_{elem_type}.inp")

        gmsh.fltk.run()

    except Exception as e:
        print("Error:", e)
    finally:
        gmsh.finalize()

if __name__ == "__main__":
    dim = "3D"      # "2D" ou "3D"
    fam = "b3D"     # par défaut B3D
    elem_type = "tri"  # "tri" ou "quad"

    if len(sys.argv) > 1:
        dim = sys.argv[1]
    if len(sys.argv) > 2:
        fam = sys.argv[2]
    if len(sys.argv) > 3:
        elem_type = sys.argv[3]

    main(model_dim=dim, specimen_family=fam, elem_type=elem_type)

