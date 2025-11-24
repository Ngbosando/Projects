import math
import gmsh

# =============================================================================
# 1. PARAMÈTRES (VOTRE CODE ORIGINAL - INCHANGÉ)
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
        "n_spec_1mm": 1,
        "n_spec_9mm": 0,
        "doutZ": 4.0,
        "eps": 0.15,
        "seq_card": [1, 0, 1, 0, 1, 0, 1, 0, 1],
        "seq_diag": [0, 1, 0, 1, 0, 1, 0, 1, 0],
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
    }

# =============================================================================
# 2. GÉOMÉTRIQUES 2D (MODIF JUSTE SUR ANCHOR)
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
    # INCHANGÉ
    length = (r_outer - r_inner) + 2 * eps
    start_r = r_inner - eps
    br = gmsh.model.occ.addRectangle(start_r, -width / 2.0, 0, length, width)
    gmsh.model.occ.rotate([(2, br)], 0, 0, 0, 0, 0, 1, math.radians(angle_deg))
    gmsh.model.occ.translate([(2, br)], cx, cy, 0)
    return br

def create_anchor_surface(cx, cy, angle_deg, params):
    # --- MODIFICATION DE L'ANCRE (TANGENTE & LONGUE) ---
    rad = math.radians(angle_deg)
    
    # 1. Point de tangence sur le cercle externe (R_start)
    # C'est là que l'ancre touche l'anneau
    R = params["R_start"]
    tan_x = cx + R * math.cos(rad)
    tan_y = cy + R * math.sin(rad)
    
    # 2. Orientation : Perpendiculaire au rayon (Tangente)
    # L'angle de l'ancre est l'angle du rayon + 90 degrés
    anchor_angle = rad + math.pi/2
    
    # 3. Dimensions : Longue pour atteindre les poutres
    # On la fait assez longue pour qu'elle rentre franchement dans la poutre externe
    # Distance approx centre -> bord = 25. Ancre part de R=23.75.
    length = 40.0 
    width = params["tb"]
    
    # Création du rectangle
    # On le centre en largeur (-width/2)
    # On le décale en longueur (-length/2) pour qu'il soit centré sur le point de tangence
    b = gmsh.model.occ.addRectangle(-length/2.0, -width/2.0, 0, length, width)
    
    # Rotation et Translation au point de tangence
    gmsh.model.occ.rotate([(2, b)], 0, 0, 0, 0, 0, 1, anchor_angle)
    gmsh.model.occ.translate([(2, b)], tan_x, tan_y, 0)
    
    return b

def create_plate_surface(x, y, width, height):
    return gmsh.model.occ.addRectangle(x, y, 0, width, height)

# =============================================================================
# 3. COMPOSANTS 2D (INCHANGÉ)
# =============================================================================

def generate_unit_cell_parts(cx, cy, params):
    parts = []
    # Ancrages diagonaux
    for ang in [45, 135, 225, 315]:
        parts.append(create_anchor_surface(cx, cy, ang, params))

    # Anneaux + poutre connecteur
    current_R = params["R_start"]
    tb = params["tb"]
    for k in range(params["num_void_layers"]):
        parts.append(create_ring_surface(cx, cy, current_R, tb))
        current_R -= tb

        void_outer = current_R
        void_inner = current_R - tb

        if params["seq_card"][k]:
            for ang in [0, 90, 180, 270]:
                parts.append(create_bridge_surface(cx, cy, void_outer, void_inner,
                                                   tb, ang, params["eps"]))
        if params["seq_diag"][k]:
            for ang in [45, 135, 225, 315]:
                parts.append(create_bridge_surface(cx, cy, void_outer, void_inner,
                                                   tb, ang, params["eps"]))
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
# 4. ASSEMBLAGE (MODIFIÉ POUR SOUDER ANCRES ET POUTRES)
# =============================================================================

def build_fused_specimen_2d(array_parts, params):
    """
    MODIFICATION : 
    On inclut les poutres externes DANS le fragment pour que les ancres
    soient soudées proprement aux poutres.
    """
    print("Building Body 1 (specimen + beams) in 2D with Fragment...")
    parts_to_fuse = [(2, p) for p in array_parts]

    tfoot = params["tfoot"]
    eps = params["eps"]
    
    # -- GESTION DES PIEDS ET POUTRES EXTERNES --
    # On ajoute tout ici pour que le fragment soude tout ensemble
    
    Lplate = params["Lplate"]
    tplate = params["tplate"]
    beam_x = (params["array_width"] - Lplate) / 2.0

    if tfoot > 0.0:
        # Pieds internes
        foot_bot = create_plate_surface(0.0, -tfoot, params["array_width"], tfoot + eps)
        foot_top_y = params["array_height"] - eps
        foot_top = create_plate_surface(0.0, foot_top_y, params["array_width"], tfoot + eps)
        parts_to_fuse.append((2, foot_bot))
        parts_to_fuse.append((2, foot_top))
        
        # Poutres externes
        top_y = params["array_height"] + tfoot
        bot_y = -tfoot - tplate
        
        y_min = -tfoot # Pour la boite de coupe
        total_h_box = (top_y + tplate) - bot_y # Pas exactement, on veut garder les poutres
    else:
        # Poutres externes directes
        top_y = params["array_height"]
        bot_y = -tplate
        y_min = 0.0

    # Création des poutres externes
    top_beam = create_plate_surface(beam_x, top_y, Lplate, tplate)
    bot_beam = create_plate_surface(beam_x, bot_y, Lplate, tplate)
    parts_to_fuse.append((2, top_beam))
    parts_to_fuse.append((2, bot_beam))

    gmsh.model.occ.synchronize()

    # --- FRAGMENTATION GLOBALE ---
    # C'est ce qui soude l'ancre à la poutre
    fragmented_shapes, _ = gmsh.model.occ.fragment(parts_to_fuse, [])
    
    # Fusionner le puzzle en un seul morceau si possible (pour nettoyer les lignes internes)
    if len(fragmented_shapes) > 1:
        fused_body_list, _ = gmsh.model.occ.fuse([fragmented_shapes[0]], fragmented_shapes[1:])
    else:
        fused_body_list = fragmented_shapes
        
    # --- NETTOYAGE FIN ---
    # On coupe ce qui dépasse des poutres (les bouts d'ancres inutiles)
    # On définit une boite qui englobe tout de la poutre du bas à la poutre du haut
    total_y_min = bot_y
    total_y_max = top_y + tplate
    total_h = total_y_max - total_y_min
    
    # On prend large en X pour ne rien couper sur les côtés, juste en haut/bas
    bbox = gmsh.model.occ.addRectangle(beam_x, total_y_min, 0, Lplate, total_h)
    
    final_body, _ = gmsh.model.occ.intersect(fused_body_list, [(2, bbox)])
    gmsh.model.occ.remove([(2, bbox)])

    # Soudure finale
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()

    return final_body

def build_external_beams_2d(params):
    # CETTE FONCTION NE SERT PLUS RIEN car les poutres sont intégrées au dessus
    # Mais on la garde vide ou on renvoie une liste vide pour ne pas casser le main
    return [], []

# =============================================================================
# 5. EMPILAGE 3D (ADAPTÉ POUR CORPS UNIQUE)
# =============================================================================

def extrude_and_stack_all_PMMA(body1_tags, body2_tags, body3_tags, params):
    # body1_tags contient maintenant TOUT (labyrinthe + poutres)
    # body2 et body3 sont vides
    gmsh.model.occ.synchronize()
    all_volumes = []

    h_spec = params["t"]
    g_small = params["din1"]
    g_large = params["din3"]
    n1 = params["n_spec_1mm"]
    nB_add = params["n_spec_9mm"]
    n_blocks = 1 + max(0, nB_add)

    if n1 == 1: H_block = h_spec
    else: H_block = n1 * h_spec + (n1 - 1) * g_small

    z_bottoms = []
    for b in range(n_blocks):
        z_block_bottom = b * (H_block + g_large)
        for i in range(n1):
            z_spec_bottom = z_block_bottom + i * (h_spec + g_small)
            z_bottoms.append(z_spec_bottom)

    # Extrusion simple de tout le bloc soudé
    for zb in z_bottoms:
        copy_list = gmsh.model.occ.copy(body1_tags)
        gmsh.model.occ.translate(copy_list, 0, 0, zb)
        ext = gmsh.model.occ.extrude(copy_list, 0, 0, h_spec)
        vols = [tag for dim, tag in ext if dim == 3]
        if vols: all_volumes.extend(vols)
        gmsh.model.occ.remove(copy_list)

    # Note: Dans cette version soudée, on ne gère pas doutZ séparément pour les poutres
    # car elles font partie du même corps 2D. 
    # Si c'est critique, il faut revenir à une séparation complexe.
    # Ici, on priorise la géométrie valide.

    gmsh.model.occ.remove(body1_tags)
    gmsh.model.occ.synchronize()
    return all_volumes, [] # Pas de poutres séparées

def extrude_single_specimen_with_beams(body1_tags, body2_tags, body3_tags, params):
    # body1_tags contient TOUT
    gmsh.model.occ.synchronize()
    h_spec = params["t"]
    
    # Extrusion unique
    ext = gmsh.model.occ.extrude(body1_tags, 0, 0, h_spec)
    all_volumes = [tag for dim, tag in ext if dim == 3]

    gmsh.model.occ.remove(body1_tags)
    gmsh.model.occ.synchronize()

    return all_volumes, [] # Pas de poutres séparées

# =============================================================================
# 6. MAIN (ADAPTÉ)
# =============================================================================

def main(model_dim="3D", specimen_family="b3D"):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("Labyrinth_Fixed_Anchors")

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
        
        # 2) Assemblage COMPLET (Labyrinthe + Poutres soudés)
        final_body = build_fused_specimen_2d(array_parts, params)
        
        # Poutres vides car intégrées
        body2_2d, body3_2d = [], []

        dim_flag = model_dim.upper()

        if dim_flag == "2D":
            print("Meshing 2D model...")
            specimen_surfaces = [tag for dim, tag in final_body]
            gmsh.model.addPhysicalGroup(2, specimen_surfaces, name="Whole_Structure")

            gmsh.option.setNumber("Mesh.MeshSizeMin", params["a"] / 60.0)
            gmsh.option.setNumber("Mesh.MeshSizeMax", params["a"] / 15.0)
            gmsh.model.mesh.generate(2)
            gmsh.write(f"output_2D_{specimen_family}.inp")

        else:
            print("Meshing 3D model...")
            if use_stack:
                vol_specimens, _ = extrude_and_stack_all_PMMA(final_body, [], [], params)
            else:
                vol_specimens, _ = extrude_single_specimen_with_beams(final_body, [], [], params)

            if vol_specimens:
                gmsh.model.addPhysicalGroup(3, vol_specimens, name="Whole_Structure_3D")

            gmsh.option.setNumber("Mesh.MeshSizeMin", params["a"] / 60.0)
            gmsh.option.setNumber("Mesh.MeshSizeMax", params["a"] / 15.0)
            gmsh.model.mesh.generate(3)
            gmsh.write(f"output_3D_{specimen_family}.inp")

        gmsh.fltk.run()

    except Exception as e:
        print("Error:", e)
    finally:
        gmsh.finalize()

if __name__ == "__main__":
    main(model_dim="3D", specimen_family="b3D")