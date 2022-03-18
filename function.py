import numpy as np
from qgis.core import QgsPointXY

def ortho_line(inters_prec, inters_suiv, inters_courant, largeur_vallee):
    x1 = inters_prec.x()
    y1 = inters_prec.y()
    x2 = inters_suiv.x()
    y2 = inters_suiv.y()
    x3 = inters_courant.x()
    y3 = inters_courant.y()
    # y = coeff_A * x + coeff_B est l'équation de la droite perpendiculaire à X1-X2 passant par X3
    coeff_A = (x2 - x1) / (y1 - y2)
    coeff_B = y3 - x3 * ((x2 - x1) / (y1 - y2))
    # linestart et lineend sont les 2 extrémités de la polyligne perpendiculaire à X1-X2 passant par X3
    if coeff_A < -1:
        y_start = y3 + largeur_vallee
        x_start = (y_start - coeff_B) / coeff_A
        y_end = y3 - largeur_vallee
        x_end = (y_end - coeff_B) / coeff_A
    elif coeff_A <= 0 and coeff_A >= -1:
        x_start = x3 + largeur_vallee
        y_start = coeff_A *  x_start + coeff_B
        x_end = x3 - largeur_vallee
        y_end = coeff_A *  x_end + coeff_B
    elif coeff_A > 0 and coeff_A <= 1:
        x_start = x3 + largeur_vallee
        y_start = coeff_A *  x_start + coeff_B
        x_end = x3 - largeur_vallee
        y_end = coeff_A *  x_end + coeff_B
    elif coeff_A > 1:
        y_start = y3 + largeur_vallee
        x_start = (y_start - coeff_B) / coeff_A
        y_end = y3 - largeur_vallee
        x_end = (y_end - coeff_B) / coeff_A
    return x_start, y_start, x_end, y_end

def vertex_add(geom, couche, feat_id, x, y, tol=0.01):
    p1, at, b1, after, d1 = geom.closestVertex(QgsPointXY(x, y))
    dist, p2, to, _ = geom.closestSegmentWithContext(QgsPointXY(x, y))
    if at == 0:
        if dist < tol:
            # insert into first segment
            couche.insertVertex(x, y, feat_id, after)
            geom.insertVertex(x, y, after)
        else:
            # insert before first vertex
            couche.insertVertex(x, y, feat_id, 0)
            geom.insertVertex(x, y, 0)
    elif after == -1:
        if dist < tol:
            # insert after last vertex
            couche.insertVertex(x, y, feat_id, at)
            geom.insertVertex(x, y, at)
        else:
            # insert into last segment
            couche.insertVertex(x, y, feat_id, at - 1)
            geom.insertVertex(x, y, at - 1)
    return geom

def tronque_profil(x_str, z_str, dist_rg, dist_rd):
    flag_rg = False
    flag_rd = False
    x_tronque = []
    z_tronque = []
    for i, xx in enumerate(x_str):
        x = float(xx)
        if x < dist_rg:
            pass
        elif x >= dist_rg and x <= dist_rd:
            if x == dist_rg or flag_rg:
                x_tronque.append(x)
                z_tronque.append(float(z_str[i]))
                flag_rg = True
            elif x == dist_rd:
                x_tronque.append(x)
                z_tronque.append(float(z_str[i]))
                flag_rd = True
            else:
                x_tronque.append(dist_rg)
                z_tronque.append((dist_rg - float(x_str[i - 1])) * (float(z_str[i]) - float(z_str[i - 1])) / (x - float(x_str[i - 1])) + float(z_str[i - 1]))
                flag_rg = True
        elif x > dist_rd and not flag_rd:
            x_tronque.append(dist_rd)
            z_tronque.append((dist_rd - float(x_str[i - 1])) * (float(z_str[i]) - float(z_str[i - 1])) / (x - float(x_str[i - 1])) + float(z_str[i - 1]))
            flag_rd = True
    return x_tronque, z_tronque

def calcul_planim_sh(ref_plani, x_z):
    # Il faut d'abord s'assurer que x_z est trié dans l'ordre croissant
    ref_sh = np.zeros(np.shape(ref_plani)[0])
    if np.shape(ref_plani)[0] > 1:
        pas_plani = ref_plani[1] - ref_plani[0]
        for i in range(len(x_z) - 1):
            d1 = x_z[i][0] # distance
            z1 = x_z[i][1] # cote
            d2 = x_z[i+1][0]
            z2 = x_z[i+1][1]
            en_eau = False
            for ip, p in enumerate(ref_plani):
                sh_cur = 0
                if en_eau:
                    sh_cur = (d2 - d1) * pas_plani
                elif (p + pas_plani) <= z1 and (p + pas_plani) <= z2:
                    sh_cur = 0
                else:
                    sh_cur = (p + pas_plani - max(z1, z2)) * (d2 - d1) + (math.fabs(z2 - z1) * (d2 - d1) * 0.5)
                    en_eau = True
                ref_sh[ip] += round(sh_cur * 1000) / 1000
    ref_sh2 = np.around(ref_sh, decimals=3)
    return ref_sh2

def tri_pt_profils(x_z_t, ordre_croissant):
    left = []
    equal = []
    right = []
    if ordre_croissant:
        if len(x_z_t) > 1:
            pivot = x_z_t[0][0]
            for i, x in enumerate(x_z_t):
                if x[0] < pivot:
                    left.append([x_z_t[i][0], x_z_t[i][1]])
                elif x[0] == pivot:
                    equal.append([x_z_t[i][0], x_z_t[i][1]])
                elif x[0] > pivot:
                    right.append([x_z_t[i][0], x_z_t[i][1]])
            return tri_pt_profils(left, ordre_croissant) + equal + tri_pt_profils(right, ordre_croissant)
        else:
            return x_z_t
    else: # Ordre décroissant
        if len(x_z_t) > 1:
            pivot = x_z_t[0][0]
            for i, x in enumerate(x_z_t):
                if x[0] > pivot:
                    left.append([x_z_t[i][0], x_z_t[i][1]])
                elif x[0] == pivot:
                    equal.append([x_z_t[i][0], x_z_t[i][1]])
                elif x[0] < pivot:
                    right.append([x_z_t[i][0], x_z_t[i][1]])
            return tri_pt_profils(left, ordre_croissant) + equal + tri_pt_profils(right, ordre_croissant)
        else:
            return x_z_t

def creation_profils(ref_plani, planim, dist_centre_profil, max_rg, max_rd): 
    # ref_plani est le tableau des altitudes planimétrées
    # planim est un tableau numpy avec le planimétrage du profil à créer
    # dist_centre_profil donne la valeur sur laquelle centrer le profil à créer
    x_gauche = ""
    z_gauche = ""
    x_droite = ""
    z_droite = ""
    x_z_gauche = []
    x_z_droite = []
    nbr_pas = len(ref_plani)
    pas = ref_plani[1] - ref_plani[0]
    miroir_p = planim[nbr_pas - 1] / pas
    for indice in range(nbr_pas):
        if indice == 0:
            esp = ""
        else:
            esp = " "
        i = nbr_pas - indice - 1
        zg = min(ref_plani[i], max_rg)
        zd = min(ref_plani[i], max_rd)
        miroir = (2 * planim[i]) / (3 * pas) + (miroir_p / 3)
        miroir_p = miroir
        if miroir > 0 and planim[i] > 0:
            if ref_plani[i] <= max_rg:
                x_gauche += str(dist_centre_profil - (miroir / 2)) + " "
                z_gauche += str(zg) + " "
                x_z_gauche.append([dist_centre_profil - (miroir / 2), zg])
            if ref_plani[i] <= max_rd:
                x_droite = str(dist_centre_profil + (miroir / 2)) + esp + x_droite
                z_droite = str(zd) + esp + z_droite
                x_z_droite = [[dist_centre_profil + (miroir / 2), zd]] + x_z_droite
        else:
            break
    return x_gauche + x_droite, z_gauche + z_droite, x_z_gauche + x_z_droite

def tri_profils(profils_source):
    left = []
    equal = []
    right = []
    if len(profils_source) > 1:
        pivot = profils_source[0]['absc']
        for profil in profils_source:
            if profil['absc'] < pivot:
                left.append(profil)
            elif profil['absc'] == pivot:
                equal.append(profil)
            elif profil['absc'] > pivot:
                right.append(profil)
        return tri_profils(left) + equal + tri_profils(right)
    else:
        return profils_source

def profils_amont_aval(abscisse, profils_source):
    profils_source_tries = tri_profils(profils_source)
    profil_amont = profils_source[0]
    profil_aval = profils_source[1]
    trouve = False
    for indice in range(len(profils_source_tries) - 1):
        if abscisse >= profils_source_tries[indice]['absc'] and abscisse <= profils_source_tries[indice + 1]['absc']:
            profil_amont = profils_source_tries[indice]
            profil_aval = profils_source_tries[indice + 1]
            trouve = True
    return trouve, profil_amont, profil_aval