import numpy as np


class CutTable():
    """docstring for CutTable"""

    def __init__(self):
        self.n_flags = 6

        self.vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
        self.tet = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
        self.edges = np.array([[0, 1], [1, 2], [2, 0], [0, 3], [1, 3], [2, 3]])

        self.table = [None] * (2**self.n_flags)
        self.track_faces = [None] * (2**self.n_flags)
        self.original_faces = [None] * (2**self.n_flags)
        self.edges_table = [None] * (2**self.n_flags)
        self.new_points = [None] * (2**self.n_flags)

        self.new_points[0] = []
        self.add_tet(0, 0, [0, 1, 2, 3], [], [], [])
        self.add_edge(0, [])

        self.create_single_vertex_cut([0, 2, 3])
        self.create_single_vertex_cut([0, 1, 4])
        self.create_single_vertex_cut([1, 2, 5])
        self.create_single_vertex_cut([3, 4, 5])


        self.create_two_vertex_cut([0, 1, 3, 5], False)
        self.create_two_vertex_cut([1, 2, 3, 4], True)
        self.create_two_vertex_cut([0, 2, 4, 5], False)


        self.create_single_vertex_on_plane_cut([0, 1])
        self.create_single_vertex_on_plane_cut([0, 2])
        self.create_single_vertex_on_plane_cut([0, 3])
        self.create_single_vertex_on_plane_cut([0, 4])

        self.create_single_vertex_on_plane_cut([1, 2])
        self.create_single_vertex_on_plane_cut([1, 4])
        self.create_single_vertex_on_plane_cut([1, 5])

        self.create_single_vertex_on_plane_cut([2, 3])
        self.create_single_vertex_on_plane_cut([2, 5])

        self.create_single_vertex_on_plane_cut([3, 4])
        self.create_single_vertex_on_plane_cut([3, 5])

        self.create_single_vertex_on_plane_cut([4, 5])


        self.create_two_vertex_on_plane_cut(0)
        self.create_two_vertex_on_plane_cut(1)
        self.create_two_vertex_on_plane_cut(2)
        self.create_two_vertex_on_plane_cut(3)
        self.create_two_vertex_on_plane_cut(4)
        self.create_two_vertex_on_plane_cut(5)

        # [0, 1, 3, 5]
        self.create_three_edges_cut1([1, 3, 5])
        self.create_three_edges_cut2([0, 3, 5])
        self.create_three_edges_cut2([0, 1, 5])
        self.create_three_edges_cut1([0, 1, 3])

        # [1, 2, 3, 4]
        self.create_three_edges_cut1([2, 3, 4])
        self.create_three_edges_cut1([1, 3, 4])
        self.create_three_edges_cut1([1, 2, 4])
        self.create_three_edges_cut1([1, 2, 3])


        # [0, 2, 4, 5]
        self.create_three_edges_cut2([2, 4, 5])
        self.create_three_edges_cut2([0, 4, 5])
        self.create_three_edges_cut2([0, 2, 5])
        self.create_three_edges_cut2([0, 2, 4])




        self.create_two_edges_cut([0, 5])
        self.create_two_edges_cut([1, 3])
        self.create_two_edges_cut([2, 4])



        for i in range(len(self.edges_table)):
            if self.edges_table[i] is None:
                continue
            for j in range(len(self.edges_table[i])):
                for k in range(len(self.edges_table[i][j])):
                    e = self.edges_table[i][j][k]
                    assert(len(e) == 2)
                    if e[1] < e[0]:
                        tmp = e[0]
                        self.edges_table[i][j][k][0] = e[1]
                        self.edges_table[i][j][k][1] = tmp
                self.edges_table[i][j].sort(key=lambda tup: (tup[0],tup[1]))


    def index_for_edges(self, edges):
        flags = [0] * self.n_flags

        for e in edges:
            flags[5-e] = 1

        return self.flags_to_int(flags)

    def __str__(self):
        res = ""
        index = 0
        for i in range(len(self.table)):
            confs = self.table[i]
            edges = self.edges_table[i]
            track = self.track_faces[i]
            orig = self.original_faces[i]
            if confs:
                res += "---------------\n"
                res += self.int_to_flags(i) + "\n"
                for i in range(len(confs)):
                    res += "\tconfiguration {}".format(i+1)
                    res += "\t"+str(confs[i]) + "\t"
                    res += str(track[i])+"\t"
                    res += str(orig[i])+"\n"
                res += "\tedges\t\t\t"+str(edges)+"\n\n"
                res += "\n"

                index += 1

        res += "\n\nn: {}".format(index)
        return res

    def int_to_flags(self, int_flag):
        tmp = self.int_to_flags_rec(int_flag)

        if len(tmp) < self.n_flags:
            tmp = "0"*(self.n_flags-len(tmp)) + tmp

        return tmp

    def int_to_flags_rec(self, int_flag):
        return str(int_flag) if int_flag<=1 else self.int_to_flags_rec(int_flag >> 1) + str(int_flag & 1)

    def flags_to_int(self, flags):
        int_flag = 0
        mask = 1

        n_args = len(flags)

        assert(n_args == self.n_flags)

        for i in range(self.n_flags):
            c = flags[self.n_flags-i-1]
            if c == 1:
                int_flag = int_flag | mask

            mask = mask << 1

        return int_flag


    def add_edges(self, index, indices1, indices2):
        indices = []

        n_found = 0
        found_index = -1

        for e1 in indices1:
            # indices.append(e1)
            for i in range(len(indices2)):
                e2 = indices2[i]
                if e1[0] == e2[0] and e1[1] == e2[1]:
                    found_index = i
                    n_found += 1
                    continue
                if e1[1] == e2[0] and e1[0] == e2[1]:
                    found_index = i
                    n_found += 1
                    continue

        assert(n_found == 1)

        for i in range(len(indices1)):
            if i == found_index:
                continue
            indices.append(indices1[i])

        for i in range(len(indices2)):
            if i == found_index:
                continue
            indices.append(indices2[i])

        assert(len(indices) == 4)

        self.add_edge(index, indices)


    def add_edge(self, index, indices):
        self.edges_table[index].append(indices)

    def add_tet(self, index, conf_index, indices, points, new_face_vertices, edges):
        if self.table[index] is None:
            assert(self.edges_table[index] is None)
            assert(self.track_faces[index] is None)
            assert(self.original_faces[index] is None)
            self.table[index] = []
            self.track_faces[index] = []
            self.original_faces[index] = []
            self.edges_table[index] = []

        for i in range(len(self.table[index]), conf_index+1):
            self.table[index].append([])
            self.track_faces[index].append([])
            self.original_faces[index].append([])


        ###########################################
        ##############reorient
        vertices = self.vertices
        if len(points) > 0:
            vertices = np.append(vertices, points, axis=0)
        if len(self.new_points[index]) > 0:
            vertices = np.append(vertices, self.new_points[index], axis=0)


        a = vertices[indices[3], :] - vertices[indices[0], :]
        b = vertices[indices[1], :] - vertices[indices[0], :]
        c = vertices[indices[2], :] - vertices[indices[0], :]

        area = np.dot(a, np.cross(b, c))
        # print("area",area, indices, a, b, c)

        if area < 0:
            tmp = indices[1]
            indices[1] = indices[3]
            indices[3] = tmp

        a = vertices[indices[3], :] - vertices[indices[0], :]
        b = vertices[indices[1], :] - vertices[indices[0], :]
        c = vertices[indices[2], :] - vertices[indices[0], :]

        area = np.dot(a, np.cross(b, c))
        assert(area > 0)
        #######################################################

        #######################################################
        #########################new face tracking
        tmp = [False, False, False, False]
        for i in range(len(indices)):
            if indices[i] in new_face_vertices:
                tmp[i] = True

        other = -1
        count = 0

        for t in range(len(tmp)):
            if tmp[t]:
                count += 1
            else:
                other = t

        tracking = [False, False, False, False]
        if count == 3:
            assert(other >= 0)
            tracking[other] = True
        #######################################################

        #######################################################
        #########################old face tracking
        count = 0
        original = [-1, -1, -1, -1]

        for v in new_face_vertices:
            if v >= 4:
                count += 1
        assert(len(edges) == count)
        # print(tracking, indices)
        for k in range(len(self.tet)):
            t = self.tet[k]
            tmp = [False, False, False, False]
            for j in t:
                tmp[j] = True
            real_index = -1
            for j in range(len(tmp)):
                if not tmp[j]:
                    real_index = j
                    break
            assert(real_index>=0)

            # print(real_index, t)
            loc_face = [indices[j] for j in t]
            # print("before",loc_face)

            for i in range(len(loc_face)):
                if loc_face[i] - 4 < len(edges) and loc_face[i] >=4:
                    ee = self.edges[edges[loc_face[i]-4]]

                    loc_face[i] = ee[0]
                    loc_face.append(ee[1])

            # print("after",loc_face)
            tmp = [False, False, False, False]
            for i in loc_face:
                if i < 4:
                    tmp[i] = True

            # print("tmp",tmp)
            count = 0
            for t in range(len(tmp)):
                if tmp[t]:
                    count += 1
                else:
                    other = t
            if count == 3 and max(loc_face) < 4:
                original[real_index] = other
        #######################################################


        self.table[index][conf_index].append(indices)
        self.track_faces[index][conf_index].append(tracking)
        self.original_faces[index][conf_index].append(original)


    def prism_cut(self, nodes, pts, max_idx):
        res = []
        edges = []

        tmp = []
        tmp.append([nodes[3], nodes[4], nodes[5], nodes[1]])
        tmp.append([nodes[5], nodes[1], nodes[2], nodes[0]])
        tmp.append([nodes[3], nodes[0], nodes[5], nodes[1]])
        edges.append([[nodes[0], nodes[5]], [nodes[5], nodes[1]], [nodes[3], nodes[1]]])
        res.append(tmp)


        ####################################################################
        tmp = []
        tmp.append([nodes[3], nodes[4], nodes[5], nodes[1]])
        tmp.append([nodes[5], nodes[1], nodes[2], nodes[3]])
        tmp.append([nodes[3], nodes[0], nodes[2], nodes[1]])
        edges.append([[nodes[3], nodes[2]], [nodes[5], nodes[1]], [nodes[3], nodes[1]]])
        res.append(tmp)




        tmp = []
        tmp.append([nodes[0], nodes[4], nodes[5], nodes[1]])
        tmp.append([nodes[5], nodes[1], nodes[2], nodes[0]])
        tmp.append([nodes[3], nodes[4], nodes[5], nodes[0]])
        edges.append([[nodes[0], nodes[5]], [nodes[5], nodes[1]], [nodes[0], nodes[4]]])
        res.append(tmp)

        ####################################################################
        tmp = []
        tmp.append([nodes[3], nodes[4], nodes[2], nodes[1]])
        tmp.append([nodes[5], nodes[4], nodes[2], nodes[3]])
        tmp.append([nodes[3], nodes[0], nodes[2], nodes[1]])
        edges.append([[nodes[3], nodes[2]], [nodes[4], nodes[2]], [nodes[3], nodes[1]]])
        res.append(tmp)

        tmp = []
        tmp.append([nodes[0], nodes[4], nodes[2], nodes[1]])
        tmp.append([nodes[0], nodes[2], nodes[5], nodes[4]])
        tmp.append([nodes[3], nodes[0], nodes[5], nodes[4]])
        edges.append([[nodes[0], nodes[5]], [nodes[4], nodes[2]], [nodes[0], nodes[4]]])
        res.append(tmp)



        ####################################################################
        tmp = []
        tmp.append([nodes[0], nodes[4], nodes[2], nodes[1]])
        tmp.append([nodes[3], nodes[2], nodes[5], nodes[4]])
        tmp.append([nodes[3], nodes[0], nodes[2], nodes[4]])
        edges.append([[nodes[2], nodes[3]], [nodes[4], nodes[2]], [nodes[0], nodes[4]]])
        res.append(tmp)



        ################invalid
        tmp = []
        tmp.append([nodes[0], nodes[3], nodes[2], max_idx])
        tmp.append([nodes[2], nodes[5], nodes[3], max_idx])
        tmp.append([nodes[4], nodes[5], nodes[1], max_idx])
        tmp.append([nodes[2], nodes[5], nodes[1], max_idx])
        tmp.append([nodes[0], nodes[4], nodes[1], max_idx])
        tmp.append([nodes[3], nodes[4], nodes[0], max_idx])
        tmp.append([nodes[0], nodes[1], nodes[2], max_idx])
        tmp.append([nodes[3], nodes[4], nodes[5], max_idx])
        edges.append([[nodes[2], nodes[3]], [nodes[5], nodes[1]], [nodes[0], nodes[4]]])
        res.append(tmp)

        ################invalid
        tmp = []
        tmp.append([nodes[0], nodes[3], nodes[5], max_idx])
        tmp.append([nodes[0], nodes[5], nodes[2], max_idx])
        tmp.append([nodes[4], nodes[5], nodes[2], max_idx])
        tmp.append([nodes[2], nodes[4], nodes[1], max_idx])
        tmp.append([nodes[3], nodes[1], nodes[0], max_idx])
        tmp.append([nodes[3], nodes[4], nodes[1], max_idx])
        tmp.append([nodes[0], nodes[1], nodes[2], max_idx])
        tmp.append([nodes[3], nodes[4], nodes[5], max_idx])
        edges.append([[nodes[0], nodes[5]], [nodes[4], nodes[2]], [nodes[3], nodes[1]]])
        res.append(tmp)





        tmp = np.append(self.vertices, pts, axis=0)
        p = np.sum(tmp[nodes, :], axis=0)/len(nodes)

        # pt = self.vertices

        return res, edges, p


    def create_two_vertex_cut(self, edges, flag):
        index = self.index_for_edges(edges)
        assert(self.table[index] is None)

        indices = [True]*6
        for e in edges:
            indices[e] = False

        missing_edge = -1
        for i in range(len(indices)):
            if indices[i]:
                missing_edge = i
                break

        assert(missing_edge >= 0)

        index1 = min(self.edges[missing_edge][0],self.edges[missing_edge][1])
        index2 = max(self.edges[missing_edge][0],self.edges[missing_edge][1])

        tip1 = -1
        tip2 = -1
        if index1 == 0:
            if index2 == 1:
                tip1 = 2
                tip2 = 3
            elif index2 == 2:
                tip1 = 1
                tip2 = 3
            elif index2 == 3:
                tip1 = 1
                tip2 = 2
        elif index1 == 1:
            if index2 == 2:
                tip1 = 0
                tip2 = 3
            elif index2 == 3:
                tip1 = 0
                tip2 = 2
        if index1 == 2:
            if index2 == 3:
                tip1 = 0
                tip2 = 1

        assert(tip1 >= 0)
        assert(tip2 >= 0)

        assert(tip1 != tip2)
        assert(index1 != index2)

        assert(index1 != tip2)
        assert(tip1 != index2)
        assert(index2 != tip2)
        assert(tip1 != index1)


        e0 = self.edges[edges[0], :]
        e1 = self.edges[edges[1], :]
        e2 = self.edges[edges[2], :]
        e3 = self.edges[edges[3], :]

        p0 = (self.vertices[e0[0], :]+self.vertices[e0[1], :]) / 2
        p1 = (self.vertices[e1[0], :]+self.vertices[e1[1], :]) / 2
        p2 = (self.vertices[e2[0], :]+self.vertices[e2[1], :]) / 2
        p3 = (self.vertices[e3[0], :]+self.vertices[e3[1], :]) / 2

        # print(tip1,tip2, index1, index2)
        # print(p0, p1, p2, p3)

        if flag:
            nodes1 = [index1, 5, 6, index2, 4, 7]
            nodes2 = [tip1, 5, 4, tip2, 6, 7]
        else:
            nodes1 = [index1, 4, 6, index2, 5, 7]
            nodes2 = [tip1, 4, 5, tip2, 6, 7]

        # print(nodes2)
        cuts1, edges1, pp1 = self.prism_cut(nodes1, [p0, p1, p2, p3], 8)
        cuts2, edges2, pp2 = self.prism_cut(nodes2, [p0, p1, p2, p3], 9)
        self.new_points[index] = []
        self.new_points[index].append(pp1)
        self.new_points[index].append(pp2)

        conf_index = 0

        matching = [
            [0, 1, 2, 6], [0, 1, 2, 6], [0, 1, 2, 6],
            [3, 4, 5, 7], [3, 4, 5, 7], [3, 4, 5, 7],

            [0, 1, 2, 6],
            [3, 4, 5, 7]
        ]
        for i in range(len(cuts1)):
            cut1 = cuts1[i]

            for j in matching[i]:
                if j > len(cuts2):
                    continue
                cut2 = cuts2[j]

                for tet in cut1:
                    self.add_tet(index, conf_index, tet, [p0, p1, p2, p3], [4, 5, 6, 7], edges)

                for tet in cut2:
                    self.add_tet(index, conf_index, tet, [p0, p1, p2, p3], [4, 5, 6, 7], edges)

                self.add_edges(index, edges1[i], edges2[j])
                conf_index += 1


    def create_single_vertex_cut(self, edges):
        index = self.index_for_edges(edges)
        assert(self.table[index] is None)

        tip = -1

        e0 = self.edges[edges[0], :]
        e1 = self.edges[edges[1], :]
        e2 = self.edges[edges[2], :]

        if e0[0] == e1[0] and e0[0] == e2[0]:
            tip = e0[0]
        elif e0[0] == e1[1] and e0[0] == e2[0]:
            tip = e0[0]
        elif e0[0] == e1[0] and e0[0] == e2[1]:
            tip = e0[0]
        elif e0[0] == e1[1] and e0[0] == e2[1]:
            tip = e0[0]
        elif e0[1] == e1[0] and e0[1] == e2[0]:
            tip = e0[1]
        elif e0[1] == e1[1] and e0[1] == e2[0]:
            tip = e0[1]
        elif e0[1] == e1[0] and e0[1] == e2[1]:
            tip = e0[1]
        elif e0[1] == e1[1] and e0[1] == e2[1]:
            tip = e0[1]

        assert(tip >= 0)

        p0 = (self.vertices[e0[0], :]+self.vertices[e0[1], :]) / 2
        p1 = (self.vertices[e1[0], :]+self.vertices[e1[1], :]) / 2
        p2 = (self.vertices[e2[0], :]+self.vertices[e2[1], :]) / 2


        index0 = e0[0] if e0[1] == tip else e0[1]
        index1 = e1[0] if e1[1] == tip else e1[1]
        index2 = e2[0] if e2[1] == tip else e2[1]


        nodes = [index0, index1, index2, 4, 5, 6]
        cuts, edges1, pp = self.prism_cut(nodes, [p0, p1, p2], 7)
        self.new_points[index] = []
        self.new_points[index].append(pp)

        conf_index = 0
        for i in range(len(cuts)):
            cut = cuts[i]
            self.add_tet(index, conf_index, [4, 5, 6, tip], [p0, p1, p2], [4, 5, 6], edges)
            for tet in cut:
                self.add_tet(index, conf_index, tet, [p0, p1, p2], [4, 5, 6], edges)
            self.add_edge(index, edges1[i])
            conf_index += 1


    def create_single_vertex_on_plane_cut(self, edges):
        index = self.index_for_edges(edges)
        assert(self.table[index] is None)

        # print(flags, self.int_to_flags(index))

        e0 = self.edges[edges[0], :]
        e1 = self.edges[edges[1], :]

        tip = -1

        if e0[0] == e1[0] or e0[0] == e1[1]:
            tip = e0[0]
        elif e0[1] == e1[0] or e0[1] == e1[1]:
            tip = e0[1]

        assert(tip >= 0)

        other_tip = -1

        # Maybe flags[other_tip] = 1 (or) flags[3-other_tip] = 1

        tmp = [False]*4
        tmp[e0[0]] = True
        tmp[e0[1]] = True
        tmp[e1[0]] = True
        tmp[e1[1]] = True

        for f in range(4):
            if not tmp[f]:
                other_tip = f

        assert(other_tip >= 0)

        index0 = e0[0] if tip == e0[1] else e0[1]
        index1 = e1[0] if tip == e1[1] else e1[1]

        p0 = (self.vertices[e0[0], :]+self.vertices[e0[1], :]) / 2
        p1 = (self.vertices[e1[0], :]+self.vertices[e1[1], :]) / 2
        self.new_points[index] = []


        self.add_tet(index, 0, [tip, 4, 5, other_tip], [p0, p1], [4, 5, other_tip], edges)
        self.add_tet(index, 0, [4, 5, other_tip, index0], [p0, p1], [4, 5, other_tip], edges)
        self.add_tet(index, 0, [5, index1, other_tip, index0], [p0, p1], [4, 5, other_tip], edges)
        self.add_edge(index, [[index0, 5]])


        self.add_tet(index, 1, [tip, 4, 5, other_tip], [p0, p1], [4, 5, other_tip], edges)
        self.add_tet(index, 1, [4, 5, other_tip, index1], [p0, p1], [4, 5, other_tip], edges)
        self.add_tet(index, 1, [4, index0, other_tip, index1], [p0, p1], [4, 5, other_tip], edges)
        self.add_edge(index, [[index1, 4]])


    def create_two_vertex_on_plane_cut(self, edge):
        index = self.index_for_edges([edge])
        assert(self.table[index] is None)

        # print(flags, self.int_to_flags(index))

        ce = self.edges[edge, :]

        other_edge = []

        for e in self.edges:
            if ce[0] == e[0] or ce[0] == e[1]:
                continue
            if ce[1] == e[0] or ce[1] == e[1]:
                continue

            other_edge = e
            break

        assert(len(other_edge) == 2)

        p = (self.vertices[ce[0], :]+self.vertices[ce[1], :]) / 2
        self.new_points[index] = []

        self.add_tet(index, 0, [4, ce[0], other_edge[0], other_edge[1]], [p], [4, other_edge[0], other_edge[1]], [edge])
        self.add_tet(index, 0, [4, ce[1], other_edge[0], other_edge[1]], [p], [4, other_edge[0], other_edge[1]], [edge])
        self.add_edge(index, [])


    def get_quad(self, e0, e1, i0, i1):
        tmp = set(e0).intersection(set(e1))
        quad = []
        tri = []
        face = []
        if len(tmp) > 0:
            tmp = tmp.pop()

            if tmp == e0[0]:
                quad.append(e0[1])
                face.append(e0[1])
            else:
                quad.append(e0[0])
                face.append(e0[0])

            quad.append(i0)
            quad.append(i1)

            if tmp == e1[0]:
                quad.append(e1[1])
                face.append(e1[1])
            else:
                quad.append(e1[0])
                face.append(e1[0])

            tri = [tmp, i0, i1]
            face.append(tmp)

        face.sort()
        tmp = -1

        if len(face) > 0:
            for i in range(len(self.tet)):
                if self.tet[i][0] == face[0] and self.tet[i][1] == face[1] and self.tet[i][2] == face[2]:
                    tmp = i
                    break

        return quad, tri, tmp


    def add_simple(self, index, conf, edges, tri0, tri1, tri4, tri5, pts):
        edge = set(tri4).intersection(tri5)
        tip4 = set(tri4).difference(edge).pop()
        tip5 = set(tri5).difference(edge).pop()

        ee1 = edge.pop()
        ee2 = edge.pop()

        index4 = -1
        index5 = -1

        for i in range(len(edges)):
            e = self.edges[edges[i]]
            count = 0
            for v in tri4:
                if v == e[0] or v == e[1]:
                    count += 1
            if count == 2:
                index4 = i
                break
        assert(index4 >= 0)


        for i in range(len(edges)):
            e = self.edges[edges[i]]
            count = 0
            for v in tri5:
                if v == e[0] or v == e[1]:
                    count += 1
            if count == 2:
                index5 = i
                break
        assert(index5 >= 0)


        e0 = self.edges[edges[index4], :]
        # e1 = self.edges[edges[1], :]
        e2 = self.edges[edges[index5], :]

        # print(index5, index4, tri4, tri5, tip4, ee2)

        self.add_tet(index, conf, [tri0[0], tri0[1], tri0[2], 7], pts, [4, 5, 6], edges)
        self.add_tet(index, conf, [tri1[0], tri1[1], tri1[2], 7], pts, [4, 5, 6], edges)

        self.add_tet(index, conf, [ee1, ee2, 4+index4, 7], pts, [4, 5, 6], edges)
        if min(e0[0],e0[1]) == min(tip4, ee2) and max(e0[0], e0[1]) == max(ee2, tip4):
            # print([ee1, tip4, 4, 7])
            self.add_tet(index, conf, [ee1, tip4, 4+index4, 7], pts, [4, 5, 6], edges)
        else:
            # print([ee2, tip4, 4, 7])
            self.add_tet(index, conf, [ee2, tip4, 4+index4, 7], pts, [4, 5, 6], edges)

        self.add_tet(index, conf, [ee1, ee2, 4+index5, 7], pts, [4, 5, 6], edges)
        if min(e2[0], e2[1]) == min(tip5, ee2) and max(e2[0], e2[1]) == max(ee2, tip5):
            self.add_tet(index, conf, [ee1, tip5, 4+index5, 7], pts, [4, 5, 6], edges)
        else:
            self.add_tet(index, conf, [ee2, tip5, 4+index5, 7], pts, [4, 5, 6], edges)



    def create_three_edges_cut1(self, edges):
        index = self.index_for_edges(edges)
        assert(self.table[index] is None)

        e0 = self.edges[edges[0], :]
        e1 = self.edges[edges[1], :]
        e2 = self.edges[edges[2], :]

        quad0 = []
        tri0 = []
        quad1 = []
        tri1 = []

        tri4 = []
        tri5 = []

        face0 = -1
        face1 = -1

        quad, tri, face = self.get_quad(e0, e1, 4, 5)
        if len(quad) > 0:
            quad0 = quad
            tri0 = tri
            face0 = face

        quad, tri, face = self.get_quad(e0, e2, 4, 6)
        if len(quad) > 0:
            if len(quad0) > 0:
                assert(len(quad1) <= 0)
                quad1 = quad
                tri1 = tri
                face1 = face
            else:
                quad0 = quad
                tri0 = tri
                face0 = face

        quad, tri, face = self.get_quad(e1, e2, 5, 6)
        if len(quad) > 0:
            assert(len(quad0) > 0)
            assert(len(quad1) <= 0)
            quad1 = quad
            tri1 = tri
            face1 = face


        # print(quad0, quad1)
        # print(face0, face1)

        tmp = [False, False, False, False]
        tmp[face0] = True
        tmp[face1] = True

        for i in range(len(tmp)):
            if not tmp[i]:
                if len(tri4) <= 0:
                    tri4 = self.tet[i].copy()
                    assert(len(tri5) <= 0)
                elif len(tri5) <= 0:
                    assert(len(tri4) > 0)
                    assert(len(tri5) <= 0)
                    tri5 = self.tet[i].copy()

        p0 = (self.vertices[e0[0], :]+self.vertices[e0[1], :]) / 2
        p1 = (self.vertices[e1[0], :]+self.vertices[e1[1], :]) / 2
        p2 = (self.vertices[e2[0], :]+self.vertices[e2[1], :]) / 2


        self.new_points[index] = []


        q00 = list(set(quad0).intersection(tri4).intersection(tri5))[0]
        q10 = list(set(quad1).intersection(tri4).intersection(tri5))[0]

        # fix me tip
        tip = list(set(tri5).intersection(tri4))
        tip = tip[1]

        tri44 = list(tri4)
        tri55 = list(tri5)
        tri44.remove(tip)
        tri55.remove(tip)

        v1 = -1
        v2 = -1
        v3 = -1
        if tip in e0 and (tri44[0] in e0 or tri44[1] in e0):
            v1 = 4
            if tri44[0] in e0:
                v2 = tri44[0]
                v3 = tri44[1]
            else:
                v2 = tri44[1]
                v3 = tri44[0]
        elif tip in e1 and (tri44[0] in e1 or tri44[1] in e1):
            v1 = 5
            if tri44[0] in e1:
                v2 = tri44[0]
                v3 = tri44[1]
            else:
                v2 = tri44[1]
                v3 = tri44[0]
        elif tip in e2 and (tri44[0] in e2 or tri44[1] in e2):
            v1 = 6
            if tri44[0] in e2:
                v2 = tri44[0]
                v3 = tri44[1]
            else:
                v2 = tri44[1]
                v3 = tri44[0]
        else:
            assert(False)

        v5 = -1
        if v3 == tri55[0]:
            v5 = tri55[1]
        elif v3 == tri55[1]:
            v5 = tri55[0]


        v4 = -1

        if v3 in e0 and v5 in e0:
            v4 = 4
        elif v3 in e1 and v5 in e1:
            v4 = 5
        elif v3 in e2 and v5 in e2:
            v4 = 6
        else:
            assert(False)

        v6 = -1

        if v1 + v4 == 9:
            v6 = 6
        elif v1 + v4 == 10:
            v6 = 5
        elif v1 + v4 == 11:
            v6 = 4
        else:
            assert(False)

        self.add_tet(index, 0, [v3, v1, v6, v2], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 0, [v3, tip, v1, v6], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 0, [v6, tip, v4, v5], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 0, [v4, tip, v3, v6], [p0, p1, p2], [4, 5, 6], edges)

        self.add_edge(index, [[tip, v6], [v6, v3]])




        self.add_tet(index, 1, [v3, v1, v6, v2], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 1, [v1, v6, v5, v4], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 1, [v1, v5, tip, v4], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 1, [v3, v1, v6, v4], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 1, [tip, v3, v1, v4], [p0, p1, p2], [4, 5, 6], edges)

        self.add_edge(index, [[v1, v5], [v6, v3]])




        self.add_tet(index, 2, [v1, v6, v5, v4], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 2, [v1, v5, tip, v4], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 2, [v1, v6, v2, v4], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 2, [v2, v3, v4, tip], [p0, p1, p2], [4, 5, 6], edges)

        self.add_edge(index, [[v1, v5], [v4, v2]])



        self.add_tet(index, 3, [v1, v6, v2, v4], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 3, [v2, v3, v4, tip], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 3, [v1, v6, tip, v4], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 3, [tip, v6, v4, v5], [p0, p1, p2], [4, 5, 6], edges)

        self.add_edge(index, [[tip, v6], [v4, v2]])


    def create_three_edges_cut2(self, edges):
        index = self.index_for_edges(edges)
        assert(self.table[index] is None)


        e0 = self.edges[edges[0], :]
        e1 = self.edges[edges[1], :]
        e2 = self.edges[edges[2], :]

        quad0 = []
        tri0 = []
        quad1 = []
        tri1 = []

        tri4 = []
        tri5 = []

        face0 = -1
        face1 = -1

        quad, tri, face = self.get_quad(e0, e1, 4, 5)
        if len(quad) > 0:
            quad0 = quad
            tri0 = tri
            face0 = face

        quad, tri, face = self.get_quad(e0, e2, 4, 6)
        if len(quad) > 0:
            if len(quad0) > 0:
                assert(len(quad1) <= 0)
                quad1 = quad
                tri1 = tri
                face1 = face
            else:
                quad0 = quad
                tri0 = tri
                face0 = face

        quad, tri, face = self.get_quad(e1, e2, 5, 6)
        if len(quad) > 0:
            assert(len(quad0) > 0)
            assert(len(quad1) <= 0)
            quad1 = quad
            tri1 = tri
            face1 = face

        tmp = [False, False, False, False]
        tmp[face0] = True
        tmp[face1] = True

        for i in range(len(tmp)):
            if not tmp[i]:
                if len(tri4) <= 0:
                    tri4 = self.tet[i].copy()
                    assert(len(tri5) <= 0)
                elif len(tri5) <= 0:
                    assert(len(tri4) > 0)
                    assert(len(tri5) <= 0)
                    tri5 = self.tet[i].copy()

        p0 = (self.vertices[e0[0], :]+self.vertices[e0[1], :]) / 2
        p1 = (self.vertices[e1[0], :]+self.vertices[e1[1], :]) / 2
        p2 = (self.vertices[e2[0], :]+self.vertices[e2[1], :]) / 2


        self.new_points[index] = []


        q00 = list(set(quad0).intersection(tri4).intersection(tri5))[0]
        q10 = list(set(quad1).intersection(tri4).intersection(tri5))[0]

        # fix me tip
        tip = list(set(tri5).intersection(tri4))
        tip = tip[0]

        tri44 = list(tri4)
        tri55 = list(tri5)
        tri44.remove(tip)
        tri55.remove(tip)

        v1 = -1
        v2 = -1
        v3 = -1
        if tip in e0 and (tri44[0] in e0 or tri44[1] in e0):
            v1 = 4
            if tri44[0] in e0:
                v2 = tri44[0]
                v3 = tri44[1]
            else:
                v2 = tri44[1]
                v3 = tri44[0]
        elif tip in e1 and (tri44[0] in e1 or tri44[1] in e1):
            v1 = 5
            if tri44[0] in e1:
                v2 = tri44[0]
                v3 = tri44[1]
            else:
                v2 = tri44[1]
                v3 = tri44[0]
        elif tip in e2 and (tri44[0] in e2 or tri44[1] in e2):
            v1 = 6
            if tri44[0] in e2:
                v2 = tri44[0]
                v3 = tri44[1]
            else:
                v2 = tri44[1]
                v3 = tri44[0]
        else:
            assert(False)

        v5 = -1
        if v3 == tri55[0]:
            v5 = tri55[1]
        elif v3 == tri55[1]:
            v5 = tri55[0]


        v4 = -1

        if v3 in e0 and v5 in e0:
            v4 = 4
        elif v3 in e1 and v5 in e1:
            v4 = 5
        elif v3 in e2 and v5 in e2:
            v4 = 6
        else:
            assert(False)

        v6 = -1

        if v1 + v4 == 9:
            v6 = 6
        elif v1 + v4 == 10:
            v6 = 5
        elif v1 + v4 == 11:
            v6 = 4
        else:
            assert(False)

        self.add_tet(index, 0, [v3, v1, v6, v2], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 0, [v3, tip, v1, v6], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 0, [v6, tip, v4, v5], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 0, [v4, tip, v3, v6], [p0, p1, p2], [4, 5, 6], edges)

        self.add_edge(index, [[tip, v6], [v6, v3]])




        self.add_tet(index, 1, [v3, v1, v6, v2], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 1, [v1, v6, v5, v4], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 1, [v1, v5, tip, v4], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 1, [v3, v1, v6, v4], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 1, [tip, v3, v1, v4], [p0, p1, p2], [4, 5, 6], edges)

        self.add_edge(index, [[v1, v5], [v6, v3]])




        self.add_tet(index, 2, [v1, v6, v5, v4], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 2, [v1, v5, tip, v4], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 2, [v1, v6, v2, v4], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 2, [v2, v3, v4, tip], [p0, p1, p2], [4, 5, 6], edges)

        self.add_edge(index, [[v1, v5], [v4, v2]])



        self.add_tet(index, 3, [v1, v6, v2, v4], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 3, [v2, v3, v4, tip], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 3, [v1, v6, tip, v4], [p0, p1, p2], [4, 5, 6], edges)
        self.add_tet(index, 3, [tip, v6, v4, v5], [p0, p1, p2], [4, 5, 6], edges)

        self.add_edge(index, [[tip, v6], [v4, v2]])


    def create_two_edges_cut(self, edges):
        index = self.index_for_edges(edges)
        assert(self.table[index] is None)

        e0 = self.edges[edges[0], :]
        e1 = self.edges[edges[1], :]

        p0 = (self.vertices[e0[0], :]+self.vertices[e0[1], :]) / 2
        p1 = (self.vertices[e1[0], :]+self.vertices[e1[1], :]) / 2



        self.new_points[index] = []
        self.add_tet(index, 0, [4, e0[0], e1[0], 5], [p0, p1], [4, 5], edges)
        self.add_tet(index, 0, [4, e0[1], e1[0], 5], [p0, p1], [4, 5], edges)
        self.add_tet(index, 0, [4, e0[0], e1[1], 5], [p0, p1], [4, 5], edges)
        self.add_tet(index, 0, [4, e0[1], e1[1], 5], [p0, p1], [4, 5], edges)

        self.add_edge(index, [])




if __name__ == '__main__':
    table = CutTable()
    print(table)

    for i in range(len(table.new_points)):
        p = table.new_points[i]
        if p is None:
            continue
        if len(p) > 0:
            print("iiii", i, p)

    for i in range(len(table.table)):
        confs = table.table[i]
        if confs:
            print(i, table.int_to_flags_rec(i), table.int_to_flags_rec(i).count("1"))
    print(".................")
    for i in range(len(table.table)):
        confs = table.table[i]
        if not confs:
            print(i, table.int_to_flags_rec(i), table.int_to_flags_rec(i).count("1"))
