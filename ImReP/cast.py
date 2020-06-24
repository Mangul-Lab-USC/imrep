from numpy import mean
import networkx as nx
from jellyfish import levenshtein_distance as edit_distance
import sys
from ImReP.utils import *

if sys.version_info.major == 2:
    str = unicode


class Cast(object):

    def __init__(self, raw_cdr3_dict):
        """
        :param :
        :type :
        :return:
        :rtype:
        """
        self.graph = nx.Graph()

        for x, y in raw_cdr3_dict.items():
            self.graph.add_node(x, weight=int(y))

        self._dist_dict = {}
        keys = list(raw_cdr3_dict.keys())

        for i in range(len(keys)):
            for j in range(i, len(keys)):
                x, u = keys[i], keys[j]
                if x != u:
                    w = 1.0 / int(edit_distance(str(x), str(u)))
                else:
                    w = 2.0
                self._dist_dict[(x, u)] = w
                self._dist_dict[(u, x)] = w
                if w >= 0.5:
                    self.graph.add_edge(x, u, weight=w)

    def __close(self, nodes, cluster, threshold):
        """
        :param :
        :type :
        :return:
        :rtype:
        """
        dists = []
        for node1 in nodes:
            d = []
            for node2 in cluster:
                d.append(self._dist_dict[(node1, node2)])
            mean_dist = mean(d)
            if mean_dist > threshold:
                dists.append((node1, mean_dist))
        if dists:
            return max(dists, key=lambda z: z[1])[0]
        return None

    def __distant(self, cluster, threshold):
        """
        :param :
        :type :
        :return:
        :rtype:
        """
        dists = []
        for node1 in cluster:
            d = []
            for node2 in cluster:
                d.append(self._dist_dict[(node1, node2)])
            mean_dist = mean(d)
            if mean_dist < threshold:
                dists.append((node1, mean_dist))
        if dists:
            return min(dists, key=lambda z: z[1])[0]
        return None

    def __cast(self, nodes, threshold):
        """
        :param :
        :type :
        :return:
        :rtype:
        """
        partition = []
        while nodes:
            degrees = [(node, self.graph.degree(node))
                       for node in self.graph.nodes()]
            max_deg_vert = max(degrees, key=lambda z: z[1])[0]
            cluster = set([max_deg_vert])

            nodes.remove(max_deg_vert)
            cl = self.__close(nodes, cluster, threshold)
            dist = self.__distant(cluster, threshold)
            while cl or dist:
                if cl:
                    cluster.add(cl)
                    nodes.remove(cl)
                if dist:
                    cluster.remove(dist)
                    nodes.add(dist)
                cl = self.__close(nodes, cluster, threshold)
                dist = self.__distant(cluster, threshold)
            new_cluster = []
            cluster_sorted = sorted(cluster)
            for node in cluster_sorted:
                new_cluster.append((node, self.graph.nodes[node]["weight"]))
            partition.append(new_cluster)
            for vert in cluster_sorted:
                self.graph.remove_node(vert)
        return sorted(partition)

    def doCast(self, threshold):
        """
        :param :
        :type :
        :return:
        :rtype:
        """
        nodes = set(self.graph.nodes())

        partition = self.__cast(nodes, threshold)

        cdr3s = []
        for part in partition:
            representative = max(part, key=lambda z: z[1])[0]
            count_cluster = sum([x[1] for x in part])
            members = list(map(lambda x: x[0], part))
            cdr3s.append([representative, count_cluster, members])
        return cdr3s

    def doCastFromFile(self, inputfile, outputfile ):
        castThreshold = {'IGH': 0.2, 'IGK': 0.2, 'IGL': 0.2, 'TRA': 0.3,
                         'TRB': 0.2, 'TRD': 0.2, 'TRG': 0.2}
        castdict = {}
        typedict = {}
        result = []
        for chain in ["IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"]:
            castdict[chain] = {}

        with open(inputfile) as f:
            content = f.readlines()
        content = map(lambda x: x.strip().split(), content)
        for line in content:
            cdr3, chtype, count = line[:3]
            typedict[cdr3] = line[3:]
            if cdr3 not in castdict[chtype]:
                castdict[chtype][cdr3] = int(count)
            else:
                castdict[chtype][cdr3] += int(count)
        for chtype, cdict in castdict.items():
            cast = Cast(cdict)
            clustered = cast.doCast(castThreshold[chtype])
            clustered = [cclone for cclone in clustered if cclone[1] > 1]
            for x, y, z in clustered:
                line = [x, chtype, y] + typedict[x]
                result.append("%s\n" % ",".join(map(str, line)))
        dumpClones2(result, outputfile)
