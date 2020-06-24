
class SuffixTree(object):

    def __init__(self, string):
        self.root = None
        self.firstChild = None
        self.edgeLabel = None
        self.next = None
        self.index = None



class IgorSuffixTree(object):

    def __init__(self, string):
        """
        :param :
        :type :
        :return:
        :rtype:
        """
        self.__stree = SuffixTree(string)

    def search_stree(self, string):
        """
        :param :
        :type :
        :return:
        :rtype:
        """
        current_index = -1
        terminal = False
        string = list(string)
        str_index = 0
        current = self.__stree.root
        while True:
            current_child = current.firstChild
            found = False
            while True:
                if current_child is None:
                    break
                edgeLabel = current_child.edgeLabel
                if string[str_index] == edgeLabel[0]:
                    found = True
                    break
                else:
                    current_child = current_child.next
            if found:
                len_matched = 0
                current_index = current_child.index
                while (len_matched < len(edgeLabel)
                       and str_index < len(string)
                       and string[str_index] == edgeLabel[len_matched]):
                    str_index += 1
                    len_matched += 1
                if len_matched < len(edgeLabel):
                    if edgeLabel[len_matched] in ["|", "$"]:
                        terminal = True
                    break
                else:
                    if str_index == len(string):
                        terminal = True
                        break
                    current = current_child
            else:
                break
        return str_index, current_index, terminal
