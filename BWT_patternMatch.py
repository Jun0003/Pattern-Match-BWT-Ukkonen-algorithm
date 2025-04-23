#%%
import sys
class Node:
    def __init__(self, start: int, end: int) -> None:
        n_askii = 127
        # start and end represents the range between parent node and this node
        self.start = start
        self.end = end
        self.child = [None] * n_askii
        self.suffix_link = None  # if node has suffix link store the node that this node link point to.
        self.suffix_index = -1 # negative means its not leaf
    
class SuffixTree:
    def __init__(self, text: str) -> None:
        self.text = text + "$"
        self.unsolved_n = None
        self.new_internal_node = None
        self.existing_node = None # This is for resolve the new node. Connect it with unsolved node by link
        self.showstopper = False
        self.lastj = 0 # The last extension j that was performed using rule 2
        self.activeN = None # Active node
        self.rem = None # Remainder store start point and end point eg, text = abcde and rem = de then rem = (3,4) [index start from 0]
        self.globalEndP = len(self.text) - 1 # Global end point 
        self.pointer = None

    # base step of Ukkonen
    def phase_one_construction(self):
        self.rootN = Node(-1, -1)
        self.rootN.suffix_link = self.rootN # suffix link point root itself
        self.rootN.suffix_index = -1
        newN = Node(0, self.globalEndP)
        self.rootN.child[ord(self.text[0])] = newN

        self.activeN = self.rootN
        newN.suffix_index = 0


    # return Active node after extension
    def extension(self, j: int, i: int):
        # childN = self.activeN.child[ord(self.text[i])].start

        if self.rem == None:
            if self.activeN.child[ord(self.text[i])] == None:
                # when we read self.rem from active node and stop on the next node
                # we update the self.rem and the active node also apply rule 2 case 1               
                leaf_node = Node(i,self.globalEndP)
                leaf_node.suffix_index = j
                self.existing_node = self.activeN
                self.activeN.child[ord(self.text[i])] = leaf_node
                self.lastj = j

            # rule 3 when rem is empty
            elif self.activeN.child[ord(self.text[i])] != None:
                self.showstopper = True
                self.rem = (i, i)
        
        elif self.rem != None:
            # Rule 3 when rem is not empty

            childN = self.activeN.child[ord(self.text[self.rem[0]])]
            rem_len = self.rem[1] - self.rem[0] + 1
            x = childN.start + rem_len
            if self.text[i] == self.text[x]: # if str[i + 1] and x is same apply rule 3
                self.showstopper = True
                self.rem = (self.rem[0], self.rem[1] + 1)
            if self.text[i] != self.text[x]:
                # phase char is not match with the pos so split the edge?

                # Split the node and edge
                internal_node = Node(self.rem[0], self.rem[1])
                new_leaf_node = Node(i,i)
                new_leaf_node.suffix_index = j
                childN.start = x
                existed_leaf_node = childN

                # connect internal node with active node 
                self.activeN.child[ord(self.text[self.rem[0]])] = internal_node

                # splits here
                internal_node.child[ord(self.text[i])] = new_leaf_node
                internal_node.child[ord(self.text[x])] = existed_leaf_node
                self.new_internal_node = internal_node # save as unsolved node, this node will be solved(ie linked) at next iteration
                
                self.lastj += 1


    def resolveSuffixLink(self):
        
        if self.new_internal_node != None and self.unsolved_n == None: # at first time new internal node is created
            self.unsolved_n = self.new_internal_node # treat new internal node as unsolved node
        else:
            if self.new_internal_node != None: # When we have unsolved node and we got new internal node 
                self.unsolved_n.suffix_link = self.new_internal_node # unsolved is connected with new internal node
                self.unsolved_n = self.new_internal_node # now new internal node will be the unsolved node and this will solved in j + 1
                
            # when we have unsolved node and we need to connect it with existing node with one branch below extending via character x  
            elif  self.unsolved_n != None and self.new_internal_node == None and self.existing_node != None: 
                self.unsolved_n.suffix_link = self.existing_node
                self.unsolved_n = None # now we don't have pending new internal node that has to be connected with link
        self.new_internal_node = None

    def skipCountDown(self):
        if self.rem != None:
            len_rem = self.rem[1] - self.rem[0] + 1 
            total_lenght_edge = 0
            skip_index = 0

            while True:
                rem_ith_char = self.text[self.rem[skip_index]]
                childNode = self.activeN.child[ord(rem_ith_char)]
                if childNode == None:
                    break
                new_edge_length = childNode.end - childNode.start + 1
                if total_lenght_edge + new_edge_length > len_rem:
                    break
                else:
                    total_lenght_edge += new_edge_length
                    self.activeN = childNode
                    self.rem = (self.rem[0] + self.activeN.end, self.rem[1])
                    if total_lenght_edge == len_rem:
                        self.rem = None
                        break

                skip_index += new_edge_length

    # if move true else false(show stop)
    def moveToNextExtension(self) -> bool:
        if self.showstopper == True:
            return False
        else:
            # Traverse the suffix link
            # when active node is root we need to substract our rem
            if self.activeN.suffix_link == self.rootN and self.activeN == self.rootN:
                # decrease the rem when link is out of the root
                if self.rem == None or self.rem[0] + 1 > self.rem[1]:  # if self.rem start is smaller than end index set as none ie does not exist
                    self.rem = None
                else:

                    self.rem = (self.rem[0] + 1, self.rem[1]) 
            else:
                self.activeN = self.activeN.suffix_link
            return True


    def ukkonen(self):
        n = len(self.text) 
        self.phase_one_construction()
        
        # i represents current phase 
        for i in range(1, n): # Iteration for Phase 
            # j represents current extention
            for j in range(self.lastj + 1, i + 1):
                # execute skip count down, this method only execute when length of rem is longer than the active node edge
                self.skipCountDown()

                # extension method here (Rule 1 or Rule 2)
                self.extension(j, i)

                # resolved unresolved node
                self.resolveSuffixLink()

                # Move to next extension or skip the phase 
                if self.moveToNextExtension() == False:
                    self.showstopper = False
                    break
                
        return self.rootN

    # Traveesr and return suuffix array
    def traverse_dfs(self):
        self.paths = []
        self.dfs_rec(self.rootN, "", self.paths)
        return(self.paths)

    def dfs_rec(self, cuur_node: Node, path, paths):
        path  = path + self.text[cuur_node.start: cuur_node.end + 1]  
        if cuur_node.suffix_index >= 0:
            paths.append( cuur_node.suffix_index)
        else:

            for next_nod in cuur_node.child:
                if next_nod != None:
                    self.dfs_rec(next_nod, path, paths)

# 1, rank table
def build_rank_table(bwt: str):
    n_askii = 127
    char_list = [(None, 0)] * n_askii
    rank_table = [None] * n_askii # The position where char apper in text

    # count the number occurence of each unique cahracter
    for char in bwt:
        char_list[ord(char)] = (char, char_list[ord(char)][1] + 1)
    
    rank = 0
    for charTuple in char_list:
        char = charTuple[0]
        nChar = charTuple[1]
        if nChar > 0:
            rank_table[ord(char)] = rank
        rank += nChar
    return rank_table

# 2, number of occurence table
# number of times x(char) appears in L[1...i] (Note, this is inclusive [] Not [) sice it is easy to get exclusive occurence by using inclusive occurence)
def build_nOcurrences_table(bwt: str):
    i = 0
    n_askii = 127
    occurenceTable = [None] * len(bwt)
    occurence_at_ith = [0] * n_askii

    for i in range(len(bwt)):
        char_at_ith = bwt[i]
        occurence_at_ith[ord(char_at_ith)] = occurence_at_ith[ord(char_at_ith)] + 1
        occurenceTable[i] = occurence_at_ith 

        occurence_at_ith = copy_askii(occurence_at_ith)


    return occurenceTable

def copy_askii(askii_lst: list):
    n_askii = 127
    copy_askii_list = [None] * n_askii

    i = 0
    for elem in askii_lst:
        copy_askii_list[i] = elem
        i += 1
    
    return copy_askii_list

## matching algorithm
def wild_card_pat_match_bwt(bwt: str, pat: str):
    rank_table = build_rank_table(bwt)
    n_occurence_matrix = build_nOcurrences_table(bwt)

    # from lecture slide... sp and ep means
    # sp = rank(pat[i]) + nOccurrences(pat[i], L[1...sp))
    # ep = rank(pat[i]) + nOccurrences(pat[i], L[1...ep]) - 1
    sp = 0
    ep = len(bwt) - 1
    # matching from back
    m = len(pat) - 1 # pointer of pat
    matching_chr = pat[m]
    result = recursive_match(m, matching_chr , bwt, pat, sp, ep, rank_table, n_occurence_matrix)

    return result

# when wild card "!" is matched, we need to branch the matching for all matched char
def recursive_match(m: int, matching_chr: chr , bwt: str, pat: str, sp: int, ep: int, rank_table: list, n_occurence_matrix: list):
    result_match_pos = []
    while m >= 0 and sp <= ep:
        if matching_chr == "!": # recursive func
            askii_id = 0
            for rank in rank_table: # substitute all char we saw in bwt into "!"
                if rank != None:
                    matching_chr = chr(askii_id)
                    if matching_chr != "$":
                        result = recursive_match(m, matching_chr, bwt, pat,sp, ep, rank_table,n_occurence_matrix)
                        if result is not None:
                            for r in result:
                                result_match_pos.append(r)
                askii_id += 1
            # return the result that brought from bottom of the recursion by returning here it prevent to go top function anymore
            return result_match_pos  # Return the collected result 
        else:
            if rank_table[ord(matching_chr)] == None:
                sp, ep = 1, 0
            else:
                rank_of_char = rank_table[ord(matching_chr)] # get rank of that char
                occurence = n_occurence_matrix[ep][ord(matching_chr)]

                # calc ep
                occurence = n_occurence_matrix[ep][ord(matching_chr)]
                ep = rank_of_char + occurence - 1

                # calc sp 
                occurence = n_occurence_matrix[sp][ord(matching_chr)]
                if bwt[sp] == matching_chr: 
                    occurence = occurence - 1
                sp = rank_of_char + occurence
        m = m - 1
        matching_chr = pat[m]

    if sp <= ep and m==-1:
        result_match_pos.append((sp, ep))
        return result_match_pos
    else: 
        result_match_pos



# I could omit this function since we can retrive what char we have in bwt string by inspecting rank table
def collect_unique_char(bwt: str):
    n_askii = 127
    unique_char_lst = [None] * n_askii
    for char in bwt:
        if unique_char_lst[ord(char)] == None:
            unique_char_lst





# build BWT function
def gen_bwt(original_str: str) -> str:
    bwt = ""

    # construct suffix tree and get array by using dfs traverse
    tree = SuffixTree(original_str)
    tree.ukkonen()
    suffix_array = tree.traverse_dfs()


    # build bwt 
    original_str = original_str + "$"
    for i in suffix_array:
        if i == 0:
            bwt = bwt + original_str[len(original_str)-1]
        else:
            bwt = bwt + original_str[i - 1]
    
    return bwt, suffix_array
            

def a2q1_function(text: str, pat: str):
    bwt, suffix_array = gen_bwt(text)
    r = wild_card_pat_match_bwt(bwt, pat) # position 


    match_pos = []
    # match the position with suffix array
    if r != None:
        for range_of_match in r:
            s, e = range_of_match[0], range_of_match[1]
            for i in range(s ,e + 1):
                match_pos.append(suffix_array[i])
        return match_pos


def read_file(file_path: str) -> str:
    f = open(file_path, 'r')
    line = f.readlines()
    f.close()
    return line

if __name__ == '__main__':
    # Retrieve the file paths from the command line arguments
    text = read_file(sys.argv[1])[0]
    pat = read_file(sys.argv[2])[0]

    result = a2q1_function(text, pat)

    with open("output a2q1.txt", 'w') as output_file:
        for i, match in enumerate(result):
            if i == len(result) - 1:  
                output_file.write(f"{match + 1}")  # Write without newline
            else:
                output_file.write(f"{match + 1}\n")  # Write with newline

