### SageMath implementation of Parking Functions!
import copy
import sys
import math

                                ######################
                                #   CLASS PARKFUNC   #
                                ######################

class ParkFunc():

    def __init__(self, n, m, func = [], w_area = [], w_label = [], specific_opt=[], check=False):
        r'''
            Initialization of the class ParkFunc. The arguments are:
            -   n : number of rows of the parking function
            -   m : number of columns of the parking function
            -   func : a list with the images of values 1,...,n
            -   w_area : a list with the direct path information
            -   w_label : a list with the direct label information
            -   specific_opt : specific options for the parking function
            -   check : if True, it check whether the parameters satisfy the condition for being a parking function

            When called, it initialize the following data:
            -   N : number of rows
            -   M : number of columns
            -   path : the path shape of the corresponding labelled Dyck path (counting all cells of the NxM grid to the east)
            -   label : the labels read from south to north
        '''
        self.N = n                          # Number of rows
        self.M = m                          # Number of columns
        self.ratio = n/m                    # Ratio for the diagonal slope

        # Recognize the specific parking function type
        if specific_opt == []:
            if m%n == 0:                # K TYPE
                specific_opt = ['k-TYPE', floor(m/n)]   # type: ignore
            else:
                specific_opt = ['rational']
        self.specific_opt = specific_opt    # Specific options for the parking function

        # Compute the path and label information
        if func != []:
            self.label = sorted([i+1 for i in range(self.N)], key = lambda k: k + m*func[k-1])
            self.path = [m - func[i-1] + 1 for i in self.label]
        else:
            self.label = w_label
            self.path = w_area

        if check:           # Check if the parking condition is satisfied
            for i in range(n):
                if w_area[i] < ceil(m - self.ratio*i):          # type: ignore
                    raise ValueError("The path does not satisfy parking function condition...")
    
    def __main_diag__(self):
        r'''
            Returns the main diagonal shape in the right format.
        '''
        return [ceil(self.M - i/self.ratio) for i in range(self.N)]  # type: ignore

    def __eq__(self,other):
        r'''
            Redefines equality between parking functions not to look at metadata.
        '''
        if self.path == other.path and self.label == other.label and self.N == other.N and self.M == other.M:
            return True
        else:
            return False

    def label_word(self):
        return self.label

    def area_word(self):
        diag = self.__main_diag__()
        return [self.path[i] - diag[i] for i in range(self.N)]

    def area(self):
        r'''
            Compute the area statistic.
        '''
        diag = self.__main_diag__()
        return sum([self.path[i] - diag[i] for i in range(self.N)])
    
    def pdinv(self):                        # PDINV: computes pathdinv statistic
        r'''
            Compute the path diagonal inversion statistic.
        '''
        if self.specific_opt[0] != 'k-TYPE':
            raise TypeError('The algorithm does not work with this shape...')
        k = self.specific_opt[1]
        dinv = 0
        # Read all cells above path
        for i in range(self.N):                  # Loop on rows
            for j in range(self.N*k-self.path[i]):       # Loop on columns
                arm = (self.N*k - self.path[i]) - (j+1)
                leg = i - max([h for h in range(self.N) if self.N*k-self.path[h]-1 < j]) - 1
                if (arm <= k*(leg+1)) and (k*leg < arm + 1):
                    dinv += 1
        return dinv

    def tdinv(self):                        # TDINV: computes tdinv statistic
        if self.specific_opt[0] != 'k-TYPE':
            raise TypeError('The algorithm does not work with this shape...')
        k = self.specific_opt[1]

        tdinv = 0
        row_cresc = [self.label.index(h+1) for h in range(self.N)]
        for i in range(self.N):
            rank_c = self.N*self.path[row_cresc[i]] + (self.N*k+1)*row_cresc[i]
            for j in range(self.N-i-1):
                rank_d = self.N*self.path[row_cresc[i+j+1]] + (self.N*k+1)*row_cresc[i+j+1]
                if rank_c < rank_d and rank_d < rank_c + (self.N*k):
                    tdinv += 1
        return tdinv

    def w_maxtdinv(self):                   # W_MAXTDINV: computes the labelling giving the max tdinv
        if self.specific_opt[0] != 'k-TYPE':
            raise TypeError('The algorithm does not work with this shape...')
        k = self.specific_opt[1]
        
        rank_list = [(self.N*k+1)*i + self.N*self.path[i] for i in range(self.N)]
        order = sorted(range(len(rank_list)), key=lambda k: rank_list[k])
        park = sorted(range(len(order)), key=lambda k: order[k])
        return [i+1 for i in park]

    def maxtdinv(self):                     # MAXTDINV: computes maxtdinv statistic
        max_pf = ParkFunc(self.N, self.M, w_area=self.path, w_label=self.w_maxtdinv(), specific_opt=self.specific_opt)
        return max_pf.tdinv()

    def dinv(self):                         # DINV: computes the dinv of a parking function
        #print("PDINV: {}\t TDINV: {}\t MAXTDINV: {}".format(self.pdinv(),self.tdinv(),self.maxtdinv()))
        return self.pdinv() + self.tdinv() - self.maxtdinv()

    def pmaj(self, infos=False):            # PMAJ: computes the pmak of a parking function
        if self.specific_opt[0] != 'k-TYPE':
            #raise Warning('The algorithm does not work with this shape...')
            return self.test_pmaj()
        k = self.specific_opt[1]
        
        pmaj = 0
        plus = 0
        current = self.N+1
        reading_word = []
        contributes = {}
        buffer = []
        for i in range(self.N*k):
            # Adding to buffer new labels with k multeplicity
            new = [self.label[j] for j in range(self.N) if self.path[j] == self.N*k-i]
            buffer = buffer + [new[floor(j/k)] for j in range(len(new)*k)]  # type: ignore
            # Get next step
            if current <= min(buffer):          # Starting new ascending chain
                plus += 1
                current = max(buffer)
            else:                               # Continuing current ascending chain
                current = max([w for w in buffer if w < current])
            
            if current not in reading_word:     # Adding pmaj if needed...
                pmaj += plus
                contributes = contributes | {current:plus}
            buffer.remove(current)                 # Removing one copy of current from buffer

            # writing down the new letter in reading word
            reading_word = reading_word + [current]
        if infos:
            return [pmaj,contributes,reading_word]
        else:
            return pmaj

    def test_pmaj(self, infos=False):
        h = math.gcd(self.N, self.M)/self.M
        k = math.gcd(self.N, self.M)/self.N
        new_path = [h*self.path[i] for i in range(self.N)]
        new_label = self.label

        pmaj = 0
        plus = 0
        current = self.N+1
        reading_word = []
        contributes = {}
        buffer = []
        for i in range(self.N*k):
            # Adding to buffer new labels with k multeplicity
            new = [new_label[j] for j in range(self.N) if new_path[j] == self.N*k-i]
            buffer = buffer + [new[floor(j/k)] for j in range(len(new)*k)]  # type: ignore
            # Get next step
            if current <= min(buffer):          # Starting new ascending chain
                plus += 1
                current = max(buffer)
            else:                               # Continuing current ascending chain
                current = max([w for w in buffer if w < current])
            
            if current not in reading_word:     # Adding pmaj if needed...
                pmaj += plus
                contributes = contributes | {current:plus}

            for i in range(h):
                if current in buffer:
                    buffer.remove(current)                 # Removing h copies of current from buffer

            # Writing down the new letter in reading word
            reading_word = reading_word + [current]
        if infos:
            return [pmaj,contributes,reading_word]
        else:
            return pmaj

    def path_word(self):
        if self.specific_opt[0] != 'k-TYPE':
            raise TypeError('The algorithm does not work with this shape...')
        k = self.specific_opt[1]
        
        rank_list = [(self.N*k+1)*i + self.N*self.path[i] for i in range(self.N)]
        path_order = sorted(range(len(rank_list)), key=lambda k: rank_list[k])      # Reading order of the columns
        path_word = [self.label[path_order[i]] for i in range(self.N)]
        return path_word

    def to_area_pmaj(self,infos = False):   # TO_AREA_PMAJ: this is the Generalized Loehr-Remmel map
        r'''
            This function returns the image of the parking function through the generalized Loehr-Remmel map.
        '''
        if self.specific_opt[0] != 'k-TYPE':
            raise TypeError('The algorithm does not work with this shape...')
        k = self.specific_opt[1]
        
        ### Create the path word, ordered by rank
        rank_list = [(self.N*k+1)*i + self.N*self.path[i] for i in range(self.N)]
        path_order = sorted(range(len(rank_list)), key=lambda k: rank_list[k])      # Reading order of the columns
        path_word = [self.label[path_order[i]] for i in range(self.N)]                      # Labels in order
        #print("The path_word: {}".format(path_word))

        ### Compute the not_tdinv
        # Compute all pairs of cars making tdinv
        tdinv_pairs = []            
        row_cresc = [self.label.index(h+1) for h in range(self.N)]
        for i in range(self.N):
            rank_c = self.N*self.path[row_cresc[i]] + (self.N*k+1)*row_cresc[i]
            for j in range(self.N-i-1):
                rank_d = self.N*self.path[row_cresc[i+j+1]] + (self.N*k+1)*row_cresc[i+j+1]
                if rank_c < rank_d and rank_d < rank_c + (self.N*k+1):
                    tdinv_pairs = tdinv_pairs + [(i+1,i+j+2)]
        if infos:
            print("INFORMATION on DINV SET")
            print("\t The pairs of tdinv are: {}".format(tdinv_pairs))

        # Compute not_tdinv with previous cars
        not_tdinv = []
        for i in range(self.N):
            temp = 0
            for j in range(i):
                if ((path_word[i],path_word[j]) not in tdinv_pairs) and ((path_word[j],path_word[i]) not in tdinv_pairs):
                    temp += 1
            not_tdinv = not_tdinv + [temp]
        #print("The not_tdinv: {}".format(not_tdinv))

        ### Compute the not_dinvcorr
        # Compute all pairs of cars making dinvcorr
        dinvcorr_pairs = []
        for i in range(self.N):
            rank_c = self.N*self.path[i] + (self.N*k+1)*i
            for j in range(self.N-i-1):
                rank_d = self.N*self.path[i+j+1] + (self.N*k+1)*(i+j+1)
                for shift in range(k-1):
                    if rank_d < rank_c + self.N*(k-1) - self.N*shift and rank_c - self.N - self.N*shift < rank_d:
                        dinvcorr_pairs = dinvcorr_pairs + [(self.label[i],self.label[i+j+1])]
        if infos:
            dinvcorr_pairs.sort()
            print("\t The pairs of dinvcorr are: {}".format(dinvcorr_pairs))

        # Compute not_dinvcorr with previous cars
        not_dinvcorr = []
        for i in range(self.N):
            temp = 0
            for j in range(i):
                temp += k - 1 - len([1 for (a,b) in dinvcorr_pairs if (a,b)==(path_word[i],path_word[j]) or (a,b)==(path_word[j],path_word[i])])
            not_dinvcorr = not_dinvcorr + [temp]

        ### Compute the not_dinv
        not_dinv = [not_tdinv[i] + not_dinvcorr[i] for i in range(self.N)]
        #print("The not_dinv: {}".format(not_dinv))

        ### Compute the w_path of the image
        not_dinv_sort = sorted(not_dinv, key=lambda k: k)
        area_word = [self.N*k-not_dinv_sort[i] for i in range(self.N)]
        
        ### Compute the w_label of the image
        label_word = [path_word[i] for i in sorted(range(self.N), key=lambda k: not_dinv[k]*self.N + path_word[k])]
        if not infos:
            return ParkFunc(self.N, self.M, w_area=area_word, w_label=label_word, specific_opt=self.specific_opt)
        else:
            return [ParkFunc(self.N, self.M, w_area=area_word, w_label=label_word, specific_opt=self.specific_opt), not_tdinv, not_dinvcorr]

    def to_dinv_area(self,infos = False):   # TO_DINV_AREA: this is the Generalized Loehr-Remmel inverse
        r'''
            This function returns the image through the generalized Loehr-Remmel inverse.
        '''
        if self.specific_opt[0] != 'k-TYPE':
            raise TypeError('The algorithm does not work with this shape...')
        k = self.specific_opt[1]

        #self.draw()

        ### Compute the pmaj contributes of each label
        [pmaj, contributes, reading_word] = self.pmaj(infos = True)
        lab_ord = [v[0] for v in contributes]       # Ordered list of the inserted labels
        area = [v[1] for v in contributes]          # Ordered list of final area contributes of each label
        #print("Label order {}".format(lab_ord))
        #print("Area {}".format(area))
        ### Compute the area contributes of each label
        #print("Path {}".format(self.path))
        codinv_dict = {self.label[i]:(self.M - self.path[i]) for i in range(self.N)}
        dinv = [k*i - codinv_dict[lab_ord[i]] for i in range(self.N)]
        #print("Dinv {}".format(dinv))

        ### Build up the path step by step
        # First step
        m = 1
        label = [lab_ord[0]]
        path = [k]
        # Consecutive steps
        for j in range(self.N-1):
            if infos:
                print("\nDIMENSION {} PARTIAL PATH {} PARTIAL LABEL {}".format(j,path,label))
            j += 1
            check = False
            #print([(m-(pos+1))*k + area[j] for pos in range(m-1)])
            #tries = [pos+1 for pos in range(m-1) if ((m+1)*k >= ((m+1)-(pos+1))*k + area[j]) and (path[pos]+k >= ((m+1)-(pos+1))*k + area[j]) and (path[pos+1] <= ((m+1)-(pos+1))*k + area[j])] + [m]
            tries = [pos+1 for pos in range(m-1) if ((pos+1)*k >= area[j]) and (path[pos]+k >= (m-pos)*k + area[j]) and (path[pos+1] <= (m-pos)*k + area[j])] + [m]
            if area[j] == 0:
                tries = [0] + tries
            ind = 0
            if infos:
                print("POSITIONS TO TRY {}".format(tries))
            while ind < len(tries) and check == False:
                temp_path1 = [path[i]+k for i in range(tries[ind])]                     # Before new step
                temp_path3 = [path[tries[ind]+i] for i in range(m-tries[ind])]        # After new step
                #print([tries[ind]+i for i in range(m-tries[ind])])        # After new step
                temp_path = temp_path1 + [(m+1-tries[ind])*k + area[j]] + temp_path3    # Provisional path
                
                temp_label1 = [label[i] for i in range(tries[ind])]                     # Before new label
                temp_label3 = [label[tries[ind]+i] for i in range(m-tries[ind])]      # After new label
                temp_label = temp_label1 + [lab_ord[j]] + temp_label3                   # Provisional label
                order = sorted(range(m+1), key = lambda i: temp_label[i])
                norm_label = sorted(range(m+1), key = lambda i: order[i])
                norm_label = [v+1 for v in norm_label]                                  # Normalized label
                #print("Temp {} \t Norm {}".format(temp_label,norm_label))
                #if infos:
                #    print("First part {}".format(temp_path1))
                #    print("Second part {}".format((m+1-tries[ind])*k + area[j]))
                #    print("Third part {}".format(temp_path3))
                #print("Trying label {}\t and path {}".format(norm_label,temp_path))
                temp = ParkFunc(m+1,k*(m+1), w_area=temp_path,w_label=norm_label)
                dinv_val = temp.dinv()
                if sum([dinv[i] for i in range(m+1)]) == dinv_val:
                    #print("\tGOOD: We have dinv {}\t with dinv_val {}\t and sum {}".format(dinv, dinv_val,sum([dinv[i] for i in range(m+1)])))
                    label = temp_label
                    path = temp_path
                    check = True
                    if infos:
                        print("Labels: {}".format(temp_label))
                        temp.draw()
                #else:
                    #print("\tBAD: We have dinv {}\t with dinv_val {}\t and sum {}".format(dinv, dinv_val,sum([dinv[i] for i in range(m+1)])))
                ind += 1
            
            if check == False:
                print("First part {}".format(temp_path1))
                print("Second part {}".format((m+1-tries[ind-1])*k + area[j]))
                print("Third part {}".format(temp_path3))
                raise ValueError("PROBLEMO!")
            m+=1
        
        return ParkFunc(self.N, self.M, w_area=path, w_label=label)


                                ##########################
                                #   PARKFUNC UTILITIES   #
                                ##########################

    def draw(self, stats=True):
        r'''
            Function that draws the path
        '''
        #print('Path {} and {}'.format(self.path,self.M))
        diag = self.__main_diag__()
        for i in range(self.N):
            i = self.N - i - 1
            row1 = '   ' * (self.M - self.path[i])
            row2 = ' {:>2}'.format(self.label[i])
            row3 = '|##' * (self.path[i] - diag[i])
            row4 = '|  ' * (diag[i]) + '|'

            buff1 = '   ' * (self.M - self.path[i]+1)
            buff2 = '+--' * (self.path[i]) + '+'
            print('\t' + buff1 + buff2)
            print('\t' + row1 + row2 + row3 + row4)
        print('\t' + '   ' + '+--' * self.M + '+')
        if stats:
            print(" Label word:\t{}".format(self.label_word()))
            print(" Area word:\t{}".format(self.area_word()))
            print(" Area: {}\t Dinv: {}\t Pmaj: {}".format(self.area(),self.dinv(),self.pmaj()))

    def labels_rises(self):                 # LABELS_RISES
        r'''
            Returns the labels corresponding to possible rises.
        '''
        return [self.label[i+1] for i in range(self.N-1) if self.path[i] == self.path[i+1]]

    def labels_valleys(self):               # LABELS_VALLEYS
        r'''
            Returns the labels corresponding to possible rises.
        '''
        return [self.label[i+1] for i in range(self.N-1) if (self.path[i]+1 == self.path[i+1] and self.label[i] < self.label[i+1]) or (self.path[i]+1 > self.path[i+1])]

def kDyckPaths(n,k):                            # kDYCKPATHS: computes the set of all (nk,n)-Dyck paths
    possible_paths = [[1]]
    for i in range(n*(k+1)-1):
        new_possible = []
        for part_path in possible_paths:
            height = sum(part_path)
            distance = sum([1-part_path[j] for j in range(i+1)])
            if height == n:                 # Reached max height, just 0's
                temp = copy.deepcopy(part_path)
                temp.append(0)
                new_possible.append(temp)
            elif k*height < distance + 1:   # Reached the "diagonal", just a 1
                temp = copy.deepcopy(part_path)
                temp.append(1)
                new_possible.append(temp)
            else:
                temp = copy.deepcopy(part_path)
                temp.append(0)
                new_possible.append(temp)
                temp = copy.deepcopy(part_path)
                temp.append(1)
                new_possible.append(temp)
        possible_paths = new_possible
    dyck_k_paths = [[sum([1-w[i+j] for j in range(len(w)-i)]) for i in range(len(w)) if w[i]==1] for w in possible_paths]
    return dyck_k_paths

def nkParkingFunctions(n,k, display=False):     # kPARKINGFUNCTIONS: computes the set of all (n,k)-parking functions
    parking_functions = []
    counter = 0
    # Compute total number of Dyck paths
    totalnum = binomial(n*(k+1)+1, n)/(n*(k+1)+1)           # type: ignore
    print("Starting the algorithm to compute all parking functions...")
    # Compute all possible set partitions
    set_partitions = {}
    for lamb in Partitions(n):                              # type: ignore
        temp = []
        for possible in OrderedSetPartitions(n,lamb):       # type: ignore
            poss = [list(part) for part in possible]
            temp.append(poss)
        lambt = tuple(lamb)
        set_partitions.update({lambt:temp})

    ### Loop on k-Dyck paths
    for path in kDyckPaths(n,k):
        if display:
            counter += 1
            perc = floor(100*counter/totalnum)              # type: ignore
            sys.stdout.write('\r')                          # Reset to start of line
            sys.stdout.write("Percentage %3d %%, computing parking functions for path %6d of %6d" % (perc, counter, totalnum))
            sys.stdout.flush()

        # Compute the partitions we want
        rec_partition = [len([1 for w in path if w==n*k-i]) for i in range(n*k) if (n*k-i) in path]
        #print(rec_partition)
        ord_partition = copy.copy(rec_partition)
        ord_partition.sort()
        ord_partition.reverse()
        #print(ord_partition)
        ord_partitiont = tuple(ord_partition)
        #print(ord_partitiont)
        part_decr_ord = sorted(range(len(rec_partition)), key=lambda i: -n*rec_partition[i])
        read_order = sorted(range(len(part_decr_ord)), key=lambda i: part_decr_ord[i])
        for set_part in set_partitions[ord_partitiont]:
            label = []
            for j in read_order:
                temp2 = list(list(set_part)[j])
                temp2.sort()
                label = label + temp2
            parking_functions.append(ParkFunc(n,n*k,w_area=path,w_label=label))
    
    print("\n...the algorithm to compute all parking functions has finished.")
    return parking_functions

def qt_Polynomial_dinv_area(set_paths,ring):           # GEN_POLY: polynomial generated by paths
    q = ring.gens()[0]
    t = ring.gens()[1]
    poly = 0*q*t
    for path in set_paths:
        poly = poly + (q**path.dinv())*(t**path.area())
    return poly

def qt_Polynomial_area_pmaj(set_paths,ring):           # GEN_POLY: polynomial generated by paths
    q = ring.gens()[0]
    t = ring.gens()[1]
    poly = 0*q*t
    for path in set_paths:
        poly = poly + (q**path.area())*(t**path.pmaj())
    return poly

def to_area_pmaj_RV(self,infos = False):   # TO_AREA_PMAJ: this is the Generalized Loehr-Remmel map
    r'''
        This function returns the image of the parking function through the generalized Loehr-Remmel map.
    '''
    if self.specific_opt[0] != 'k-TYPE':
        raise TypeError('The algorithm does not work with this shape...')
    k = self.specific_opt[1]
    
    ### Create the path word, ordered by rank
    rank_list = [(self.N*k+1)*i + self.N*self.path[i] for i in range(self.N)]
    path_order = sorted(range(len(rank_list)), key=lambda k: rank_list[k])      # Reading order of the columns
    path_word = [self.label[path_order[i]] for i in range(self.N)]                      # Labels in order
    #print("The path_word: {}".format(path_word))

    ### Compute the not_tdinv
    # Compute all pairs of cars making tdinv
    tdinv_pairs = []            
    row_cresc = [self.label.index(h+1) for h in range(self.N)]
    for i in range(self.N):
        rank_c = self.N*self.path[row_cresc[i]] + (self.N*k+1)*row_cresc[i]
        for j in range(self.N-i-1):
            rank_d = self.N*self.path[row_cresc[i+j+1]] + (self.N*k+1)*row_cresc[i+j+1]
            if rank_c < rank_d and rank_d < rank_c + (self.N*k+1):
                tdinv_pairs = tdinv_pairs + [(i+1,i+j+2)]

    # Compute not_tdinv with previous cars
    not_tdinv = []
    for i in range(self.N):
        temp = 0
        for j in range(i):
            if ((path_word[i],path_word[j]) not in tdinv_pairs) and ((path_word[j],path_word[i]) not in tdinv_pairs):
                temp += 1
        not_tdinv = not_tdinv + [temp]
    #print("The not_tdinv: {}".format(not_tdinv))

    ### Compute the not_dinvcorr
    # Compute all pairs of cars making dinvcorr
    dinvcorr_pairs = []
    for i in range(self.N):
        rank_c = self.N*self.path[i] + (self.N*k+1)*i
        for j in range(self.N-i-1):
            rank_d = self.N*self.path[i+j+1] + (self.N*k+1)*(i+j+1)
            for shift in range(k-1):
                if rank_d < rank_c + self.N*(k-1) - self.N*shift and rank_c - self.N - self.N*shift < rank_d:
                    dinvcorr_pairs = dinvcorr_pairs + [(self.label[i],self.label[i+j+1])]
    # Compute not_dinvcorr with previous cars
    not_dinvcorr = []
    for i in range(self.N):
        temp = 0
        for j in range(i):
            temp += k - 1 - len([1 for (a,b) in dinvcorr_pairs if (a,b)==(path_word[i],path_word[j]) or (a,b)==(path_word[j],path_word[i])])
        not_dinvcorr = not_dinvcorr + [temp]

    ### Compute the not_dinv
    not_dinv = [not_tdinv[i] + not_dinvcorr[i] for i in range(self.N)]
    #print("The not_dinv: {}".format(not_dinv))

    ### Compute the w_path of the image
    not_dinv_sort = sorted(not_dinv, key=lambda k: k)
    area_word = [self.N*k-not_dinv_sort[i] for i in range(self.N)]
    
    ### Compute the w_label of the image
    label_word = [path_word[i] for i in sorted(range(self.N), key=lambda k: not_dinv[k]*self.N + path_word[k])]
    if not infos:
        return ParkFunc(self.N, self.M, w_area=area_word, w_label=label_word, specific_opt=self.specific_opt)
    else:
        return [ParkFunc(self.N, self.M, w_area=area_word, w_label=label_word, specific_opt=self.specific_opt), not_tdinv, not_dinvcorr]