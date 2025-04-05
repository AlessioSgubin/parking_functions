#
#       A program to test Parking Functions and the bijection with k-Delay
#
import copy
import sys

TRY = 7

def kDyckPaths(n,k):                        # kDYCKPATHS: computes the set of all (nk,n)-Dyck paths
    possible_paths = [[1]]
    for i in range(n*(k+1)-1):
        new_possible = []
        for part_path in possible_paths:
            height = sum(part_path)
            distance = sum([1-part_path[j] for j in range(i+1)])
            if height == n:                 # Reached max height, just 0's
                temp = deepcopy(part_path)
                temp.append(0)
                new_possible.append(temp)
            elif k*height < distance + 1:   # Reached the "diagonal", just a 1
                temp = deepcopy(part_path)
                temp.append(1)
                new_possible.append(temp)
            else:
                temp = deepcopy(part_path)
                temp.append(0)
                new_possible.append(temp)
                temp = deepcopy(part_path)
                temp.append(1)
                new_possible.append(temp)
        possible_paths = new_possible
    dyck_k_paths = [[sum([1-w[i+j] for j in range(len(w)-i)]) for i in range(len(w)) if w[i]==1] for w in possible_paths]
    return dyck_k_paths

def kParkingFunctions(n,k):                 # kPARKINGFUNCTIONS: computes the set of all (n,k)-parking functions
    parking_functions = []
    counter = 0
    # Compute total number of Dyck paths
    totalnum = binomial(n*(k+1)+1, n)/(n*(k+1)+1)
    print("Starting the algorithm to compute all parking functions.")
    ### Loop on k-Dyck paths
    for path in kDyckPaths(n,k):
        counter += 1
        perc = floor(100*counter/totalnum)
        sys.stdout.write('\r')                          # Reset to start of line
        sys.stdout.write("Percentage %3d %%, computing parking functions for path %6d of %6d" % (perc, counter, totalnum))
        sys.stdout.flush()
        #print("Computing {}".format(path))
        possible_labels = [[[],set([i+1 for i in range(n)])]]
        for row in range(n*k):
            #print(possible_labels)
            new_possible_labels = []
            for [label,remaining] in possible_labels:
                s = len([j for j in path if j == n*k - row])
                for sub in Subsets(remaining,s):
                    lab = sorted(sub)
                    new_possible_labels = new_possible_labels + [[(label + lab),remaining-set(sub)]]
            possible_labels = new_possible_labels
        parking_functions = parking_functions + [[path, label] for [label,rem] in possible_labels]
    print("\nCompleted computation!\n")
    return parking_functions

def nkParkingFunctions(n,k):                # nkPARKINGFUNCTIONS: computes the set of all (n,k)-parking functions
    parking_functions = []
    counter = 0
    # Compute total number of Dyck paths
    totalnum = binomial(n*(k+1)+1, n)/(n*(k+1)+1)
    print("Starting the algorithm to compute all parking functions.")
    # Compute all possible set partitions
    set_partitions = {}
    for lamb in Partitions(n):
        temp = []
        for possible in OrderedSetPartitions(n,lamb):
            poss = [list(part) for part in possible]
            temp.append(poss)
        lambt = tuple(lamb)
        set_partitions.update({lambt:temp})

    ### Loop on k-Dyck paths
    for path in kDyckPaths(n,k):
        counter += 1
        perc = floor(100*counter/totalnum)
        sys.stdout.write('\r')                          # Reset to start of line
        sys.stdout.write("Percentage %3d %%, computing parking functions for path %6d of %6d" % (perc, counter, totalnum))
        sys.stdout.flush()
        #print("Computing {}".format(path))

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
            parking_functions.append(copy.copy([path,label]))

    print("\nCompleted computation!\n")
    return parking_functions


def area(n,k,w):                            # AREA: computes the area of a given path
    return sum([w[i]-n*k+i*k for i in range(n)])

def pdinv(n,k,w):                           # PDINV: computes the path dinv of a given path
    dinv = 0
    # Read all cells above path
    for i in range(n):                  # Loop on rows
        for j in range(n*k-w[i]):       # Loop on columns
            arm = (n*k - w[i]) - (j+1)
            leg = i - max([h for h in range(n) if n*k-w[h]-1 < j]) - 1
            if (arm <= k*(leg+1)) and (k*leg < arm + 1):
                dinv += 1
    return dinv

def tdinv(n,k,w_path,w_label):              # TDINV
    tdinv = 0
    row_cresc = [w_label.index(h+1) for h in range(n)]
    for i in range(n):
        rank_c = n*w_path[row_cresc[i]] + (n*k+1)*row_cresc[i]
        for j in range(n-i-1):
            rank_d = n*w_path[row_cresc[i+j+1]] + (n*k+1)*row_cresc[i+j+1]
            #print("Cell {} has rank {} and cell {} has rank {}".format(i+1,rank_c,i+j+2,rank_d))
            if rank_c < rank_d and rank_d < rank_c + (n*k):
                tdinv += 1
    return tdinv

def maxtdinv(n,k,w_path):                   # MAXTDINV
    rank_list = [(n*k+1)*i + n*w_path[i] for i in range(n)]
    order = sorted(range(len(rank_list)), key=lambda k: rank_list[k])
    park = sorted(range(len(order)), key=lambda k: order[k])
    park = [i+1 for i in park]
    return tdinv(n,k,w_path,park)

def max_w_read(n,k,w_path):                 # MAX_W_READ: returns the labelling maximising tdinv
    rank_list = [(n*k+1)*i + n*w_path[i] for i in range(n)]
    order = sorted(range(len(rank_list)), key=lambda k: rank_list[k])
    park = sorted(range(len(order)), key=lambda k: order[k])
    return [i+1 for i in park]

def dinv(n,k,w_path,w_label):               # DINV: computes the dinv of a parking function
    #print("La statistica pdinv è {}".format(pdinv(n,k,w_path)))
    #print("La statistica maxtdinv è {}".format(maxtdinv(n,k,w_path)))
    #print("La statistica tdinv è {}".format(tdinv(n,k,w_path,w_label)))
    return pdinv(n,k,w_path) + tdinv(n,k,w_path,w_label) - maxtdinv(n,k,w_path)

def delay(n,k,w_path,w_label,contr=False):              # DELAY: computes the delay of a parking function
    delay = 0
    plus = 0
    current = n+1
    reading_word = []
    contributes = []
    buffer = []
    for i in range(n*k):
        # Adding to buffer new labels with k multeplicity
        new = [w_label[j] for j in range(n) if w_path[j] == n*k-i]
        buffer = buffer + [new[floor(j/k)] for j in range(len(new)*k)]
        # Get next step
        if current <= min(buffer):          # Starting new ascending chain
            plus += 1
            current = max(buffer)
        else:                               # Continuing current ascending chain
            current = max([w for w in buffer if w < current])
        
        if current not in reading_word:     # Adding delay if needed...
            delay += plus
            contributes = contributes + [[current,plus]]
        buffer.remove(current)                 # Removing one copy of current from buffer

        # writing down the new letter in reading word
        reading_word = reading_word + [current]
    if contr:
        return [delay,contributes,reading_word]
    else:
        return delay

def gen_poly(n,k,set_paths,ring):           # GEN_POLY: polynomial generated by paths
    q = ring.gens()[0]
    t = ring.gens()[1]
    poly = 0*q*t
    for path in set_paths:
        poly = poly + (q**pdinv(n,k,path))*(t**area(n,k,path))
    return poly

def path_word(n,k,w_path,w_label):
    rank_list = [(n*k+1)*i + n*w_path[i] for i in range(n)]
    path_order = sorted(range(len(rank_list)), key=lambda k: rank_list[k])      # Reading order of the columns
    path_word = [w_label[path_order[i]] for i in range(n)]
    return path_word

def Psi(n,k,w_path,w_label,infos = False):                # PSI: this is the bijection
    # This function tests the bijection between (dinv,area) -> (area,delay)

    ### Create the path word, ordered by rank
    rank_list = [(n*k+1)*i + n*w_path[i] for i in range(n)]
    path_order = sorted(range(len(rank_list)), key=lambda k: rank_list[k])      # Reading order of the columns
    path_word = [w_label[path_order[i]] for i in range(n)]                      # Labels in order
    #print("The path_word: {}".format(path_word))

    ### Compute the not_tdinv
    # Compute all pairs of cars making tdinv
    tdinv_pairs = []            
    row_cresc = [w_label.index(h+1) for h in range(n)]
    for i in range(n):
        rank_c = n*w_path[row_cresc[i]] + (n*k+1)*row_cresc[i]
        for j in range(n-i-1):
            rank_d = n*w_path[row_cresc[i+j+1]] + (n*k+1)*row_cresc[i+j+1]
            if rank_c < rank_d and rank_d < rank_c + (n*k+1):
                tdinv_pairs = tdinv_pairs + [(i+1,i+j+2)]

    # Compute not_tdinv with previous cars
    not_tdinv = []
    for i in range(n):
        temp = 0
        for j in range(i):
            if ((path_word[i],path_word[j]) not in tdinv_pairs) and ((path_word[j],path_word[i]) not in tdinv_pairs):
                temp += 1
        not_tdinv = not_tdinv + [temp]
    #print("The not_tdinv: {}".format(not_tdinv))

    ### Compute the not_dinvcorr
    # Compute all pairs of cars making dinvcorr
    dinvcorr_pairs = []
    for i in range(n):
        rank_c = n*w_path[i] + (n*k+1)*i
        for j in range(n-i-1):
            rank_d = n*w_path[i+j+1] + (n*k+1)*(i+j+1)
            for shift in range(k-1):
                if rank_d < rank_c + n*(k-1) - n*shift and rank_c - n - n*shift < rank_d:
                    dinvcorr_pairs = dinvcorr_pairs + [(w_label[i],w_label[i+j+1])]
    # Compute not_dinvcorr with previous cars
    not_dinvcorr = []
    for i in range(n):
        temp = 0
        for j in range(i):
            temp += k - 1 - len([1 for (a,b) in dinvcorr_pairs if (a,b)==(path_word[i],path_word[j]) or (a,b)==(path_word[j],path_word[i])])
        not_dinvcorr = not_dinvcorr + [temp]
    #print("The not_dinvcorr: {}".format(not_dinvcorr))
    '''
    if not all(a <= b for a, b in zip(not_dinvcorr, not_dinvcorr[1:])):
        print("NO\nNO\nNO")
    '''

    ### Compute the not_dinv
    not_dinv = [not_tdinv[i] + not_dinvcorr[i] for i in range(n)]
    #print("The not_dinv: {}".format(not_dinv))

    ### Compute the w_path of the image
    not_dinv_sort = sorted(not_dinv, key=lambda k: k)
    area_word = [n*k-not_dinv_sort[i] for i in range(n)]
    
    ### Compute the w_label of the image
    label_word = [path_word[i] for i in sorted(range(n), key=lambda k: not_dinv[k]*n + path_word[k])]
    if not infos:
        return [area_word, label_word]
    else:
        return [area_word, label_word, not_tdinv, not_dinvcorr]

def sweep_map(n,k,w_path):
    # Write down the sequence of north (N)/east (E) steps
    steps = ['N']
    for i in range(n-1):
        steps = steps + ['E' for j in range(w_path[i]-w_path[i+1])] + ['N']
    steps = steps + ['E' for j in range(w_path[n-1]+1)]
    # Write the rank of each "wand"
    rank_steps = [0]
    current = 0
    for i in range(n*(k+1)):
        if steps[i] == 'N':
            current += n*k+1
        else:
            current -= n
        rank_steps = rank_steps + [current]
    # Write down sequence of steps reordered by rank
    new_order = sorted(range(len(rank_steps)), key=lambda i: rank_steps[i])
    sweep_steps = [steps[i] for i in new_order]
    # Compute new w_path
    w_sweep_path = []
    current = n*k
    for w in sweep_steps:
        if w == 'N':
            w_sweep_path = w_sweep_path + [current]
        else:
            current -= 1
    return w_sweep_path

def inverse_sweep_map(n,k,w_path):
    M = n*k+1
    N = n
    # Write down the sequence of north (N)/east (E) steps
    steps = ['N']
    for i in range(n-1):
        steps = steps + ['E' for j in range(w_path[i]-w_path[i+1])] + ['N']
    steps = steps + ['E' for j in range(w_path[n-1]+1)]
    # Write down the arrows in the grid with starting/ending heights
    grid = []
    min_rank = 0
    for st in steps:
        if st == 'N':
            point_arr = min_rank + M
            grid.append([min_rank,point_arr])
            min_rank += 1
        else:
            if min_rank - N < 0:
                min_rank = N
            point_arr = min_rank - N
            grid.append([min_rank,point_arr])
            min_rank += 1
    # Move all the arrows according to the algorithm
    num_rows = max([max(a,b) for [a,b] in grid]) + 1
    row_counts = [sum([sign(b-a) for [a,b] in grid if min([a,b]) <= i < max([a,b])]) for i in range(num_rows)]
    while row_counts != [0 for i in range(num_rows)]:
        #print(row_counts)
        #print(grid)
        lev = min([i for i in range(num_rows) if row_counts[i] > 0])
        moved_arr = min([i for i in range(len(grid)) if grid[i][0] == lev])
        to_shift = moved_arr
        to_lev = lev
        while to_shift < len(grid):
            if grid[to_shift][0] == to_lev:
                grid[to_shift][0] += 1
                grid[to_shift][1] += 1
                to_shift += 1
                to_lev += 1
            else:
                to_shift = len(grid)
        num_rows = max([max(a,b) for [a,b] in grid]) + 1
        row_counts = [sum([sign(b-a) for [a,b] in grid if min([a,b]) <= i < max([a,b])]) for i in range(num_rows)]
    # Reorder the path
    new_path = []
    rank = min([min([a,b]) for [a,b] in grid])
    #print("Finished: {}".format(grid))
    while grid != []:
        ind = min([i for i in range(len(grid)) if grid[i][0] == rank])
        step = grid[ind]
        rank = step[1]
        new_path.append(step)
        grid.remove(step)
    # Compute new w_path
    sweep_steps = [sign(b - a) for [a,b] in new_path]
    w_sweep_path = []
    current = n*k
    for w in sweep_steps:
        if w == 1:
            w_sweep_path = w_sweep_path + [current]
        else:
            current -= 1
    return w_sweep_path


def area_word(n,k,w_path):
    return [w_path[i] - (n-i)*k for i in range(n)]

def tdinv_word(n,k,w_path,w_label,path_word = []):
    ### Create the path word, ordered by rank
    rank_list = [(n*k+1)*i + n*w_path[i] for i in range(n)]
    path_order = sorted(range(len(rank_list)), key=lambda k: rank_list[k])      # Reading order of the columns
    if path_word == []:
        path_word = [w_label[path_order[i]] for i in range(n)]                      # Labels in reading order
    #print("The path_word: {}".format(path_word))

    ### Compute the not_tdinv
    # Compute all pairs of cars making tdinv
    tdinv_pairs = []            
    row_cresc = [w_label.index(h+1) for h in range(n)]
    for i in range(n):
        rank_c = n*w_path[row_cresc[i]] + (n*k+1)*row_cresc[i]
        for j in range(n-i-1):
            rank_d = n*w_path[row_cresc[i+j+1]] + (n*k+1)*row_cresc[i+j+1]
            if rank_c < rank_d and rank_d < rank_c + (n*k+1):
                tdinv_pairs = tdinv_pairs + [(i+1,i+j+2)]

    # Compute tdinv with previous cars
    tdinv = []
    for i in range(n):
        temp = 0
        for j in range(i):
            if ((path_word[i],path_word[j]) in tdinv_pairs) or ((path_word[j],path_word[i]) in tdinv_pairs):
                temp += 1
        tdinv = tdinv + [temp]
    return tdinv

def not_tdinv_word(n,k,w_path,w_label,path_word = []):
    ### Create the path word, ordered by rank
    rank_list = [(n*k+1)*i + n*w_path[i] for i in range(n)]
    path_order = sorted(range(len(rank_list)), key=lambda k: rank_list[k])      # Reading order of the columns
    path_word = [w_label[path_order[i]] for i in range(n)]                      # Labels in order
    
    ### Compute the not_tdinv
    # Compute all pairs of cars making tdinv
    tdinv_pairs = []            
    row_cresc = [w_label.index(h+1) for h in range(n)]
    for i in range(n):
        rank_c = n*w_path[row_cresc[i]] + (n*k+1)*row_cresc[i]
        for j in range(n-i-1):
            rank_d = n*w_path[row_cresc[i+j+1]] + (n*k+1)*row_cresc[i+j+1]
            if rank_c < rank_d and rank_d < rank_c + (n*k+1):
                tdinv_pairs = tdinv_pairs + [(i+1,i+j+2)]
    # Compute not_tdinv with previous cars
    not_tdinv = []
    for i in range(n):
        temp = 0
        for j in range(i):
            if ((path_word[i],path_word[j]) not in tdinv_pairs) and ((path_word[j],path_word[i]) not in tdinv_pairs):
                temp += 1
        not_tdinv = not_tdinv + [temp]
    return not_tdinv

def max_tdinv_label(n,k,w_path):
    # Compute w_label for maximum tdinv
    rank_list = [(n*k+1)*i + n*w_path[i] for i in range(n)]
    order = sorted(range(len(rank_list)), key=lambda k: rank_list[k])
    park = sorted(range(len(order)), key=lambda k: order[k])
    w_label_max = [i+1 for i in park]
    return w_label_max

def max_tdinv_word(n,k,w_path):
    # Compute w_label for maximum tdinv
    rank_list = [(n*k+1)*i + n*w_path[i] for i in range(n)]
    order = sorted(range(len(rank_list)), key=lambda k: rank_list[k])
    park = sorted(range(len(order)), key=lambda k: order[k])
    w_label_max = [i+1 for i in park]
    tdinv = tdinv_word(n,k,w_path,w_label_max)
    return tdinv

def not_pdinv_word(n,k,w_path):
    rank_list = [(n*k+1)*i + n*w_path[i] for i in range(n)]
    path_word = sorted(range(len(rank_list)), key=lambda k: rank_list[k])
    not_dinv_word = []
    # Computing the correlation E steps and N steps via diagonal shadows
    bijection_E_steps = []      # Stores the corresponding N step to each E step
    for i in range(n*k):
        height_E = max([j+1 for j in range(n) if n*k-w_path[j] <= i])
        corresp = max([j for j in range(height_E) if ((i - n*k + w_path[j]) >= k*(height_E - j - 1)) and ((i - n*k + w_path[j]) < k*(height_E - j - 1) + k)])
        bijection_E_steps.append(corresp + 1)
    #print(w_path)
    #print(bijection_E_steps)
    
    # Read all cells above path, searching for pdinv pairs
    pdinv_pairs = []
    for i in range(n):                  # Loop on rows
        for j in range(n*k-w_path[i]):       # Loop on columns
            arm = (n*k - w_path[i]) - (j+1)
            leg = i - max([h for h in range(n) if n*k-w_path[h]-1 < j]) - 1
            if (arm <= k*(leg+1)) and (k*leg < arm + 1):
                pdinv_pairs.append((bijection_E_steps[j],i+1))
    #print(pdinv_pairs)

    # Construct the pdinv_word using the reading order
    not_pdinv = []
    for i in range(n):
        temp = 0
        for j in range(i):
            temp += k - len([1 for (a,b) in pdinv_pairs if (a,b)==(path_word[i]+1,path_word[j]+1) or (a,b)==(path_word[j]+1,path_word[i]+1)])
        not_pdinv.append(temp)
    return not_pdinv


def dinvcorr_word(n,k,w_path,w_label,path_word = []):
    ### Create the path word, ordered by rank
    if path_word == []:
        rank_list = [(n*k+1)*i + n*w_path[i] for i in range(n)]
        path_order = sorted(range(len(rank_list)), key=lambda k: rank_list[k])      # Reading order of the columns
        path_word = [w_label[path_order[i]] for i in range(n)]                      # Labels in reading order
    #print("The path_word: {}".format(path_word))
    ### Compute the not_dinvcorr
    # Compute all pairs of cars making dinvcorr
    dinvcorr_pairs = []
    for i in range(n):
        rank_c = n*w_path[i] + (n*k+1)*i
        for j in range(n-i-1):
            rank_d = n*w_path[i+j+1] + (n*k+1)*(i+j+1)
            for shift in range(k-1):
                if rank_d < rank_c + n*(k-1) - n*shift and rank_c - n - n*shift < rank_d:
                    dinvcorr_pairs = dinvcorr_pairs + [(w_label[i],w_label[i+j+1])]
    # Compute dinvcorr with previous cars
    dinvcorr = []
    for i in range(n):
        temp = 0
        for j in range(i):
            temp += len([1 for (a,b) in dinvcorr_pairs if (a,b)==(path_word[i],path_word[j]) or (a,b)==(path_word[j],path_word[i])])
        dinvcorr = dinvcorr + [temp]
    return dinvcorr

def not_dinvcorr_word(n,k,w_path,w_label,path_word = []):
    ### Create the path word, ordered by rank
    if path_word == []:
        rank_list = [(n*k+1)*i + n*w_path[i] for i in range(n)]
        path_order = sorted(range(len(rank_list)), key=lambda k: rank_list[k])      # Reading order of the columns
        path_word = [w_label[path_order[i]] for i in range(n)]                      # Labels in reading order
    ### Compute the not_dinvcorr
    # Compute all pairs of cars making dinvcorr
    dinvcorr_pairs = []
    for i in range(n):
        rank_c = n*w_path[i] + (n*k+1)*i
        for j in range(n-i-1):
            rank_d = n*w_path[i+j+1] + (n*k+1)*(i+j+1)
            for shift in range(k-1):
                if rank_d < rank_c + n*(k-1) - n*shift and rank_c - n - n*shift < rank_d:
                    dinvcorr_pairs = dinvcorr_pairs + [(w_label[i],w_label[i+j+1])]
    # Compute not_dinvcorr with previous cars
    not_dinvcorr = []
    for i in range(n):
        temp = 0
        for j in range(i):
            temp += k - 1 - len([1 for (a,b) in dinvcorr_pairs if (a,b)==(path_word[i],path_word[j]) or (a,b)==(path_word[j],path_word[i])])
        not_dinvcorr = not_dinvcorr + [temp]
    return not_dinvcorr

def tries_w_reading(n,k,w_path,w_label):
    ### For now just try to get the reading word...
    [d,contrib] = delay(n,k,w_path,w_label,contr=True)
    rank_list = [(n*k+1)*i + n*w_path[i] for i in range(n)]
    rank_order = [w_label[i] for i in sorted(range(n), key=lambda j: rank_list[j])]
    #order = sorted(range(n), key=lambda i: n*contrib[i][1] + rank_order.index(contrib[i][0]))
    order = sorted(range(n), key=lambda i: n*k*contrib[i][1] + w_path[w_label.index(contrib[i][0])])
    w_reading = [contrib[i][0] for i in order]
    return w_reading

def w_reading(n,k,w_path,w_label):
    rank_list = [(n*k+1)*i + n*w_path[i] for i in range(n)]
    path_order = sorted(range(len(rank_list)), key=lambda k: rank_list[k])
    path_word = [w_label[path_order[i]] for i in range(n)]
    return path_word


def w_delay_by_contributes(n,k,contrib):
    # contrib = [area_label_1, area_label_2, ...]
    rounds = max(contrib) + k
    word = []
    buffer = []
    for n_loop in range(rounds):
        current = len(contrib)
        new_labels = [j for j in range(len(contrib)) if contrib[j] == n_loop for l in range(k)]
        buffer = buffer + new_labels
        if buffer == []:
            check = False
        else:
            check = (current > min(buffer))
        while check:
            temp = max([j for j in buffer if j < current])
            buffer.remove(temp)
            word = word + [temp]
            current = temp
            if buffer == []:
                check = False
            else:
                check = (current > min(buffer))
    return [w+1 for w in word]
    

def Phi(n,k,w_path,w_label,infos = False):
    ### Compute delay and delay_word
    dela,info_delay,w_delay = delay(n,k,w_path,w_label,contr=True)
    
    ### Compute the path of the image
    path_max_tdinv = []
    areas = [info_delay[i][1] for i in range(len(info_delay))]
    areas.sort()
    w_delay_max = w_delay_by_contributes(n,k,areas)
    current_label = 0
    weak_des = [0] + [i+1 for i in range(n*k-1) if w_delay_max[i]<=w_delay_max[i+1]]
    #print("\nArea {}".format(areas))
    #print("Delay_max {}".format(w_delay_max))
    for a in range(max(areas) + 1):
        num_areas = len([1 for w in areas if w == a])
        added_labels = [current_label + (num_areas-i) for i in range(num_areas)]
        for j in range(len(added_labels)):
            '''
            if a-1 not in areas:        # No labels with a-1 area
                path_max_tdinv.append(n*k - w_delay_max.index(added_labels[0]))
            else:                       # Some labels with a-1 area
                if added_labels[j] > max(previous_areas):           # label greater than those from previous
                    path_max_tdinv.append(n*k - w_delay_max.index(max(previous_areas)) - 1)
                    #path_max_tdinv.append(n*k - len(path_max_tdinv))
                else:                                               # otherwise
                    path_max_tdinv.append(n*k - w_delay_max.index(min([u for u in previous_areas if u > added_labels[j]])) - 1)
            '''
            first_occ = w_delay_max.index(added_labels[0])          # Position of the first occurrence of added_labels[j]
            if a > 0:                                           # We are NOT in the first descending run
                # Compute the first position from which all following labels are less than added_labels[j]
                lowest_pos = min([ l for l in range(first_occ) if max(w_delay_max[l:first_occ]) < added_labels[j] ]) + 1
                # Compute the lowest position possible in the previous descending run
                lowest_last_des = weak_des[a-1] + 1
                posit = max([lowest_pos,lowest_last_des])
                path_max_tdinv.append(n*k - posit)
            else:                                               # We are in the first descending run
                path_max_tdinv.append(n*k)
        current_label += num_areas
    #print("Delay word {} and max_path {}".format(w_delay_max,path_max_tdinv))
    #print("Path_max_dinv {}".format(path_max_tdinv))

    path_image = inverse_sweep_map(n,k,path_max_tdinv)              # Obtain the image path by inverting the sweep map

    return path_max_tdinv,path_image,w_delay_max
    
def visual_degree(p):
    q_max = p.degrees()[0]
    t_max = p.degrees()[1]
    max_deg = max([q_max,t_max])
    # Compute the matrix
    coeff_matrix = [[p.coefficient([q_deg, t_deg]) for t_deg in range(max_deg + 1)] for q_deg in range(max_deg + 1)]
    # Plot the matrix
    return matrix_plot(matrix(coeff_matrix))



### Now let's test the bijection...
n = 6
k = 2

if TRY == 1:
    for [path1,label1] in nkParkingFunctions(n,k):
        area1 = area(n,k,path1)
        dinv1 = dinv(n,k,path1,label1)

        path2,label2,not_tdinv1,not_dinvcorr1 = Psi(n,k,path1,label1, infos=True)
        area2 = area(n,k,path2)

        area1_word = [- k*(n-i) + path1[i] for i in range(n)]
        [delay2,contrib] = delay(n,k,path2,label2,contr=True)
        delay_ord = sorted(range(n), key=lambda l: label1.index(contrib[l][0]))
        delay_word = [contrib[i][1] for i in delay_ord]
        check = (area1_word == delay_word)

        w_read1 = w_reading(n,k,path1,label1)
        w_read2 = tries_w_reading(n,k,path2,label2)
        check2 = (w_read1 == w_read2)

        path3 = sweep_map(n,k,path1)
        difference = [path3[i]-path2[i] for i in range(n)]
        check3 = (sum(difference) == maxtdinv(n,k,path1) - tdinv(n,k,path1,label1))

        difference2 = [not_dinvcorr1[i]-not_tdinv1[i] for i in range(n)]
        check4 = (difference == difference2)

        maxtdinv1 = max_tdinv_word(n,k,path1)
        tdinv1 = tdinv_word(n,k,path1,label1)
        difference3 = [- tdinv1[i] + maxtdinv1[i] for i in range(n)]
        check5 = (difference == difference3)

        '''
        if (area1 == delay2) and (dinv1 == area2):
            print("Funziona: \t[{},{}] \t---> \t[{},{}]".format(path1,label1,path2,label2))
            print("|\tdinv1 = area2\t{}".format(dinv1))
            print("|\tarea1 = delay2\t{}".format(area1))
            print("|\tw_area1 =?= w_delay2: {}\t {}".format(check,area1_word))
            print("|\tRight w_reading: {}\t {}\t{}".format(check2,w_read1,w_read2))
            print("")
        else:
            print("\nNO")
            print("Path1 {} con Label1 {} ha dinv {} e area {}".format(path1,label1,dinv1,area1))
            print("Path2 {} con Label2 {} ha area {} e delay {}\n".format(path2,label2,area2,delay2))
        '''
        '''
        if not check5:
            print("Funziona: \t[{},{}] \t---> \t[{},{}]".format(path1,label1,path2,label2))
            #print("|\tdinv1 = area2\t{}".format(dinv1))
            #print("|\tarea1 = delay2\t{}".format(area1))
            #print("|\tw_area1 =?= w_delay2: {}\t {}".format(check,area1_word))
            #print("|\tRight w_reading: {}\t {}\t{}".format(check2,w_read1,w_read2))
            #print("|\tPaths {}\t {}\t Difference {}".format(path2,path3,[path3[i]-path2[i] for i in range(n)]))
            print("|\tDifferences\t{}\t {}".format(difference, difference3))
            print("|\t\t\t{}\t {} with {}\t{}".format(tdinv1,maxtdinv1,tdinv(n,k,path1,label1),maxtdinv(n,k,path1)))
            print("")
        '''
        print("Funziona: \t[{},{}] \t---> \t[{},{}]".format(path1,label1,path2,label2))
        
        reading_new = Phi_bruteforce(n,k,path2,label2)
        reading_p1 = w_reading(n,k,path1,label1)
        
        if reading_new != reading_p1:
            print("NO! {}\t{}".format(reading_p1,reading_new))
            print(reading_new)
        else:
            print("YES!")

if TRY == 2:
    nkPF = nkParkingFunctions(n,k)
    nkPF2 = copy.copy(nkPF)
    print("There are {} parking functions".format(len(nkPF)))
    for [path1,label1] in nkPF:
        print("Doing {},{}".format(path1,label1))
        path2,label2 = Psi(n,k,path1,label1)

        if [path2,label2] in nkPF2:
            nkPF2.remove([path2,label2])
        else:
            print("Not bijective...")
    if nkPF2 == []:
        print("The map is bijective!")
    else:
        print("The map is NOT bijective!")

if TRY == 3:
    nkPF = nkParkingFunctions(n,k)
    for [path1,label1] in nkPF:
        # Compute area word
        area1 = area_word(n,k,path1)
        # Compute path word
        path_word1 = path_word(n,k,path1,label1)
        # Compute dinvcorr word (order: path word)
        not_dinvcorr1 = not_dinvcorr_word(n,k,path1,label1)
        ordering = [path_word1.index(label1[i]) for i in range(n)]
        not_dinvcorr1 = [not_dinvcorr_temp[i] for i in ordering]
        # Compute tdinv word
        not_tdinv1 = not_tdinv_word(n,k,path1,label1)
        not_tdinv1 = [not_tdinv_temp[i] for i in ordering]

        diff = [not_dinvcorr1[i] - area1[i] for i in range(n)]
        diff2 = [i + diff[i] - not_tdinv_temp[i] for i in range(n)]
        print("Difference word {} \tis {} \t {}".format([path1,label1],diff, diff2))

if TRY == 4:
    nkPF = nkParkingFunctions(n,k)
    for [path1,label1] in nkPF:
        # Compute max_tdinv_label and images, delays
        label1_max = max_tdinv_label(n,k,path1)
        [path2,label2, not_tdinv1, not_dinvcorr1] = Psi(n,k,path1,label1,infos=True)
        not_pdinv1 = not_pdinv_word(n,k,path1)
        not_pdinv1.sort()
        [path3,label3] = Psi(n,k,path1,label1_max)
        del2,contr2,w_delay2 = delay(n,k,path2,label2,contr=True)
        del3,contr3,w_delay3 = delay(n,k,path3,label3,contr=True)
        
        obtained_path = sweep_map(n,k,path1)
        diff = [path3[i]+not_pdinv1[i] for i in range(n)] 

        if path3 == obtained_path:
            print("Original path is {}\t{}".format(path1,label1))
            print("Original not_tdinv \t{}".format(not_tdinv1))
            print("Original not_dinvc \t{}".format(not_dinvcorr1))
            print("Original not_pdinv \t{}".format(not_pdinv1))
            print("Maximum  path is {}\t{}".format(path1,label1_max))
            print("Path obtained is \t{}\t{} with delay word {}".format(path2,label2,w_delay2))
            print("Path maximize is \t{}\t{} with delay word {}".format(path3,label3,w_delay3))
            print("Sum {}\t\t THIS WORKS!".format(diff))
            print("")
            print("Sweep map is {}".format(obtained_path))

if TRY == 5:
    # Testing inverse sweep map
    kDyck = kDyckPaths(n,k)
    for path in kDyck:
        path2 = sweep_map(n,k,path)
        path3 = inverse_sweep_map(n,k,path2)
        if path != path3:
            print("Path {} \t path2 {} \t path3 {}".format(path,path2,path3))

if TRY == 6:
    # Testing new Phi map
    nkPF = nkParkingFunctions(n,k)
    for [path1,label1] in nkPF:
        [path2,label2] = Psi(n,k,path1,label1)
        delay2,contrs,w_del2 = delay(n,k,path2,label2,contr=True)

        label1_max = max_tdinv_label(n,k,path1)
        [path3,label3, not_tdinv3,not_dinvcorr3] = Psi(n,k,path1,label1_max,infos=True)

        path3_c,path1_c,w_del_max = Phi(n,k,path2,label2)
        '''
        if path3 != path3_c:
            print("Path3 {} \tand Path_max_tdinv {}".format(path3,path3_c))
            #print(path1)
        '''
        if path1 != path1_c and max(area_word(n,k,path1)) <= 1:
            tdinv1 = tdinv_word(n,k,path1,label1)
            maxtdinv1 = max_tdinv_word(n,k,path1)
            pdinv1 = not_pdinv_word(n,k,path1)
            tdinv_max = tdinv_word(n,k,path1,label1_max)
            maxtdinv_max = max_tdinv_word(n,k,path1)
            pdinv_max = not_pdinv_word(n,k,path1)
            print("")
            print("Path1 {} \tand label1 {}".format(path1,label1))
            print("\t Tdinv {}\n\t Mdinv {}\n\t Pdinv {}".format(tdinv1,maxtdinv1,pdinv1))
            print("Path1 {} \tand label1_max {}".format(path1,label1_max))
            print("\t Tdinv {}\n\t Mdinv {}\n\t Pdinv {}".format(tdinv_max,maxtdinv_max,pdinv_max))
            print("Path2 {} \tand label2 {}".format(path2,label2))
            print("Path3 {} \tand Path_max_tdinv {}".format(path3,path3_c))
            print("Now word delay {}".format(w_del2))
            print("Max word delay {}".format(w_del_max))

if TRY == 7:
    nPF = nkParkingFunctions(n,1)
    for [path1,label1] in nPF:
        # Compute image
        [path2,label2] = Psi(n,1,path1,label1)
        # Compute statistics
        area1 = area(n,1,path1)
        dinv1 = dinv(n,1,path1,label1)
        area2 = area(n,1,path2)
        dinv2 = dinv(n,1,path2,max_w_read(n,1,path2))
        dela2 = delay(n,1,path2,label2)
        if dinv2 != dela2:
            print("Path 1 is:\t{}\t and has dinv: {}\t area: {}".format([path1,label1],dinv1,area1))
            print("Path 2 is:\t{}\t and has area: {}\t dinv: {} \t delay:{}\n".format([path2,label2],area2,dinv2,dela2))

    