### TEST FILE FOR PARKFUNC CLASS

def visual_degree(p):
    q_max = p.degrees()[0]
    t_max = p.degrees()[1]
    max_deg = max([q_max,t_max])
    # Compute the matrix
    coeff_matrix = [[p.coefficient([q_deg, t_deg]) for t_deg in range(max_deg + 1)] for q_deg in range(max_deg + 1)]
    # Plot the matrix
    return matrix_plot(matrix(coeff_matrix))


# Parameters
n = 6
k = 2
TRY = 3

if TRY == 1:
    nkPF = nkParkingFunctions(n,k)
    counter = 0
    totalnum = len(nkPF)
    for pf in nkPF:
        counter += 1
        perc = floor(100*counter/totalnum)              # type: ignore
        sys.stdout.write('\r')                          # Reset to start of line
        sys.stdout.write("Percentage %3d %%, computing parking functions for path %6d of %6d" % (perc, counter, totalnum))
        sys.stdout.flush()
        img = pf.to_area_pmaj()
        img2 = img.to_dinv_area(infos=False)
        if pf != img2:
            print("\nPROBLEM!")
            pf.draw(stats=False)
            img2.draw(stats=False)
            img.to_dinv_area(infos=True)
        #else:
        #    print("\n\n ALL GOOD!")
        #    pf.draw(stats=True)
        #    img.draw(stats=True)
            #img2.draw(stats=True)
    print("")

if TRY == 2:
    nkPF = nkParkingFunctions(n,k, display=True)
    size = k*(n*(n-1))/2 + 1
    counter = 0
    totalnum = len(nkPF)
    distr = [[0 for i in range(size)] for j in range(size)]

    for pf in nkPF:
        counter += 1
        perc = floor(100*counter/totalnum)              # type: ignore
        sys.stdout.write('\r')                          # Reset to start of line
        sys.stdout.write("Percentage %3d %%, computing parking functions for path %6d of %6d" % (perc, counter, totalnum))
        sys.stdout.flush()

        x = pf.dinv()
        y = pf.pmaj()
        distr[x][y] += 1
    X = matrix_plot(matrix(distr))

if TRY == 3:
    pf = ParkFunc(7,21,func=[4,1,7,8,1,7,16])
    pf.draw()
    img = pf.to_dinv_area(infos=True)
    img.draw()