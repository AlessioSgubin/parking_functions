### TEST FILE FOR PARKFUNC CLASS

# Parameters
n = 4
k = 3
TRY = 1

if TRY == 1:
    nkPF = nkParkingFunctions(n,k)
    for pf in nkPF:
        img = pf.to_area_pmaj()
        img2 = img.to_dinv_area()
        if pf != img2:
            print("\nPROBLEM!")
            pf.draw(stats=False)
            img2.draw(stats=False)
            img.to_dinv_area(infos=True)
