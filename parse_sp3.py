# parse_sp3.py

filename = ''
sat_ephem = {}
time = []

with open(filename,'r') as sp3file:
    lines = sp3file.readlines()
    num_sat = int(lines[2].split()[1])

    for s in range(num_sat):
        sat_ephem['{:02}'.format(s+1)] = []

    for j in range(22,len(lines)-1,num_sat+1):

        t = lines[j][2:-6]
        time.append(dt.datetime.strptime(t, ' %Y %m %d %H %M %S.%f'))

        for i in range(num_sat):
            s = lines[j+i+1].split()
            sat_ephem[s[0][2:]].append([float(s[1]),float(s[2]),float(s[3])])

for sat in sat_ephem:
    sat_ephem[sat] = np.array(sat_ephem[sat]).T*1000.
