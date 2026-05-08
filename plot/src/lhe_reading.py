import numpy as np

class Data:
    momenta_l_minus = []
    momenta_l_plus = []
    momenta_boson = []
    masses_boson = []

def remove_whitespace(line):
    line = line.replace('\n', '')
    line = line.strip()

    while '  ' in line:
        line = line.replace('  ', ' ')

    return line

def extract_momentum(line):
    cols = line.split(' ')

    e = float(cols[9])
    px = float(cols[6])
    py = float(cols[7])
    pz = float(cols[8])

    return np.array([e, px, py, pz])

def read_lhe(filepath):
    in_event = False
    event_line = 1

    data = Data()
    with open(filepath, mode='r') as file:
        for line in file.readlines():
            line = remove_whitespace(line)

            if line.startswith('#'): continue

            if in_event:
                if event_line == '</event>':
                    in_event = False
                    continue
                elif event_line != 1:
                    if line.startswith('11'):
                        data.momenta_l_minus.append(extract_momentum(line))
                    elif line.startswith('-11'):
                        data.momenta_l_plus.append(extract_momentum(line))
                    elif line.startswith('23'):
                        data.momenta_boson.append(extract_momentum(line))
                        data.masses_boson.append(float(line.split(' ')[10]))
                event_line += 1
            elif line == '<event>':
                in_event = True
                event_line = 1

    return data
