# Read and parse a xvg files
# xvg file example:
# # This is a comment line
# @ This is a comment line
#    0.0000000    0.3644627
#   10.0000000    0.3536768
#   20.0000000    0.3509805

# 2 columns
def xvg_parse (filename : str):
    times = []
    values = []
    # Read the specified file line per line
    with open(filename, 'r') as file:
        lines = list(file)
        for line in lines:
            # Skip comment lines
            first_character = line[0]
            if first_character in ['#','@']:
                continue
            # Useful lines are splitted by spaces
            # Splits are saved in 2 columns: times and values
            [time, value] = line.split()
            times.append(float(time))
            values.append(float(value))
    return { 'times': times, 'values': values }

# 3 columns
def xvg_parse_3c (filename : str):
    c1 = []
    c2 = []
    c3 = []
    # Read the specified file line per line
    with open(filename, 'r') as file:
        lines = list(file)
        for line in lines:
            # Skip comment lines
            first_character = line[0]
            if first_character in ['#','@']:
                continue
            # Useful lines are splitted by spaces
            # Splits are saved in 2 columns: times and values
            [a, b, c] = line.split()
            c1.append(float(a))
            c2.append(float(b))
            c3.append(float(c))
    return { 'c1': c1, 'c2': c2, 'c3': c3 }