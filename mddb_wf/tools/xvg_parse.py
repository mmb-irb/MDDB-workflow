# Read and parse a xvg files
# xvg file example:
# # This is a comment line
# @ This is a comment line
#    0.0000000    0.3644627
#   10.0000000    0.3536768
#   20.0000000    0.3509805

# Columns is a list with the name of each column
# It returns a dict with as many entries as columns specifided
# Each entry has the column name as key and the mined data (list) as value
# WARNING: If the xvg file has less columns than specified then the leftover entries will have empty lists
# WARNING: If the xvg file has more columns than specified then it will fail
def xvg_parse (filename : str, columns : list) -> dict:
    data = [ [] for column in columns ]
    # Read the specified file line per line
    with open(filename, 'r') as file:
        lines = list(file)
        for line in lines:
            # Skip comment lines
            first_character = line[0]
            if first_character in ['#','@']:
                continue
            # Useful lines are splitted by spaces
            # Splits are saved in as many columns as specified
            splits = line.split()
            for s, split in enumerate(splits):
                data[s].append(float(split))
    # Create the formatted dict to be returned
    results = {}
    for c, column in enumerate(columns):
        results[column] = data[c]
    return results