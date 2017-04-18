import string

__all__ = ['dna2nM_convert', 'echo_round', 'floatify',
           'process_column_argument']

def dna2nM_convert(dnaconc, dnalength):
    '''
    Converts DNA ng/ul to nM

    DNA_conc(ng/ul) * 1e6 (uL/L) * (bp*mol)/660g * 1/dna_length(bp = DNA (nM)
    Double stranded DNA is 660g/(bp*mol), single stranded DNA is 330g/(bp*mol)
    '''
    return (dnaconc*1e6)/(660*dnalength)

def echo_round(x, base = 25):
    '''
    Round a volume to an exact number the Echo can pipette

    Arguments:
        x - Desired volume.
        base - Size of an Echo drop, in units of the volume (default nL).

    Returns: x rounded to the nearest multiple of base.
    '''
    return int(base * round(float(x)/base))

def floatify(element):
    '''
    Convert a string to a float, if possible; otherwise, return the string
    right back. Code ripped off code from Jacob Gabrielson's answer on
    StackOverflow thread 736043.

    Empty strings are converted to "0" to make vector addition work on relevant
    ranges.
    '''
    if not type(element) == str:
        raise ValueError("Why are you trying to floatify a " + type(element) + \
                         "?")
    if element == "":
        return 0.0
    partition=element.partition('.')
    if element.isdigit() or \
       (partition[0].isdigit() and partition[1]=='.' and \
        partition[2].isdigit()) \
        or \
        (partition[0]=='' and partition[1]=='.' and partition[2].isdigit()) \
        or \
        (partition[0].isdigit() and partition[1]=='.' and partition[2]==''):
        return float(element)
    else:
        return element

def process_column_argument(col):
    '''
    Convenience function for processing a column-name argument. Converts
    from either a string ("C") or a zero-indexed int (2) to a zero-indexed
    int. Nonetype arguments are returned as None.
    '''
    if col == None:
        return None
    if type(col) == str:
        if col == "":
            raise ValueError("Column argument can't be an empty string.")
        upper_col = col.upper()
        col_num = 0
        while len(upper_col) > 0:
            col_num *= 26
            col_num += string.ascii_uppercase.find(upper_col[-1]) + 1
            upper_col = upper_col[:-1]
        #Remember, it's zero-indexed!
        col_num -= 1
        return col_num
    elif type(col) == int:
        return col
    else:
        raise TypeError("Column argument must be a string ('C') or a " +\
                        "zero-indexed int (2)")
