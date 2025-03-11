import random

def generate_column_by_number_of_sections(L, S, s_in_c):
    data = [0]
    while max(data) == 0:
        sections = sorted(random.choices(S, k=s_in_c))
        level_breakpoints = [0] + sorted(random.choices(L, k = s_in_c+1)) #+1 because there can be 0 on top and bottom
        data = [0 * level_breakpoints[0]]
        for i in range(s_in_c):
            data += [sections[i]] * (level_breakpoints[i+1]-level_breakpoints[i])
        data += [0] * (len(L) - len(data))
    return data

def generate_column_by_levels_per_section(L,S,min_l_per_s, max_l_per_s):
    data = [0]
    while max(data) == 0:
        level_lengths_extra = [random.randint(min_l_per_s,max_l_per_s) for i in range(int(L/min_l_per_s))]
        level_lengths = []
        i = 0
        s = 0
        while i < len(level_lengths_extra) and s < L:
            level_lengths.append(level_lengths_extra[i])
            s += level_lengths_extra[i]
            i += 1
        
        first_level_zero = random.randint(0,1)
        last_level_zero = random.randint(0,1)
        
        sections = sorted(random.choices(S,k = len(level_lengths) - first_level_zero - last_level_zero))

        data = [0] * level_lengths[0]*first_level_zero

        for i in range(len(sections)):
            data += [sections[i]] * level_lengths[i]

        data += [0] * (L-len(data))*last_level_zero
        data = data[:L]

    return data

def generate_instance(S, C, L, opt = "num_sections", min_s_in_c = 0, max_s_in_c = 0, min_l_per_s = 0, max_l_per_s = 0):
    #S = number of sections
    #C = number of columns
    #L = number of levels
    #opt = option if by sections or by levels
    #min/max_s_in_c: min/max number of sections in a column
    #min/max_l_per_s: min/max number of levels per section

    sections = list(range(1,S+1))
    levels = list(range(L))

    data = []
    if opt == "num_sections":
        for c in range(C):
            s_in_c = random.randint(min_s_in_c, max_s_in_c)
            data.append(generate_column_by_number_of_sections(levels, sections, s_in_c))
    elif opt == "num_levels":
        for c in range(C):
            data.append(generate_column_by_levels_per_section(L, sections, min_l_per_s, max_l_per_s))
    
    return data

def generate_section_costs(S,min_cost, max_cost):
    return sorted([0] + [random.randint(min_cost,max_cost) for i in range(S)])
    