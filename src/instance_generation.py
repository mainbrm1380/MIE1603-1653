import random
import json

def generate_column_by_number_of_sections(L, S, s_in_c):
    data = [0]
    while max(data) == 0:
        sections = sorted(random.choices(S, k=s_in_c))
        level_breakpoints = [0] + sorted(random.choices(L, k = s_in_c+1)) #+1 because there can be 0 on top and bottom
        data = [0 * level_breakpoints[0]]
        for i in range(s_in_c):
            data += [sections[i]] * (level_breakpoints[i+1]-level_breakpoints[i])
        data += [0] * (len(L) - len(data))

    # if len(data) != len(L):
    #     pass

    return data

def generate_column_by_levels_per_section(L,S,min_l_per_s, max_l_per_s):
    data = [0]
    while max(data) == 0:
        level_lengths_extra = [random.randint(min_l_per_s,max_l_per_s) for i in range(int(L/min_l_per_s))]
        level_lengths = []
        i = 0
        s = 0
        while s < L:
            level_lengths.append(level_lengths_extra[i])
            s += level_lengths_extra[i]
            i += 1
        
        first_level_zero = random.randint(0,1)
        last_level_zero = random.randint(0,1)
        
        sections = sorted(random.choices(S,k = len(level_lengths) - first_level_zero - last_level_zero))

        data = [0] * level_lengths[0]*first_level_zero

        for i in range(first_level_zero,len(sections)+first_level_zero):
            data += [sections[i-first_level_zero]] * level_lengths[i]

        data += [0] * (L-len(data))*last_level_zero
        data = data[:L]

    # if len(data) != L:
    #     pass

    return data

def generate_instance(S, C, L, opt, min_s_in_c = 0, max_s_in_c = 0, min_l_per_s = 0, max_l_per_s = 0):
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


def generate_data_set(classes, num_per_gen_type):
    data = {}
    instance_num = 1
    for cl in classes.keys():
        for gen_type in ["num_sections", "num_levels"]:
            for i in range(num_per_gen_type):
                instance = {"name":"","class": cl, "generation_type": gen_type, "columns": [], "section_costs": []}
                args = {"opt":gen_type}
                args.update(classes[cl])
                instance["name"] = cl + "_" + gen_type.split("_")[-1] + "_" + str(i+1)
                instance["columns"] = generate_instance(**args)
                instance["section_costs"] = generate_section_costs(classes[cl]["S"], 1000, 5000)
                data[instance_num] = instance
                instance_num += 1
    return data
                
def main():
    classes = {"toy": {"S":5,"C":5,"L":5,"min_s_in_c":1,"max_s_in_c":5,"min_l_per_s":1,"max_l_per_s":5}, 
                "industrial": {"S":10,"C":25,"L":5,"min_s_in_c":1,"max_s_in_c":5,"min_l_per_s":1,"max_l_per_s":5}, 
                "mid-rise": {"S":20,"C":30,"L":10,"min_s_in_c":2,"max_s_in_c":5,"min_l_per_s":2,"max_l_per_s":5}, 
                "ltc": {"S":30,"C":50,"L":6,"min_s_in_c":2,"max_s_in_c":6,"min_l_per_s":1,"max_l_per_s":3},
                "condo": {"S":30,"C":80,"L":60,"min_s_in_c":6,"max_s_in_c":15,"min_l_per_s":4,"max_l_per_s":10},
                "hospital": {"S":30,"C":200,"L":25,"min_s_in_c":5,"max_s_in_c":13,"min_l_per_s":2,"max_l_per_s":5}}
    data = generate_data_set(classes, 5)
    with open("data.json","w") as f:
        json.dump(data,f)
    
main()