import os
import json
import math

def isfloat(item):

    if item == '0\n':
        return True

    try:
        return float(item)
    except:
        return False

folderpath = os.getcwd()
results_folder = os.path.join(folderpath,"logs")
instances = {}
j = 0
# instance = "cp-test.txt"

for instance in os.listdir(results_folder):
    j += 1
    fname = os.path.join(results_folder, instance)

    opt = fname[-5]

    f = open(fname,"r")
    cont = False
    empty_counter = 0

    instances[j] = {"instance":"","opt":opt,"dual":{"bound":[],"time":[]},"primal":{"bound":[],"time":[]},"gap":{"gap":[],"time":[]},"hit_time_limit":False,"hit_memory_limit":False,"best_dual":0,"best_primal":0,"explored":0,"iterations":0,"infeasible":False,"optimal":False,"time":0,"groups":[],"sections":[], "presolve_time":0,"presolve_remove_rows":0,"presolve_remove_cols":0,"presolved_rows":0,"presolved_cols":0,"presolved_nonzeros":0,"num_cont_vars":0,"num_int_vars":0,"num_bin_vars":0,"original_rows":0,"original_cols":0,"original_nonzeros":0,"presolved_nonzeros":0,"original_cont_vars":0,"original_int_vars":0,"original_bin_vars":0,"memory":0,"best_gap":0,"num_groups":0}
    
    for l in f:
        if "Instance Name: " in l:
            instances[j]["instance"] = l.split(" ")[-1].strip()

        elif "Optimize a model with" in l:
            line = l.strip().split(" ")
            instances[j]["original_rows"] = int(line[4])
            instances[j]["original_cols"] = int(line[6])
            instances[j]["original_nonzeros"] = int(line[9])

        elif "Variable types: " in l:
            line = l.strip().split(" ")
            if instances[j]["original_bin_vars"] == 0:
              instances[j]["original_cont_vars"] = int(line[2])
              instances[j]["original_int_vars"] = int(line[4])
              instances[j]["original_bin_vars"] = int(line[6][1:])
            else:
              instances[j]["num_cont_vars"] = int(line[2])
              instances[j]["num_int_vars"] = int(line[4])
              instances[j]["num_bin_vars"] = int(line[6][1:]) 


        elif "Presolve removed" in l:
            line = l.split(" ")
            instances[j]["presolve_remove_rows"] = int(line[2])
            instances[j]["presolve_remove_cols"] = int(line[5])

        elif "Presolve time" in l:
            line = l.split(" ")
            instances[j]["presolve_time"] = float(line[-1][:-2])

        elif "Presolved: " in l:
            line = l.split(" ")
            instances[j]["presolved_rows"] = int(line[1])
            instances[j]["presolved_cols"] = int(line[3])
            instances[j]["presolved_nonzeros"] = int(line[5])

        elif "Expl Unexpl" in l:
            cont = True

        elif cont:

            if l == '\n':
                empty_counter += 1
                if empty_counter == 2:
                    cont = False
            else:
                line = [i for i in l.split(" ") if i != '']
                time = float(line[-1].strip()[:-1])
                try:
                    gap_index = b = [y for y,x in enumerate(line) if "%" in x][0]
                    new_gap = float(line[gap_index][:-1])
                    new_dual = float(line[gap_index-1])
                    new_primal = float(line[gap_index-2])
                    instances[j]["primal"]["bound"].append(new_primal)
                    instances[j]["primal"]["time"].append(time)
                    instances[j]["dual"]["bound"].append(new_dual)
                    instances[j]["dual"]["time"].append(time)     
                    instances[j]["gap"]["gap"].append(new_gap)
                    instances[j]["gap"]["time"].append(time)
                except:
                    new_dual = float(line[-4])
                    instances[j]["dual"]["bound"].append(new_dual)
                    instances[j]["dual"]["time"].append(time)     
                
                

        elif "Model is infeasible" in l:
            instances[j]["infeasible"] = True
        elif "Memory Used" in l:
            instances[j]["memory"] = float(l.split(" ")[-1].strip())
        elif "Explored" in l:
            line = l.split(" ")
            explored = int(line[1])
            iterations = int(line[3][1:])
            instances[j]["explored"] = explored
            instances[j]["iterations"] = iterations

        elif "Optimal" in l:
            instances[j]["optimal"] = True
        elif "Memory limit reached" in l:
            instances[j]["memory_limit_reached"] = True
        elif "best bound" in l:
            dual = l.split(" ")[5]
            if dual == "-,":
                dual = 0
            else:
                dual = float(dual[:-1])
            instances[j]["best_dual"] = dual
            instances[j]["best_gap"] = l.split(" ")[-1].strip()[:-1] 
        elif "Obj:" in l:
            primal = l.split(" ")[1]
            if primal == "inf":
                primal == 0
            else:
                primal = float(primal)
            
            if math.isinf(primal):
                primal = 0
            instances[j]["best_primal"] = primal
        elif "Groups: " in l:
            instances[j]["num_groups"] = int(l.split(" ")[-1].strip())
        elif "Time:" in l:
            timeused = float(l.split(" ")[-1])
            instances[j]["time"] = timeused
        elif "Grouped elements" in l:
            groups = eval(l.split("Grouped elements:  ")[-1])
            instances[j]["groups"] = groups
        elif "Element sections" in l:
            sections = eval(l.split("Element sections:  ")[-1])
            instances[j]["sections"] = sections
            
results_fname = os.path.join(folderpath,"results","MIP-results.json")

with open(results_fname,'w') as res:
    json.dump(instances,res)