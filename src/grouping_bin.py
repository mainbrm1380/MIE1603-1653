import gurobipy as gp
from gurobipy import GRB
import json
import psutil
import os
import sys
import numpy as np
import networkx as nx

import matplotlib.pyplot as plt
# import matplotlib.pyplot as plt
# import matplotlib.colors as mcolors
# import seaborn as sns

process = psutil.Process()

#1 ADD CONSTAINTS THAT CERTAIN ELEMENTS CAN'T BE IN SAME GROUP

#2 REMOVE SYMMETY

def exlusive_groups(col_exists, ncol, nlvl):
  infeasible_dict = dict()
  for (col,lvl) in col_exists:
    if col_exists[col,lvl] == 1:
      infeasible = set()
      for (col2,lvl2) in col_exists:
         if col2 != col and lvl2 != lvl and col_exists[col2, lvl2] == 1:
            if col_exists[col,lvl2] == 0 or col_exists[col2, lvl] == 0:
                infeasible.add((col2,lvl2))
      infeasible_dict[col,lvl] = infeasible
  
  return infeasible_dict

def exlusive_groups_by_col(col_exists, ncol, nlvl):
  infeasible_dict = dict()
  for (col,lvl) in col_exists:
    if col_exists[col,lvl] == 1:
      infeasible = set()
      for col2 in range(ncol):
        if col_exists[col2,lvl] == 0:
           infeasible.add(col2)
      infeasible_dict[col,lvl] = infeasible
  return infeasible_dict

def break_partial_symmetry(infeasible_dict):
  group_dict = dict()
  ind_by_len = sorted(infeasible_dict.keys(), key= lambda x: len(infeasible_dict[x]), reverse = True)
  print(ind_by_len)

  sim_sets = set()
  list_of_sets = []
  list_of_indices = []
  list_of_largest_sets = []
  for ind,st in infeasible_dict.items():
    if st not in list_of_sets and len(st) > 0:
      list_of_sets.append(st)
      list_of_indices.append([ind])
    elif len(st) > 0:
      list_of_indices[[i for i,j in enumerate(list_of_sets) if j == st][0]].append(ind)

  #pick largest set
  #go through indices that point to that set



def get_heuristic_group_ub(columns):
  col_ub = 0
  level_ub = 0
  for col in columns:
    col_ub += len(set([elem for elem in col if elem > 0]))
  for lvl in range(len(columns[0])):
    level_ub += len(set([col[lvl] for col in columns if col[lvl] > 0]))
  return min(col_ub,level_ub), max(col_ub,level_ub)


def lazy_callback(model, where):
    if where == GRB.Callback.MIPSOL:  # Check if a new solution is found
        # Retrieve the values of decision variables at the current solution
        column_vals = model.cbGetSolution(column_in_group)
        level_vals = model.cbGetSolution(level_in_group)
        x_vals = model.cbGetSolution(x)

        # Iterate through indices and add violated constraints
        for g in range(max_min_groups):
            for l in range(n_levels):
                for c in range(n_cols):
                    lhs = column_vals[g, c] + level_vals[g, l]  # Left-hand side
                    rhs = 1 + x_vals[g, c, l]  # Right-hand side
                    if lhs > rhs + 1e-6:  # Constraint is violated
                        model.cbLazy(column_in_group[g, c] + level_in_group[g, l] <= 1 + x[g, c, l])

if __name__ == "__main__":
 
  print("===INSTANCE START")

  
  # script, instance, opt = sys.argv
  instance = "2"

  opt = int(2)
  # instance = "4"
  # opt = "7"

  #opt: 1 - only int, 2 - cont, 4 - lazy,  6 - cont+lazy
  #opt: 3 - only bin, 5 - bin+cont, 7: bin+lazy, 8:bin+lazy+cont


  folderpath = os.getcwd()
  data_path = os.path.join(folderpath,"data.json")

  with open(data_path, 'r') as file:
    data = json.load(file)

  i = instance
  i_name = data[i]["name"]
  print(f"Instance Name: {instance}-{i_name}")
  test_data = data[i]["columns"]
  SectionCost = data[i]["section_costs"]
  n_cols = len(test_data)
  n_levels = len(test_data[0])

  max_min_groups, max_max_groups = get_heuristic_group_ub(test_data) #min(n_cols,n_levels) 
  
  min_cost = sum([SectionCost[i] for col in test_data for i in col])
  
  GroupCost = round(min_cost/max_max_groups)
  max_min_groups = min(100, max_min_groups)
  n_cols = len(test_data)
  n_levels = len(test_data[0])
  
  M = n_levels + 2
  M_sections = max([max(i) for i in test_data]) + 1

  Xgcl = [(g,c,l) for l in range(n_levels) for c in range(n_cols) for g in range(max_min_groups)]
  Gs = [(g,s) for g in range(max_min_groups) for s in range(M_sections)]
  Es = [(c,l,s) for s in range(M_sections) for l in range(n_levels) for c in range(n_cols)]

  #c,l
  S = {(i,j): test_data[i][j] for j in range(n_levels) for i in range(n_cols)}
  S_bin = {(i,j,s): 1 if s >= test_data[i][j] else 0 for s in range(M_sections) for j in range(n_levels) for i in range(n_cols)}
  E = {(i,j): 1 if test_data[i][j] > 0 else 0 for j in range(n_levels) for i in range(n_cols)}

  infeasible_dict = exlusive_groups(E,n_cols,n_levels)
  infeasible_dict_cols = exlusive_groups_by_col(E,n_cols,n_levels)

  #COLORING FOR INITIAL GROUPS

  G = nx.Graph(infeasible_dict)
  
  nx.draw(G,with_labels=True,node_size=100)

  plt.show()

  # break_partial_symmetry(infeasible_dict)

  # # Create a new model
  # model = gp.Model("Grouping-Optimization")

  # # # Create variables
  
  # # print("add vars")
  # x = model.addVars(Xgcl, vtype = GRB.BINARY, name="x")

  # # SectionCost = [0,100,400,800,1600,3200,6400,12800,14000]
  # sections = set(i for i in range(len(SectionCost)))

  # group_exists = model.addVars(max_min_groups, vtype = GRB.BINARY, name="group_exists")
  # column_in_group = model.addVars(max_min_groups, n_cols, vtype = GRB.BINARY, name="col_in_group")
  # level_in_group = model.addVars(max_min_groups, n_levels, vtype = GRB.BINARY, name="level_in_group")

  # # can these be continuous? or must be int?

  # #section size variables
  # group_section = model.addVars(Gs, vtype = GRB.BINARY, name="group_section")
  # element_section = model.addVars(Es, vtype=GRB.BINARY, name="element_section")
  # # Zs = model.addVars(Xgcl, vtype = GRB.BINARY, name = "Zs") #if group section = element section

  # #each group is same section
  # #group section has to be >= section of element
  # #take min of cost_of_section * num of elements with that section
  # #^either count # of sections in each group, and linearize cost_of_section var * # sections
  # #that is for integer version, binary version has binary var for which section a grou phas, times its cost, times number of elements
  # # or has bin var for section of each element, and multiple by cost of that section.
  # # 4 variations?? 

  # if opt in [2,5,6,8]:
  #   group_lower_bound = model.addVars(max_min_groups, vtype = GRB.CONTINUOUS, name="group_lb", lb = 0, ub = n_levels)
  #   group_upper_bound = model.addVars(max_min_groups, vtype = GRB.CONTINUOUS, name="group_ub", lb = 0, ub = n_levels)
  #   group_level_range = model.addVars(max_min_groups, vtype = GRB.CONTINUOUS, name="group_range", lb = 0, ub = n_levels)
  #   Zu = model.addVars(max_min_groups, n_levels, vtype = GRB.CONTINUOUS, name="Zu")
  #   Zl = model.addVars(max_min_groups, n_levels, vtype = GRB.CONTINUOUS, name="Zl")
  # else:
  #   group_lower_bound = model.addVars(max_min_groups, vtype = GRB.INTEGER, name="group_lb", lb = 0, ub = n_levels)
  #   group_upper_bound = model.addVars(max_min_groups, vtype = GRB.INTEGER, name="group_ub", lb = 0, ub = n_levels)
  #   group_level_range = model.addVars(max_min_groups, vtype = GRB.INTEGER, name="group_range", lb = 0, ub = n_levels)
  #   Zu = model.addVars(max_min_groups, n_levels, vtype = GRB.BINARY, name="Zu")
  #   Zl = model.addVars(max_min_groups, n_levels, vtype = GRB.BINARY, name="Zl")

  # if opt > 0: #CHANGE
  #   for [c,l], inf_set in infeasible_dict.items():
  #     if inf_set != {}:
  #       M_grp = len(inf_set)
  #       model.addConstrs(M_grp*x[g,c,l] + gp.quicksum(x[g,c2,l2] for (c2,l2) in inf_set) <= M_grp for g in range(max_min_groups))

  # if opt > 0: #CHANGE
  #   for [c,l], inf_set in infeasible_dict_cols.items():
  #     if inf_set != {}:
  #       M_grp = len(inf_set)
  #       model.addConstrs(M_grp*x[g,c,l] + gp.quicksum(column_in_group[g,c2] for c2 in inf_set) <= M_grp for g in range(max_min_groups))

  # # print("add elem constrs")
  # for c in range(n_cols):
  #   for l in range(n_levels):

  #     model.addConstr(gp.quicksum(element_section[c,l,s] for s in range(M_sections)) == E[c,l])

  #     #sum of Xgcl over all groups must = if col at that level exists
  #     model.addConstr(gp.quicksum(x[g,c,l] for g in range(max_min_groups)) == E[c,l])
  #     model.addConstr(gp.quicksum(s * element_section[c,l,s] for s in range(M_sections)) >= S[c,l])
  #     # model.addConstr(gp.quicksum(element_section[c,l,s] for s in range(M_sections)) == 1)
      
  #   for l in range(1,n_levels):
  #     if E[c,l] == 1 and E[c,l-1] == 1:
  #       model.addConstr(gp.quicksum(s*element_section[c,l,s] for s in range(M_sections)) <= gp.quicksum((s*element_section[c,l-1,s] for s in range(M_sections))))

  # # print("add grp constrs")
  # for g in range(max_min_groups):
  #   print(g, g/max_min_groups) 
  #   #range of group g = upper bound - lower bound 
  #   model.addConstr(group_level_range[g] == group_upper_bound[g] - group_lower_bound[g])

  #   #1 section per group
  #   model.addConstr(gp.quicksum(group_section[g,s] for s in range(M_sections)) == 1)

  #   for c in range(n_cols):
  #     for l in range(n_levels):

  #       #if element is in a group, it's section is at least the group's
  #       #if not in the group, it's greater than (at most) 0
  #       #TODO: these are problems constraints.
  #       model.addConstr(gp.quicksum(s * element_section[c,l,s] for s in range(M_sections)) >= gp.quicksum(s * group_section[g,s] for s in range(M_sections)) - (M_sections+1)*(1-x[g,c,l]))
  #       model.addConstr(gp.quicksum(s * element_section[c,l,s] for s in range(M_sections)) <= gp.quicksum(s * group_section[g,s] for s in range(M_sections)) + (M_sections+1)*(1-x[g,c,l]))
        
  #       #if element is in group, that column is in the group 
  #       model.addConstr(column_in_group[g,c] >= x[g,c,l])

  #       #if element is in group, that level is in the group
  #       model.addConstr(level_in_group[g,l] >= x[g,c,l]) 

  #       #if element is in that group, that group exists
  #       model.addConstr(group_exists[g] >= x[g,c,l])

  # # print("add lvl constrs")
  # for g in range(max_min_groups):
  #   for l in range(n_levels):

  #     #calculate upper, lower level
  #     #NOTE: need to min(Ug) and max(Lg) for this to work
  #     #NOTE: in julia, everything is 1 based, might need to adjust formulation to make everything work as 0-based
  #     model.addConstr(l*level_in_group[g,l] <= group_upper_bound[g])
  #     model.addConstr(l*level_in_group[g,l] + M*(1-level_in_group[g,l]) >= group_lower_bound[g])     

  #     #set Zl
  #     model.addConstr(M*Zl[g,l] >= l-group_lower_bound[g]+1)
  #     model.addConstr(M*(1-Zl[g,l]) >= group_lower_bound[g]-l)

  #     #set Zu
  #     model.addConstr(M*Zu[g,l] >= group_upper_bound[g]-l+1)
  #     model.addConstr(M*(1-Zu[g,l]) >= l-group_upper_bound[g])

  #     #if a level is within lower/upper bound, it's in the group 
  #     model.addConstr(1+level_in_group[g,l] >= Zu[g,l]+Zl[g,l])

  #     # for c in range(n_cols):
  #     #   #if column and level are in the group, so is the element
  #     #   #todo: lazily constraint??
  #     #   model.addConstr(column_in_group[g,c]+level_in_group[g,l] <= 1 + x[g,c,l])

  # model.setObjective(gp.quicksum(element_section[c,l,s] * SectionCost[s] for s in range(M_sections) for c in range(n_cols) for l in range(n_levels)) + GroupCost*gp.quicksum(group_exists) - gp.quicksum(group_upper_bound) + gp.quicksum(group_lower_bound),GRB.MINIMIZE)

  # model.setParam('TimeLimit', 1800)
  # model.setParam('SoftMemLimit', 10)

  # if opt in [4,6,7,8]:
  #   model.Params.LazyConstraints = 1
  #   model.optimize(lazy_callback)
  # else:
  #   for g in range(max_min_groups):
  #     for l in range(n_levels):
  #       for c in range(n_cols):
  #       #if column and level are in the group, so is the element
  #         model.addConstr(column_in_group[g,c]+level_in_group[g,l] <= 1 + x[g,c,l])
  #   model.optimize()

  # # return model, x, element_section, GroupCost, column_in_group, level_in_group, max_min_groups

  # # model.write("group-optim-toy.lp")

  # # Optimize model

  # ns = 0
  # for v in model.getVars():
  #     if "group_exists" in v.VarName and v.X > 0.5:
  #       ns += 1

  # print(f"Obj: {model.ObjVal:g}")
  # print(f"Time: {model.Runtime:g}")
  # print("Memory Used (MiB): {}".format(round(process.memory_info().rss / 1024 ** 2,2)))
  # print("Groups: ", ns)

  # grouped_elements = np.full((n_cols, n_levels), 0)  # -1 as default (if element doesn't exist)

  # for g, i, j in x.keys():
  #     if x[g, i, j].X > 0.5:  # Check if x[g, i, j] is active
  #         grouped_elements[i, j] = g+1

  # print("Original columns: ", test_data)
  # print("Section costs: ", SectionCost)
  # print("Group cost: ", GroupCost)

  # print("Grouped elements: ", grouped_elements.tolist())

  # section_of_elements = np.full((n_cols, n_levels), 0)

  # for c, l, s in element_section.keys():
  #     if element_section[c, l, s].X > 0.5:
  #       section_of_elements[c, l] = s

  # print("Element sections: ", section_of_elements.tolist())

  # heatmap(grouped_elements, section_of_elements)

  # print("---ALGORITHM END")



