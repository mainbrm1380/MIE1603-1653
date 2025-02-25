import gurobipy as gp
from gurobipy import GRB

test_data = [[2,1,1,0],[1,1,1,0]] #e_ij = 1 if col i exists at level j, bottom to top


n_cols = len(test_data)
n_levels = len(test_data[0])
max_min_groups = 5 #min(n_cols,n_levels)
M = n_levels + 2

Xgcl = [(g,c,l) for l in range(n_levels) for c in range(n_cols) for g in range(max_min_groups)]

S = {(i,j): test_data[i][j] for j in range(n_levels) for i in range(n_cols)}
E = {(i,j): 1 if test_data[i][j] > 0 else 0 for j in range(n_levels) for i in range(n_cols)}
sections = set(l for l in S.values())

# Create a new model
model = gp.Model("Grouping-Optimization")

# # Create variables

x = model.addVars(Xgcl, vtype = GRB.BINARY, name="x")

group_exists = model.addVars(max_min_groups, vtype = GRB.BINARY, name="group")
column_in_group = model.addVars(max_min_groups, n_cols, vtype = GRB.BINARY, name="col_in_group")
level_in_group = model.addVars(max_min_groups, n_levels, vtype = GRB.BINARY, name="level_in_group")

# can these be continuous? or must be int?
group_lower_bound = model.addVars(max_min_groups, vtype = GRB.INTEGER, name="group_lb", lb = 0, ub = n_levels)
group_upper_bound = model.addVars(max_min_groups, vtype = GRB.INTEGER, name="group_ub", lb = 0, ub = n_levels)
group_level_range = model.addVars(max_min_groups, vtype = GRB.INTEGER, name="group_range", lb = 0, ub = n_levels)

#section size variables
#done as binary? or as integer? try both: use worse/better for ablation
group_section = model.addVars(max_min_groups, vtype = GRB.INTEGER, name="group_section", lb = 0, ub = max(sections))
element_section = model.addVars(S.keys(), vtype=GRB.INTEGER, name="element_section", lb=0, ub = max(sections))

#each group is same section
#group section has to be >= section of element
#take min of cost_of_section * num of elements with that section
#^either count # of sections in each group, and linearize cost_of_section var * # sections
#that is for integer version, binary version has binary var for which section a grou phas, times its cost, times number of elements
# or has bin var for section of each element, and multiple by cost of that section.
# 4 variations?? 

Zu = model.addVars(max_min_groups, n_levels, vtype = GRB.BINARY, name="Zu")
Zl = model.addVars(max_min_groups, n_levels, vtype = GRB.BINARY, name="Zl")

for c in range(n_cols):
  for l in range(n_levels):
    #sum of Xgcl over all groups must = if col at that level exists
    model.addConstr(gp.quicksum(x[g,c,l] for g in range(max_min_groups)) == E[c,l])

for g in range(max_min_groups):
  #range of group g = upper bound - lower bound 
  model.addConstr(group_level_range[g] == group_upper_bound[g] - group_lower_bound[g])

  for c in range(n_cols):
    for l in range(n_levels):
      #if element is in group, that column is in the group 
      model.addConstr(column_in_group[g,c] >= x[g,c,l])

      #if element is in group, that level is in the group
      model.addConstr(level_in_group[g,l] >= x[g,c,l]) 

      #if element is in that group, that group exists
      model.addConstr(group_exists[g] >= x[g,c,l])
  
  for g in range(max_min_groups):
    for l in range(n_levels):

      #calculate upper, lower level
      #NOTE: need to min(Ug) and max(Lg) for this to work
      #NOTE: in julia, everything is 1 based, might need to adjust formulation to make everything work as 0-based
      model.addConstr(l*level_in_group[g,l] <= group_upper_bound[g])
      model.addConstr(l*level_in_group[g,l] + M*(1-level_in_group[g,l]) >= group_lower_bound[g])     

      #set Zl
      model.addConstr(M*Zl[g,l] >= (l+1)-group_lower_bound[g]+1)
      model.addConstr(M*(1-Zl[g,l]) >= group_lower_bound[g]-(l+1))

      #set Zu
      model.addConstr(M*Zu[g,l] >= group_upper_bound[g]-(l+1)+1)
      model.addConstr(M*(1-Zu[g,l]) >= (l+1)-group_upper_bound[g])

			#if a level is within lower/upper bound, it's in the group 
      model.addConstr(1+level_in_group[g,l] >= Zu[g,l]+Zl[g,l])

      for c in range(n_cols):
			  #if column and level are in the group, so is the element
        #todo: lazily constraint??
        model.addConstr(column_in_group[g,c]+level_in_group[g,l] <= 1 + x[g,c,l])

model.setObjective(gp.quicksum(group_lower_bound) - gp.quicksum(group_upper_bound) + 100*gp.quicksum(group_exists),GRB.MINIMIZE)

model.write("group-optim-toy.lp")

# Optimize model
model.optimize()

# for v in model.getVars():
#     if v.X > 0.5:
#         print(f"{v.VarName}") #{v.X:g}")

print(f"Obj: {model.ObjVal:g}")
print(f"Time: {model.Runtime:g}")