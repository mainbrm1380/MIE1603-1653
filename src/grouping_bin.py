import gurobipy as gp
from gurobipy import GRB
from instance_generation import generate_instance, generate_section_costs

# test_data = [[4,2,1,1,0],[1,1,1,1,0], [3,2,1,1,1],[4,3,3,3,1],[0,0,2,2,0],[8,8,7,6,3]] #e_ij = 1 if col i exists at level j, bottom to top
#INSTANCE GENERATION
num_Sections = 20
num_Cols = 10
num_Levels = 10
min_section_cost = 100
max_section_cost = 10000

test_data = generate_instance(num_Sections,num_Cols,num_Levels,min_s_in_c=1,max_s_in_c=4)
SectionCost = generate_section_costs(num_Sections,min_section_cost,max_section_cost)


n_cols = len(test_data)
n_levels = len(test_data[0])
max_min_groups = 15 #min(n_cols,n_levels) #TODO: NEED NEW HEURISTIC FOR THIS
M = n_levels + 2
M_sections = max([max(i) for i in test_data]) + 1

Xgcl = [(g,c,l) for l in range(n_levels) for c in range(n_cols) for g in range(max_min_groups)]
Gs = [(g,s) for g in range(max_min_groups) for s in range(M_sections)]
Es = [(c,l,s) for s in range(M_sections) for l in range(n_levels) for c in range(n_cols)]

#c,l
S = {(i,j): test_data[i][j] for j in range(n_levels) for i in range(n_cols)}
S_bin = {(i,j,s): 1 if s >= test_data[i][j] else 0 for s in range(M_sections) for j in range(n_levels) for i in range(n_cols)}
E = {(i,j): 1 if test_data[i][j] > 0 else 0 for j in range(n_levels) for i in range(n_cols)}


# Create a new model
model = gp.Model("Grouping-Optimization")

# # Create variables

x = model.addVars(Xgcl, vtype = GRB.BINARY, name="x")

GroupCost = 0
# SectionCost = [0,100,400,800,1600,3200,6400,12800,14000]
sections = set(i for i in range(len(SectionCost)))

group_exists = model.addVars(max_min_groups, vtype = GRB.BINARY, name="group_exists")
column_in_group = model.addVars(max_min_groups, n_cols, vtype = GRB.BINARY, name="col_in_group")
level_in_group = model.addVars(max_min_groups, n_levels, vtype = GRB.BINARY, name="level_in_group")

# can these be continuous? or must be int?
group_lower_bound = model.addVars(max_min_groups, vtype = GRB.INTEGER, name="group_lb", lb = 0, ub = n_levels)
group_upper_bound = model.addVars(max_min_groups, vtype = GRB.INTEGER, name="group_ub", lb = 0, ub = n_levels)
group_level_range = model.addVars(max_min_groups, vtype = GRB.INTEGER, name="group_range", lb = 0, ub = n_levels)

#section size variables
group_section = model.addVars(Gs, vtype = GRB.BINARY, name="group_section")
element_section = model.addVars(Es, vtype=GRB.BINARY, name="element_section")
# Zs = model.addVars(Xgcl, vtype = GRB.BINARY, name = "Zs") #if group section = element section

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
    model.addConstr(gp.quicksum(s * element_section[c,l,s] for s in range(M_sections)) >= S[c,l])
    # model.addConstr(gp.quicksum(element_section[c,l,s] for s in range(M_sections)) == 1)
    
  for l in range(1,n_levels):
    if E[c,l] == 1:
      model.addConstr(gp.quicksum(s*element_section[c,l,s] for s in range(M_sections)) >= gp.quicksum((s*element_section[c,l-1,s] for s in range(M_sections))))

for g in range(max_min_groups):
  #range of group g = upper bound - lower bound 
  model.addConstr(group_level_range[g] == group_upper_bound[g] - group_lower_bound[g])

  #1 section per group
  model.addConstr(gp.quicksum(group_section[g,s] for s in range(M_sections)) == 1)

  for c in range(n_cols):
    for l in range(n_levels):

      #if element is in a group, it's section is at least the group's
      #if not in the group, it's greater than (at most) 0
      #TODO: these are problems constraints.
      model.addConstr(gp.quicksum(s * element_section[c,l,s] for s in range(M_sections)) >= gp.quicksum(s * group_section[g,s] for s in range(M_sections)) - (M_sections+1)*(1-x[g,c,l]))
      model.addConstr(gp.quicksum(s * element_section[c,l,s] for s in range(M_sections)) <= gp.quicksum(s * group_section[g,s] for s in range(M_sections)) + (M_sections+1)*(1-x[g,c,l]))
      
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
    model.addConstr(M*Zl[g,l] >= l-group_lower_bound[g]+1)
    model.addConstr(M*(1-Zl[g,l]) >= group_lower_bound[g]-l)

    #set Zu
    model.addConstr(M*Zu[g,l] >= group_upper_bound[g]-l+1)
    model.addConstr(M*(1-Zu[g,l]) >= l-group_upper_bound[g])

    #if a level is within lower/upper bound, it's in the group 
    model.addConstr(1+level_in_group[g,l] >= Zu[g,l]+Zl[g,l])

    for c in range(n_cols):
      #if column and level are in the group, so is the element
      #todo: lazily constraint??
      model.addConstr(column_in_group[g,c]+level_in_group[g,l] <= 1 + x[g,c,l])

model.setObjective(gp.quicksum(element_section[c,l,s] * SectionCost[s] for s in range(M_sections) for c in range(n_cols) for l in range(n_levels)) + GroupCost*gp.quicksum(group_exists),GRB.MINIMIZE)

# model.write("group-optim-toy.lp")

# Optimize model
model.optimize()

ns = 0
for v in model.getVars():
    if "group_exists" in v.VarName and v.X > 0.5:
      ns += 1
        # print(f"{v.VarName}, {v.x:g}") #{v.X:g}")

print(f"Obj: {model.ObjVal:g}")
print(f"Time: {model.Runtime:g}")

print("Groups: ", ns)

print(test_data)
print(SectionCost)