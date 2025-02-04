using JuMP #for optimization
import HiGHS #for optimization
using Dates 
using ODBC, DataFrames #to connect to Access

global PHI_C = 0.65
global PHI_S = 0.85

Base.@kwdef mutable struct Subgroup

	id::Int
	sgName::String
	sgID::String
	columns::Union{String, Vector{Any}}
	stories::Union{String, Vector{Any}}
	mat::String = "Concrete"
	shp::String
	b::Float64
	h::Float64
	fc::Float64
	Ec::Float64
	Gc::Float64
	fy::Float64 = 400.0
	Lu::Float64 #TODO
	Cm::Float64 = 0.0
	k::Float64 = 0.9
	Pf::Float64 #TODO
	Mfx::Float64 = 0.275
	Mfy::Float64 = 0.275
	rho::Float64 = 0.01
	DL::Float64 = 100.0
	LL::Float64 = 40.0
	SNOW::Float64 = 0.0
	SeismicTie::Bool = false
	IsValid::Bool = true
	DateCreated::String
	LastUpdated::String
	updSCon::Bool = false
	updspCol::Bool = false
	updSBC::Bool = false
	minRhoPct::Float64 = 0.01
	UtilPct::Float64 = 1.0
	MultilineComment::String = ""
	VerticalSpliceType::String = ""

end

function readColumns(columnTable, storyDict)

	typeDict = Dict{String, Int64}()
	groupDict = Dict{String, Dict{String, Any}}()
	colDict = Dict{Any,Dict{String, Any}}()

	#get types 
	#convert colID to colNum
	#convert storyID to storyNum
	#get area of each col 
	#get rho of each col 

	#columnTable = columnTable[2:end, begin:end] #skip first row

	i = 1
	for col in eachrow(columnTable)
		typeVal = col["mat"] * "," * col["shp"] * "," * string(col["b"]) * "," * string(col["h"]) * "," * string(col["fc"])
		if typeVal ∉ keys(typeDict)
			typeDict[typeVal] = i
			i += 1
		end

		colDict[col["colID"]*"-"*col["storID"]] = Dict{Any,Any}("type" => typeDict[typeVal], "area" => 0, "story" => storyDict[col["storID"]]["order"])

		#TODO: check if can do by typeval or only typeDict[typeVal]
		if typeVal ∉ keys(groupDict)
			groupDict[typeVal] = Dict("columns" => [col["colID"]], "levels" => [storyDict[col["storID"]]["order"]])
		else
			push!(groupDict[typeVal]["columns"], col["colID"])
			push!(groupDict[typeVal]["levels"], storyDict[col["storID"]]["order"])
		end

	end

	return groupDict

end

function readAccessDB(fpath)
	
	connection_string = "Driver={Microsoft Access Driver (*.mdb, *.accdb)};Dbq="*fpath

	conn = ODBC.Connection(connection_string)

	dfProps = DBInterface.execute(conn, "SELECT colID,storID,mat,shp,b,h,fc FROM tblProperties") |> DataFrame

	dfStories = DBInterface.execute(conn, "SELECT storID,Elevation,Order FROM tblStories") |> DataFrame

	dfLoads = DBInterface.execute(conn, "SELECT colID, storID, FACT FROM tblLoads") |> DataFrame

	dfColumns = leftjoin(dfProps,dfLoads, on = [:colID, :storID])

	sort!(dfStories, [:Order])

	storyDict = Dict()
	reverseStoryDict = Dict()
	prev_elev = 0
	for row in eachrow(dfStories)
		storyDict[row["storID"]] = Dict("order" => row["Order"], "height" => round(row["Elevation"]-prev_elev,sigdigits=4))
		prev_elev = row["Elevation"]
		reverseStoryDict[row["Order"]] = row["storID"]
	end

	groupDict = readColumns(dfColumns, storyDict)

	return storyDict, groupDict, reverseStoryDict

end

function writeAccessDb(fpath, resultsDict)

	for group in values(resultsDict)
		if length(group.columns) > 1
			group.columns = join(group.columns,",")
		else
			group.columns = group.columns[1]
		end

		if length(group.stories) > 1
			group.stories = join(group.stories,",")
		else
			group.stories = group.stories[1]
		end
	end

	propDf = DataFrame(Any[collect(keys(resultsDict)), collect(values(resultsDict))])

	println(profDf)

end

function calcArea(b, h)
	#if h = 0, it's a circle
	if h == 0
		return pi*(b/2)^2
	elseif b == 0
		return pi*(h/2)^2
	else
		return b*h
	end
end

function getMin(b,h)
	if h == 0
		return b
	elseif b == 0 
		return h
	else
		return minimum(b,h)
	end
end

calcAlpha(fc) = maximum(0.85 - 0.0015 * fc, 0.67)

function calcRho(Pf, fc, fy, b, h)

	Ag = calcArea(b, h)

	hmin = getMin(b,h)
	c = minimum(0.2+0.002*hmin, 0.8)

	#TODO: check units are ok 
	c = minimum(0.01, (Pf / c / Ag - fcRed) / (fyRed - fcRed))
	fcRed = PHI_C * calcAlpha(fc) * fc
	fyRed = PHI_S * fy

	return maximum(0.01, (Pf / c / Ag - fcRed) / (fyRed - fcRed))
end

function reduce_columns(column_list, level_list)
	reducedColumns = []
	reducedLevels = []
	reducedColDict=  Dict()
	lvlColDict = Dict()
	uniqueCols = unique(column_list)

	for col in uniqueCols
		lvlSet = Set([level_list[l] for l in eachindex(level_list) if column_list[l]==col])
		if lvlSet ∉ keys(lvlColDict)
			lvlColDict[lvlSet] = Set([col])
		else
			push!(lvlColDict[lvlSet],col)
		end
	end

	for (i, (levels, cols)) in enumerate(lvlColDict)

		for lvl in levels 
			push!(reducedColumns, i)
			push!(reducedLevels, lvl)
		end
		reducedColDict[i] = cols
	end

	return reducedColumns, reducedLevels, reducedColDict	
end

function reduce_levels(column_list, level_list)
	#TODO: awful first implementation, much to clean up/speed up if it becomes speed issue
	reducedColumns = []
	reducedLevels = []
	reducedColDict=  Dict()
	lvlColDict = Dict()
	reverseLevelDict = Dict()
	uniqueLevels = sort!(unique(level_list))
	
	global i = 1
	
	while length(uniqueLevels) > 1 
	
	  #pop when done 
	  curLevel = popfirst!(uniqueLevels)
	  reverseLevelDict[curLevel] = i
	  levelList = [curLevel]
	
	  colSet_lvl = Set([column_list[l] for l in eachindex(column_list) if level_list[l]==curLevel])#popfirst!(uniqueLevels)])
	  colSet_next_lvl = Set([column_list[l] for l in eachindex(column_list) if level_list[l]==uniqueLevels[1]])
	
	  while colSet_lvl == colSet_next_lvl && length(uniqueLevels)>0
		push!(levelList, uniqueLevels[1])
		reverseLevelDict[uniqueLevels[1]] = i
		deleteat!(uniqueLevels, 1)
		if length(uniqueLevels) > 0
			colSet_next_lvl = Set([column_list[l] for l in eachindex(column_list) if level_list[l]==uniqueLevels[1]])
		end
	  end
	
	  lvlColDict[i] = levelList
	  global i += 1
	
	end
	
	if length(uniqueLevels) == 1 
		lvlColDict[i] = [uniqueLevels[1]]
		reverseLevelDict[uniqueLevels[1]] = i
	end
	
	reducedSet = Set()
	
	for i in eachindex(level_list)
		push!(reducedSet, (column_list[i], reverseLevelDict[level_list[i]]))
	end
	
	reducedCols = [i[1] for i in reducedSet]
	reducedLevels = [i[2] for i in reducedSet]
	
	return reducedCols, reducedLevels, lvlColDict, reverseLevelDict

end

function runModel(groupCols, groupLevels, typeVal, preference, timeLimit = 180)

	model = Model(HiGHS.Optimizer)
	set_silent(model)

	#whether or not 
	Bcl = zeros(Int8, maximum(groupCols), maximum(groupLevels))

	for c in 1:lastindex(groupCols)
		Bcl[groupCols[c], groupLevels[c]] = 1
	end

	#CONSTANTS
	numCols = maximum(groupCols)
	numLevels = maximum(groupLevels)

	#min max groups = min(# cols,# levels)
	numGroups = min(numCols, numLevels)

	M_Level = numLevels + 1

	# VARIABLES

	#whether elemnent in column c at level l is in group g 
	@variable(model, Xgcl[1:numGroups, 1:numCols, 1:numLevels], Bin, start=0)

	#whether group g exists 
	@variable(model, Xg[1:numGroups], Bin, start=0)

	#whether column c is in group g 
	@variable(model, Xc[1:numGroups, 1:numCols], Bin, start=0)

	#whether level l is in group g 
	@variable(model, Xl[1:numGroups, 1:numLevels], Bin, start=0)

	#upper bound level for group g
	@variable(model, 0<=Ug[1:numGroups]<=maximum(groupLevels), Int, start=0)

	#lower bound level for group g
	@variable(model, 0<=Lg[1:numGroups]<=maximum(groupLevels), Int, start=0)

	#range of levels for group g 
	@variable(model, 0<=Rg[1:numGroups]<=maximum(groupLevels), Int, start=0)

	#aux variables for making sure a level is in the gruop if it's within upper/lower bound
	@variable(model, Zu[1:numGroups, 1:numLevels], Bin, start=0)
	@variable(model, Zl[1:numGroups, 1:numLevels], Bin, start=0)

	# CONSTRAINTS 
	for c = 1:numCols
			for l = 1:numLevels
					#if an element exists, it is only in 1 group 
					#if an element doesn't exist, it is not in any groups
					@constraint(model, sum(Xgcl[1:numGroups, c, l]) == Bcl[c,l])
			end
	end

	for g = 1:numGroups 
		
		#range of group g = upper bound - lower bound 
		@constraint(model, Rg[g] == Ug[g]-Lg[g])

		for c = 1:numCols 
			for l = 1:numLevels
				#if element is in group, that column is in the group 
				@constraint(model, Xc[g,c] >= Xgcl[g,c,l])

				#if element is in group, that level is in the group 
				@constraint(model, Xl[g,l] >= Xgcl[g,c,l])

				#if element is in that group, that group exists
				@constraint(model, Xg[g] >= Xgcl[g,c,l])
			end
		end
	end

	for g = 1:numGroups
		for l = 1:numLevels
			#calculate upper, lower level
			#note, need to min(Ug) and max(Lg) for this to work

			@constraint(model, l*Xl[g,l] <= Ug[g] )
			@constraint(model, l*Xl[g,l]+(1-Xl[g,l])*M_Level >= Lg[g])

			@constraint(model, M_Level*Zl[g,l] >= l-Lg[g] + 1)
			@constraint(model, M_Level*(1-Zl[g,l]) >= Lg[g]-l)

			@constraint(model, M_Level*Zu[g,l] >= Ug[g]-l + 1)
			@constraint(model, M_Level*(1-Zu[g,l]) >= l-Ug[g])

			#if a level is within lower/upper bound, it's in the group 
			@constraint(model, 1+Xl[g,l] >= Zu[g,l]+Zl[g,l])

			for c = 1:numCols 
				#if column and level are in the group, so is the element 
				@constraint(model, Xc[g,c]+Xl[g,l] <= 1 + Xgcl[g,c,l])
			end
		end
	end

	# WARM START 
	## take min of num cols and num levels, group based on that 

	# if numCols <= numGroups 
	# 		#prioritize grouping by columns 
	# 	for c in 1:numCols
		
	# 		col_floors = [col_levels[i] for i in 1:length(cols) if cols[i] == c]
	# 		col_floors = Set(col_floors)
	# 		println(c," column: ",col_floors)

	# 		set_start_value(Xg[c], 1)
	# 		set_start_value(Xc[c,c], 1)
	# 		set_start_value(Ug[c], maximum(col_floors))
	# 		set_start_value(Lg[c], minimum(col_floors))
	# 		set_start_value(Rg[c], maximum(col_floors)-minimum(col_floors))

	# 		for l in col_floors
	# 			set_start_value(Xgcl[c,c,l], 1)
	# 			set_start_value(Xl[c,l], 1)
	# 		end

	# 		for l in 1:numLevels
	# 			if l >= minimum(col_floors)
	# 				set_start_value(Zl[c,l], 1)
	# 			end
	# 			if l <= maximum(col_floors)
	# 				set_start_value(Zu[c,l], 1)
	# 			end
	# 		end

	# 	end
	# else
	# 	for l in 1:numLevels

	# 		if l in col_levels
		
	# 			level_cols = [i for i in cols if col_levels[i] == l]
	# 			println(l, " level: ",level_cols)

	# 			set_start_value(Xg[l], 1)
	# 			set_start_value(Ug[l], l)
	# 			set_start_value(Lg[l], l)
	# 			set_start_value(Rg[l], 0)
	# 			set_start_value(Xl[l,l], 1)

	# 			for i in level_cols
	# 				set_start_value(Xgcl[l,i,l], 1)
	# 				set_start_value(Xc[l,i], 1)
	# 			end

	# 			set_start_value(Zl[l,l], 1)
	# 			set_start_value(Zu[l,l], 1)

	# 		end
	# 	end
	# end

	set_time_limit_sec(model, timeLimit)

	if preference == "MIN_COLUMNS" #taller groups 
		@objective(model, Min, sum(Lg) - sum(Ug) + 10000*sum(Xg) + 1000*sum(Xc))
	elseif preference == "MIN_LEVELS" #wider groups
		@objective(model, Min, sum(Lg) - sum(Ug) + 10000*sum(Xg) + 1000*sum(Xl))
	elseif preference == "MINIMAL"
		@objective(model, Min, sum(Xg))
	else preference == "DEFAULT"
		@objective(model, Min, sum(Lg) - sum(Ug) + 10000*sum(Xg))
	end

	optimize!(model)

	#println(solution_summary(model; verbose = false))

	groupedDict = Dict(i => Dict() for i in 1:numGroups)

	for i in 1:numGroups
		if value(Xg[i])>0.9

			cols_in_group = []

			levels_in_group = [round(j) for j in max(1,value(Lg[i])):value(Ug[i])]

			for c in 1:numCols
				if value(Xc[i,c]) > 0.9
					push!(cols_in_group, c)
				end
			end

			groupedDict[i] = Dict("columns" => cols_in_group, "levels" => levels_in_group)

		end
	end

	return groupedDict

end

function runModelWithTonnage()

	## FOR TONNAGE 

	levelHeights = [2.8,2.8,2.8,4.5,19.7]
	#colRho = 4*[rand(Float64) for i = 1:lastindex(groupCols)].+0.5
	println(colRho)
	colRho = [4.0922748737587105, 1.5672860888866786, 2.9148734381097077, 1.5836455134279879, 3.193023112010257, 1.057838126427391, 4.267899003287026, 1.7983520174702745, 2.8940732441884642, 1.0142432405086925, 3.7188721173026766, 2.7916000841710593, 4.219759061367123, 3.681817159647598, 3.332197643320552, 0.8259810530989795, 1.629898044233554, 1.0325221092266608, 2.122190828271761, 1.80178499815801, 4.257195164488189, 1.0325068397668917, 1.8484087902361548, 1.3271047191344363, 4.06345530598249, 1.769734946331773, 2.4067497877608073, 2.676946346154035, 3.557374763580838, 1.0190248991779387, 2.5103849768144424]
	area = 1*0.3

	#whether or not col exits
	Bcl = zeros(Int8, maximum(groupCols), maximum(groupLevels))

	for c in 1:lastindex(groupCols)
		Bcl[groupCols[c], groupLevels[c]] = 1
	end

	#rho of col 
	Rcl = zeros(Float64, maximum(groupCols), maximum(groupLevels))

	for c in 1:lastindex(groupCols)
		Rcl[groupCols[c], groupLevels[c]] = colRho[c]
	end

	#CONSTANTS
	numCols = maximum(groupCols)
	numLevels = maximum(groupLevels)

	#min max groups = min(# cols,# levels)
	numGroups = min(numCols, numLevels)+1

	M_Level = numLevels + 1
	M_Rho = 10

	# VARIABLES

	#group rho
	@variable(model, GR[1:numGroups] >= 0)

	#grouped column rho 
	@variable(model, GCR[1:numCols, 1:numLevels] >= 0)


	#whether elemnent in column c at level l is in group g 
	@variable(model, Xgcl[1:numGroups, 1:numCols, 1:numLevels], Bin, start=0)

	#whether group g exists 
	@variable(model, Xg[1:numGroups], Bin, start=0)

	#whether column c is in group g 
	@variable(model, Xc[1:numGroups, 1:numCols], Bin, start=0)

	#whether level l is in group g 
	@variable(model, Xl[1:numGroups, 1:numLevels], Bin, start=0)

	#upper bound level for group g
	@variable(model, 0<=Ug[1:numGroups]<=maximum(groupLevels), Int, start=0)

	#lower bound level for group g
	@variable(model, 0<=Lg[1:numGroups]<=maximum(groupLevels), Int, start=0)

	#range of levels for group g 
	@variable(model, 0<=Rg[1:numGroups]<=maximum(groupLevels), Int, start=0)

	#aux variables for making sure a level is in the gruop if it's within upper/lower bound
	@variable(model, Zu[1:numGroups, 1:numLevels], Bin, start=0)
	@variable(model, Zl[1:numGroups, 1:numLevels], Bin, start=0)

	# CONSTRAINTS 

	for g = 1:numGroups
		for c = 1:numCols 
			for l = 1:numLevels
				@constraint(model, GR[g] >= Xgcl[g,c,l]*Rcl[c,l])
				@constraint(model, GCR[c,l] >= GR[g] - M_Rho*(1-Xgcl[g,c,l]))
			end
		end
	end
end

function parseTypeVals(typeVal)

	typeVals = split(typeVal,",") #1 = Material, 2 = Shape, 3 = b, 4 = h, 5 = fc

	typeValDict = Dict{String, Any}()
	typeValDict["mat"] = typeVals[1]
	typeValDict["shp"] = typeVals[2]
	typeValDict["b"] = parse(Float64,typeVals[3])
	typeValDict["h"] = parse(Float64,typeVals[4])
	typeValDict["fc"] = parse(Float64,typeVals[5])
	gamma = 2396 #based on unit weight of 23.5 kN/m3
	poisson = 0.2
	typeValDict["Ec"] = (3300*sqrt(typeValDict["fc"])+6900)*(gamma/2300)^1.5
	typeValDict["Gc"] = typeValDict["Ec"]/(2*(1+poisson))

	if typeValDict["shp"] == "Rectangular"
		typeValDict["sgID_suffix"] = "COL"*typeVals[3]*"X"*typeVals[4]*"C"*typeVals[5]
	elseif typeValDict["shp"] == "Circular"
		typeValDict["sgID_suffix"] = "COL"*typeVals[3]*"DIA_"*"C"*typeVals[5]
	end
	
	return typeValDict

end

function groupZeroSparsity(resultsDict, i, typeVal, groupLevelDict, groupColDict, reverseStoryDict)

	#reduced levels already splits level groups which aren't adjacent!
	typeValDict = parseTypeVals(typeVal)
	
	cols = []

	for colGroup in values(groupColDict)
		for col in colGroup
			push!(cols, col)
		end

		sort!(cols)

	end
	
	for (levelGroup, lvls) in groupLevelDict
		i += 1
		
		levels = []
		for lvlOrder in lvls
			push!(levels, reverseStoryDict[lvlOrder])
		end

		sgName = cols[begin] * " " * levels[begin] * "-" * levels[end]
		sgID = sgName * typeValDict["sgID_suffix"]
		curTime = Dates.format(now(), "yyyy-mm-dd II:MM:SS p")

		#TODO: LU, PF
		newGroup = Subgroup(id = i, sgName = sgName, sgID = sgID, columns = cols, stories = levels, mat = typeValDict["mat"], shp = typeValDict["shp"], b = typeValDict["b"], h = typeValDict["h"], fc = typeValDict["fc"], Ec = typeValDict["Ec"], Gc = typeValDict["Gc"], Lu = 1.0, Pf = 1.0, DateCreated = curTime, LastUpdated = curTime)		

		resultsDict[i] = newGroup

	end

	return resultsDict, i

end

function groupByColumns(resultsDict, i, typeVal, groupLevelDict, groupColDict, reverseStoryDict, groupLevels, groupCols)

	#reduced levels already splits level groups which aren't adjacent!
	typeValDict = parseTypeVals(typeVal)

	for (colNum, colGroup) in groupColDict

		i += 1
		cols = []

		for col in colGroup
			push!(cols, col)
		end

		sort!(cols)

		lvlNums = Set([groupLevels[i] for i in eachindex(groupLevels) if groupCols[i] == colNum])
		
		levels = []

		for lvlNum in lvlNums
			for lvlOrder in groupLevelDict[lvlNum]
				push!(levels, reverseStoryDict[lvlOrder])
			end
		end

		sgName = cols[begin] * " " * levels[begin] * "-" * levels[end]
		sgID = sgName * typeValDict["sgID_suffix"]
		curTime = Dates.format(now(), "yyyy-mm-dd II:MM:SS p")

		#TODO: LU, PF
		newGroup = Subgroup(id = i, sgName = sgName, sgID = sgID, columns = cols, stories = levels, mat = typeValDict["mat"], shp = typeValDict["shp"], b = typeValDict["b"], h = typeValDict["h"], fc = typeValDict["fc"], Ec = typeValDict["Ec"], Gc = typeValDict["Gc"], Lu = 1.0, Pf = 1.0, DateCreated = curTime, LastUpdated = curTime)		

		resultsDict[i] = newGroup

	end
	
	return resultsDict, i

end

function groupByLevels(resultsDict, i, typeVal, groupLevelDict, groupColDict, reverseStoryDict, groupLevels, groupCols)

	#reduced levels already splits level groups which aren't adjacent!
	typeValDict = parseTypeVals(typeVal)

	for (lvlNum, lvlGroup) in groupLevelDict

		i += 1

		levels = []

		for lvlOrder in lvlGroup
			push!(levels, reverseStoryDict[lvlOrder])
		end

		cols = []

		colGroups = Set([groupCols[i] for i in eachindex(groupCols) if groupLevels[i] == lvlNum])

		for colGroup in colGroups
			for col in groupColDict[colGroup]
				push!(cols, col)
			end
		end	

		sort!(cols)

		sgName = cols[begin] * " " * levels[begin] * "-" * levels[end]
		sgID = sgName * typeValDict["sgID_suffix"]
		curTime = Dates.format(now(), "yyyy-mm-dd II:MM:SS p")

		#TODO: LU, PF
		newGroup = Subgroup(id = i, sgName = sgName, sgID = sgID, columns = cols, stories = levels, mat = typeValDict["mat"], shp = typeValDict["shp"], b = typeValDict["b"], h = typeValDict["h"], fc = typeValDict["fc"], Ec = typeValDict["Ec"], Gc = typeValDict["Gc"], Lu = 1.0, Pf = 1.0, DateCreated = curTime, LastUpdated = curTime)		

		resultsDict[i] = newGroup

	end
	
	return resultsDict, i

end

function addToResults(resultsDict, i, typeVal, groupedDict, groupLevelDict, groupColDict, reverseStoryDict)

	typeValDict = parseTypeVals(typeVal)
	
	for group in values(groupedDict)
		i += 1
		
		cols = []
		
		for colNum in group["columns"]
			for col in groupColDict[colNum]
				push!(cols, col)
			end
		end

		sort!(cols)

		levels = []
		for lvlNum in group["levels"]
			for lvlOrder in groupLevelDict[lvlNum]
				push!(levels, reverseStoryDict[lvlOrder])
			end
		end

		sgName = cols[begin] * " " * levels[begin] * "-" * levels[end]
		sgID = sgName * typeValDict["sgID_suffix"]
		curTime = Dates.format(now(), "yyyy-mm-dd II:MM:SS p")

		#TODO: LU, PF
		newGroup = Subgroup(id = i, sgName = sgName, sgID = sgID, columns = cols, stories = levels, mat = typeValDict["mat"], shp = typeValDict["shp"], b = typeValDict["b"], h = typeValDict["h"], fc = typeValDict["fc"], Ec = typeValDict["Ec"], Gc = typeValDict["Gc"], Lu = 1.0, Pf = 1.0, DateCreated = curTime, LastUpdated = curTime)		

		resultsDict[i] = newGroup

	end

	return resultsDict, i

end

function parseTo_tblSubGroups(group)

	return 0

end

function julia_main(fpath, model_preference = "MIN_COLUMNS")::Cint

	start = now()
	println("===START===")
	n = 0

	resultsDict = Dict()

	storyDict, groupDict, reverseStoryDict = readAccessDB(fpath)

	for (typeVal, group) in groupDict

		groupCols = group["columns"]
		groupLevels = group["levels"]
		
		groupCols, groupLevels, groupColDict = reduce_columns(group["columns"], group["levels"])

		groupCols, groupLevels, groupLevelDict, reverseLevelDict = reduce_levels(groupCols, groupLevels)

		uniqueColumns = Set(groupCols)
		uniqueLevels = Set(groupLevels)

		numColumns = length(groupCols)
		numUniqueColumns = length(uniqueColumns)
		numUniqueLevels = length(uniqueLevels)

		#how many 0s are in the array
		sparsity = 1 - numColumns/(numUniqueColumns*numUniqueLevels)

		typeCharacteristics = split(typeVal,",")
		
		#resultsDict[typeVal] = Dict("Material" => typeCharacteristics[1], "Shape" => typeCharacteristics[2], "b" => typeCharacteristics[3], "h" => typeCharacteristics[4], "Groups" => Dict())

		#TODO: adjust so each group is it's own dict, and add TYPE (zero sparse, model, 2 unique col, 2 unique row)

		#if array is "full" then can just make it all one big happy family
		if sparsity == 0
			println("zero sparse")
			#check for any level discontinuities, group by level intervals
			
			resultsDict, n = groupZeroSparsity(resultsDict, n, typeVal,  groupLevelDict, groupColDict, reverseStoryDict)

		else

			if numUniqueColumns == 2
				println("unique col")
				#group into 2 groups by columns

				resultsDict, n = groupByColumns(resultsDict,n, typeVal,  groupLevelDict, groupColDict, reverseStoryDict, groupLevels, groupCols)

			elseif numUniqueLevels == 2
				println("unique level")
				#group into 2 groups by level

				resultsDict, n = groupByLevelsresultsDict,(n, typeVal,  groupLevelDict, groupColDict, reverseStoryDict, groupLevels, groupCols)

			else
				println("---MODEL---")

				groupedDict = runModel(groupCols, groupLevels, typeVal, model_preference)

				resultsDict, n = addToResults(resultsDict, n, typeVal, groupedDict, groupLevelDict, groupColDict, reverseStoryDict)

			end
		end

		#

	end

	println("=== GROUP RESULTS ===")

	# for (groupNum, vals) in resultsDict

	# 	#parseTo_tblSubGroups(vals)

	# 	# println("\tGROUP: ", groupNum)
	# 	# println(vals)
	# 	#println("\t\tLEVELS: ", val["Levels"])
	# 	#println("\t\tCOLUMNS: ", val["Columns"])

	# end
					
	writeAccessDb(fpath, resultsDict)

	println("COMPLETED IN: ", now()-start)

	return 0

end


#julia_main("U:\\lyin\\For people\\For Daniel Pekar\\20241025\\TOR.129991.0001-20 CALEDONIA-20241023-Building A TOWER Column Workflow.accdb")

julia_main("H:\\Column Hub\\Nicole Butkovic\\Bayview - 20230818\\TOR.116564-20230818-NCB-TOWER A PODIUM Column Design Workflow.accdb")