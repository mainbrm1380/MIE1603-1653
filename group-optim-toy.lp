\ Model Grouping-Optimization-int
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  4048 group_exists[0] + 4048 group_exists[1] + 4048 group_exists[2]
   + 4048 group_exists[3] + 4048 group_exists[4] + cost_section[0,0]
   + cost_section[1,0] + cost_section[2,0] + cost_section[3,0]
   + cost_section[4,0] + cost_section[0,1] + cost_section[1,1]
   + cost_section[2,1] + cost_section[3,1] + cost_section[4,1]
   + cost_section[0,2] + cost_section[1,2] + cost_section[2,2]
   + cost_section[3,2] + cost_section[4,2] + cost_section[0,3]
   + cost_section[1,3] + cost_section[2,3] + cost_section[3,3]
   + cost_section[4,3] + cost_section[0,4] + cost_section[1,4]
   + cost_section[2,4] + cost_section[3,4] + cost_section[4,4]
   + group_lb[0] + group_lb[1] + group_lb[2] + group_lb[3] + group_lb[4]
   - group_ub[0] - group_ub[1] - group_ub[2] - group_ub[3] - group_ub[4]
Subject To
 R0: y[0,0,0] + y[0,0,1] + y[0,0,2] + y[0,0,3] + y[0,0,4] + y[0,0,5] = 1
 R1: y[1,0,0] + y[1,0,1] + y[1,0,2] + y[1,0,3] + y[1,0,4] + y[1,0,5] = 1
 R2: y[2,0,0] + y[2,0,1] + y[2,0,2] + y[2,0,3] + y[2,0,4] + y[2,0,5] = 1
 R3: y[3,0,0] + y[3,0,1] + y[3,0,2] + y[3,0,3] + y[3,0,4] + y[3,0,5] = 1
 R4: y[4,0,0] + y[4,0,1] + y[4,0,2] + y[4,0,3] + y[4,0,4] + y[4,0,5] = 1
 R5: y[0,1,0] + y[0,1,1] + y[0,1,2] + y[0,1,3] + y[0,1,4] + y[0,1,5] = 1
 R6: y[1,1,0] + y[1,1,1] + y[1,1,2] + y[1,1,3] + y[1,1,4] + y[1,1,5] = 1
 R7: y[2,1,0] + y[2,1,1] + y[2,1,2] + y[2,1,3] + y[2,1,4] + y[2,1,5] = 1
 R8: y[3,1,0] + y[3,1,1] + y[3,1,2] + y[3,1,3] + y[3,1,4] + y[3,1,5] = 1
 R9: y[4,1,0] + y[4,1,1] + y[4,1,2] + y[4,1,3] + y[4,1,4] + y[4,1,5] = 1
 R10: y[0,2,0] + y[0,2,1] + y[0,2,2] + y[0,2,3] + y[0,2,4] + y[0,2,5] = 1
 R11: y[1,2,0] + y[1,2,1] + y[1,2,2] + y[1,2,3] + y[1,2,4] + y[1,2,5] = 1
 R12: y[2,2,0] + y[2,2,1] + y[2,2,2] + y[2,2,3] + y[2,2,4] + y[2,2,5] = 1
 R13: y[3,2,0] + y[3,2,1] + y[3,2,2] + y[3,2,3] + y[3,2,4] + y[3,2,5] = 1
 R14: y[4,2,0] + y[4,2,1] + y[4,2,2] + y[4,2,3] + y[4,2,4] + y[4,2,5] = 1
 R15: y[0,3,0] + y[0,3,1] + y[0,3,2] + y[0,3,3] + y[0,3,4] + y[0,3,5] = 1
 R16: y[1,3,0] + y[1,3,1] + y[1,3,2] + y[1,3,3] + y[1,3,4] + y[1,3,5] = 1
 R17: y[2,3,0] + y[2,3,1] + y[2,3,2] + y[2,3,3] + y[2,3,4] + y[2,3,5] = 1
 R18: y[3,3,0] + y[3,3,1] + y[3,3,2] + y[3,3,3] + y[3,3,4] + y[3,3,5] = 1
 R19: y[4,3,0] + y[4,3,1] + y[4,3,2] + y[4,3,3] + y[4,3,4] + y[4,3,5] = 1
 R20: y[0,4,0] + y[0,4,1] + y[0,4,2] + y[0,4,3] + y[0,4,4] + y[0,4,5] = 1
 R21: y[1,4,0] + y[1,4,1] + y[1,4,2] + y[1,4,3] + y[1,4,4] + y[1,4,5] = 1
 R22: y[2,4,0] + y[2,4,1] + y[2,4,2] + y[2,4,3] + y[2,4,4] + y[2,4,5] = 1
 R23: y[3,4,0] + y[3,4,1] + y[3,4,2] + y[3,4,3] + y[3,4,4] + y[3,4,5] = 1
 R24: y[4,4,0] + y[4,4,1] + y[4,4,2] + y[4,4,3] + y[4,4,4] + y[4,4,5] = 1
 R25: x[0,0,0] + x[1,0,0] + x[2,0,0] + x[3,0,0] + x[4,0,0] = 0
 R26: element_section[0,0] >= 0
 R27: x[0,0,1] + x[1,0,1] + x[2,0,1] + x[3,0,1] + x[4,0,1] = 1
 R28: element_section[0,1] >= 2
 R29: x[0,0,2] + x[1,0,2] + x[2,0,2] + x[3,0,2] + x[4,0,2] = 0
 R30: element_section[0,2] >= 0
 R31: x[0,0,3] + x[1,0,3] + x[2,0,3] + x[3,0,3] + x[4,0,3] = 0
 R32: element_section[0,3] >= 0
 R33: x[0,0,4] + x[1,0,4] + x[2,0,4] + x[3,0,4] + x[4,0,4] = 0
 R34: element_section[0,4] >= 0
 R35: x[0,1,0] + x[1,1,0] + x[2,1,0] + x[3,1,0] + x[4,1,0] = 1
 R36: element_section[1,0] >= 2
 R37: x[0,1,1] + x[1,1,1] + x[2,1,1] + x[3,1,1] + x[4,1,1] = 1
 R38: element_section[1,1] >= 2
 R39: x[0,1,2] + x[1,1,2] + x[2,1,2] + x[3,1,2] + x[4,1,2] = 1
 R40: element_section[1,2] >= 2
 R41: x[0,1,3] + x[1,1,3] + x[2,1,3] + x[3,1,3] + x[4,1,3] = 1
 R42: element_section[1,3] >= 2
 R43: x[0,1,4] + x[1,1,4] + x[2,1,4] + x[3,1,4] + x[4,1,4] = 0
 R44: element_section[1,4] >= 0
 R45: - element_section[1,0] + element_section[1,1] <= 0
 R46: - element_section[1,1] + element_section[1,2] <= 0
 R47: - element_section[1,2] + element_section[1,3] <= 0
 R48: x[0,2,0] + x[1,2,0] + x[2,2,0] + x[3,2,0] + x[4,2,0] = 1
 R49: element_section[2,0] >= 2
 R50: x[0,2,1] + x[1,2,1] + x[2,2,1] + x[3,2,1] + x[4,2,1] = 1
 R51: element_section[2,1] >= 2
 R52: x[0,2,2] + x[1,2,2] + x[2,2,2] + x[3,2,2] + x[4,2,2] = 1
 R53: element_section[2,2] >= 2
 R54: x[0,2,3] + x[1,2,3] + x[2,2,3] + x[3,2,3] + x[4,2,3] = 1
 R55: element_section[2,3] >= 2
 R56: x[0,2,4] + x[1,2,4] + x[2,2,4] + x[3,2,4] + x[4,2,4] = 1
 R57: element_section[2,4] >= 2
 R58: - element_section[2,0] + element_section[2,1] <= 0
 R59: - element_section[2,1] + element_section[2,2] <= 0
 R60: - element_section[2,2] + element_section[2,3] <= 0
 R61: - element_section[2,3] + element_section[2,4] <= 0
 R62: x[0,3,0] + x[1,3,0] + x[2,3,0] + x[3,3,0] + x[4,3,0] = 1
 R63: element_section[3,0] >= 4
 R64: x[0,3,1] + x[1,3,1] + x[2,3,1] + x[3,3,1] + x[4,3,1] = 1
 R65: element_section[3,1] >= 4
 R66: x[0,3,2] + x[1,3,2] + x[2,3,2] + x[3,3,2] + x[4,3,2] = 1
 R67: element_section[3,2] >= 4
 R68: x[0,3,3] + x[1,3,3] + x[2,3,3] + x[3,3,3] + x[4,3,3] = 1
 R69: element_section[3,3] >= 4
 R70: x[0,3,4] + x[1,3,4] + x[2,3,4] + x[3,3,4] + x[4,3,4] = 0
 R71: element_section[3,4] >= 0
 R72: - element_section[3,0] + element_section[3,1] <= 0
 R73: - element_section[3,1] + element_section[3,2] <= 0
 R74: - element_section[3,2] + element_section[3,3] <= 0
 R75: x[0,4,0] + x[1,4,0] + x[2,4,0] + x[3,4,0] + x[4,4,0] = 1
 R76: element_section[4,0] >= 3
 R77: x[0,4,1] + x[1,4,1] + x[2,4,1] + x[3,4,1] + x[4,4,1] = 1
 R78: element_section[4,1] >= 3
 R79: x[0,4,2] + x[1,4,2] + x[2,4,2] + x[3,4,2] + x[4,4,2] = 1
 R80: element_section[4,2] >= 3
 R81: x[0,4,3] + x[1,4,3] + x[2,4,3] + x[3,4,3] + x[4,4,3] = 1
 R82: element_section[4,3] >= 3
 R83: x[0,4,4] + x[1,4,4] + x[2,4,4] + x[3,4,4] + x[4,4,4] = 1
 R84: element_section[4,4] >= 3
 R85: - element_section[4,0] + element_section[4,1] <= 0
 R86: - element_section[4,1] + element_section[4,2] <= 0
 R87: - element_section[4,2] + element_section[4,3] <= 0
 R88: - element_section[4,3] + element_section[4,4] <= 0
 R89: group_lb[0] - group_ub[0] + group_range[0] = 0
 R90: - 5 x[0,0,0] - group_section[0] + element_section[0,0] >= -5
 R91: 5 x[0,0,0] - group_section[0] + element_section[0,0] <= 5
 R92: - x[0,0,0] + col_in_group[0,0] >= 0
 R93: - x[0,0,0] + level_in_group[0,0] >= 0
 R94: - x[0,0,0] + group_exists[0] >= 0
 R95: - 5 x[0,0,1] - group_section[0] + element_section[0,1] >= -5
 R96: 5 x[0,0,1] - group_section[0] + element_section[0,1] <= 5
 R97: - x[0,0,1] + col_in_group[0,0] >= 0
 R98: - x[0,0,1] + level_in_group[0,1] >= 0
 R99: - x[0,0,1] + group_exists[0] >= 0
 R100: - 5 x[0,0,2] - group_section[0] + element_section[0,2] >= -5
 R101: 5 x[0,0,2] - group_section[0] + element_section[0,2] <= 5
 R102: - x[0,0,2] + col_in_group[0,0] >= 0
 R103: - x[0,0,2] + level_in_group[0,2] >= 0
 R104: - x[0,0,2] + group_exists[0] >= 0
 R105: - 5 x[0,0,3] - group_section[0] + element_section[0,3] >= -5
 R106: 5 x[0,0,3] - group_section[0] + element_section[0,3] <= 5
 R107: - x[0,0,3] + col_in_group[0,0] >= 0
 R108: - x[0,0,3] + level_in_group[0,3] >= 0
 R109: - x[0,0,3] + group_exists[0] >= 0
 R110: - 5 x[0,0,4] - group_section[0] + element_section[0,4] >= -5
 R111: 5 x[0,0,4] - group_section[0] + element_section[0,4] <= 5
 R112: - x[0,0,4] + col_in_group[0,0] >= 0
 R113: - x[0,0,4] + level_in_group[0,4] >= 0
 R114: - x[0,0,4] + group_exists[0] >= 0
 R115: - 5 x[0,1,0] - group_section[0] + element_section[1,0] >= -5
 R116: 5 x[0,1,0] - group_section[0] + element_section[1,0] <= 5
 R117: - x[0,1,0] + col_in_group[0,1] >= 0
 R118: - x[0,1,0] + level_in_group[0,0] >= 0
 R119: - x[0,1,0] + group_exists[0] >= 0
 R120: - 5 x[0,1,1] - group_section[0] + element_section[1,1] >= -5
 R121: 5 x[0,1,1] - group_section[0] + element_section[1,1] <= 5
 R122: - x[0,1,1] + col_in_group[0,1] >= 0
 R123: - x[0,1,1] + level_in_group[0,1] >= 0
 R124: - x[0,1,1] + group_exists[0] >= 0
 R125: - 5 x[0,1,2] - group_section[0] + element_section[1,2] >= -5
 R126: 5 x[0,1,2] - group_section[0] + element_section[1,2] <= 5
 R127: - x[0,1,2] + col_in_group[0,1] >= 0
 R128: - x[0,1,2] + level_in_group[0,2] >= 0
 R129: - x[0,1,2] + group_exists[0] >= 0
 R130: - 5 x[0,1,3] - group_section[0] + element_section[1,3] >= -5
 R131: 5 x[0,1,3] - group_section[0] + element_section[1,3] <= 5
 R132: - x[0,1,3] + col_in_group[0,1] >= 0
 R133: - x[0,1,3] + level_in_group[0,3] >= 0
 R134: - x[0,1,3] + group_exists[0] >= 0
 R135: - 5 x[0,1,4] - group_section[0] + element_section[1,4] >= -5
 R136: 5 x[0,1,4] - group_section[0] + element_section[1,4] <= 5
 R137: - x[0,1,4] + col_in_group[0,1] >= 0
 R138: - x[0,1,4] + level_in_group[0,4] >= 0
 R139: - x[0,1,4] + group_exists[0] >= 0
 R140: - 5 x[0,2,0] - group_section[0] + element_section[2,0] >= -5
 R141: 5 x[0,2,0] - group_section[0] + element_section[2,0] <= 5
 R142: - x[0,2,0] + col_in_group[0,2] >= 0
 R143: - x[0,2,0] + level_in_group[0,0] >= 0
 R144: - x[0,2,0] + group_exists[0] >= 0
 R145: - 5 x[0,2,1] - group_section[0] + element_section[2,1] >= -5
 R146: 5 x[0,2,1] - group_section[0] + element_section[2,1] <= 5
 R147: - x[0,2,1] + col_in_group[0,2] >= 0
 R148: - x[0,2,1] + level_in_group[0,1] >= 0
 R149: - x[0,2,1] + group_exists[0] >= 0
 R150: - 5 x[0,2,2] - group_section[0] + element_section[2,2] >= -5
 R151: 5 x[0,2,2] - group_section[0] + element_section[2,2] <= 5
 R152: - x[0,2,2] + col_in_group[0,2] >= 0
 R153: - x[0,2,2] + level_in_group[0,2] >= 0
 R154: - x[0,2,2] + group_exists[0] >= 0
 R155: - 5 x[0,2,3] - group_section[0] + element_section[2,3] >= -5
 R156: 5 x[0,2,3] - group_section[0] + element_section[2,3] <= 5
 R157: - x[0,2,3] + col_in_group[0,2] >= 0
 R158: - x[0,2,3] + level_in_group[0,3] >= 0
 R159: - x[0,2,3] + group_exists[0] >= 0
 R160: - 5 x[0,2,4] - group_section[0] + element_section[2,4] >= -5
 R161: 5 x[0,2,4] - group_section[0] + element_section[2,4] <= 5
 R162: - x[0,2,4] + col_in_group[0,2] >= 0
 R163: - x[0,2,4] + level_in_group[0,4] >= 0
 R164: - x[0,2,4] + group_exists[0] >= 0
 R165: - 5 x[0,3,0] - group_section[0] + element_section[3,0] >= -5
 R166: 5 x[0,3,0] - group_section[0] + element_section[3,0] <= 5
 R167: - x[0,3,0] + col_in_group[0,3] >= 0
 R168: - x[0,3,0] + level_in_group[0,0] >= 0
 R169: - x[0,3,0] + group_exists[0] >= 0
 R170: - 5 x[0,3,1] - group_section[0] + element_section[3,1] >= -5
 R171: 5 x[0,3,1] - group_section[0] + element_section[3,1] <= 5
 R172: - x[0,3,1] + col_in_group[0,3] >= 0
 R173: - x[0,3,1] + level_in_group[0,1] >= 0
 R174: - x[0,3,1] + group_exists[0] >= 0
 R175: - 5 x[0,3,2] - group_section[0] + element_section[3,2] >= -5
 R176: 5 x[0,3,2] - group_section[0] + element_section[3,2] <= 5
 R177: - x[0,3,2] + col_in_group[0,3] >= 0
 R178: - x[0,3,2] + level_in_group[0,2] >= 0
 R179: - x[0,3,2] + group_exists[0] >= 0
 R180: - 5 x[0,3,3] - group_section[0] + element_section[3,3] >= -5
 R181: 5 x[0,3,3] - group_section[0] + element_section[3,3] <= 5
 R182: - x[0,3,3] + col_in_group[0,3] >= 0
 R183: - x[0,3,3] + level_in_group[0,3] >= 0
 R184: - x[0,3,3] + group_exists[0] >= 0
 R185: - 5 x[0,3,4] - group_section[0] + element_section[3,4] >= -5
 R186: 5 x[0,3,4] - group_section[0] + element_section[3,4] <= 5
 R187: - x[0,3,4] + col_in_group[0,3] >= 0
 R188: - x[0,3,4] + level_in_group[0,4] >= 0
 R189: - x[0,3,4] + group_exists[0] >= 0
 R190: - 5 x[0,4,0] - group_section[0] + element_section[4,0] >= -5
 R191: 5 x[0,4,0] - group_section[0] + element_section[4,0] <= 5
 R192: - x[0,4,0] + col_in_group[0,4] >= 0
 R193: - x[0,4,0] + level_in_group[0,0] >= 0
 R194: - x[0,4,0] + group_exists[0] >= 0
 R195: - 5 x[0,4,1] - group_section[0] + element_section[4,1] >= -5
 R196: 5 x[0,4,1] - group_section[0] + element_section[4,1] <= 5
 R197: - x[0,4,1] + col_in_group[0,4] >= 0
 R198: - x[0,4,1] + level_in_group[0,1] >= 0
 R199: - x[0,4,1] + group_exists[0] >= 0
 R200: - 5 x[0,4,2] - group_section[0] + element_section[4,2] >= -5
 R201: 5 x[0,4,2] - group_section[0] + element_section[4,2] <= 5
 R202: - x[0,4,2] + col_in_group[0,4] >= 0
 R203: - x[0,4,2] + level_in_group[0,2] >= 0
 R204: - x[0,4,2] + group_exists[0] >= 0
 R205: - 5 x[0,4,3] - group_section[0] + element_section[4,3] >= -5
 R206: 5 x[0,4,3] - group_section[0] + element_section[4,3] <= 5
 R207: - x[0,4,3] + col_in_group[0,4] >= 0
 R208: - x[0,4,3] + level_in_group[0,3] >= 0
 R209: - x[0,4,3] + group_exists[0] >= 0
 R210: - 5 x[0,4,4] - group_section[0] + element_section[4,4] >= -5
 R211: 5 x[0,4,4] - group_section[0] + element_section[4,4] <= 5
 R212: - x[0,4,4] + col_in_group[0,4] >= 0
 R213: - x[0,4,4] + level_in_group[0,4] >= 0
 R214: - x[0,4,4] + group_exists[0] >= 0
 R215: group_lb[1] - group_ub[1] + group_range[1] = 0
 R216: - 5 x[1,0,0] - group_section[1] + element_section[0,0] >= -5
 R217: 5 x[1,0,0] - group_section[1] + element_section[0,0] <= 5
 R218: - x[1,0,0] + col_in_group[1,0] >= 0
 R219: - x[1,0,0] + level_in_group[1,0] >= 0
 R220: - x[1,0,0] + group_exists[1] >= 0
 R221: - 5 x[1,0,1] - group_section[1] + element_section[0,1] >= -5
 R222: 5 x[1,0,1] - group_section[1] + element_section[0,1] <= 5
 R223: - x[1,0,1] + col_in_group[1,0] >= 0
 R224: - x[1,0,1] + level_in_group[1,1] >= 0
 R225: - x[1,0,1] + group_exists[1] >= 0
 R226: - 5 x[1,0,2] - group_section[1] + element_section[0,2] >= -5
 R227: 5 x[1,0,2] - group_section[1] + element_section[0,2] <= 5
 R228: - x[1,0,2] + col_in_group[1,0] >= 0
 R229: - x[1,0,2] + level_in_group[1,2] >= 0
 R230: - x[1,0,2] + group_exists[1] >= 0
 R231: - 5 x[1,0,3] - group_section[1] + element_section[0,3] >= -5
 R232: 5 x[1,0,3] - group_section[1] + element_section[0,3] <= 5
 R233: - x[1,0,3] + col_in_group[1,0] >= 0
 R234: - x[1,0,3] + level_in_group[1,3] >= 0
 R235: - x[1,0,3] + group_exists[1] >= 0
 R236: - 5 x[1,0,4] - group_section[1] + element_section[0,4] >= -5
 R237: 5 x[1,0,4] - group_section[1] + element_section[0,4] <= 5
 R238: - x[1,0,4] + col_in_group[1,0] >= 0
 R239: - x[1,0,4] + level_in_group[1,4] >= 0
 R240: - x[1,0,4] + group_exists[1] >= 0
 R241: - 5 x[1,1,0] - group_section[1] + element_section[1,0] >= -5
 R242: 5 x[1,1,0] - group_section[1] + element_section[1,0] <= 5
 R243: - x[1,1,0] + col_in_group[1,1] >= 0
 R244: - x[1,1,0] + level_in_group[1,0] >= 0
 R245: - x[1,1,0] + group_exists[1] >= 0
 R246: - 5 x[1,1,1] - group_section[1] + element_section[1,1] >= -5
 R247: 5 x[1,1,1] - group_section[1] + element_section[1,1] <= 5
 R248: - x[1,1,1] + col_in_group[1,1] >= 0
 R249: - x[1,1,1] + level_in_group[1,1] >= 0
 R250: - x[1,1,1] + group_exists[1] >= 0
 R251: - 5 x[1,1,2] - group_section[1] + element_section[1,2] >= -5
 R252: 5 x[1,1,2] - group_section[1] + element_section[1,2] <= 5
 R253: - x[1,1,2] + col_in_group[1,1] >= 0
 R254: - x[1,1,2] + level_in_group[1,2] >= 0
 R255: - x[1,1,2] + group_exists[1] >= 0
 R256: - 5 x[1,1,3] - group_section[1] + element_section[1,3] >= -5
 R257: 5 x[1,1,3] - group_section[1] + element_section[1,3] <= 5
 R258: - x[1,1,3] + col_in_group[1,1] >= 0
 R259: - x[1,1,3] + level_in_group[1,3] >= 0
 R260: - x[1,1,3] + group_exists[1] >= 0
 R261: - 5 x[1,1,4] - group_section[1] + element_section[1,4] >= -5
 R262: 5 x[1,1,4] - group_section[1] + element_section[1,4] <= 5
 R263: - x[1,1,4] + col_in_group[1,1] >= 0
 R264: - x[1,1,4] + level_in_group[1,4] >= 0
 R265: - x[1,1,4] + group_exists[1] >= 0
 R266: - 5 x[1,2,0] - group_section[1] + element_section[2,0] >= -5
 R267: 5 x[1,2,0] - group_section[1] + element_section[2,0] <= 5
 R268: - x[1,2,0] + col_in_group[1,2] >= 0
 R269: - x[1,2,0] + level_in_group[1,0] >= 0
 R270: - x[1,2,0] + group_exists[1] >= 0
 R271: - 5 x[1,2,1] - group_section[1] + element_section[2,1] >= -5
 R272: 5 x[1,2,1] - group_section[1] + element_section[2,1] <= 5
 R273: - x[1,2,1] + col_in_group[1,2] >= 0
 R274: - x[1,2,1] + level_in_group[1,1] >= 0
 R275: - x[1,2,1] + group_exists[1] >= 0
 R276: - 5 x[1,2,2] - group_section[1] + element_section[2,2] >= -5
 R277: 5 x[1,2,2] - group_section[1] + element_section[2,2] <= 5
 R278: - x[1,2,2] + col_in_group[1,2] >= 0
 R279: - x[1,2,2] + level_in_group[1,2] >= 0
 R280: - x[1,2,2] + group_exists[1] >= 0
 R281: - 5 x[1,2,3] - group_section[1] + element_section[2,3] >= -5
 R282: 5 x[1,2,3] - group_section[1] + element_section[2,3] <= 5
 R283: - x[1,2,3] + col_in_group[1,2] >= 0
 R284: - x[1,2,3] + level_in_group[1,3] >= 0
 R285: - x[1,2,3] + group_exists[1] >= 0
 R286: - 5 x[1,2,4] - group_section[1] + element_section[2,4] >= -5
 R287: 5 x[1,2,4] - group_section[1] + element_section[2,4] <= 5
 R288: - x[1,2,4] + col_in_group[1,2] >= 0
 R289: - x[1,2,4] + level_in_group[1,4] >= 0
 R290: - x[1,2,4] + group_exists[1] >= 0
 R291: - 5 x[1,3,0] - group_section[1] + element_section[3,0] >= -5
 R292: 5 x[1,3,0] - group_section[1] + element_section[3,0] <= 5
 R293: - x[1,3,0] + col_in_group[1,3] >= 0
 R294: - x[1,3,0] + level_in_group[1,0] >= 0
 R295: - x[1,3,0] + group_exists[1] >= 0
 R296: - 5 x[1,3,1] - group_section[1] + element_section[3,1] >= -5
 R297: 5 x[1,3,1] - group_section[1] + element_section[3,1] <= 5
 R298: - x[1,3,1] + col_in_group[1,3] >= 0
 R299: - x[1,3,1] + level_in_group[1,1] >= 0
 R300: - x[1,3,1] + group_exists[1] >= 0
 R301: - 5 x[1,3,2] - group_section[1] + element_section[3,2] >= -5
 R302: 5 x[1,3,2] - group_section[1] + element_section[3,2] <= 5
 R303: - x[1,3,2] + col_in_group[1,3] >= 0
 R304: - x[1,3,2] + level_in_group[1,2] >= 0
 R305: - x[1,3,2] + group_exists[1] >= 0
 R306: - 5 x[1,3,3] - group_section[1] + element_section[3,3] >= -5
 R307: 5 x[1,3,3] - group_section[1] + element_section[3,3] <= 5
 R308: - x[1,3,3] + col_in_group[1,3] >= 0
 R309: - x[1,3,3] + level_in_group[1,3] >= 0
 R310: - x[1,3,3] + group_exists[1] >= 0
 R311: - 5 x[1,3,4] - group_section[1] + element_section[3,4] >= -5
 R312: 5 x[1,3,4] - group_section[1] + element_section[3,4] <= 5
 R313: - x[1,3,4] + col_in_group[1,3] >= 0
 R314: - x[1,3,4] + level_in_group[1,4] >= 0
 R315: - x[1,3,4] + group_exists[1] >= 0
 R316: - 5 x[1,4,0] - group_section[1] + element_section[4,0] >= -5
 R317: 5 x[1,4,0] - group_section[1] + element_section[4,0] <= 5
 R318: - x[1,4,0] + col_in_group[1,4] >= 0
 R319: - x[1,4,0] + level_in_group[1,0] >= 0
 R320: - x[1,4,0] + group_exists[1] >= 0
 R321: - 5 x[1,4,1] - group_section[1] + element_section[4,1] >= -5
 R322: 5 x[1,4,1] - group_section[1] + element_section[4,1] <= 5
 R323: - x[1,4,1] + col_in_group[1,4] >= 0
 R324: - x[1,4,1] + level_in_group[1,1] >= 0
 R325: - x[1,4,1] + group_exists[1] >= 0
 R326: - 5 x[1,4,2] - group_section[1] + element_section[4,2] >= -5
 R327: 5 x[1,4,2] - group_section[1] + element_section[4,2] <= 5
 R328: - x[1,4,2] + col_in_group[1,4] >= 0
 R329: - x[1,4,2] + level_in_group[1,2] >= 0
 R330: - x[1,4,2] + group_exists[1] >= 0
 R331: - 5 x[1,4,3] - group_section[1] + element_section[4,3] >= -5
 R332: 5 x[1,4,3] - group_section[1] + element_section[4,3] <= 5
 R333: - x[1,4,3] + col_in_group[1,4] >= 0
 R334: - x[1,4,3] + level_in_group[1,3] >= 0
 R335: - x[1,4,3] + group_exists[1] >= 0
 R336: - 5 x[1,4,4] - group_section[1] + element_section[4,4] >= -5
 R337: 5 x[1,4,4] - group_section[1] + element_section[4,4] <= 5
 R338: - x[1,4,4] + col_in_group[1,4] >= 0
 R339: - x[1,4,4] + level_in_group[1,4] >= 0
 R340: - x[1,4,4] + group_exists[1] >= 0
 R341: group_lb[2] - group_ub[2] + group_range[2] = 0
 R342: - 5 x[2,0,0] - group_section[2] + element_section[0,0] >= -5
 R343: 5 x[2,0,0] - group_section[2] + element_section[0,0] <= 5
 R344: - x[2,0,0] + col_in_group[2,0] >= 0
 R345: - x[2,0,0] + level_in_group[2,0] >= 0
 R346: - x[2,0,0] + group_exists[2] >= 0
 R347: - 5 x[2,0,1] - group_section[2] + element_section[0,1] >= -5
 R348: 5 x[2,0,1] - group_section[2] + element_section[0,1] <= 5
 R349: - x[2,0,1] + col_in_group[2,0] >= 0
 R350: - x[2,0,1] + level_in_group[2,1] >= 0
 R351: - x[2,0,1] + group_exists[2] >= 0
 R352: - 5 x[2,0,2] - group_section[2] + element_section[0,2] >= -5
 R353: 5 x[2,0,2] - group_section[2] + element_section[0,2] <= 5
 R354: - x[2,0,2] + col_in_group[2,0] >= 0
 R355: - x[2,0,2] + level_in_group[2,2] >= 0
 R356: - x[2,0,2] + group_exists[2] >= 0
 R357: - 5 x[2,0,3] - group_section[2] + element_section[0,3] >= -5
 R358: 5 x[2,0,3] - group_section[2] + element_section[0,3] <= 5
 R359: - x[2,0,3] + col_in_group[2,0] >= 0
 R360: - x[2,0,3] + level_in_group[2,3] >= 0
 R361: - x[2,0,3] + group_exists[2] >= 0
 R362: - 5 x[2,0,4] - group_section[2] + element_section[0,4] >= -5
 R363: 5 x[2,0,4] - group_section[2] + element_section[0,4] <= 5
 R364: - x[2,0,4] + col_in_group[2,0] >= 0
 R365: - x[2,0,4] + level_in_group[2,4] >= 0
 R366: - x[2,0,4] + group_exists[2] >= 0
 R367: - 5 x[2,1,0] - group_section[2] + element_section[1,0] >= -5
 R368: 5 x[2,1,0] - group_section[2] + element_section[1,0] <= 5
 R369: - x[2,1,0] + col_in_group[2,1] >= 0
 R370: - x[2,1,0] + level_in_group[2,0] >= 0
 R371: - x[2,1,0] + group_exists[2] >= 0
 R372: - 5 x[2,1,1] - group_section[2] + element_section[1,1] >= -5
 R373: 5 x[2,1,1] - group_section[2] + element_section[1,1] <= 5
 R374: - x[2,1,1] + col_in_group[2,1] >= 0
 R375: - x[2,1,1] + level_in_group[2,1] >= 0
 R376: - x[2,1,1] + group_exists[2] >= 0
 R377: - 5 x[2,1,2] - group_section[2] + element_section[1,2] >= -5
 R378: 5 x[2,1,2] - group_section[2] + element_section[1,2] <= 5
 R379: - x[2,1,2] + col_in_group[2,1] >= 0
 R380: - x[2,1,2] + level_in_group[2,2] >= 0
 R381: - x[2,1,2] + group_exists[2] >= 0
 R382: - 5 x[2,1,3] - group_section[2] + element_section[1,3] >= -5
 R383: 5 x[2,1,3] - group_section[2] + element_section[1,3] <= 5
 R384: - x[2,1,3] + col_in_group[2,1] >= 0
 R385: - x[2,1,3] + level_in_group[2,3] >= 0
 R386: - x[2,1,3] + group_exists[2] >= 0
 R387: - 5 x[2,1,4] - group_section[2] + element_section[1,4] >= -5
 R388: 5 x[2,1,4] - group_section[2] + element_section[1,4] <= 5
 R389: - x[2,1,4] + col_in_group[2,1] >= 0
 R390: - x[2,1,4] + level_in_group[2,4] >= 0
 R391: - x[2,1,4] + group_exists[2] >= 0
 R392: - 5 x[2,2,0] - group_section[2] + element_section[2,0] >= -5
 R393: 5 x[2,2,0] - group_section[2] + element_section[2,0] <= 5
 R394: - x[2,2,0] + col_in_group[2,2] >= 0
 R395: - x[2,2,0] + level_in_group[2,0] >= 0
 R396: - x[2,2,0] + group_exists[2] >= 0
 R397: - 5 x[2,2,1] - group_section[2] + element_section[2,1] >= -5
 R398: 5 x[2,2,1] - group_section[2] + element_section[2,1] <= 5
 R399: - x[2,2,1] + col_in_group[2,2] >= 0
 R400: - x[2,2,1] + level_in_group[2,1] >= 0
 R401: - x[2,2,1] + group_exists[2] >= 0
 R402: - 5 x[2,2,2] - group_section[2] + element_section[2,2] >= -5
 R403: 5 x[2,2,2] - group_section[2] + element_section[2,2] <= 5
 R404: - x[2,2,2] + col_in_group[2,2] >= 0
 R405: - x[2,2,2] + level_in_group[2,2] >= 0
 R406: - x[2,2,2] + group_exists[2] >= 0
 R407: - 5 x[2,2,3] - group_section[2] + element_section[2,3] >= -5
 R408: 5 x[2,2,3] - group_section[2] + element_section[2,3] <= 5
 R409: - x[2,2,3] + col_in_group[2,2] >= 0
 R410: - x[2,2,3] + level_in_group[2,3] >= 0
 R411: - x[2,2,3] + group_exists[2] >= 0
 R412: - 5 x[2,2,4] - group_section[2] + element_section[2,4] >= -5
 R413: 5 x[2,2,4] - group_section[2] + element_section[2,4] <= 5
 R414: - x[2,2,4] + col_in_group[2,2] >= 0
 R415: - x[2,2,4] + level_in_group[2,4] >= 0
 R416: - x[2,2,4] + group_exists[2] >= 0
 R417: - 5 x[2,3,0] - group_section[2] + element_section[3,0] >= -5
 R418: 5 x[2,3,0] - group_section[2] + element_section[3,0] <= 5
 R419: - x[2,3,0] + col_in_group[2,3] >= 0
 R420: - x[2,3,0] + level_in_group[2,0] >= 0
 R421: - x[2,3,0] + group_exists[2] >= 0
 R422: - 5 x[2,3,1] - group_section[2] + element_section[3,1] >= -5
 R423: 5 x[2,3,1] - group_section[2] + element_section[3,1] <= 5
 R424: - x[2,3,1] + col_in_group[2,3] >= 0
 R425: - x[2,3,1] + level_in_group[2,1] >= 0
 R426: - x[2,3,1] + group_exists[2] >= 0
 R427: - 5 x[2,3,2] - group_section[2] + element_section[3,2] >= -5
 R428: 5 x[2,3,2] - group_section[2] + element_section[3,2] <= 5
 R429: - x[2,3,2] + col_in_group[2,3] >= 0
 R430: - x[2,3,2] + level_in_group[2,2] >= 0
 R431: - x[2,3,2] + group_exists[2] >= 0
 R432: - 5 x[2,3,3] - group_section[2] + element_section[3,3] >= -5
 R433: 5 x[2,3,3] - group_section[2] + element_section[3,3] <= 5
 R434: - x[2,3,3] + col_in_group[2,3] >= 0
 R435: - x[2,3,3] + level_in_group[2,3] >= 0
 R436: - x[2,3,3] + group_exists[2] >= 0
 R437: - 5 x[2,3,4] - group_section[2] + element_section[3,4] >= -5
 R438: 5 x[2,3,4] - group_section[2] + element_section[3,4] <= 5
 R439: - x[2,3,4] + col_in_group[2,3] >= 0
 R440: - x[2,3,4] + level_in_group[2,4] >= 0
 R441: - x[2,3,4] + group_exists[2] >= 0
 R442: - 5 x[2,4,0] - group_section[2] + element_section[4,0] >= -5
 R443: 5 x[2,4,0] - group_section[2] + element_section[4,0] <= 5
 R444: - x[2,4,0] + col_in_group[2,4] >= 0
 R445: - x[2,4,0] + level_in_group[2,0] >= 0
 R446: - x[2,4,0] + group_exists[2] >= 0
 R447: - 5 x[2,4,1] - group_section[2] + element_section[4,1] >= -5
 R448: 5 x[2,4,1] - group_section[2] + element_section[4,1] <= 5
 R449: - x[2,4,1] + col_in_group[2,4] >= 0
 R450: - x[2,4,1] + level_in_group[2,1] >= 0
 R451: - x[2,4,1] + group_exists[2] >= 0
 R452: - 5 x[2,4,2] - group_section[2] + element_section[4,2] >= -5
 R453: 5 x[2,4,2] - group_section[2] + element_section[4,2] <= 5
 R454: - x[2,4,2] + col_in_group[2,4] >= 0
 R455: - x[2,4,2] + level_in_group[2,2] >= 0
 R456: - x[2,4,2] + group_exists[2] >= 0
 R457: - 5 x[2,4,3] - group_section[2] + element_section[4,3] >= -5
 R458: 5 x[2,4,3] - group_section[2] + element_section[4,3] <= 5
 R459: - x[2,4,3] + col_in_group[2,4] >= 0
 R460: - x[2,4,3] + level_in_group[2,3] >= 0
 R461: - x[2,4,3] + group_exists[2] >= 0
 R462: - 5 x[2,4,4] - group_section[2] + element_section[4,4] >= -5
 R463: 5 x[2,4,4] - group_section[2] + element_section[4,4] <= 5
 R464: - x[2,4,4] + col_in_group[2,4] >= 0
 R465: - x[2,4,4] + level_in_group[2,4] >= 0
 R466: - x[2,4,4] + group_exists[2] >= 0
 R467: group_lb[3] - group_ub[3] + group_range[3] = 0
 R468: - 5 x[3,0,0] - group_section[3] + element_section[0,0] >= -5
 R469: 5 x[3,0,0] - group_section[3] + element_section[0,0] <= 5
 R470: - x[3,0,0] + col_in_group[3,0] >= 0
 R471: - x[3,0,0] + level_in_group[3,0] >= 0
 R472: - x[3,0,0] + group_exists[3] >= 0
 R473: - 5 x[3,0,1] - group_section[3] + element_section[0,1] >= -5
 R474: 5 x[3,0,1] - group_section[3] + element_section[0,1] <= 5
 R475: - x[3,0,1] + col_in_group[3,0] >= 0
 R476: - x[3,0,1] + level_in_group[3,1] >= 0
 R477: - x[3,0,1] + group_exists[3] >= 0
 R478: - 5 x[3,0,2] - group_section[3] + element_section[0,2] >= -5
 R479: 5 x[3,0,2] - group_section[3] + element_section[0,2] <= 5
 R480: - x[3,0,2] + col_in_group[3,0] >= 0
 R481: - x[3,0,2] + level_in_group[3,2] >= 0
 R482: - x[3,0,2] + group_exists[3] >= 0
 R483: - 5 x[3,0,3] - group_section[3] + element_section[0,3] >= -5
 R484: 5 x[3,0,3] - group_section[3] + element_section[0,3] <= 5
 R485: - x[3,0,3] + col_in_group[3,0] >= 0
 R486: - x[3,0,3] + level_in_group[3,3] >= 0
 R487: - x[3,0,3] + group_exists[3] >= 0
 R488: - 5 x[3,0,4] - group_section[3] + element_section[0,4] >= -5
 R489: 5 x[3,0,4] - group_section[3] + element_section[0,4] <= 5
 R490: - x[3,0,4] + col_in_group[3,0] >= 0
 R491: - x[3,0,4] + level_in_group[3,4] >= 0
 R492: - x[3,0,4] + group_exists[3] >= 0
 R493: - 5 x[3,1,0] - group_section[3] + element_section[1,0] >= -5
 R494: 5 x[3,1,0] - group_section[3] + element_section[1,0] <= 5
 R495: - x[3,1,0] + col_in_group[3,1] >= 0
 R496: - x[3,1,0] + level_in_group[3,0] >= 0
 R497: - x[3,1,0] + group_exists[3] >= 0
 R498: - 5 x[3,1,1] - group_section[3] + element_section[1,1] >= -5
 R499: 5 x[3,1,1] - group_section[3] + element_section[1,1] <= 5
 R500: - x[3,1,1] + col_in_group[3,1] >= 0
 R501: - x[3,1,1] + level_in_group[3,1] >= 0
 R502: - x[3,1,1] + group_exists[3] >= 0
 R503: - 5 x[3,1,2] - group_section[3] + element_section[1,2] >= -5
 R504: 5 x[3,1,2] - group_section[3] + element_section[1,2] <= 5
 R505: - x[3,1,2] + col_in_group[3,1] >= 0
 R506: - x[3,1,2] + level_in_group[3,2] >= 0
 R507: - x[3,1,2] + group_exists[3] >= 0
 R508: - 5 x[3,1,3] - group_section[3] + element_section[1,3] >= -5
 R509: 5 x[3,1,3] - group_section[3] + element_section[1,3] <= 5
 R510: - x[3,1,3] + col_in_group[3,1] >= 0
 R511: - x[3,1,3] + level_in_group[3,3] >= 0
 R512: - x[3,1,3] + group_exists[3] >= 0
 R513: - 5 x[3,1,4] - group_section[3] + element_section[1,4] >= -5
 R514: 5 x[3,1,4] - group_section[3] + element_section[1,4] <= 5
 R515: - x[3,1,4] + col_in_group[3,1] >= 0
 R516: - x[3,1,4] + level_in_group[3,4] >= 0
 R517: - x[3,1,4] + group_exists[3] >= 0
 R518: - 5 x[3,2,0] - group_section[3] + element_section[2,0] >= -5
 R519: 5 x[3,2,0] - group_section[3] + element_section[2,0] <= 5
 R520: - x[3,2,0] + col_in_group[3,2] >= 0
 R521: - x[3,2,0] + level_in_group[3,0] >= 0
 R522: - x[3,2,0] + group_exists[3] >= 0
 R523: - 5 x[3,2,1] - group_section[3] + element_section[2,1] >= -5
 R524: 5 x[3,2,1] - group_section[3] + element_section[2,1] <= 5
 R525: - x[3,2,1] + col_in_group[3,2] >= 0
 R526: - x[3,2,1] + level_in_group[3,1] >= 0
 R527: - x[3,2,1] + group_exists[3] >= 0
 R528: - 5 x[3,2,2] - group_section[3] + element_section[2,2] >= -5
 R529: 5 x[3,2,2] - group_section[3] + element_section[2,2] <= 5
 R530: - x[3,2,2] + col_in_group[3,2] >= 0
 R531: - x[3,2,2] + level_in_group[3,2] >= 0
 R532: - x[3,2,2] + group_exists[3] >= 0
 R533: - 5 x[3,2,3] - group_section[3] + element_section[2,3] >= -5
 R534: 5 x[3,2,3] - group_section[3] + element_section[2,3] <= 5
 R535: - x[3,2,3] + col_in_group[3,2] >= 0
 R536: - x[3,2,3] + level_in_group[3,3] >= 0
 R537: - x[3,2,3] + group_exists[3] >= 0
 R538: - 5 x[3,2,4] - group_section[3] + element_section[2,4] >= -5
 R539: 5 x[3,2,4] - group_section[3] + element_section[2,4] <= 5
 R540: - x[3,2,4] + col_in_group[3,2] >= 0
 R541: - x[3,2,4] + level_in_group[3,4] >= 0
 R542: - x[3,2,4] + group_exists[3] >= 0
 R543: - 5 x[3,3,0] - group_section[3] + element_section[3,0] >= -5
 R544: 5 x[3,3,0] - group_section[3] + element_section[3,0] <= 5
 R545: - x[3,3,0] + col_in_group[3,3] >= 0
 R546: - x[3,3,0] + level_in_group[3,0] >= 0
 R547: - x[3,3,0] + group_exists[3] >= 0
 R548: - 5 x[3,3,1] - group_section[3] + element_section[3,1] >= -5
 R549: 5 x[3,3,1] - group_section[3] + element_section[3,1] <= 5
 R550: - x[3,3,1] + col_in_group[3,3] >= 0
 R551: - x[3,3,1] + level_in_group[3,1] >= 0
 R552: - x[3,3,1] + group_exists[3] >= 0
 R553: - 5 x[3,3,2] - group_section[3] + element_section[3,2] >= -5
 R554: 5 x[3,3,2] - group_section[3] + element_section[3,2] <= 5
 R555: - x[3,3,2] + col_in_group[3,3] >= 0
 R556: - x[3,3,2] + level_in_group[3,2] >= 0
 R557: - x[3,3,2] + group_exists[3] >= 0
 R558: - 5 x[3,3,3] - group_section[3] + element_section[3,3] >= -5
 R559: 5 x[3,3,3] - group_section[3] + element_section[3,3] <= 5
 R560: - x[3,3,3] + col_in_group[3,3] >= 0
 R561: - x[3,3,3] + level_in_group[3,3] >= 0
 R562: - x[3,3,3] + group_exists[3] >= 0
 R563: - 5 x[3,3,4] - group_section[3] + element_section[3,4] >= -5
 R564: 5 x[3,3,4] - group_section[3] + element_section[3,4] <= 5
 R565: - x[3,3,4] + col_in_group[3,3] >= 0
 R566: - x[3,3,4] + level_in_group[3,4] >= 0
 R567: - x[3,3,4] + group_exists[3] >= 0
 R568: - 5 x[3,4,0] - group_section[3] + element_section[4,0] >= -5
 R569: 5 x[3,4,0] - group_section[3] + element_section[4,0] <= 5
 R570: - x[3,4,0] + col_in_group[3,4] >= 0
 R571: - x[3,4,0] + level_in_group[3,0] >= 0
 R572: - x[3,4,0] + group_exists[3] >= 0
 R573: - 5 x[3,4,1] - group_section[3] + element_section[4,1] >= -5
 R574: 5 x[3,4,1] - group_section[3] + element_section[4,1] <= 5
 R575: - x[3,4,1] + col_in_group[3,4] >= 0
 R576: - x[3,4,1] + level_in_group[3,1] >= 0
 R577: - x[3,4,1] + group_exists[3] >= 0
 R578: - 5 x[3,4,2] - group_section[3] + element_section[4,2] >= -5
 R579: 5 x[3,4,2] - group_section[3] + element_section[4,2] <= 5
 R580: - x[3,4,2] + col_in_group[3,4] >= 0
 R581: - x[3,4,2] + level_in_group[3,2] >= 0
 R582: - x[3,4,2] + group_exists[3] >= 0
 R583: - 5 x[3,4,3] - group_section[3] + element_section[4,3] >= -5
 R584: 5 x[3,4,3] - group_section[3] + element_section[4,3] <= 5
 R585: - x[3,4,3] + col_in_group[3,4] >= 0
 R586: - x[3,4,3] + level_in_group[3,3] >= 0
 R587: - x[3,4,3] + group_exists[3] >= 0
 R588: - 5 x[3,4,4] - group_section[3] + element_section[4,4] >= -5
 R589: 5 x[3,4,4] - group_section[3] + element_section[4,4] <= 5
 R590: - x[3,4,4] + col_in_group[3,4] >= 0
 R591: - x[3,4,4] + level_in_group[3,4] >= 0
 R592: - x[3,4,4] + group_exists[3] >= 0
 R593: group_lb[4] - group_ub[4] + group_range[4] = 0
 R594: - 5 x[4,0,0] - group_section[4] + element_section[0,0] >= -5
 R595: 5 x[4,0,0] - group_section[4] + element_section[0,0] <= 5
 R596: - x[4,0,0] + col_in_group[4,0] >= 0
 R597: - x[4,0,0] + level_in_group[4,0] >= 0
 R598: - x[4,0,0] + group_exists[4] >= 0
 R599: - 5 x[4,0,1] - group_section[4] + element_section[0,1] >= -5
 R600: 5 x[4,0,1] - group_section[4] + element_section[0,1] <= 5
 R601: - x[4,0,1] + col_in_group[4,0] >= 0
 R602: - x[4,0,1] + level_in_group[4,1] >= 0
 R603: - x[4,0,1] + group_exists[4] >= 0
 R604: - 5 x[4,0,2] - group_section[4] + element_section[0,2] >= -5
 R605: 5 x[4,0,2] - group_section[4] + element_section[0,2] <= 5
 R606: - x[4,0,2] + col_in_group[4,0] >= 0
 R607: - x[4,0,2] + level_in_group[4,2] >= 0
 R608: - x[4,0,2] + group_exists[4] >= 0
 R609: - 5 x[4,0,3] - group_section[4] + element_section[0,3] >= -5
 R610: 5 x[4,0,3] - group_section[4] + element_section[0,3] <= 5
 R611: - x[4,0,3] + col_in_group[4,0] >= 0
 R612: - x[4,0,3] + level_in_group[4,3] >= 0
 R613: - x[4,0,3] + group_exists[4] >= 0
 R614: - 5 x[4,0,4] - group_section[4] + element_section[0,4] >= -5
 R615: 5 x[4,0,4] - group_section[4] + element_section[0,4] <= 5
 R616: - x[4,0,4] + col_in_group[4,0] >= 0
 R617: - x[4,0,4] + level_in_group[4,4] >= 0
 R618: - x[4,0,4] + group_exists[4] >= 0
 R619: - 5 x[4,1,0] - group_section[4] + element_section[1,0] >= -5
 R620: 5 x[4,1,0] - group_section[4] + element_section[1,0] <= 5
 R621: - x[4,1,0] + col_in_group[4,1] >= 0
 R622: - x[4,1,0] + level_in_group[4,0] >= 0
 R623: - x[4,1,0] + group_exists[4] >= 0
 R624: - 5 x[4,1,1] - group_section[4] + element_section[1,1] >= -5
 R625: 5 x[4,1,1] - group_section[4] + element_section[1,1] <= 5
 R626: - x[4,1,1] + col_in_group[4,1] >= 0
 R627: - x[4,1,1] + level_in_group[4,1] >= 0
 R628: - x[4,1,1] + group_exists[4] >= 0
 R629: - 5 x[4,1,2] - group_section[4] + element_section[1,2] >= -5
 R630: 5 x[4,1,2] - group_section[4] + element_section[1,2] <= 5
 R631: - x[4,1,2] + col_in_group[4,1] >= 0
 R632: - x[4,1,2] + level_in_group[4,2] >= 0
 R633: - x[4,1,2] + group_exists[4] >= 0
 R634: - 5 x[4,1,3] - group_section[4] + element_section[1,3] >= -5
 R635: 5 x[4,1,3] - group_section[4] + element_section[1,3] <= 5
 R636: - x[4,1,3] + col_in_group[4,1] >= 0
 R637: - x[4,1,3] + level_in_group[4,3] >= 0
 R638: - x[4,1,3] + group_exists[4] >= 0
 R639: - 5 x[4,1,4] - group_section[4] + element_section[1,4] >= -5
 R640: 5 x[4,1,4] - group_section[4] + element_section[1,4] <= 5
 R641: - x[4,1,4] + col_in_group[4,1] >= 0
 R642: - x[4,1,4] + level_in_group[4,4] >= 0
 R643: - x[4,1,4] + group_exists[4] >= 0
 R644: - 5 x[4,2,0] - group_section[4] + element_section[2,0] >= -5
 R645: 5 x[4,2,0] - group_section[4] + element_section[2,0] <= 5
 R646: - x[4,2,0] + col_in_group[4,2] >= 0
 R647: - x[4,2,0] + level_in_group[4,0] >= 0
 R648: - x[4,2,0] + group_exists[4] >= 0
 R649: - 5 x[4,2,1] - group_section[4] + element_section[2,1] >= -5
 R650: 5 x[4,2,1] - group_section[4] + element_section[2,1] <= 5
 R651: - x[4,2,1] + col_in_group[4,2] >= 0
 R652: - x[4,2,1] + level_in_group[4,1] >= 0
 R653: - x[4,2,1] + group_exists[4] >= 0
 R654: - 5 x[4,2,2] - group_section[4] + element_section[2,2] >= -5
 R655: 5 x[4,2,2] - group_section[4] + element_section[2,2] <= 5
 R656: - x[4,2,2] + col_in_group[4,2] >= 0
 R657: - x[4,2,2] + level_in_group[4,2] >= 0
 R658: - x[4,2,2] + group_exists[4] >= 0
 R659: - 5 x[4,2,3] - group_section[4] + element_section[2,3] >= -5
 R660: 5 x[4,2,3] - group_section[4] + element_section[2,3] <= 5
 R661: - x[4,2,3] + col_in_group[4,2] >= 0
 R662: - x[4,2,3] + level_in_group[4,3] >= 0
 R663: - x[4,2,3] + group_exists[4] >= 0
 R664: - 5 x[4,2,4] - group_section[4] + element_section[2,4] >= -5
 R665: 5 x[4,2,4] - group_section[4] + element_section[2,4] <= 5
 R666: - x[4,2,4] + col_in_group[4,2] >= 0
 R667: - x[4,2,4] + level_in_group[4,4] >= 0
 R668: - x[4,2,4] + group_exists[4] >= 0
 R669: - 5 x[4,3,0] - group_section[4] + element_section[3,0] >= -5
 R670: 5 x[4,3,0] - group_section[4] + element_section[3,0] <= 5
 R671: - x[4,3,0] + col_in_group[4,3] >= 0
 R672: - x[4,3,0] + level_in_group[4,0] >= 0
 R673: - x[4,3,0] + group_exists[4] >= 0
 R674: - 5 x[4,3,1] - group_section[4] + element_section[3,1] >= -5
 R675: 5 x[4,3,1] - group_section[4] + element_section[3,1] <= 5
 R676: - x[4,3,1] + col_in_group[4,3] >= 0
 R677: - x[4,3,1] + level_in_group[4,1] >= 0
 R678: - x[4,3,1] + group_exists[4] >= 0
 R679: - 5 x[4,3,2] - group_section[4] + element_section[3,2] >= -5
 R680: 5 x[4,3,2] - group_section[4] + element_section[3,2] <= 5
 R681: - x[4,3,2] + col_in_group[4,3] >= 0
 R682: - x[4,3,2] + level_in_group[4,2] >= 0
 R683: - x[4,3,2] + group_exists[4] >= 0
 R684: - 5 x[4,3,3] - group_section[4] + element_section[3,3] >= -5
 R685: 5 x[4,3,3] - group_section[4] + element_section[3,3] <= 5
 R686: - x[4,3,3] + col_in_group[4,3] >= 0
 R687: - x[4,3,3] + level_in_group[4,3] >= 0
 R688: - x[4,3,3] + group_exists[4] >= 0
 R689: - 5 x[4,3,4] - group_section[4] + element_section[3,4] >= -5
 R690: 5 x[4,3,4] - group_section[4] + element_section[3,4] <= 5
 R691: - x[4,3,4] + col_in_group[4,3] >= 0
 R692: - x[4,3,4] + level_in_group[4,4] >= 0
 R693: - x[4,3,4] + group_exists[4] >= 0
 R694: - 5 x[4,4,0] - group_section[4] + element_section[4,0] >= -5
 R695: 5 x[4,4,0] - group_section[4] + element_section[4,0] <= 5
 R696: - x[4,4,0] + col_in_group[4,4] >= 0
 R697: - x[4,4,0] + level_in_group[4,0] >= 0
 R698: - x[4,4,0] + group_exists[4] >= 0
 R699: - 5 x[4,4,1] - group_section[4] + element_section[4,1] >= -5
 R700: 5 x[4,4,1] - group_section[4] + element_section[4,1] <= 5
 R701: - x[4,4,1] + col_in_group[4,4] >= 0
 R702: - x[4,4,1] + level_in_group[4,1] >= 0
 R703: - x[4,4,1] + group_exists[4] >= 0
 R704: - 5 x[4,4,2] - group_section[4] + element_section[4,2] >= -5
 R705: 5 x[4,4,2] - group_section[4] + element_section[4,2] <= 5
 R706: - x[4,4,2] + col_in_group[4,4] >= 0
 R707: - x[4,4,2] + level_in_group[4,2] >= 0
 R708: - x[4,4,2] + group_exists[4] >= 0
 R709: - 5 x[4,4,3] - group_section[4] + element_section[4,3] >= -5
 R710: 5 x[4,4,3] - group_section[4] + element_section[4,3] <= 5
 R711: - x[4,4,3] + col_in_group[4,4] >= 0
 R712: - x[4,4,3] + level_in_group[4,3] >= 0
 R713: - x[4,4,3] + group_exists[4] >= 0
 R714: - 5 x[4,4,4] - group_section[4] + element_section[4,4] >= -5
 R715: 5 x[4,4,4] - group_section[4] + element_section[4,4] <= 5
 R716: - x[4,4,4] + col_in_group[4,4] >= 0
 R717: - x[4,4,4] + level_in_group[4,4] >= 0
 R718: - x[4,4,4] + group_exists[4] >= 0
 R719: - group_ub[0] <= 0
 R720: - 7 level_in_group[0,0] - group_lb[0] >= -7
 R721: group_lb[0] + 7 Zl[0,0] >= 1
 R722: - group_lb[0] - 7 Zl[0,0] >= -7
 R723: - group_ub[0] + 7 Zu[0,0] >= 1
 R724: group_ub[0] - 7 Zu[0,0] >= -7
 R725: level_in_group[0,0] - Zu[0,0] - Zl[0,0] >= -1
 R726: - x[0,0,0] + col_in_group[0,0] + level_in_group[0,0] <= 1
 R727: - x[0,1,0] + col_in_group[0,1] + level_in_group[0,0] <= 1
 R728: - x[0,2,0] + col_in_group[0,2] + level_in_group[0,0] <= 1
 R729: - x[0,3,0] + col_in_group[0,3] + level_in_group[0,0] <= 1
 R730: - x[0,4,0] + col_in_group[0,4] + level_in_group[0,0] <= 1
 R731: level_in_group[0,1] - group_ub[0] <= 0
 R732: - 6 level_in_group[0,1] - group_lb[0] >= -7
 R733: group_lb[0] + 7 Zl[0,1] >= 2
 R734: - group_lb[0] - 7 Zl[0,1] >= -8
 R735: - group_ub[0] + 7 Zu[0,1] >= 0
 R736: group_ub[0] - 7 Zu[0,1] >= -6
 R737: level_in_group[0,1] - Zu[0,1] - Zl[0,1] >= -1
 R738: - x[0,0,1] + col_in_group[0,0] + level_in_group[0,1] <= 1
 R739: - x[0,1,1] + col_in_group[0,1] + level_in_group[0,1] <= 1
 R740: - x[0,2,1] + col_in_group[0,2] + level_in_group[0,1] <= 1
 R741: - x[0,3,1] + col_in_group[0,3] + level_in_group[0,1] <= 1
 R742: - x[0,4,1] + col_in_group[0,4] + level_in_group[0,1] <= 1
 R743: 2 level_in_group[0,2] - group_ub[0] <= 0
 R744: - 5 level_in_group[0,2] - group_lb[0] >= -7
 R745: group_lb[0] + 7 Zl[0,2] >= 3
 R746: - group_lb[0] - 7 Zl[0,2] >= -9
 R747: - group_ub[0] + 7 Zu[0,2] >= -1
 R748: group_ub[0] - 7 Zu[0,2] >= -5
 R749: level_in_group[0,2] - Zu[0,2] - Zl[0,2] >= -1
 R750: - x[0,0,2] + col_in_group[0,0] + level_in_group[0,2] <= 1
 R751: - x[0,1,2] + col_in_group[0,1] + level_in_group[0,2] <= 1
 R752: - x[0,2,2] + col_in_group[0,2] + level_in_group[0,2] <= 1
 R753: - x[0,3,2] + col_in_group[0,3] + level_in_group[0,2] <= 1
 R754: - x[0,4,2] + col_in_group[0,4] + level_in_group[0,2] <= 1
 R755: 3 level_in_group[0,3] - group_ub[0] <= 0
 R756: - 4 level_in_group[0,3] - group_lb[0] >= -7
 R757: group_lb[0] + 7 Zl[0,3] >= 4
 R758: - group_lb[0] - 7 Zl[0,3] >= -10
 R759: - group_ub[0] + 7 Zu[0,3] >= -2
 R760: group_ub[0] - 7 Zu[0,3] >= -4
 R761: level_in_group[0,3] - Zu[0,3] - Zl[0,3] >= -1
 R762: - x[0,0,3] + col_in_group[0,0] + level_in_group[0,3] <= 1
 R763: - x[0,1,3] + col_in_group[0,1] + level_in_group[0,3] <= 1
 R764: - x[0,2,3] + col_in_group[0,2] + level_in_group[0,3] <= 1
 R765: - x[0,3,3] + col_in_group[0,3] + level_in_group[0,3] <= 1
 R766: - x[0,4,3] + col_in_group[0,4] + level_in_group[0,3] <= 1
 R767: 4 level_in_group[0,4] - group_ub[0] <= 0
 R768: - 3 level_in_group[0,4] - group_lb[0] >= -7
 R769: group_lb[0] + 7 Zl[0,4] >= 5
 R770: - group_lb[0] - 7 Zl[0,4] >= -11
 R771: - group_ub[0] + 7 Zu[0,4] >= -3
 R772: group_ub[0] - 7 Zu[0,4] >= -3
 R773: level_in_group[0,4] - Zu[0,4] - Zl[0,4] >= -1
 R774: - x[0,0,4] + col_in_group[0,0] + level_in_group[0,4] <= 1
 R775: - x[0,1,4] + col_in_group[0,1] + level_in_group[0,4] <= 1
 R776: - x[0,2,4] + col_in_group[0,2] + level_in_group[0,4] <= 1
 R777: - x[0,3,4] + col_in_group[0,3] + level_in_group[0,4] <= 1
 R778: - x[0,4,4] + col_in_group[0,4] + level_in_group[0,4] <= 1
 R779: - group_ub[1] <= 0
 R780: - 7 level_in_group[1,0] - group_lb[1] >= -7
 R781: group_lb[1] + 7 Zl[1,0] >= 1
 R782: - group_lb[1] - 7 Zl[1,0] >= -7
 R783: - group_ub[1] + 7 Zu[1,0] >= 1
 R784: group_ub[1] - 7 Zu[1,0] >= -7
 R785: level_in_group[1,0] - Zu[1,0] - Zl[1,0] >= -1
 R786: - x[1,0,0] + col_in_group[1,0] + level_in_group[1,0] <= 1
 R787: - x[1,1,0] + col_in_group[1,1] + level_in_group[1,0] <= 1
 R788: - x[1,2,0] + col_in_group[1,2] + level_in_group[1,0] <= 1
 R789: - x[1,3,0] + col_in_group[1,3] + level_in_group[1,0] <= 1
 R790: - x[1,4,0] + col_in_group[1,4] + level_in_group[1,0] <= 1
 R791: level_in_group[1,1] - group_ub[1] <= 0
 R792: - 6 level_in_group[1,1] - group_lb[1] >= -7
 R793: group_lb[1] + 7 Zl[1,1] >= 2
 R794: - group_lb[1] - 7 Zl[1,1] >= -8
 R795: - group_ub[1] + 7 Zu[1,1] >= 0
 R796: group_ub[1] - 7 Zu[1,1] >= -6
 R797: level_in_group[1,1] - Zu[1,1] - Zl[1,1] >= -1
 R798: - x[1,0,1] + col_in_group[1,0] + level_in_group[1,1] <= 1
 R799: - x[1,1,1] + col_in_group[1,1] + level_in_group[1,1] <= 1
 R800: - x[1,2,1] + col_in_group[1,2] + level_in_group[1,1] <= 1
 R801: - x[1,3,1] + col_in_group[1,3] + level_in_group[1,1] <= 1
 R802: - x[1,4,1] + col_in_group[1,4] + level_in_group[1,1] <= 1
 R803: 2 level_in_group[1,2] - group_ub[1] <= 0
 R804: - 5 level_in_group[1,2] - group_lb[1] >= -7
 R805: group_lb[1] + 7 Zl[1,2] >= 3
 R806: - group_lb[1] - 7 Zl[1,2] >= -9
 R807: - group_ub[1] + 7 Zu[1,2] >= -1
 R808: group_ub[1] - 7 Zu[1,2] >= -5
 R809: level_in_group[1,2] - Zu[1,2] - Zl[1,2] >= -1
 R810: - x[1,0,2] + col_in_group[1,0] + level_in_group[1,2] <= 1
 R811: - x[1,1,2] + col_in_group[1,1] + level_in_group[1,2] <= 1
 R812: - x[1,2,2] + col_in_group[1,2] + level_in_group[1,2] <= 1
 R813: - x[1,3,2] + col_in_group[1,3] + level_in_group[1,2] <= 1
 R814: - x[1,4,2] + col_in_group[1,4] + level_in_group[1,2] <= 1
 R815: 3 level_in_group[1,3] - group_ub[1] <= 0
 R816: - 4 level_in_group[1,3] - group_lb[1] >= -7
 R817: group_lb[1] + 7 Zl[1,3] >= 4
 R818: - group_lb[1] - 7 Zl[1,3] >= -10
 R819: - group_ub[1] + 7 Zu[1,3] >= -2
 R820: group_ub[1] - 7 Zu[1,3] >= -4
 R821: level_in_group[1,3] - Zu[1,3] - Zl[1,3] >= -1
 R822: - x[1,0,3] + col_in_group[1,0] + level_in_group[1,3] <= 1
 R823: - x[1,1,3] + col_in_group[1,1] + level_in_group[1,3] <= 1
 R824: - x[1,2,3] + col_in_group[1,2] + level_in_group[1,3] <= 1
 R825: - x[1,3,3] + col_in_group[1,3] + level_in_group[1,3] <= 1
 R826: - x[1,4,3] + col_in_group[1,4] + level_in_group[1,3] <= 1
 R827: 4 level_in_group[1,4] - group_ub[1] <= 0
 R828: - 3 level_in_group[1,4] - group_lb[1] >= -7
 R829: group_lb[1] + 7 Zl[1,4] >= 5
 R830: - group_lb[1] - 7 Zl[1,4] >= -11
 R831: - group_ub[1] + 7 Zu[1,4] >= -3
 R832: group_ub[1] - 7 Zu[1,4] >= -3
 R833: level_in_group[1,4] - Zu[1,4] - Zl[1,4] >= -1
 R834: - x[1,0,4] + col_in_group[1,0] + level_in_group[1,4] <= 1
 R835: - x[1,1,4] + col_in_group[1,1] + level_in_group[1,4] <= 1
 R836: - x[1,2,4] + col_in_group[1,2] + level_in_group[1,4] <= 1
 R837: - x[1,3,4] + col_in_group[1,3] + level_in_group[1,4] <= 1
 R838: - x[1,4,4] + col_in_group[1,4] + level_in_group[1,4] <= 1
 R839: - group_ub[2] <= 0
 R840: - 7 level_in_group[2,0] - group_lb[2] >= -7
 R841: group_lb[2] + 7 Zl[2,0] >= 1
 R842: - group_lb[2] - 7 Zl[2,0] >= -7
 R843: - group_ub[2] + 7 Zu[2,0] >= 1
 R844: group_ub[2] - 7 Zu[2,0] >= -7
 R845: level_in_group[2,0] - Zu[2,0] - Zl[2,0] >= -1
 R846: - x[2,0,0] + col_in_group[2,0] + level_in_group[2,0] <= 1
 R847: - x[2,1,0] + col_in_group[2,1] + level_in_group[2,0] <= 1
 R848: - x[2,2,0] + col_in_group[2,2] + level_in_group[2,0] <= 1
 R849: - x[2,3,0] + col_in_group[2,3] + level_in_group[2,0] <= 1
 R850: - x[2,4,0] + col_in_group[2,4] + level_in_group[2,0] <= 1
 R851: level_in_group[2,1] - group_ub[2] <= 0
 R852: - 6 level_in_group[2,1] - group_lb[2] >= -7
 R853: group_lb[2] + 7 Zl[2,1] >= 2
 R854: - group_lb[2] - 7 Zl[2,1] >= -8
 R855: - group_ub[2] + 7 Zu[2,1] >= 0
 R856: group_ub[2] - 7 Zu[2,1] >= -6
 R857: level_in_group[2,1] - Zu[2,1] - Zl[2,1] >= -1
 R858: - x[2,0,1] + col_in_group[2,0] + level_in_group[2,1] <= 1
 R859: - x[2,1,1] + col_in_group[2,1] + level_in_group[2,1] <= 1
 R860: - x[2,2,1] + col_in_group[2,2] + level_in_group[2,1] <= 1
 R861: - x[2,3,1] + col_in_group[2,3] + level_in_group[2,1] <= 1
 R862: - x[2,4,1] + col_in_group[2,4] + level_in_group[2,1] <= 1
 R863: 2 level_in_group[2,2] - group_ub[2] <= 0
 R864: - 5 level_in_group[2,2] - group_lb[2] >= -7
 R865: group_lb[2] + 7 Zl[2,2] >= 3
 R866: - group_lb[2] - 7 Zl[2,2] >= -9
 R867: - group_ub[2] + 7 Zu[2,2] >= -1
 R868: group_ub[2] - 7 Zu[2,2] >= -5
 R869: level_in_group[2,2] - Zu[2,2] - Zl[2,2] >= -1
 R870: - x[2,0,2] + col_in_group[2,0] + level_in_group[2,2] <= 1
 R871: - x[2,1,2] + col_in_group[2,1] + level_in_group[2,2] <= 1
 R872: - x[2,2,2] + col_in_group[2,2] + level_in_group[2,2] <= 1
 R873: - x[2,3,2] + col_in_group[2,3] + level_in_group[2,2] <= 1
 R874: - x[2,4,2] + col_in_group[2,4] + level_in_group[2,2] <= 1
 R875: 3 level_in_group[2,3] - group_ub[2] <= 0
 R876: - 4 level_in_group[2,3] - group_lb[2] >= -7
 R877: group_lb[2] + 7 Zl[2,3] >= 4
 R878: - group_lb[2] - 7 Zl[2,3] >= -10
 R879: - group_ub[2] + 7 Zu[2,3] >= -2
 R880: group_ub[2] - 7 Zu[2,3] >= -4
 R881: level_in_group[2,3] - Zu[2,3] - Zl[2,3] >= -1
 R882: - x[2,0,3] + col_in_group[2,0] + level_in_group[2,3] <= 1
 R883: - x[2,1,3] + col_in_group[2,1] + level_in_group[2,3] <= 1
 R884: - x[2,2,3] + col_in_group[2,2] + level_in_group[2,3] <= 1
 R885: - x[2,3,3] + col_in_group[2,3] + level_in_group[2,3] <= 1
 R886: - x[2,4,3] + col_in_group[2,4] + level_in_group[2,3] <= 1
 R887: 4 level_in_group[2,4] - group_ub[2] <= 0
 R888: - 3 level_in_group[2,4] - group_lb[2] >= -7
 R889: group_lb[2] + 7 Zl[2,4] >= 5
 R890: - group_lb[2] - 7 Zl[2,4] >= -11
 R891: - group_ub[2] + 7 Zu[2,4] >= -3
 R892: group_ub[2] - 7 Zu[2,4] >= -3
 R893: level_in_group[2,4] - Zu[2,4] - Zl[2,4] >= -1
 R894: - x[2,0,4] + col_in_group[2,0] + level_in_group[2,4] <= 1
 R895: - x[2,1,4] + col_in_group[2,1] + level_in_group[2,4] <= 1
 R896: - x[2,2,4] + col_in_group[2,2] + level_in_group[2,4] <= 1
 R897: - x[2,3,4] + col_in_group[2,3] + level_in_group[2,4] <= 1
 R898: - x[2,4,4] + col_in_group[2,4] + level_in_group[2,4] <= 1
 R899: - group_ub[3] <= 0
 R900: - 7 level_in_group[3,0] - group_lb[3] >= -7
 R901: group_lb[3] + 7 Zl[3,0] >= 1
 R902: - group_lb[3] - 7 Zl[3,0] >= -7
 R903: - group_ub[3] + 7 Zu[3,0] >= 1
 R904: group_ub[3] - 7 Zu[3,0] >= -7
 R905: level_in_group[3,0] - Zu[3,0] - Zl[3,0] >= -1
 R906: - x[3,0,0] + col_in_group[3,0] + level_in_group[3,0] <= 1
 R907: - x[3,1,0] + col_in_group[3,1] + level_in_group[3,0] <= 1
 R908: - x[3,2,0] + col_in_group[3,2] + level_in_group[3,0] <= 1
 R909: - x[3,3,0] + col_in_group[3,3] + level_in_group[3,0] <= 1
 R910: - x[3,4,0] + col_in_group[3,4] + level_in_group[3,0] <= 1
 R911: level_in_group[3,1] - group_ub[3] <= 0
 R912: - 6 level_in_group[3,1] - group_lb[3] >= -7
 R913: group_lb[3] + 7 Zl[3,1] >= 2
 R914: - group_lb[3] - 7 Zl[3,1] >= -8
 R915: - group_ub[3] + 7 Zu[3,1] >= 0
 R916: group_ub[3] - 7 Zu[3,1] >= -6
 R917: level_in_group[3,1] - Zu[3,1] - Zl[3,1] >= -1
 R918: - x[3,0,1] + col_in_group[3,0] + level_in_group[3,1] <= 1
 R919: - x[3,1,1] + col_in_group[3,1] + level_in_group[3,1] <= 1
 R920: - x[3,2,1] + col_in_group[3,2] + level_in_group[3,1] <= 1
 R921: - x[3,3,1] + col_in_group[3,3] + level_in_group[3,1] <= 1
 R922: - x[3,4,1] + col_in_group[3,4] + level_in_group[3,1] <= 1
 R923: 2 level_in_group[3,2] - group_ub[3] <= 0
 R924: - 5 level_in_group[3,2] - group_lb[3] >= -7
 R925: group_lb[3] + 7 Zl[3,2] >= 3
 R926: - group_lb[3] - 7 Zl[3,2] >= -9
 R927: - group_ub[3] + 7 Zu[3,2] >= -1
 R928: group_ub[3] - 7 Zu[3,2] >= -5
 R929: level_in_group[3,2] - Zu[3,2] - Zl[3,2] >= -1
 R930: - x[3,0,2] + col_in_group[3,0] + level_in_group[3,2] <= 1
 R931: - x[3,1,2] + col_in_group[3,1] + level_in_group[3,2] <= 1
 R932: - x[3,2,2] + col_in_group[3,2] + level_in_group[3,2] <= 1
 R933: - x[3,3,2] + col_in_group[3,3] + level_in_group[3,2] <= 1
 R934: - x[3,4,2] + col_in_group[3,4] + level_in_group[3,2] <= 1
 R935: 3 level_in_group[3,3] - group_ub[3] <= 0
 R936: - 4 level_in_group[3,3] - group_lb[3] >= -7
 R937: group_lb[3] + 7 Zl[3,3] >= 4
 R938: - group_lb[3] - 7 Zl[3,3] >= -10
 R939: - group_ub[3] + 7 Zu[3,3] >= -2
 R940: group_ub[3] - 7 Zu[3,3] >= -4
 R941: level_in_group[3,3] - Zu[3,3] - Zl[3,3] >= -1
 R942: - x[3,0,3] + col_in_group[3,0] + level_in_group[3,3] <= 1
 R943: - x[3,1,3] + col_in_group[3,1] + level_in_group[3,3] <= 1
 R944: - x[3,2,3] + col_in_group[3,2] + level_in_group[3,3] <= 1
 R945: - x[3,3,3] + col_in_group[3,3] + level_in_group[3,3] <= 1
 R946: - x[3,4,3] + col_in_group[3,4] + level_in_group[3,3] <= 1
 R947: 4 level_in_group[3,4] - group_ub[3] <= 0
 R948: - 3 level_in_group[3,4] - group_lb[3] >= -7
 R949: group_lb[3] + 7 Zl[3,4] >= 5
 R950: - group_lb[3] - 7 Zl[3,4] >= -11
 R951: - group_ub[3] + 7 Zu[3,4] >= -3
 R952: group_ub[3] - 7 Zu[3,4] >= -3
 R953: level_in_group[3,4] - Zu[3,4] - Zl[3,4] >= -1
 R954: - x[3,0,4] + col_in_group[3,0] + level_in_group[3,4] <= 1
 R955: - x[3,1,4] + col_in_group[3,1] + level_in_group[3,4] <= 1
 R956: - x[3,2,4] + col_in_group[3,2] + level_in_group[3,4] <= 1
 R957: - x[3,3,4] + col_in_group[3,3] + level_in_group[3,4] <= 1
 R958: - x[3,4,4] + col_in_group[3,4] + level_in_group[3,4] <= 1
 R959: - group_ub[4] <= 0
 R960: - 7 level_in_group[4,0] - group_lb[4] >= -7
 R961: group_lb[4] + 7 Zl[4,0] >= 1
 R962: - group_lb[4] - 7 Zl[4,0] >= -7
 R963: - group_ub[4] + 7 Zu[4,0] >= 1
 R964: group_ub[4] - 7 Zu[4,0] >= -7
 R965: level_in_group[4,0] - Zu[4,0] - Zl[4,0] >= -1
 R966: - x[4,0,0] + col_in_group[4,0] + level_in_group[4,0] <= 1
 R967: - x[4,1,0] + col_in_group[4,1] + level_in_group[4,0] <= 1
 R968: - x[4,2,0] + col_in_group[4,2] + level_in_group[4,0] <= 1
 R969: - x[4,3,0] + col_in_group[4,3] + level_in_group[4,0] <= 1
 R970: - x[4,4,0] + col_in_group[4,4] + level_in_group[4,0] <= 1
 R971: level_in_group[4,1] - group_ub[4] <= 0
 R972: - 6 level_in_group[4,1] - group_lb[4] >= -7
 R973: group_lb[4] + 7 Zl[4,1] >= 2
 R974: - group_lb[4] - 7 Zl[4,1] >= -8
 R975: - group_ub[4] + 7 Zu[4,1] >= 0
 R976: group_ub[4] - 7 Zu[4,1] >= -6
 R977: level_in_group[4,1] - Zu[4,1] - Zl[4,1] >= -1
 R978: - x[4,0,1] + col_in_group[4,0] + level_in_group[4,1] <= 1
 R979: - x[4,1,1] + col_in_group[4,1] + level_in_group[4,1] <= 1
 R980: - x[4,2,1] + col_in_group[4,2] + level_in_group[4,1] <= 1
 R981: - x[4,3,1] + col_in_group[4,3] + level_in_group[4,1] <= 1
 R982: - x[4,4,1] + col_in_group[4,4] + level_in_group[4,1] <= 1
 R983: 2 level_in_group[4,2] - group_ub[4] <= 0
 R984: - 5 level_in_group[4,2] - group_lb[4] >= -7
 R985: group_lb[4] + 7 Zl[4,2] >= 3
 R986: - group_lb[4] - 7 Zl[4,2] >= -9
 R987: - group_ub[4] + 7 Zu[4,2] >= -1
 R988: group_ub[4] - 7 Zu[4,2] >= -5
 R989: level_in_group[4,2] - Zu[4,2] - Zl[4,2] >= -1
 R990: - x[4,0,2] + col_in_group[4,0] + level_in_group[4,2] <= 1
 R991: - x[4,1,2] + col_in_group[4,1] + level_in_group[4,2] <= 1
 R992: - x[4,2,2] + col_in_group[4,2] + level_in_group[4,2] <= 1
 R993: - x[4,3,2] + col_in_group[4,3] + level_in_group[4,2] <= 1
 R994: - x[4,4,2] + col_in_group[4,4] + level_in_group[4,2] <= 1
 R995: 3 level_in_group[4,3] - group_ub[4] <= 0
 R996: - 4 level_in_group[4,3] - group_lb[4] >= -7
 R997: group_lb[4] + 7 Zl[4,3] >= 4
 R998: - group_lb[4] - 7 Zl[4,3] >= -10
 R999: - group_ub[4] + 7 Zu[4,3] >= -2
 R1000: group_ub[4] - 7 Zu[4,3] >= -4
 R1001: level_in_group[4,3] - Zu[4,3] - Zl[4,3] >= -1
 R1002: - x[4,0,3] + col_in_group[4,0] + level_in_group[4,3] <= 1
 R1003: - x[4,1,3] + col_in_group[4,1] + level_in_group[4,3] <= 1
 R1004: - x[4,2,3] + col_in_group[4,2] + level_in_group[4,3] <= 1
 R1005: - x[4,3,3] + col_in_group[4,3] + level_in_group[4,3] <= 1
 R1006: - x[4,4,3] + col_in_group[4,4] + level_in_group[4,3] <= 1
 R1007: 4 level_in_group[4,4] - group_ub[4] <= 0
 R1008: - 3 level_in_group[4,4] - group_lb[4] >= -7
 R1009: group_lb[4] + 7 Zl[4,4] >= 5
 R1010: - group_lb[4] - 7 Zl[4,4] >= -11
 R1011: - group_ub[4] + 7 Zu[4,4] >= -3
 R1012: group_ub[4] - 7 Zu[4,4] >= -3
 R1013: level_in_group[4,4] - Zu[4,4] - Zl[4,4] >= -1
 R1014: - x[4,0,4] + col_in_group[4,0] + level_in_group[4,4] <= 1
 R1015: - x[4,1,4] + col_in_group[4,1] + level_in_group[4,4] <= 1
 R1016: - x[4,2,4] + col_in_group[4,2] + level_in_group[4,4] <= 1
 R1017: - x[4,3,4] + col_in_group[4,3] + level_in_group[4,4] <= 1
 R1018: - x[4,4,4] + col_in_group[4,4] + level_in_group[4,4] <= 1
 GC0: y[0,0,0] = 1 -> element_section[0,0] = 0
 GC1: y[0,0,0] = 1 -> cost_section[0,0] = 0
 GC2: y[0,0,1] = 1 -> element_section[0,0] = 1
 GC3: y[0,0,1] = 1 -> cost_section[0,0] = 1131
 GC4: y[0,0,2] = 1 -> element_section[0,0] = 2
 GC5: y[0,0,2] = 1 -> cost_section[0,0] = 1876
 GC6: y[0,0,3] = 1 -> element_section[0,0] = 3
 GC7: y[0,0,3] = 1 -> cost_section[0,0] = 4140
 GC8: y[0,0,4] = 1 -> element_section[0,0] = 4
 GC9: y[0,0,4] = 1 -> cost_section[0,0] = 4304
 GC10: y[0,0,5] = 1 -> element_section[0,0] = 5
 GC11: y[0,0,5] = 1 -> cost_section[0,0] = 4944
 GC12: y[1,0,0] = 1 -> element_section[1,0] = 0
 GC13: y[1,0,0] = 1 -> cost_section[1,0] = 0
 GC14: y[1,0,1] = 1 -> element_section[1,0] = 1
 GC15: y[1,0,1] = 1 -> cost_section[1,0] = 1131
 GC16: y[1,0,2] = 1 -> element_section[1,0] = 2
 GC17: y[1,0,2] = 1 -> cost_section[1,0] = 1876
 GC18: y[1,0,3] = 1 -> element_section[1,0] = 3
 GC19: y[1,0,3] = 1 -> cost_section[1,0] = 4140
 GC20: y[1,0,4] = 1 -> element_section[1,0] = 4
 GC21: y[1,0,4] = 1 -> cost_section[1,0] = 4304
 GC22: y[1,0,5] = 1 -> element_section[1,0] = 5
 GC23: y[1,0,5] = 1 -> cost_section[1,0] = 4944
 GC24: y[2,0,0] = 1 -> element_section[2,0] = 0
 GC25: y[2,0,0] = 1 -> cost_section[2,0] = 0
 GC26: y[2,0,1] = 1 -> element_section[2,0] = 1
 GC27: y[2,0,1] = 1 -> cost_section[2,0] = 1131
 GC28: y[2,0,2] = 1 -> element_section[2,0] = 2
 GC29: y[2,0,2] = 1 -> cost_section[2,0] = 1876
 GC30: y[2,0,3] = 1 -> element_section[2,0] = 3
 GC31: y[2,0,3] = 1 -> cost_section[2,0] = 4140
 GC32: y[2,0,4] = 1 -> element_section[2,0] = 4
 GC33: y[2,0,4] = 1 -> cost_section[2,0] = 4304
 GC34: y[2,0,5] = 1 -> element_section[2,0] = 5
 GC35: y[2,0,5] = 1 -> cost_section[2,0] = 4944
 GC36: y[3,0,0] = 1 -> element_section[3,0] = 0
 GC37: y[3,0,0] = 1 -> cost_section[3,0] = 0
 GC38: y[3,0,1] = 1 -> element_section[3,0] = 1
 GC39: y[3,0,1] = 1 -> cost_section[3,0] = 1131
 GC40: y[3,0,2] = 1 -> element_section[3,0] = 2
 GC41: y[3,0,2] = 1 -> cost_section[3,0] = 1876
 GC42: y[3,0,3] = 1 -> element_section[3,0] = 3
 GC43: y[3,0,3] = 1 -> cost_section[3,0] = 4140
 GC44: y[3,0,4] = 1 -> element_section[3,0] = 4
 GC45: y[3,0,4] = 1 -> cost_section[3,0] = 4304
 GC46: y[3,0,5] = 1 -> element_section[3,0] = 5
 GC47: y[3,0,5] = 1 -> cost_section[3,0] = 4944
 GC48: y[4,0,0] = 1 -> element_section[4,0] = 0
 GC49: y[4,0,0] = 1 -> cost_section[4,0] = 0
 GC50: y[4,0,1] = 1 -> element_section[4,0] = 1
 GC51: y[4,0,1] = 1 -> cost_section[4,0] = 1131
 GC52: y[4,0,2] = 1 -> element_section[4,0] = 2
 GC53: y[4,0,2] = 1 -> cost_section[4,0] = 1876
 GC54: y[4,0,3] = 1 -> element_section[4,0] = 3
 GC55: y[4,0,3] = 1 -> cost_section[4,0] = 4140
 GC56: y[4,0,4] = 1 -> element_section[4,0] = 4
 GC57: y[4,0,4] = 1 -> cost_section[4,0] = 4304
 GC58: y[4,0,5] = 1 -> element_section[4,0] = 5
 GC59: y[4,0,5] = 1 -> cost_section[4,0] = 4944
 GC60: y[0,1,0] = 1 -> element_section[0,1] = 0
 GC61: y[0,1,0] = 1 -> cost_section[0,1] = 0
 GC62: y[0,1,1] = 1 -> element_section[0,1] = 1
 GC63: y[0,1,1] = 1 -> cost_section[0,1] = 1131
 GC64: y[0,1,2] = 1 -> element_section[0,1] = 2
 GC65: y[0,1,2] = 1 -> cost_section[0,1] = 1876
 GC66: y[0,1,3] = 1 -> element_section[0,1] = 3
 GC67: y[0,1,3] = 1 -> cost_section[0,1] = 4140
 GC68: y[0,1,4] = 1 -> element_section[0,1] = 4
 GC69: y[0,1,4] = 1 -> cost_section[0,1] = 4304
 GC70: y[0,1,5] = 1 -> element_section[0,1] = 5
 GC71: y[0,1,5] = 1 -> cost_section[0,1] = 4944
 GC72: y[1,1,0] = 1 -> element_section[1,1] = 0
 GC73: y[1,1,0] = 1 -> cost_section[1,1] = 0
 GC74: y[1,1,1] = 1 -> element_section[1,1] = 1
 GC75: y[1,1,1] = 1 -> cost_section[1,1] = 1131
 GC76: y[1,1,2] = 1 -> element_section[1,1] = 2
 GC77: y[1,1,2] = 1 -> cost_section[1,1] = 1876
 GC78: y[1,1,3] = 1 -> element_section[1,1] = 3
 GC79: y[1,1,3] = 1 -> cost_section[1,1] = 4140
 GC80: y[1,1,4] = 1 -> element_section[1,1] = 4
 GC81: y[1,1,4] = 1 -> cost_section[1,1] = 4304
 GC82: y[1,1,5] = 1 -> element_section[1,1] = 5
 GC83: y[1,1,5] = 1 -> cost_section[1,1] = 4944
 GC84: y[2,1,0] = 1 -> element_section[2,1] = 0
 GC85: y[2,1,0] = 1 -> cost_section[2,1] = 0
 GC86: y[2,1,1] = 1 -> element_section[2,1] = 1
 GC87: y[2,1,1] = 1 -> cost_section[2,1] = 1131
 GC88: y[2,1,2] = 1 -> element_section[2,1] = 2
 GC89: y[2,1,2] = 1 -> cost_section[2,1] = 1876
 GC90: y[2,1,3] = 1 -> element_section[2,1] = 3
 GC91: y[2,1,3] = 1 -> cost_section[2,1] = 4140
 GC92: y[2,1,4] = 1 -> element_section[2,1] = 4
 GC93: y[2,1,4] = 1 -> cost_section[2,1] = 4304
 GC94: y[2,1,5] = 1 -> element_section[2,1] = 5
 GC95: y[2,1,5] = 1 -> cost_section[2,1] = 4944
 GC96: y[3,1,0] = 1 -> element_section[3,1] = 0
 GC97: y[3,1,0] = 1 -> cost_section[3,1] = 0
 GC98: y[3,1,1] = 1 -> element_section[3,1] = 1
 GC99: y[3,1,1] = 1 -> cost_section[3,1] = 1131
 GC100: y[3,1,2] = 1 -> element_section[3,1] = 2
 GC101: y[3,1,2] = 1 -> cost_section[3,1] = 1876
 GC102: y[3,1,3] = 1 -> element_section[3,1] = 3
 GC103: y[3,1,3] = 1 -> cost_section[3,1] = 4140
 GC104: y[3,1,4] = 1 -> element_section[3,1] = 4
 GC105: y[3,1,4] = 1 -> cost_section[3,1] = 4304
 GC106: y[3,1,5] = 1 -> element_section[3,1] = 5
 GC107: y[3,1,5] = 1 -> cost_section[3,1] = 4944
 GC108: y[4,1,0] = 1 -> element_section[4,1] = 0
 GC109: y[4,1,0] = 1 -> cost_section[4,1] = 0
 GC110: y[4,1,1] = 1 -> element_section[4,1] = 1
 GC111: y[4,1,1] = 1 -> cost_section[4,1] = 1131
 GC112: y[4,1,2] = 1 -> element_section[4,1] = 2
 GC113: y[4,1,2] = 1 -> cost_section[4,1] = 1876
 GC114: y[4,1,3] = 1 -> element_section[4,1] = 3
 GC115: y[4,1,3] = 1 -> cost_section[4,1] = 4140
 GC116: y[4,1,4] = 1 -> element_section[4,1] = 4
 GC117: y[4,1,4] = 1 -> cost_section[4,1] = 4304
 GC118: y[4,1,5] = 1 -> element_section[4,1] = 5
 GC119: y[4,1,5] = 1 -> cost_section[4,1] = 4944
 GC120: y[0,2,0] = 1 -> element_section[0,2] = 0
 GC121: y[0,2,0] = 1 -> cost_section[0,2] = 0
 GC122: y[0,2,1] = 1 -> element_section[0,2] = 1
 GC123: y[0,2,1] = 1 -> cost_section[0,2] = 1131
 GC124: y[0,2,2] = 1 -> element_section[0,2] = 2
 GC125: y[0,2,2] = 1 -> cost_section[0,2] = 1876
 GC126: y[0,2,3] = 1 -> element_section[0,2] = 3
 GC127: y[0,2,3] = 1 -> cost_section[0,2] = 4140
 GC128: y[0,2,4] = 1 -> element_section[0,2] = 4
 GC129: y[0,2,4] = 1 -> cost_section[0,2] = 4304
 GC130: y[0,2,5] = 1 -> element_section[0,2] = 5
 GC131: y[0,2,5] = 1 -> cost_section[0,2] = 4944
 GC132: y[1,2,0] = 1 -> element_section[1,2] = 0
 GC133: y[1,2,0] = 1 -> cost_section[1,2] = 0
 GC134: y[1,2,1] = 1 -> element_section[1,2] = 1
 GC135: y[1,2,1] = 1 -> cost_section[1,2] = 1131
 GC136: y[1,2,2] = 1 -> element_section[1,2] = 2
 GC137: y[1,2,2] = 1 -> cost_section[1,2] = 1876
 GC138: y[1,2,3] = 1 -> element_section[1,2] = 3
 GC139: y[1,2,3] = 1 -> cost_section[1,2] = 4140
 GC140: y[1,2,4] = 1 -> element_section[1,2] = 4
 GC141: y[1,2,4] = 1 -> cost_section[1,2] = 4304
 GC142: y[1,2,5] = 1 -> element_section[1,2] = 5
 GC143: y[1,2,5] = 1 -> cost_section[1,2] = 4944
 GC144: y[2,2,0] = 1 -> element_section[2,2] = 0
 GC145: y[2,2,0] = 1 -> cost_section[2,2] = 0
 GC146: y[2,2,1] = 1 -> element_section[2,2] = 1
 GC147: y[2,2,1] = 1 -> cost_section[2,2] = 1131
 GC148: y[2,2,2] = 1 -> element_section[2,2] = 2
 GC149: y[2,2,2] = 1 -> cost_section[2,2] = 1876
 GC150: y[2,2,3] = 1 -> element_section[2,2] = 3
 GC151: y[2,2,3] = 1 -> cost_section[2,2] = 4140
 GC152: y[2,2,4] = 1 -> element_section[2,2] = 4
 GC153: y[2,2,4] = 1 -> cost_section[2,2] = 4304
 GC154: y[2,2,5] = 1 -> element_section[2,2] = 5
 GC155: y[2,2,5] = 1 -> cost_section[2,2] = 4944
 GC156: y[3,2,0] = 1 -> element_section[3,2] = 0
 GC157: y[3,2,0] = 1 -> cost_section[3,2] = 0
 GC158: y[3,2,1] = 1 -> element_section[3,2] = 1
 GC159: y[3,2,1] = 1 -> cost_section[3,2] = 1131
 GC160: y[3,2,2] = 1 -> element_section[3,2] = 2
 GC161: y[3,2,2] = 1 -> cost_section[3,2] = 1876
 GC162: y[3,2,3] = 1 -> element_section[3,2] = 3
 GC163: y[3,2,3] = 1 -> cost_section[3,2] = 4140
 GC164: y[3,2,4] = 1 -> element_section[3,2] = 4
 GC165: y[3,2,4] = 1 -> cost_section[3,2] = 4304
 GC166: y[3,2,5] = 1 -> element_section[3,2] = 5
 GC167: y[3,2,5] = 1 -> cost_section[3,2] = 4944
 GC168: y[4,2,0] = 1 -> element_section[4,2] = 0
 GC169: y[4,2,0] = 1 -> cost_section[4,2] = 0
 GC170: y[4,2,1] = 1 -> element_section[4,2] = 1
 GC171: y[4,2,1] = 1 -> cost_section[4,2] = 1131
 GC172: y[4,2,2] = 1 -> element_section[4,2] = 2
 GC173: y[4,2,2] = 1 -> cost_section[4,2] = 1876
 GC174: y[4,2,3] = 1 -> element_section[4,2] = 3
 GC175: y[4,2,3] = 1 -> cost_section[4,2] = 4140
 GC176: y[4,2,4] = 1 -> element_section[4,2] = 4
 GC177: y[4,2,4] = 1 -> cost_section[4,2] = 4304
 GC178: y[4,2,5] = 1 -> element_section[4,2] = 5
 GC179: y[4,2,5] = 1 -> cost_section[4,2] = 4944
 GC180: y[0,3,0] = 1 -> element_section[0,3] = 0
 GC181: y[0,3,0] = 1 -> cost_section[0,3] = 0
 GC182: y[0,3,1] = 1 -> element_section[0,3] = 1
 GC183: y[0,3,1] = 1 -> cost_section[0,3] = 1131
 GC184: y[0,3,2] = 1 -> element_section[0,3] = 2
 GC185: y[0,3,2] = 1 -> cost_section[0,3] = 1876
 GC186: y[0,3,3] = 1 -> element_section[0,3] = 3
 GC187: y[0,3,3] = 1 -> cost_section[0,3] = 4140
 GC188: y[0,3,4] = 1 -> element_section[0,3] = 4
 GC189: y[0,3,4] = 1 -> cost_section[0,3] = 4304
 GC190: y[0,3,5] = 1 -> element_section[0,3] = 5
 GC191: y[0,3,5] = 1 -> cost_section[0,3] = 4944
 GC192: y[1,3,0] = 1 -> element_section[1,3] = 0
 GC193: y[1,3,0] = 1 -> cost_section[1,3] = 0
 GC194: y[1,3,1] = 1 -> element_section[1,3] = 1
 GC195: y[1,3,1] = 1 -> cost_section[1,3] = 1131
 GC196: y[1,3,2] = 1 -> element_section[1,3] = 2
 GC197: y[1,3,2] = 1 -> cost_section[1,3] = 1876
 GC198: y[1,3,3] = 1 -> element_section[1,3] = 3
 GC199: y[1,3,3] = 1 -> cost_section[1,3] = 4140
 GC200: y[1,3,4] = 1 -> element_section[1,3] = 4
 GC201: y[1,3,4] = 1 -> cost_section[1,3] = 4304
 GC202: y[1,3,5] = 1 -> element_section[1,3] = 5
 GC203: y[1,3,5] = 1 -> cost_section[1,3] = 4944
 GC204: y[2,3,0] = 1 -> element_section[2,3] = 0
 GC205: y[2,3,0] = 1 -> cost_section[2,3] = 0
 GC206: y[2,3,1] = 1 -> element_section[2,3] = 1
 GC207: y[2,3,1] = 1 -> cost_section[2,3] = 1131
 GC208: y[2,3,2] = 1 -> element_section[2,3] = 2
 GC209: y[2,3,2] = 1 -> cost_section[2,3] = 1876
 GC210: y[2,3,3] = 1 -> element_section[2,3] = 3
 GC211: y[2,3,3] = 1 -> cost_section[2,3] = 4140
 GC212: y[2,3,4] = 1 -> element_section[2,3] = 4
 GC213: y[2,3,4] = 1 -> cost_section[2,3] = 4304
 GC214: y[2,3,5] = 1 -> element_section[2,3] = 5
 GC215: y[2,3,5] = 1 -> cost_section[2,3] = 4944
 GC216: y[3,3,0] = 1 -> element_section[3,3] = 0
 GC217: y[3,3,0] = 1 -> cost_section[3,3] = 0
 GC218: y[3,3,1] = 1 -> element_section[3,3] = 1
 GC219: y[3,3,1] = 1 -> cost_section[3,3] = 1131
 GC220: y[3,3,2] = 1 -> element_section[3,3] = 2
 GC221: y[3,3,2] = 1 -> cost_section[3,3] = 1876
 GC222: y[3,3,3] = 1 -> element_section[3,3] = 3
 GC223: y[3,3,3] = 1 -> cost_section[3,3] = 4140
 GC224: y[3,3,4] = 1 -> element_section[3,3] = 4
 GC225: y[3,3,4] = 1 -> cost_section[3,3] = 4304
 GC226: y[3,3,5] = 1 -> element_section[3,3] = 5
 GC227: y[3,3,5] = 1 -> cost_section[3,3] = 4944
 GC228: y[4,3,0] = 1 -> element_section[4,3] = 0
 GC229: y[4,3,0] = 1 -> cost_section[4,3] = 0
 GC230: y[4,3,1] = 1 -> element_section[4,3] = 1
 GC231: y[4,3,1] = 1 -> cost_section[4,3] = 1131
 GC232: y[4,3,2] = 1 -> element_section[4,3] = 2
 GC233: y[4,3,2] = 1 -> cost_section[4,3] = 1876
 GC234: y[4,3,3] = 1 -> element_section[4,3] = 3
 GC235: y[4,3,3] = 1 -> cost_section[4,3] = 4140
 GC236: y[4,3,4] = 1 -> element_section[4,3] = 4
 GC237: y[4,3,4] = 1 -> cost_section[4,3] = 4304
 GC238: y[4,3,5] = 1 -> element_section[4,3] = 5
 GC239: y[4,3,5] = 1 -> cost_section[4,3] = 4944
 GC240: y[0,4,0] = 1 -> element_section[0,4] = 0
 GC241: y[0,4,0] = 1 -> cost_section[0,4] = 0
 GC242: y[0,4,1] = 1 -> element_section[0,4] = 1
 GC243: y[0,4,1] = 1 -> cost_section[0,4] = 1131
 GC244: y[0,4,2] = 1 -> element_section[0,4] = 2
 GC245: y[0,4,2] = 1 -> cost_section[0,4] = 1876
 GC246: y[0,4,3] = 1 -> element_section[0,4] = 3
 GC247: y[0,4,3] = 1 -> cost_section[0,4] = 4140
 GC248: y[0,4,4] = 1 -> element_section[0,4] = 4
 GC249: y[0,4,4] = 1 -> cost_section[0,4] = 4304
 GC250: y[0,4,5] = 1 -> element_section[0,4] = 5
 GC251: y[0,4,5] = 1 -> cost_section[0,4] = 4944
 GC252: y[1,4,0] = 1 -> element_section[1,4] = 0
 GC253: y[1,4,0] = 1 -> cost_section[1,4] = 0
 GC254: y[1,4,1] = 1 -> element_section[1,4] = 1
 GC255: y[1,4,1] = 1 -> cost_section[1,4] = 1131
 GC256: y[1,4,2] = 1 -> element_section[1,4] = 2
 GC257: y[1,4,2] = 1 -> cost_section[1,4] = 1876
 GC258: y[1,4,3] = 1 -> element_section[1,4] = 3
 GC259: y[1,4,3] = 1 -> cost_section[1,4] = 4140
 GC260: y[1,4,4] = 1 -> element_section[1,4] = 4
 GC261: y[1,4,4] = 1 -> cost_section[1,4] = 4304
 GC262: y[1,4,5] = 1 -> element_section[1,4] = 5
 GC263: y[1,4,5] = 1 -> cost_section[1,4] = 4944
 GC264: y[2,4,0] = 1 -> element_section[2,4] = 0
 GC265: y[2,4,0] = 1 -> cost_section[2,4] = 0
 GC266: y[2,4,1] = 1 -> element_section[2,4] = 1
 GC267: y[2,4,1] = 1 -> cost_section[2,4] = 1131
 GC268: y[2,4,2] = 1 -> element_section[2,4] = 2
 GC269: y[2,4,2] = 1 -> cost_section[2,4] = 1876
 GC270: y[2,4,3] = 1 -> element_section[2,4] = 3
 GC271: y[2,4,3] = 1 -> cost_section[2,4] = 4140
 GC272: y[2,4,4] = 1 -> element_section[2,4] = 4
 GC273: y[2,4,4] = 1 -> cost_section[2,4] = 4304
 GC274: y[2,4,5] = 1 -> element_section[2,4] = 5
 GC275: y[2,4,5] = 1 -> cost_section[2,4] = 4944
 GC276: y[3,4,0] = 1 -> element_section[3,4] = 0
 GC277: y[3,4,0] = 1 -> cost_section[3,4] = 0
 GC278: y[3,4,1] = 1 -> element_section[3,4] = 1
 GC279: y[3,4,1] = 1 -> cost_section[3,4] = 1131
 GC280: y[3,4,2] = 1 -> element_section[3,4] = 2
 GC281: y[3,4,2] = 1 -> cost_section[3,4] = 1876
 GC282: y[3,4,3] = 1 -> element_section[3,4] = 3
 GC283: y[3,4,3] = 1 -> cost_section[3,4] = 4140
 GC284: y[3,4,4] = 1 -> element_section[3,4] = 4
 GC285: y[3,4,4] = 1 -> cost_section[3,4] = 4304
 GC286: y[3,4,5] = 1 -> element_section[3,4] = 5
 GC287: y[3,4,5] = 1 -> cost_section[3,4] = 4944
 GC288: y[4,4,0] = 1 -> element_section[4,4] = 0
 GC289: y[4,4,0] = 1 -> cost_section[4,4] = 0
 GC290: y[4,4,1] = 1 -> element_section[4,4] = 1
 GC291: y[4,4,1] = 1 -> cost_section[4,4] = 1131
 GC292: y[4,4,2] = 1 -> element_section[4,4] = 2
 GC293: y[4,4,2] = 1 -> cost_section[4,4] = 1876
 GC294: y[4,4,3] = 1 -> element_section[4,4] = 3
 GC295: y[4,4,3] = 1 -> cost_section[4,4] = 4140
 GC296: y[4,4,4] = 1 -> element_section[4,4] = 4
 GC297: y[4,4,4] = 1 -> cost_section[4,4] = 4304
 GC298: y[4,4,5] = 1 -> element_section[4,4] = 5
 GC299: y[4,4,5] = 1 -> cost_section[4,4] = 4944
Bounds
 group_section[0] <= 5
 group_section[1] <= 5
 group_section[2] <= 5
 group_section[3] <= 5
 group_section[4] <= 5
 element_section[0,0] <= 5
 element_section[1,0] <= 5
 element_section[2,0] <= 5
 element_section[3,0] <= 5
 element_section[4,0] <= 5
 element_section[0,1] <= 5
 element_section[1,1] <= 5
 element_section[2,1] <= 5
 element_section[3,1] <= 5
 element_section[4,1] <= 5
 element_section[0,2] <= 5
 element_section[1,2] <= 5
 element_section[2,2] <= 5
 element_section[3,2] <= 5
 element_section[4,2] <= 5
 element_section[0,3] <= 5
 element_section[1,3] <= 5
 element_section[2,3] <= 5
 element_section[3,3] <= 5
 element_section[4,3] <= 5
 element_section[0,4] <= 5
 element_section[1,4] <= 5
 element_section[2,4] <= 5
 element_section[3,4] <= 5
 element_section[4,4] <= 5
 group_lb[0] <= 5
 group_lb[1] <= 5
 group_lb[2] <= 5
 group_lb[3] <= 5
 group_lb[4] <= 5
 group_ub[0] <= 5
 group_ub[1] <= 5
 group_ub[2] <= 5
 group_ub[3] <= 5
 group_ub[4] <= 5
 group_range[0] <= 5
 group_range[1] <= 5
 group_range[2] <= 5
 group_range[3] <= 5
 group_range[4] <= 5
Binaries
 x[0,0,0] x[1,0,0] x[2,0,0] x[3,0,0] x[4,0,0] x[0,1,0] x[1,1,0] x[2,1,0]
 x[3,1,0] x[4,1,0] x[0,2,0] x[1,2,0] x[2,2,0] x[3,2,0] x[4,2,0] x[0,3,0]
 x[1,3,0] x[2,3,0] x[3,3,0] x[4,3,0] x[0,4,0] x[1,4,0] x[2,4,0] x[3,4,0]
 x[4,4,0] x[0,0,1] x[1,0,1] x[2,0,1] x[3,0,1] x[4,0,1] x[0,1,1] x[1,1,1]
 x[2,1,1] x[3,1,1] x[4,1,1] x[0,2,1] x[1,2,1] x[2,2,1] x[3,2,1] x[4,2,1]
 x[0,3,1] x[1,3,1] x[2,3,1] x[3,3,1] x[4,3,1] x[0,4,1] x[1,4,1] x[2,4,1]
 x[3,4,1] x[4,4,1] x[0,0,2] x[1,0,2] x[2,0,2] x[3,0,2] x[4,0,2] x[0,1,2]
 x[1,1,2] x[2,1,2] x[3,1,2] x[4,1,2] x[0,2,2] x[1,2,2] x[2,2,2] x[3,2,2]
 x[4,2,2] x[0,3,2] x[1,3,2] x[2,3,2] x[3,3,2] x[4,3,2] x[0,4,2] x[1,4,2]
 x[2,4,2] x[3,4,2] x[4,4,2] x[0,0,3] x[1,0,3] x[2,0,3] x[3,0,3] x[4,0,3]
 x[0,1,3] x[1,1,3] x[2,1,3] x[3,1,3] x[4,1,3] x[0,2,3] x[1,2,3] x[2,2,3]
 x[3,2,3] x[4,2,3] x[0,3,3] x[1,3,3] x[2,3,3] x[3,3,3] x[4,3,3] x[0,4,3]
 x[1,4,3] x[2,4,3] x[3,4,3] x[4,4,3] x[0,0,4] x[1,0,4] x[2,0,4] x[3,0,4]
 x[4,0,4] x[0,1,4] x[1,1,4] x[2,1,4] x[3,1,4] x[4,1,4] x[0,2,4] x[1,2,4]
 x[2,2,4] x[3,2,4] x[4,2,4] x[0,3,4] x[1,3,4] x[2,3,4] x[3,3,4] x[4,3,4]
 x[0,4,4] x[1,4,4] x[2,4,4] x[3,4,4] x[4,4,4] group_exists[0]
 group_exists[1] group_exists[2] group_exists[3] group_exists[4]
 col_in_group[0,0] col_in_group[0,1] col_in_group[0,2] col_in_group[0,3]
 col_in_group[0,4] col_in_group[1,0] col_in_group[1,1] col_in_group[1,2]
 col_in_group[1,3] col_in_group[1,4] col_in_group[2,0] col_in_group[2,1]
 col_in_group[2,2] col_in_group[2,3] col_in_group[2,4] col_in_group[3,0]
 col_in_group[3,1] col_in_group[3,2] col_in_group[3,3] col_in_group[3,4]
 col_in_group[4,0] col_in_group[4,1] col_in_group[4,2] col_in_group[4,3]
 col_in_group[4,4] level_in_group[0,0] level_in_group[0,1]
 level_in_group[0,2] level_in_group[0,3] level_in_group[0,4]
 level_in_group[1,0] level_in_group[1,1] level_in_group[1,2]
 level_in_group[1,3] level_in_group[1,4] level_in_group[2,0]
 level_in_group[2,1] level_in_group[2,2] level_in_group[2,3]
 level_in_group[2,4] level_in_group[3,0] level_in_group[3,1]
 level_in_group[3,2] level_in_group[3,3] level_in_group[3,4]
 level_in_group[4,0] level_in_group[4,1] level_in_group[4,2]
 level_in_group[4,3] level_in_group[4,4] y[0,0,0] y[0,0,1] y[0,0,2]
 y[0,0,3] y[0,0,4] y[0,0,5] y[1,0,0] y[1,0,1] y[1,0,2] y[1,0,3] y[1,0,4]
 y[1,0,5] y[2,0,0] y[2,0,1] y[2,0,2] y[2,0,3] y[2,0,4] y[2,0,5] y[3,0,0]
 y[3,0,1] y[3,0,2] y[3,0,3] y[3,0,4] y[3,0,5] y[4,0,0] y[4,0,1] y[4,0,2]
 y[4,0,3] y[4,0,4] y[4,0,5] y[0,1,0] y[0,1,1] y[0,1,2] y[0,1,3] y[0,1,4]
 y[0,1,5] y[1,1,0] y[1,1,1] y[1,1,2] y[1,1,3] y[1,1,4] y[1,1,5] y[2,1,0]
 y[2,1,1] y[2,1,2] y[2,1,3] y[2,1,4] y[2,1,5] y[3,1,0] y[3,1,1] y[3,1,2]
 y[3,1,3] y[3,1,4] y[3,1,5] y[4,1,0] y[4,1,1] y[4,1,2] y[4,1,3] y[4,1,4]
 y[4,1,5] y[0,2,0] y[0,2,1] y[0,2,2] y[0,2,3] y[0,2,4] y[0,2,5] y[1,2,0]
 y[1,2,1] y[1,2,2] y[1,2,3] y[1,2,4] y[1,2,5] y[2,2,0] y[2,2,1] y[2,2,2]
 y[2,2,3] y[2,2,4] y[2,2,5] y[3,2,0] y[3,2,1] y[3,2,2] y[3,2,3] y[3,2,4]
 y[3,2,5] y[4,2,0] y[4,2,1] y[4,2,2] y[4,2,3] y[4,2,4] y[4,2,5] y[0,3,0]
 y[0,3,1] y[0,3,2] y[0,3,3] y[0,3,4] y[0,3,5] y[1,3,0] y[1,3,1] y[1,3,2]
 y[1,3,3] y[1,3,4] y[1,3,5] y[2,3,0] y[2,3,1] y[2,3,2] y[2,3,3] y[2,3,4]
 y[2,3,5] y[3,3,0] y[3,3,1] y[3,3,2] y[3,3,3] y[3,3,4] y[3,3,5] y[4,3,0]
 y[4,3,1] y[4,3,2] y[4,3,3] y[4,3,4] y[4,3,5] y[0,4,0] y[0,4,1] y[0,4,2]
 y[0,4,3] y[0,4,4] y[0,4,5] y[1,4,0] y[1,4,1] y[1,4,2] y[1,4,3] y[1,4,4]
 y[1,4,5] y[2,4,0] y[2,4,1] y[2,4,2] y[2,4,3] y[2,4,4] y[2,4,5] y[3,4,0]
 y[3,4,1] y[3,4,2] y[3,4,3] y[3,4,4] y[3,4,5] y[4,4,0] y[4,4,1] y[4,4,2]
 y[4,4,3] y[4,4,4] y[4,4,5] Zu[0,0] Zu[0,1] Zu[0,2] Zu[0,3] Zu[0,4] Zu[1,0]
 Zu[1,1] Zu[1,2] Zu[1,3] Zu[1,4] Zu[2,0] Zu[2,1] Zu[2,2] Zu[2,3] Zu[2,4]
 Zu[3,0] Zu[3,1] Zu[3,2] Zu[3,3] Zu[3,4] Zu[4,0] Zu[4,1] Zu[4,2] Zu[4,3]
 Zu[4,4] Zl[0,0] Zl[0,1] Zl[0,2] Zl[0,3] Zl[0,4] Zl[1,0] Zl[1,1] Zl[1,2]
 Zl[1,3] Zl[1,4] Zl[2,0] Zl[2,1] Zl[2,2] Zl[2,3] Zl[2,4] Zl[3,0] Zl[3,1]
 Zl[3,2] Zl[3,3] Zl[3,4] Zl[4,0] Zl[4,1] Zl[4,2] Zl[4,3] Zl[4,4]
Generals
 group_section[0] group_section[1] group_section[2] group_section[3]
 group_section[4] element_section[0,0] element_section[1,0]
 element_section[2,0] element_section[3,0] element_section[4,0]
 element_section[0,1] element_section[1,1] element_section[2,1]
 element_section[3,1] element_section[4,1] element_section[0,2]
 element_section[1,2] element_section[2,2] element_section[3,2]
 element_section[4,2] element_section[0,3] element_section[1,3]
 element_section[2,3] element_section[3,3] element_section[4,3]
 element_section[0,4] element_section[1,4] element_section[2,4]
 element_section[3,4] element_section[4,4] group_lb[0] group_lb[1]
 group_lb[2] group_lb[3] group_lb[4] group_ub[0] group_ub[1] group_ub[2]
 group_ub[3] group_ub[4] group_range[0] group_range[1] group_range[2]
 group_range[3] group_range[4]
End
