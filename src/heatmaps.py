import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors

def test_heatmap(section_of_elements, grouped_elements, plot_sections = True):

  unique_values = np.unique(section_of_elements)

  if plot_sections:
    unique_values = unique_values[unique_values > 0]  # Exclude zero
  else:
    unique_values = unique_values[unique_values != -1]  # Exclude -1

  palette = sns.color_palette("viridis", len(unique_values))

  # add white color for zero (-1) values
  cmap = mcolors.ListedColormap([(1, 1, 1)] + palette)  # White for zero (-1), then mapped colors

  # mapping of values to color indices
  value_to_color_idx = {val: i+1 for i, val in enumerate(unique_values)}

  if plot_sections:
    value_to_color_idx[0] = 0
  else:
    value_to_color_idx[-1] = 0

  # array values to color indices
  color_indices = np.vectorize(value_to_color_idx.get)(section_of_elements)

  if plot_sections:
    sns.heatmap(np.rot90(section_of_elements), annot=np.rot90(grouped_elements), cmap=cmap,
              linewidths=0.5, linecolor="black", cbar=True)
  else:
    sns.heatmap(np.rot90(section_of_elements), annot=np.rot90(grouped_elements), cmap=cmap,
              linewidths=0.5, linecolor="black", cbar=False)

  if plot_sections:
    plt.title("Section of elements")
  else:
    plt.title("Grouping of elements")
  plt.xlabel("Column")
  plt.ylabel("Level")
  plt.yticks(ticks=np.arange(n_levels), labels=range(n_levels - 1, -1, -1), rotation=0)

  plt.show()

grouped_elements = np.full((n_cols, n_levels), -1)  # -1 as default (if element doesn't exist)

for g, i, j in x.keys():
    if x[g, i, j].X > 0.5:  # Check if x[g, i, j] is active
        grouped_elements[i, j] = g

print(grouped_elements)

section_of_elements = np.full((n_cols, n_levels), 0)

for i, j in element_section.keys():
    section_of_elements[i, j] = element_section[i, j].X

print(section_of_elements)

unique_values = np.unique(test_data_rotated)
unique_values = unique_values[unique_values > 0]

palette = sns.color_palette("viridis", len(unique_values))

cmap = mcolors.ListedColormap([(1, 1, 1)] + palette)

value_to_color_idx = {val: i+1 for i, val in enumerate(unique_values)}
value_to_color_idx[0] = 0

color_indices = np.vectorize(value_to_color_idx.get)(test_data_rotated)

sns.heatmap(test_data_rotated, annot=True, cmap=cmap,
            linewidths=0.5, linecolor="black", cbar=True)

plt.title("Section size of each Element")
plt.xlabel("Column")
plt.ylabel("Level")
plt.yticks(ticks=np.arange(n_levels), labels=range(n_levels - 1, -1, -1), rotation=0)


plt.show()
test_heatmap(grouped_elements, section_of_elements, False)

test_heatmap(section_of_elements, grouped_elements)