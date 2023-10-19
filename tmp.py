import numpy as np
import matplotlib.pyplot as plt

# Sample symmetric matrices
matrix1 = np.array([[1, 2, 3],
                    [2, 4, 5],
                    [3, 5, 6]])

matrix2 = np.array([[7, 8, 9],
                    [8, 10, 11],
                    [9, 11, 12]])

# Combine matrices
size = matrix1.shape[0]
combined_matrix = np.zeros((size, size))
combined_matrix[np.tril_indices(size)] = matrix1[np.tril_indices(size)]
combined_matrix[np.triu_indices(size, k=1)] = matrix2[np.triu_indices(size, k=1)]

# Create heatmap plot
fig, ax = plt.subplots()
heatmap = ax.imshow(combined_matrix, cmap='hot')

# Add colorbar
cbar = plt.colorbar(heatmap)

# Set x and y labels
ax.set_xticks(np.arange(size))
ax.set_yticks(np.arange(size))
ax.set_xticklabels(np.arange(1, size + 1))
ax.set_yticklabels(np.arange(1, size + 1))
ax.set_xlabel('Column')
ax.set_ylabel('Row')

# Set legend
cbar.ax.set_ylabel('Value')

# Set title
ax.set_title('Combined Matrix Heatmap')

# Display the plot
plt.show()
