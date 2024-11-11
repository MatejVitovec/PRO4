import numpy as np

# Funkce pro načtení bodů ze souboru
def load_points(filename):
    points = []
    with open(filename, 'r') as file:
        for line in file:
            x, y = map(float, line.split())
            points.append((x, y))
    return np.array(points)

# Funkce pro normalizaci bodů
def normalize_points(points):
    min_vals = np.min(points, axis=0)
    max_vals = np.max(points, axis=0)
    return (points - min_vals) / (max_vals - min_vals), min_vals, max_vals

# Funkce pro denormalizaci bodů
def denormalize_points(points, min_vals, max_vals):
    return points * (max_vals - min_vals) + min_vals

# Funkce pro nalezení nejbližšího souseda
def find_nearest_neighbor(point, points):
    distances = np.linalg.norm(points - point, axis=1)
    nearest_index = np.argmin(distances)
    return nearest_index

# Funkce pro seřazení bodů do uzavřené hladké křivky
def sort_points(points):
    sorted_points = [points[0]]
    remaining_points = points[1:]
    
    while len(remaining_points) > 0:
        nearest_index = find_nearest_neighbor(sorted_points[-1], remaining_points)
        sorted_points.append(remaining_points[nearest_index])
        remaining_points = np.delete(remaining_points, nearest_index, axis=0)
    
    return np.array(sorted_points)

# Funkce pro uložení seřazených bodů do souboru
def save_sorted_points(filename, points):
    with open(filename, 'w') as file:
        for point in points:
            file.write(f"{point[0]} {point[1]}\n")

# Hlavní část skriptu
def main():
    input_filename = 'body.txt'  # Změňte na název vašeho vstupního souboru
    output_filename = 'serazene_body.txt'  # Název výstupního souboru

    points = load_points(input_filename)
    
    # Normalizace bodů
    normalized_points, min_vals, max_vals = normalize_points(points)
    
    # Seřazení bodů
    sorted_normalized_points = sort_points(normalized_points)
    
    # Denormalizace bodů
    sorted_points = denormalize_points(sorted_normalized_points, min_vals, max_vals)
    
    # Uložení seřazených bodů
    save_sorted_points(output_filename, sorted_points)

if __name__ == '__main__':
    main()
