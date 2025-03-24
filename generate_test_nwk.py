import random

def generate_newick(n_clusters=300, snp_distance=10):
    newick_str = ""
    clusters = []
    
    for i in range(n_clusters):
        sample1 = f"S{i*2}"
        sample2 = f"S{i*2+1}"
        cluster = f"({sample1}:{snp_distance},{sample2}:{snp_distance})"
        clusters.append(cluster)
    
    while len(clusters) > 1:
        new_clusters = []
        for j in range(0, len(clusters), 2):
            if j + 1 < len(clusters):
                branch_length = random.randint(20, 100)  # Random deeper divergence
                new_cluster = f"({clusters[j]}:{branch_length},{clusters[j+1]}:{branch_length})"
            else:
                new_cluster = clusters[j]  # Keep the last element if odd number
            new_clusters.append(new_cluster)
        clusters = new_clusters
    
    newick_str = clusters[0] + ";"
    return newick_str

newick_tree = generate_newick()
print(newick_tree)
