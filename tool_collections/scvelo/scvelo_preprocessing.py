import argparse
import scvelo as scv
import anndata
 
# Initialize parser
parser = argparse.ArgumentParser()
 
# args filter_and_normalize
parser.add_argument("-o", "--Output", help = "Output")
parser.add_argument("-i", "--Input", help = "Input file")
parser.add_argument("--min", "--mincounts", help = "int (default: None). Minimum number of counts required for a gene to pass filtering (spliced). ", type=int)
parser.add_argument("--mu", "--mincountsu", help = "int (default: None). Minimum number of counts required for a gene to pass filtering (unspliced).", type=int)
parser.add_argument("--mc", "--mincells", help = "int (default: None). Minimum number of cells expressed required to pass filtering (spliced)", type=int)
parser.add_argument("--mcu", "--mincellsu", help = "int (default: None). Minimum number of cells expressed required to pass filtering (unspliced)", type=int)
parser.add_argument("--msc", "--minsharedcounts", help = "int, optional (default: None). Minimum number of counts (both unspliced and spliced) required for a gene.", type=int)
parser.add_argument("--msce", "--minsharedcells", help = "int, optional (default: None). Minimum number of cells required to be expressed (both unspliced and spliced)", type=int)
parser.add_argument("-n", "--ntop", help = "int (default: None). Number of genes to keep", type=int)
parser.add_argument("-r", "--retain", help = "list, optional (default: None). List of gene names to be retained independent of thresholds.")
parser.add_argument("-s", "--subset", help = " bool (default: True). Whether to subset highly variable genes or to store in .var[highly_variable].", default=True, type=bool)
parser.add_argument("-f", "--flavor", help = "‘seurat’, ‘cell_ranger’, ‘svr’}, optional (default: ‘seurat’).  Choose the flavor for computing normalized dispersion. If choosing ‘seurat’, this expects non-logarithmized data ", default="seurat")
parser.add_argument("-l", "--log", help = "bool (default: True). Take logarithm.", default=True, type=bool)
parser.add_argument("-y", "--layernormalize", help = " list of str (default: None). List of layers to be normalized. If set to None, the layers {‘X’, ‘spliced’, ‘unspliced’} are considered for normalization upon testing whether they have already been normalized (by checking type of entries: int -> unprocessed, float -> processed).")
parser.add_argument("--cpca", "--cpercellafter", help = "float or None, optional (default: None). If None, after normalization, each cell has a total count equal to the median of the counts_per_cell before normalization", type=float)
parser.add_argument("--cpc", "--cpercell", help = "np.array, optional (default: None). Precomputed counts per cell")
parser.add_argument("--knc", "--keyncount", help = "str, optional (default: ‘n_counts’). Name of the field in adata.obs where the total counts per cell are stored", default="n_counts")
parser.add_argument("--mppc", "--maxproportionpercell", help = "int (default: None). Exclude genes counts that account for more than a specific proportion of cell size, e.g. 0.05", type=float) # example is float though ?
parser.add_argument("--uis", "--useinitialsize", help = "bool (default: True). Whether to use initial cell sizes oder actual cell sizes", default=True, type=bool)
#parser.add_argument("--ly", "--layer", help = "str or list (default: [‘spliced’, ‘unspliced’]). Keys for layers to be also considered for normalization", default="spliced,unspliced")

# args pp moments

parser.add_argument("--nn", "--nneighbors", help = "int (default: 30). Number of neighbors to use.", default=30, type=int)
parser.add_argument("--np", "--npcs", help = "int (default: None). Number of principal components to use. If not specified, the full space is used of a pre-computed PCA, or 30 components are used when PCA is computed internally",  type=int)
parser.add_argument("--md", "--mode", help = " ‘connectivities’ or ‘distances’ (default: ‘connectivities’). Distance metric to use for moment computation.", default="connectivities")
parser.add_argument("--mt", "--method", help = "{{‘umap’, ‘hnsw’, ‘sklearn’, None}} (default: ‘umap’). Method to compute neighbors, only differs in runtime. Connectivities are computed with adaptive kernel width as proposed in Haghverdi et al. 2016 (https://doi.org/10.1038/nmeth.3971).", default="umap")
parser.add_argument("--ur", "--userep", help = "None, ‘X’ or any key for .obsm (default: None). Use the indicated representation. If None, the representation is chosen automatically: for .n_vars < 50, .X is used, otherwise ‘X_pca’ is used.")
parser.add_argument("--uhv", "--usehighlyvariable", help = "bool (default: True). Whether to use highly variable genes only, stored in .var[‘highly_variable’].", default=True, type=bool)




 
# Read arguments from command line
args = parser.parse_args()
 
if args.Input:
    print("Displaying Input as: % s" % args.Input)
    adata = scv.read(args.Input, cache=True)
else :
    print("Input file missing")


#layerslist= args.ly.split(',')
if args.layernormalize:
    layer_norm= args.layernormalize.split(',')
else:
    layer_norm=args.layernormalize

#print(layerslist)
#scv.set_figure_params()

#adata = scv.read(args.Input, cache=True)


##compute the first- and second-order moments (means and uncentered variances) for velocity estimation

scv.pp.filter_and_normalize(adata, min_counts=args.min, min_counts_u=args.mu, min_cells=args.mc, min_cells_u=args.mcu, min_shared_counts=args.msc, min_shared_cells=args.msce, n_top_genes=args.ntop, retain_genes=args.retain, subset_highly_variable=args.subset, flavor=args.flavor, log=args.log, layers_normalize=layer_norm, counts_per_cell_after=args.cpca, counts_per_cell=args.cpc, key_n_counts=args.knc, max_proportion_per_cell=args.mppc, use_initial_size=args.uis)

scv.pp.moments(adata, n_neighbors=args.nn, n_pcs=args.np, mode=args.md, method=args.mt, use_rep=args.ur, use_highly_variable=args.uhv)


adata.write_loom(args.Output)

# estimation of velocities
#scv.tl.velocity(adata, mode='stochastic')

#scv.tl.velocity_graph(adata)