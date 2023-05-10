
# <a name="top"></a> Seeking Critical Nodes in Digraphs

* [The Repository](#repo)
* [BF](#bf)
* [CNH](#cnh)
* [Standard](#standard)
* [Iterative](#iterative)
* [Cite the Paper](#paper)


## <a name="repo"></a> The Repository ([top](#top))


The `src` directory contains 3 directories called: `bf`, `cnh` and `igraph-centralities`. Inside `bf` and `cnh` you will find the Makefile files used to compile the code of the brute force (BF) algorithm and the Critical Node Heuristic (CNH) respectivally. Inside `igraph-centralities` you will find the directories `iterative` and `standard` containing the Makefile files to compile the code of the random, standard and iterative strategies.

```
src/
 |-- bf/ 
     |-- include/
     |-- src/
     |-- Makefile
 |-- cnh
     |-- include/
     |-- src/
     |-- Makefile
 |-- igraph-centralities
     | -- common/ 
     | -- include/
     | -- iterative/
           |-- Makefile
     | -- standard/
           |-- Makefile
```

## <a name="bf"></a> BF ([top](#top))

You can compile it using the command `make all`, afterword the executable `digraphCriticalNodesBF` is created inside the directory bin and you can run it as follows:

```sh
./bin/digraphCriticalNodesBF --graph <FILE_NAME> --format <TYPE> --stop <NODES_NUMBER>
```
Below there is a short description of the parameters read by the program.

```
Usage: digraphCriticalNodesBF --graph <FILE_NAME> --format <TYPE> --stop <NODES_NUMBER>

	Starting from size 1, we remove incremental sets of critical nodes from the graph, we stop when:
	1) we reach a 0 connectivity value or 2) we reach the size defined by the <stop> parameter.
	The procedure computes the critical nodes of the graph based on the formula:

		f(G) = \sum_{i=1}^{l}inom{C_i}{2}

	Where C_i are the strongly connected components of G, and a node x is a critical node
	of G if it minimizes f(G\x), i.e x = arg min_{v \in V} f(G\v) 

	-s, --stop <NODES_NUMBER>
		Maximum size of the critical nodes set to remove, 0 is interpreted as all.
	-f, --graph <FILE_NAME>
		Read the graph from file <FILE_NAME>.
	-t, --format <TYPE>
		The format of the file containing the graph, admitted types are:
			- 0 (Edges List) The first line of the file is expected to contain the number of nodes
			    and edges of the graph in the following format: "# Nodes: <NUMBER> Edges: <NUMBER>"
```

## <a name="cnh"></a> CNH ([top](#top))

You can compile it using the command `make all`, afterword executables `cnh` and `cnh_large` are created inside the directory bin and you can run them as follows:

```sh
./bin/cnh -f <INPUT_FILE> -t <FILE_TYPE> -o <OUTPUT_FILE> [-r <ROOT>] [-k <NUMBER>] [-h]
```
Below there is a short description of the parameters read by the program.

```
Usage: cnh -f <INPUT_FILE> -t <FILE_TYPE> -o <OUTPUT_FILE> [-r <ROOT>] [-k <NUMBER>] [-h]

	-f <INPUT_FILE> 
	   Read the graph from file <INPUT_FILE>.
	-t <FILE_TYPE> 
	   The supported file types are 0 or 1; 0 for the DIMACS format and 1 for the RMAT format.
	-o <OUTPUT_FILE> 
	   Write the results in the file <OUTPUT_FILE>.
	-r <ROOT> 
	   Specify the source vertex
	-k <NUMBER>
	   Specify the number of vertex to remove. It should be 0 < NUMBER < #Vertices. The default value is 1. 
	-h prints this short help
```

## <a name="standard"></a> Standard ([top](#top))

You can compile the code using the command `make all`, afterword executables `igraph_random` and `igraph_standard` are created in the local directory. 

In order to compile the code you need the [C igraph library](https://igraph.org/c/), see the `-ligraph` directive and the variable `IGRAPH_HEADERS` in the Makefile.  `IGRAPH_HEADERS` is supposed to contain the path to the headers of the igraph library. Moreover you will probably need to set `LD_LIBRARY_PATH` and `LIBRARY_PATH` accordingly.


#### igraph_random

You can run `igraph_random` as follows:

```sh
./igraph_random --graph <FILE_NAME> --format <TYPE> --outfile <FILE_NAME> --percentage <PERCENTAGE>
```

Below there is a short description of the parameters read by the program.

```
Usage: igraph_random --graph <FILE_NAME> --format <TYPE> --outfile <FILE_NAME> --percentage <PERCENTAGE>

	-f, --graph <FILE_NAME>
		Read the graph from file <FILE_NAME>.
	-o, --outfile <FILE_NAME>
		The name of the output file.
	-p, --percentage <PERCENTAGE>.
		The percentage of nodes to remove.
	-t, --format <TYPE>
		The format of the file containing the graph, supported types are:
			- 0 (DIMACS)
			- 1 (Edges List) The first line of the file is expected to contain the number of nodes
			    and edges of the graph in the following format: "# Nodes: <NUMBER> Edges: <NUMBER>"
```


#### igraph_standard

You can run `igraph_standard ` as follows:

```sh
./igraph_standard --graph <FILE_NAME> --format <TYPE> --outfile <FILE_NAME> --metric <METRIC_TYPE> --percentage <PERCENTAGE>
```

Below there is a short description of the parameters read by the program.


```
Usage: igraph_standard --graph <FILE_NAME> --format <TYPE> --outfile <FILE_NAME> --metric <METRIC_TYPE> --percentage <PERCENTAGE>

	-f, --graph <FILE_NAME>
		Read the graph from file <FILE_NAME>.
	-m, --metric <METRIC_TYPS>
		Supported metrics are:
			- 1 Betwenness centrality
			- 2 Closeness (all)
			- 21 Closeness (in)
			- 22 Closeness (out)
			- 3 Degree (all)
			- 31 Degree (in)
			- 32 Degree (out)
			- 4 Pagerank
	-o, --outfile <FILE_NAME>
		The name of the output file.
	-p, --percentage <PERCENTAGE>.
		The percentage of nodes to remove.
	-t, --format <TYPE>
		The format of the file containing the graph, supported types are:
			- 0 (DIMACS)
			- 1 (Edges List) The first line of the file is expected to contain the number of nodes
			    and edges of the graph in the following format: "# Nodes: <NUMBER> Edges: <NUMBER>" 
```
## <a name="iterative"></a> Iterative ([top](#top))

You can compile the code using the command `make all`, afterword executable `igraph_iterative` is created in the local directory.

In order to compile the code you need the [C igraph library](https://igraph.org/c/), see the `-ligraph` directive and the variable `IGRAPH_HEADERS` in the Makefile.  `IGRAPH_HEADERS` is supposed to contain the path to the headers of the igraph library. Moreover you will probably need to set `LD_LIBRARY_PATH` and `LIBRARY_PATH` accordingly.


You can run `igraph_iterative` as follows:

```sh
./igraph_iterative --graph <FILE_NAME> --format <TYPE> --outfile <FILE_NAME> --metric <METRIC_TYPE> --percentage <PERCENTAGE>
```

Below there is a short description of the parameters read by the program.


```
Usage: igraph_iterative --graph <FILE_NAME> --format <TYPE> --outfile <FILE_NAME> --metric <METRIC_TYPE> --percentage <PERCENTAGE>

	-f, --graph <FILE_NAME>
		Read the graph from file <FILE_NAME>.
	-m, --metric <METRIC_TYPS>
		Supported metrics are:
			- 1 Betwenness centrality
			- 2 Closeness (all)
			- 21 Closeness (in)
			- 22 Closeness (out)
			- 3 Degree (all)
			- 31 Degree (in)
			- 32 Degree (out)
			- 4 Pagerank
	-o, --outfile <FILE_NAME>
		The output filename <FILE_NAME>.
	-p, --percentage <PERCENTAGE>.
		The percentage of nodes to output.
	-t, --format <TYPE>
		The format of the file containing the graph, admitted types are:
			- 0 (DIMACS)
			- 1 (Edges List) The first line of the file is expected to contain the number of nodes
			    and edges of the graph in the following format: "# Nodes: <NUMBER> Edges: <NUMBER>"  
```

## <a name="paper"></a> Cite the Paper ([top](#top))

[Seeking critical nodes in digraphs](https://www.sciencedirect.com/science/article/pii/S1877750323000728)

```
@article{BERNASCHI2023102012,
title = {Seeking critical nodes in digraphs},
journal = {Journal of Computational Science},
volume = {69},
pages = {102012},
year = {2023},
issn = {1877-7503},
doi = {https://doi.org/10.1016/j.jocs.2023.102012},
url = {https://www.sciencedirect.com/science/article/pii/S1877750323000728},
author = {Massimo Bernaschi and Alessandro Celestini and Marco Cianfriglia and Stefano Guarino and Giuseppe F. Italiano and Enrico Mastrostefano and Lena Rebecca Zastrow},
keywords = {Critical nodes, Networks connectivity, Centrality measures, Network analysis},
abstract = {The Critical Node Detection Problem (CNDP) consists in finding the set of nodes, defined critical, whose removal maximally degrades the graph. In this work we focus on finding the set of critical nodes whose removal minimizes the pairwise connectivity of a direct graph (digraph). Such problem has been proved to be NP-hard, thus we need efficient heuristics to detect critical nodes in real-world applications. We aim at understanding which is the best heuristic we can apply to identify critical nodes in practice, i.e., taking into account time constrains and real-world networks. We present an in-depth analysis of several heuristics we ran on both real-world and on synthetic graphs. We define and evaluate two different strategies for each heuristic: standard and iterative. Our main findings show that an algorithm recently proposed to solve the CNDP and that can be used as heuristic for the general case provides the best results in real-world graphs, and it is also the fastest. However, there are few exceptions that are thoroughly analyzed and discussed. We show that among the heuristics we analyzed, few of them cannot be applied to very large graphs, when the iterative strategy is used, due to their time complexity. Finally, we suggest possible directions to further improve the heuristic providing the best results.}
}
```






