# genetic-gerrymandering
## Simulating district gerrymandering using a genetic algorithm

This code simulates the practice of partisan gerrymandering using a genetic algorithm, a type of evolutionary algorithm. The purpose is to simulate and predict how an incumbent party will alter the electoral map to their advantage.

### How it works
  - State election shape files and election results come from Harvard Election Data Archive: https://projects.iq.harvard.edu/eda .
  - A partitition (redistricting) of the individual voting precincts within the state are representated as a graph/adjacency matrix, with each node (precinct) labeled with its parent district and Democrat/Republican vote distribution.
  - To initialize the population, districts are randomly partitioned and balanced by the number of voters. The random balanced partition operation is repeated on sub-graphs throughout the optimization as part of the crossover and mutation operations.
  - Crossover and mutation operations are defined on partitions:
    - **Crossover**: two existing partitions are searched for districts (an equal number from each) which are all disjoint and therefore compatible as members of a new partition. These districts are used as the basis of a new partition, with the remaining precincts randomly partitioned. See function `cx_dissolve_difference`.
    - **Mutation**: two adjacent precincts from different districts and with similar voting populations are swapped, slightly changing the dividing border between their districts. See function `mu_swap_geo`.
  - The optimization maximizes the sum of log ratios of incumbent party votes to opponent party votes of each district.

### How to run
1. Install and activate conda environment `conda env create --file geo.yml && conda activate geo`
2. Download a state file from https://projects.iq.harvard.edu/eda
3. Run the script, specifying your file and number of target districts: `python geo.py shp/pa_final.shp -n 18`
4. Can also specify optional parameters:
    - `-p` population size (default = 50)
    - `-l` lambda, the number of children to produce at each generation (default = 1/3 population)
    - `-g` number of generations (default = 100)
    - `--cxpb` probability that an offspring is produced by crossover (default = 0.8)
    - `--mutpb` probability that an offspring is produced by mutation (default = 0.2)


### Example output
An example of a random balanced partition, an individual from the initial population (no gerrymandering)

![PA random partition](output_g0.png?raw=true "PA random partition")

A partition from the 5th generation

![PA generation 10](output_g10.png?raw=true "PA generation 10")

A partition from the 100th generation

![PA generation 100](output_g100.png?raw=true "PA generation 100")

Optimization logging:
```
gen	nevals	avg       	std      	min       	max        
0  	50    	-0.0299536	0.0178801	-0.0633902	-0.00173615
1  	5     	-0.0273432	0.0158491	-0.0554915	-0.00173615
2  	5     	-0.0247061	0.0145179	-0.0483671	-0.00173615
3  	5     	-0.0223505	0.0133894	-0.0446486	-0.00173615
4  	5     	-0.0200328	0.0119169	-0.0432254	-0.00173615
5  	5     	-0.0189999	0.0108833	-0.0398585	-0.00173615
6  	5     	-0.0178442	0.00983655	-0.0348279	-0.00173615
...
```