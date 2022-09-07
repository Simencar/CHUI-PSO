# CHUI-PSO

This is a repository for CHUI-PSO, a Particle swarm optimization algorithm for mining closed high-utility itemsets

How to run:

* Set the "input" string in CHUI_PSO.java to the path of the dataset on file. The file must be in SPMF format, [here are some example datasets.](https://www.philippe-fournier-viger.com/spmf/index.php?link=datasets.php) 
* Set the "output" string to any .txt file path. The discovered patterns are written to this file during execution.
* Run main.java

The minimum utility threshold, population size, and number of iterations can be set in CHUI_PSO.java.

The algorithm can also discover high-utility itemsets by changing the variable "closed" to False in CHUI_PSO.java
