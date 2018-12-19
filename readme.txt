The results from all growth simulations (1,051,596 simulations of two organisms with two carbon sources, with and without oxygen) are contained in the MATLAB structure 'MData.' Results are also provided in CSV format in 'MData.csv.zip.'

MData contains the following fields:
	
	- aer and anaer: Data structures containing the following results for the simulations with and without oxygen, respectively:

		- allGrewAlone: Binary vectors denoting whether or not both organisms grew in the first iteration of the algorithm.

		- allGrewCross: Binary vectors denoting whether or not both organisms grew in the last iteration of the algorithm (after cross-feeding).

		- expansions: Numerical vector detailing the number of medium expansions that occured in each simulation

		- secMetCountsAlone: Numerical vector containing the number of metabolites secreted at the first iteration of each simulation.

		- secMetCountsCross: Numerical vector containing the number of metabolites secreted at the last iteration of each simulation.

		- sharedMedium: Cell array detailing the minimal medium composition.

		- sub: Sparse binary double matrix detailing phenotypic information for the first (subject) organism in each simulation. Rows are simulations, column identities are contained in MData.col (Data given as row, column, value pairs in CSV version).

		- par: Sparse binary double matrix detailing phenotypic information for the second (partner) organism in each simulation. Rows are simulations, column identities are contained in MData.col (Data given as row, column, value pairs in CSV version).

		- Motifs: Data structure with two fields:

			- MotifsGen: Cell array detailing the general interaction type in each simulation (C: commensal, M: mutualistic, N: non-interacting).

			- MotifsSp: Cell array detailing the specific interaction motifs in each simulation (detailed in main text).

	- col: Data structure defining which columns in the 'sub' and 'par' matrices correspond to phenotypic information (Species: organisms, CSources: carbon sources, absMets: absorbed metabolites corresponding to MData.fullMetList, secMets: secreted metabolites corresponding to MData.fullMetList, O2: oxygen availability, growthAlone: whether the organism grew alone in minimal medium, growthCross: whether the organism grew after cross-feeding).

	- fullMetList: Cell array of all extracellular metabolites in all organism models.

	- speciesPairCombos: Cell array of organisms used in each of the simulations.

	- CSource pair combos: Cell array of carbon sources used in each of the simulations.
