# Notes on producing results using MicroCOSM

This document provides details on how to replicate the results presented in our submission to *PLOS Medicine* (‘Money, jobs or schooling? A model-based evaluation of economic strengthening in South Africa and its impact on HIV, sexually transmitted infections and teenage births’). It is not intended as a user guide, but in Appendix B we have provided some general guidance on the model code. The code was written in C++, and was compiled and run using the gcc compiler.

The default code that we have shared is the code used to run the model up to 2025, with a sample of 100 parameter combinations drawn from the main posterior distribution, not allowing for any economic strengthening interventions (i.e. the baseline scenario). The lists that follow provide more detail on the specific code changes need to replicate the published model outputs. For convenience, we have also included two spreadsheet files: Likelihood.xlsx (which shows how the likelihood values were calculated) and Posteriors.xlsx (which shows how the posterior sample was drawn, based on the likelihood values).

## To do the calibration to the data from randomized controlled trials (RCTs):

- Set `ProjectionTerm = 27`.

- Set `BaselineStart = 2005`.

- Set `FixedUncertainty = 0`.

- Set `StructuralRCTcalib = 1`.

- In the OneSimulation function, make sure the command `if (StructuralRCTcalib == 1){ RunStructRCTs(); }` is active (not commented out).

- Set `ParamCombs` and `samplesize` to 5000 (or whatever sample size you want to use in the initial calibration step).

- Run the model. The command line to compile in gcc is `g++ -std=c++14 Microsimulation.cpp StatFunctions.cpp mersenne.cpp -o Microsimulation.exe`. After the model has compiled, the command to run the model is just `./Microsimulation`.

- After you’ve finished producing results, copy the contents of the CashTransferOut, SchoolSupportOut and VocationTrainingOut files into cell A6 of the corresponding sheets of the **Likelihood.xlsx** workbook. This calculates the likelihood values (relative to the maximum) for all parameter combinations. 

- Repeat the previous two steps after changing the random seed (in the OneYear function, make sure that the lines `int seed = 100;` and `rg.RandomInit(seed);` are active), but paste the new results to cells BC6 in each sheet (not A6). 

- Copy column G-J of the ‘Summary’ sheet in **Likelihood.xlsx** into columns AS-AV in the first sheet of the **Posteriors.xlsx**.

- The final sheet in this **Posteriors.xlsx** workbook summarizes the posterior distributions (after transforming the U(0,1) values to the relevant scale). Currently the posteriors are calculated for all the data combined, but if you wanted to calculate the posteriors specifically for the cash transfer data (say), the easiest way to do this would be to change the formulas in the first column of the RandomUniformStruct sheet to refer to column AX not BA.

## To generate the validation results, where we compare the model associations between structural drivers:

- Set `FixedUncertainty = 1`.

- Set `StructuralDriverCalib = 1` (in addition to setting `StructuralRCTcalib = 1`). 

- Set `ProjectionTerm = 41`. (You only need 40 for the purpose of validation against the survey SES associations, but if you’re also wanting to show validation against other data sources, it’s better to go slightly longer.)

- Comment out the call to `RunStructRCTs`.

- Set `ParamCombs` and `samplesize` to 100.

- Copy and paste the contents of the ‘RandomUniformStruct’ sheet in the **Posteriors.xlsx** workbook into the RandomUniformStruct.txt file.

- Run the model (see the instructions above).

- Paste the contents of the StructuralOutputs.txt file into the StructuralOutput sheet in the Validation.xlsx workbook. Similarly for GenOutputs.

## To calculate the PAFs for low socio-economic status (i.e. the proportions of incident HIV and STI cases and new teenage pregnancies attributable to low income and poor education, over the 2000-2020 period):

- Set `FixedUncertainty = 1`.

- Set `StructuralRCTcalib = 1`.

- Set `ProjectionTerm = 36`.

- Set `BaselineStart = 2000` (we only need this for the purpose of keeping the teenage fertility rates constant at the 2000 levels after 2000).

- Comment out the call to `RunStructRCTs` in the `OneSimulation` function.

- Keep `ParamCombs` and `samplesize` at 100.

- Copy and paste the contents of the ‘RandomUniformStruct’ sheet in the **Posteriors.xlsx** workbook into the RandomUniformStruct.txt file.

- See Appendix A on code changes in the counterfactual scenario. 

- Run the model.

- Copy the NewHIVexp, NewHIVexpM and NewHIVexpF output files into rows 278-327 of the relevant sheets in the PAFs workbook. (Note that the row references differ depending on the iteration.)

- Similarly copy the NewCTNGTV, NewCTNGTV_M and NewCTNGTV_F output files into rows 278-327 of the PAFsSTIs.xlsx workbook.

- And copy the TeenBirths into the PAFsTeenPreg workbook (there’s no sex stratification). 

- Undo the changes listed in Appendix A and repeat the above process, but this time save the outputs in rows 8-57 of each sheet (the baseline scenario).

- Repeat all the above steps with different seeds, pasting the results into different ranges in the relevant PAFs workbooks. For the second iteration, in the `OneYear` function, make sure that the lines `int seed = 100;` and `rg.RandomInit(seed);` are active. Then for the third iteration, set the seed to 200, for the fourth iteration set the seed to 300, and so on.

## To generate the results for the cash transfer baseline scenario (note that these are needed for the intervention impact calculations):

- Set `FixedUncertainty = 1`.

- Set `StructuralRCTcalib = 1`.

- Set `ProjectionTerm = 56`.

- Set `BaselineStart = 2025` (we only need this for the purpose of keeping the ‘base’ fertility rates in teenagers constant).

- Comment out the call to `RunStructRCTs` in the `OneSimulation` function.

- Keep `ParamCombs` and `samplesize` at 100.

- Copy and paste the contents of the ‘RandomUniform cash’ sheet in the **Posteriors.xlsx** workbook into the RandomUniformStruct.txt file.

- Run the model.

- Copy the NewHIVexp, NewHIVexpM and NewHIVexpF output files into the ‘NewHIVexp baseline’ sheet in the ‘Cash transfers’ workbook.

- Similarly copy the NewCTNGTV, NewCTNGTV_M and NewCTNGTV_F output files into the ‘NewCTNGTV baseline’ sheet in the ‘Cash transfers’ workbook.

- Copy the secondary behavioural outputs into the ‘Input’ sheet of the BehavOutput workbook. Then copy the summary outputs from the ‘Outputs’ sheet of the same workbook into the ‘Behav baseline’ sheet of ‘Cash transfers’.

- Copy GenOutputs.txt contents into the ‘GenOutput baseline’ sheet in the ‘Cash transfers’ workbook. 

- Repeat all the above steps with different seeds, pasting the results into different ranges in the ‘Cash transfers’ workbook. For the second iteration, in the `OneYear` function, make sure that the lines `int seed = 100;” and “rg.RandomInit(seed);` are active. Then for the third iteration, set the seed to 200, for the fourth iteration set the seed to 300, and so on.

Note that we go through similar steps for the <ins>school support</ins> and <ins>vocational training</ins> scenarios (as inputs, copy and paste the contents of the ‘RandomUniform XX’ sheet in **Posteriors.xlsx**).

## To calculate the effect of a cash transfer intervention, starting in 2025 (note that these are essentially the same steps as when running the baseline scenario for cash transfers, except where highlighted):

- Keep `FixedUncertainty = 1`.

- Set `StructuralRCTcalib = 1`.

- Set `ProjectionTerm = 56`.

- Comment out the call to `RunStructRCTs` in the `OneSimulation` function.

- Set <mark>`StructIntScenario = 3`</mark> and `BaselineStart = 2025`.

- Keep `ParamCombs` and `samplesize` at 100.

- Copy and paste the contents of the ‘RandomUniform cash’ sheet in the Posteriors.xlsx workbook into the RandomUniformStruct.txt file.

- Run the model.

- Copy the NewHIVexp, NewHIVexpM and NewHIVexpF output files into the <mark>‘NewHIVexp cash’ sheet</mark> in the ‘Cash transfers’ workbook. Similarly with the <mark>‘GenOutput cash’ sheet</mark>.

- Similarly copy the NewCTNGTV, NewCTNGTV_M and NewCTNGTV_F output files into the <mark>‘NewCTNGTV cash’ sheet</mark> in the ‘Cash transfers’ workbook.

- Similarly copy the TeenBirths output files into the <mark>‘TeenBirths cash’ sheet</mark> in the ‘Cash transfers’ workbook.

- Copy the secondary behavioural outputs into the ‘Input’ sheet of the BehavOutput workbook. Then copy the summary outputs from the ‘Outputs’ sheet of the same workbook into the <mark>‘Behav cash’ sheet</mark> of ‘Cash transfers’.

- Repeat all the above steps with different seeds, pasting the results into different ranges in the ‘Cash transfers’ workbook. For the second iteration, in the`OneYear` function, make sure that the lines `int seed = 100;` and `rg.RandomInit(seed);` are active. Then for the third iteration, set the seed to 200, for the fourth iteration set the seed to 300, and so on.

We go through similar steps in simulating the impact of school support and vocational training interventions (but set StructIntScenario to 4 and 5 respectively).










## Appendix A: Counterfactual code changes in the ‘high SES’ scenario

There are various ways to define the counterfactual in which everyone has high education attainment. Firstly we consider the changes in sexual risk behaviour associated with educational attainment, employment and household income:

- At the start of the `OneYear` function, if `CurrYear=2000`, set `RRcasualLogDropIncomeF = 1`, `RRdebutLogDropIncomeF = 1` and `RR_STemployedM = 1`. 

- Change `TempAdj = pow(ORcondomPerYrSchool, 5.0)` in the `UpdateCondomEffects` function, if `CurrYear` >= 2000.

- In the `SimulateSexActs` function, in the part where we apply `RR_FSWcontactEmployedM`, modify the condition so that it’s `if(Employed==1 || CurrYear >=2000){ … }`. Make similar changes in the `UpdateFSWprofile` function. (But note that you should **NOT** change the code in the `AssignBehav` function.)

- In the `UpdateHetCasual` function, where the `RRcasualEmployedM` adjustment is applied, modify the condition from `(SexInd==0 && Employed==1)` to `(SexInd==0 && (Employed==1 || CurrYear>=2000))`. 

Secondly we include effects of SES on healthcare seeking:

- Where the `RR_VMMClogIncome` adjustment is applied in the `UpdateMaleCirc` function, in the line prior to `if(r2[ic] < circprob * …` add the line `if(Register[ic].CumSelectionProb > 0.0 && Register[ic].CumSelectionProb < 1.0 && CurrYear >= 2000){ Register[ic].CumSelectionProb = 1.0; }`.

- In the `UpdateContraception` function, where `ORinjectableEdu` is used to calculate `CumOR`, add the line `if(CurrYear >=2000){ CumOR = pow(ORinjectableEdu, 3.0); }`. 

- Ditto in the `TestStartContr` function.

- In the `TestStartContr` function, where `CumOR` is first calculated using `EduEffectOR`, add the line `if(CurrYear >=2000){ CumOR = pow(EduEffectContr, 3.0); }`.

- In the `CalcHIVtesting` function, replace the line `VCTprob *= pow(VCTadjEdu, ih – 11);` with `if(CurrYear < 2000){ VCTprob *= pow(VCTadjEdu, ih – 11); }` and follow it with `else{ VCTprob *= pow(VCTadjEdu, 2.0);}`.

And finally we include the effects of current schooling on sexual debut:

- In the `BalanceSexualPartners` function, where the `RRdebutInSchool[ig]` adjustment is applied, modify the code to `if(InSchool == 1 || (CurrAge <= 21 && CurrYear >=2000))`.

Note that we are not changing rates of marriage (since we do not consider marriage to be a sexual risk behaviour). There is nevertheless an implicit change in marriage rates with the changes to sexual debut (because virgins are assumed not to enter into marriage).






## Appendix B: General notes on running the MicroCOSM model

The most important files in the MicroCOSM project are:

- **Microsimulation.h** is the main header file. This lists the main variables, classes, objects created from those functions and functions.

- **Microsimulation.cpp** is the main source code file. This code contains details on how the various functions are defined.

- **RatesByYear.txt** is an input file in which we specify all time-varying input parameters. For each row of numbers we specify parameter values for every year from 1985 to 2070 (86 years).

- **SexAssumps.txt** is an input file in which we specify all assumptions about sexual behaviour.

- **STDepidemiology.txt** is an input file in which we specify all assumptions about sexually transmitted infections (STIs or STDs), including transmission probabilities, average time spent in different disease stages, treatment probabilities, etc.

- **RandomUniformXX.txt**: There are separate input files for each group of parameters (the ‘XX’). These files store the parameters that give the best fits to the corresponding calibration data; for example, RandomUniformBV contains the parameters that give the best fit to the bacterial vaginosis (BV) prevalence data. In each file there are different rows for each parameter combination and different columns for each parameter. The parameters are all stored as uniform (0, 1) variates, and have to be transformed using the appropriate prior distributions before they can be used in the model.

- **StatFunctions.h** is a header file for special statistical functions, and StatFunctions.cpp is the corresponding source code file.

- **randomc.h** is the header file for the random number generator (the ‘Mersenne twister’), and Mersenne.cpp is the corresponding source code file.



Some of the key variables in the model:

- **`Register[]`** is a vector of all individuals in the population. Each individual in the population is represented by a different element in this vector. The ID of each individual is 1 + their position in the vector. For example, `Register[0]` corresponds to the person with ID 1. It follows that if we refer to someone with ID 0, they don’t exist (for example, if a person’s primary sexual partner has ID 0, it means they don’t have a primary sexual partner). People’s characteristics can easily be referenced using `Register[];` for example, `Register[20].Imprisoned` refers to the incarceration status of the individual with ID 21. Note that although the Register is dynamically expanded each time a birth occurs, it does not contract when people die; we keep people in the Register even after they’ve died. To check if someone is alive, refer to `Register[x].AliveInd`, which is 0 if they’re dead and 1 if they’re alive.

- **`HHregister[]`** is a vector of all households in the population. As with the Register vector, the ID of each household is 1 + its position in the vector. Each individual has a HouseholdID, which represents the ID of the household they belong to; if `HouseholdID = 0`, it means they’re homeless. As with the Register vector, the `HHregister` is dynamically expanded each time a new household is formed, but households remain in the `HHregister` even after they cease to exist. To check if a household still exists refer to `HHregister[x].Active`, which is 0 if the household has no members, and 1 if the household has living members. For simplicity, we assume a person can only be a member of one household (meaning they live there more than half their time there).

- **`RSApop`** is somewhat redundant. You can think of it as representing the South African population. 

- `**CurrYear**` is a counter, which gets updated at the start of each projection year. It starts at 1985 (the first year of the projection). Note that projection year runs from mid-year to mid-year, so when we refer to 1985, for example, we are referring to the period from mid-1985 to mid-1986. `CurrYear` only changes to 1986 in mid-1986.

- **`ProjectionTerm`** is the number of years for which you want to run the model. For example, if you specify `ProjectionTerm = 20`, the model will run from mid-1985 to mid-2005 (and outputs will be produced up to the 2004 projection year).

- **`CycleS`** determines the annual frequency at which we update sexual behaviour characteristics (e.g. finding new sexual partnerships, ending relationships). The default of 48 means that we update sexual behaviour characteristics at approximately weekly time steps (we use a multiple of 12, instead of 52, to simplify other updates that occur at monthly time steps).

- **`CycleD`** determines the annual frequency at which we update HIV and STI transmission. The default of 48 works well for HIV and HSV-2, but be aware that higher frequencies might be appropriate for the STIs that are more short-term (e.g. gonorrhoea) – consider possibly using values of 96 or 192.

- Arrays such as **`r[]`**, **`r2[]`**, **`revent[]`**, etc., are used to store random numbers (uniform (0, 1) variates).

- **`ParamCombs`** represents the number of parameter combinations you want to run. For example, if your set `ParamCombs = 5` it means you want to run the model with 5 different parameter combinations.

- **`IterationsPerPC`** represents the number of times you want to run the model per parameter combination (sometimes we want to run the model multiple times with the same parameter combination, then take the average across simulations, in order to reduce stochastic variation in the outputs). But in almost all cases we set `IterationsPerPC` to 1.

- **`samplesize`** is the product of `ParameCombs` and `IterationsPerPC`; it represents the total number of times you want to run the model.

- **`CurrSim`** is a counter of how far we have got in running through samplesize. If `CurrSim = 1`, we are busy with the first simulation.

- **`FixedUncertainty`** is either 0 or 1. We would typically set it to 0 if we’re running the model many 1000s of times, using different randomly generated parameter combination, and just want to calculate likelihood statistics for each of the randomly selected parameter combinations (from which we would subsequently select the best-fitting parameter combinations). This is often referred to as “running the model in calibration mode”. We would set `FixedUncertainty` to 1 if we want to produce detailed outputs for the best-fitting parameter combinations. For model testing purposes it is generally easier to set `FixedUncertainty = 1`.

## Running the model in calibration mode

Running the model in calibration mode means that we want to run the model many thousands of times, using different randomly generated parameter combinations, to identify the model results that give the best fit to the data. As a first step, you have to set `FixedUncertainty = 0`. Secondly, consider which datasets you want to include in defining the likelihood function. For example, if you want to fit the model to HIV prevalence data, set `HIVcalib = 1`, or if you want to fit the model to bacterial vaginosis (BV) prevalence data, set `BVcalib = 1`. Then set `ParamCombs`, `IterationsPerPC` and `samplesize` to appropriate values (consider the time it takes to run one simulation and the available computational resources). Remember that you may need to adjust ProjectionYear depending on what data you’re calibrating to. For example, if you’re calibrating your model to gonorrhoea prevalence data, but you only have gonorrhoea prevalence data up to 2012, you could set `ProjectionTerm = 28`.

## Variables that should not be changed

Although you may be tempted to change certain variables, there are some variables which you should NOT change (or the model will produce garbage results):

- `StartYear = 1985`

- `InitPop` and `Max-` variables: Although the default of 20 000 can be changed, if you increase it much above the default, you may run into problems with vector sizes exceeding the maximum allowed. The maximum vector sizes (`MaxPop`, `MaxCSWs`, `MaxCasual`) can be increased if you need to get around this problem, but generally I would not recommend this. Because the run times increase exponentially as the population size increases, it’s quicker to simulate 2 populations of 20 000 individuals (and aggregate/average the results) than it is to simulate a single population of 40 000 individuals.

- `ConstantInitialHIV = 1`

- `SetInitPrev1990 = 1`

- `HIVind = 1`

In general, if you’re unsure about what a variable does, you should keep it at its default value. Don’t change variables unless you have a good reason to!

## General format of the output files

In almost all output files, the number of rows corresponds to the specified samplesize variable (i.e. each row corresponds to a different parameter combination). In most output files, the first two columns are parameter combination IDs (or blanks) and the subsequent columns store the results for each year (for 1985, 1986, etc).

