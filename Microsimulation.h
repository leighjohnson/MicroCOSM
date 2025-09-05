#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>

//using std::cout;
//using std::vector;
//using std::ifstream;
using namespace std;

// ---------------------------------------------------------------------------------------
// General parameters and settings
// ---------------------------------------------------------------------------------------

int StartYear = 1985;
int CurrYear;
int BaselineStart = 2025; // 2005 for calibration, 2025 for intervention scenarios
int BehavCycleCount;
int STDcycleCount;
int ProjectionTerm = 41; // 27 for calibration, 56 for intervention scenarios

const int InitPop = 20000;
const int MaxPop = 70000;
const int MaxCSWs = 200;
const int MaxCasual = 4000;
const int MaxHHsize = 50; // maximum number of people in a single household

int CofactorType = 2; // Nature of STD cofactors in HIV transmission: 0 = no cofactors,
					  // 1 = multiplicative cofactors; 2 = saturation cofactors
int CalculatePAFs = 0; // 1 if you want to calculate PAFs of STDs for HIV transmission,
					   // 0 otherwise
int SyphilisImmunity = 1; // 1 if you want to allow for resistance to syphilis reinfection
						  // when recovered but still seropositive, 0 otherwise
int KeepCurrYearFixed = 0; // 1 if you want to keep the current year fixed at the start year,
						   // 0 otherwise
int NoViralTransm22 = 0; // 1 if you want to assume no transmission of viral STDs in long-
						 // term mutually monogamous relationships, 0 otherwise
int NoViralTransm12 = 0; // 1 if you want to assume no transmission of viral STDs in long-
						 // term relationships if infected partner is monogamous, 0 otherwise
int NoBacterialTransm22 = 0; // 1 if you want to assume no transmission of non-viral STDs in
							 // long-term mutually monogamous relationships, 0 otherwise
int NoBacterialTransm12 = 0; // 1 if you want to assume no transmission of non-viral STDs in
							 // long-term relationships if infected partner is monogamous

int TransitionCalc = 0; // 0 if calculating dependent probs using tradnal actuarial methods;
						// 1 if calculating dependent probs using exact method
int AllowBalancing = 0; // 0 if no explicit balancing of rate of sexual partner acquisition
						// 1 if there is balancing analagous to TSHISA model
int AllowPartnerRateAdj = 0; // 1 if there is individual variation in rates of starting and
							 // ending partnerships
							 // 0 if there is no indiv variation, analogous to TSHISA model
int AllowHIVsuscepAdj = 0; // 1 if there is individual variation in susceptibility to HIV,
						   // independent of the other sources of variation already allowed for
						   // 0 if there is no indiv variation, analogous to TSHISA model
int ConstantInitialHIV = 1; // 1 if the initial # HIV infections is to be the same across 
							// all simulations, 0 if initial HIV prevalence can vary
int SetInitPrev1990 = 1; // 1 if initial HIV prevalence is assigned in 1990
int InclMSM = 1; // 1 if we are modelling MSM, 0 if we are modelling only heterosexuals
int HCTuncertainty = 0; // 1 if allowing for uncertainty in HIV testing parameters
int IncludePopWeights = 1; // 1 if weighting by age/sex/race to SA population totals, 0 if not
int StructIntScenario = 0; // 0 = no intervention, 1 = single-session alcohol, 2 = multi-session
						   // alcohol, 3 = cash transfer, 4 = school support, 5 = vocational
						   // training, 6/7 = gender transformative (community-/indiv-level)

// Change the following indicators to specify which STDs you want to model (1 = include,
// 0 = exclude):
int HIVind = 1;
int HSVind = 1;
int TPind = 1;
int HDind = 0;
int NGind = 1;
int CTind = 1;
int TVind = 1;
int BVind = 1;
int VCind = 1;

// Change the following indicators to specify for which STDs you want to generate 
// calibration outputs (1 = include, 0 = exclude)
int HIVcalib = 0;
int HSVcalib = 0;
int TPcalib = 0;
int HDcalib = 0;
int NGcalib = 0;
int CTcalib = 0;
int TVcalib = 0;
int BVcalib = 0;
int VCcalib = 0;
int MSMcalib = 0;
int MSMcalibHIV = 0;
int StructuralDriverCalib = 0;
int StructuralRCTcalib = 1;

int CycleS = 48; // Number of sexual behaviour cycles per year. Must be a multiple of 12.
int CycleD = 48; // Number of STD cycles per year (NB: must be integer multiple of CycleS)

// ----------------------------------------------------------------------------------------
// Sexual behaviour parameters and arrays
// ----------------------------------------------------------------------------------------

double HighPropnM, HighPropnF; // % of males and females with propensity for >1 partner
double ConscientiousEffectSex; // Increase in odds high risk per SD increase conscientiousness
double HighPropn15to24[2]; // % of sexually experienced 15-24 year olds with propensity for
						   // >1 partner, by sex
double AssortativeM, AssortativeF; // Degree of assortative mixing in males and females
double GenderEquality; // Gender equality factor
double AnnNumberClients; // Average annual number of sex acts FSW has with clients
double DebutMedian[2]; // Median age of sexual debut in high-risk group, by sex
double DebutShape[2]; // Shape of log-logistic distribution for sexual debut in high-risk group
double SexualDebut[81][2]; // Continuous rate at which individuals in high risk group start
						   // sexual activity (by age and sex)
double DebutAdjLow[2]; // Factor by which the rate of sexual debut in the high risk group is
					   // multiplied in order to get the rate of debut in the low risk group
double DebutAdjRace[3]; // Factors by which baseline sexual debut rates are multiplied by race
double RRdebutLogDropIncomeF; // RR debut in females per log decrease in income below
							  // the national mean log income
double RRdebutInSchool[2]; // RR sexual debut in youth who are in school, by sex

// Short-term partner acquisition

double PartnershipFormation[2][2]; // Average annual number of new partners, by risk (first 
								   // index) & sex (2nd index)
double BasePartnerAcqH[2]; // Rate at which a single 15-19 yr old in high risk group acquires
						   // next sexual partner, by sex (base for other calcs)
double AgeEffectPartners[81][2]; // Factor by which average rate of partnership formation
								 // is multiplied to get age-specific rates (by age & sex)
double GammaMeanST[2]; // Mean of gamma density used to determine AgeEffectPartners
double GammaStdDevST[2]; // Std deviation of gamma density used to determine AgeEffectPartners
double PartnerEffectNew[2][2]; // Factor by which average rate of partnership formation is
							   // multiplied to get prob of acquiring ADDITIONAL partner
							   // (1st index is type of current rel, 2nd index is sex)
double RaceEffectNew[3]; // Adjustment to rate of secondary partner acquisition by race
double RaceMixing[3]; // Relative prob of selecting a partner of a different race
double HIVeffectPartners[6]; // Factor by which average rate of partnership formation is
							 // multiplied to get rates by HIV stage
double RR_STemployedM; // RR of ST partners in men who are employed

// Marriage

double MarriageConstant[2]; // Scale parameter for log-logistic distribution (M, F)
double MarriageTrend[2]; // Effect of birth cohort on log-logistic distribution (M, F)
double MarriageShape[3][2]; // Shape parameter for log-logistic distribution (M, F)
double MarriageMin[2]; // Minimum age at which marriage can occur (M, F)
double DivorceAdj; // Multiplier applied to empirically-derived rates of union dissolution
double DivorceTrend; // Annual change in rates of union dissolution (multiplier)
double ORremarriage[2]; // Odds of remarriage in recently divorced/widowed, relative to never-married
double ORmarriage2ndaryEdu[2]; // Odds ratio for marriage in people with completed 2ndary edu
double ORmarriageTertiary[2]; // Odds ratio for marriage in people with completed tertiary edu
double ORmarriageInSchool; // Odds ratio for marriage in people currently in school/university
double MarriageIncidence[81][2][3]; // Continuous rate at which marriages are formed, by age &
								 // sex & race (note: this is a departure from the Excel model)
double InSTrelationship[81][2][3]; // % of unmarried pop in a ST relationship, by age, sex & race
double RaceMarriage[3]; // Relative incidence of marriage by race
double BaseOddsMarriage[2]; // Base odds of marriage in 1985, by sex
double OR_MarriedAge[3][2]; // Effect of age, age^2, age^3 on odds of marriage in 1985, by sex
double OR_MarriedEmployed[2]; // Effect of employment on odds of marriage in 1985, by sex
double OR_MarriedRace[3][2]; // Effect of race on odds of marriage in 1985, by sex
double OR_MarriedEdu[5][2]; // Effect of education on odds of marriage in 1985, by sex
double OR_MarriedAgeEmployed[2]; // Interaction between age and exmployment, effect on marriage

// Sex work

double FSWcontactConstant; // Constant term used in determining rate at which men in high 
						   // risk group have contact with FSWs
double MeanFSWcontacts; // Ave annual # sex acts with FSWs, per sexually experienced male
double AgeEffectFSWcontact[81]; // Factor by which average rate of FSW contact is multiplied
								// to get age-specific rates
double GammaMeanFSW; // Mean of gamma density used to determine AgeEffectFSWcontact 
double GammaStdDevFSW; // Std deviation of gamma density used to determine AgeEffectFSWcontact
double PartnerEffectFSWcontact[5]; // Factor by which average rate of FSW contact is 
								   // multiplied to take account of current partners (0=>
								   // no partner, 1=>1ST, 2=>1LT, 3=>2ST, 4=>1ST & 1LT)
double SWurban[2]; // RR of SW contact in urban and rural areas
double RR_FSWcontactEmployedM; // RR of FSW contact for men who are employed
double InitFSWageDbn[16]; // Initial proportion of sex workers at each age
double FSWageMean[81]; // Mean age of FSWs, used to calculate InitFSWageDbn
double FSWageSD[81]; // Standard deviation of ages of FSWs, used to calculate InitFSWageDbn
double FSWentry[16][3]; // Relative rates of entry into FSW group, by age and race
double FSWexit; // Rate of exit from the FSW group
double DurFSW[81]; // Average duration of sex work (in years), used to calculate FSWexit
double HIVeffectFSWentry[6]; // Effect of HIV stage on rate of entry into FSW group
double HIVeffectFSWexit[6]; // Effect of HIV stage on rate of exit from FSW group
double RR_FSWentryDiagnosed; // RR of entry into sex work if diagnosed positive

// Relationship duration and age mixing

double MeanDurSTrel[2][2]; // Mean duration (in years) of short-term relationships, by male
						   // risk group (1st index) & female risk group (2nd index)
double RRbreakupBinge; // RR of ST relationship breakup if individual is a binge drinker
double LTseparation[16][2]; // Annual rate (cts) of separation/divorce for LT unions, by
							// age and sex
double AgePrefF[16][16]; // Proportion of male partners in each age group, for women of each
						 // age (1st index: female age, 2nd index: male age)
double AgePrefM[16][16]; // Proportion of female partners in each age group, for men of each
						 // age (1st index: male age, 2nd index: female age)
double AgePrefAdj2ndary; // >1 implies greater age disparity in 2ndary partnerships

// Coital frequencies and condom use

double FreqSexST[16][2]; // Ave # sex acts per ST relationship, per sexual behaviour cycle
						 // (by age and sex)
double FreqSexLT[16][2]; // Ave # sex acts per LT relationship, per sexual behaviour cycle
						 // (by age and sex)
double BaselineCondomUse; // % of sex acts protected among 15-19 females in ST rels in 1998
double BaselineCondomSvy; // As above, but before applying bias adjustment
double RelEffectCondom[3]; // Effect of partnership type on odds of condom use at baseline
double AgeEffectCondom[3]; // Effect of age on odds of condom use, by partnership type 
double RatioInitialTo1998[3]; // Ratio of inital odds of condom use to odds in 1998
double RatioUltTo1998[3]; // Ratio of ultimate odds of condom use to odds in 1998
double MedianToBehavChange[3]; // Median time (in years since 1985) to condom behav change
double ShapeBehavChange[3]; // Weibull shape parameter determining speed of behaviour change
double CondomUseST[16][2]; // Propn of sex acts that are protected in ST rels, by age and sex
double CondomUseLT[16][2]; // Propn of sex acts that are protected in LT rels, by age and sex
double CondomUseFSW; // Propn of sex acts that are protected in FSW-client relationships
double CondomEdu[14]; // Effect of educational attainment on odds of condom use
double CondomEduMedian; // Effect of educational attainment on median time to behaviour change
double CondomRace[3]; // Effect of race on odds of condom use
double SDpartnerRateAdj; // Std deviation of rates of partnership formation and dissolution, if
						 // AllowPartnerRateAdj = 1
double SDcondomPref; // Std deviation of innate condom preference
double BehavBiasVar[3][2]; // Sample variance of bias estimates (on logit scale), by type of
						   // behav (1st index) & sex (2nd index)
double CondomScaling; // Parameter to allow for bias in reporting of condoms (1 = no bias)
double ORcondomBingePW; // Reduction in odds of condom use per day of binge drinking, per week
double ORcondomPerYrSchool; // To determine CondomEdu

// Casual sex

double CasualEntryHet[2]; // Annual rate of entry to casual sex, high risk M & F aged 20 
double CasualAgeAdjHet[2]; // Factor by which entry to casual sex reduces per yr of age, M & F
double CasualLowAdjHet[2]; // Factor by which casual sex is reduced in low risk heterosexuals
double CasualHighAdjHet; // Factor by which entry to casual sex is reduced in high risk 
						 // heterosexuals who have a regular partner
double CasualSexFreqHet; // Monthly frequency of casual sex for heterosexuals having casual sex
double RRcasualBinge[2]; // RR entry into casual sex if people binge drink >=once per month (M, F)
double RRcasualEmployedM; // RR entry into casual sex for men who are employed
double RRcasualLogDropIncomeF; // RR entry into casual sex in women, per log decrease in income 
							   // per capita below national mean log income
double CondomUseCasualHet[16][2]; // Prob of condom use in heterosexual casual sex, by age & sex
int TotCurrCasualHet[3][2]; // # people engaging in casual heterosexual sex currently, by race & sex

int MaxNewPartnerInd; // Indicates if the maximum number of new partners has been formed (1 = yes)
double CumAgePrefM[16][16];
double CumAgePrefF[16][16]; 

// ----------------------------------------------------------------------------------------
// Arrays for balancing male and female sexual activity
// ----------------------------------------------------------------------------------------

double DesiredSTpartners[2][3]; // Desired number of new partners, by risk group (1st index)
								// and sex (2nd index - 3rd level corresponds to MSM)
double DesiredPartnerRiskM[2][2]; // Desired proportion of partners in high and low risk
								  // groups, by male risk group
double DesiredPartnerRiskF[2][2]; // Desired proportion of partners in high and low risk
								  // groups, by female risk group
double AdjSTrateM[2][2]; // Adjustment to rate at which males form partnerships with females
double AdjSTrateF[2][2]; // Adjustment to rate at which females form partnerships with males
double DesiredMarriagesM[2][2]; // Number of new marriages desired by males
double DesiredMarriagesF[2][2]; // Number of new marriages desired by females
double AdjLTrateM[2][2]; // Adjustment to rate at which males marry females
double AdjLTrateF[2][2]; // Adjustment to rate at which females marry males
double ActualPropnLTH[2][2]; // Proportion of long-term partners who are in the high-risk
							 // group, by risk group (1st index) and sex (2nd index)
double ActualPropnSTH[2][2]; // Proportion of short-term partners who are in the high-risk
							 // group, by risk group (1st index) and sex (2nd index)
double DesiredFSWcontacts[3]; // Total numbers of contacts with sex workers (per annum) by men
						      // in the high risk group, by race
double RequiredNewFSW; // Number of women becoming sex workers in current behaviour cycle,
					   // in order to meet excess male demand
int TotCurrFSW[3]; // Total female sex workers at current time, by race

// ----------------------------------------------------------------------------------------
// Demographic parameters and arrays
// ----------------------------------------------------------------------------------------

double HIVnegFert[35][3]; // Fertility rates in HIV-negative women, by age and race
double SexuallyExpFert[35][3]; // Fertility rates in sexually experienced women
double SexuallyExpFertBase[35][3];
double FertilityTable[35][81][3]; // Fertility rates in HIV-negative women by age, year, race
double InfantMort1st6mM[81][3]; // Prob of death in 1st 6 months of life (males) by year, race
double InfantMort1st6mF[81][3]; // Prob of death in 1st 6 months of life (females), by year, race
double NonAIDSmortM[91][81][3]; // Non-AIDS mortality rates in males by age, year, race
double NonAIDSmortF[91][81][3]; // Non-AIDS mortality rates in females by age, year, race
double StartPop[91][6]; // Numbers of males and females at each individual age (0, ..., 
						// 90) and race as at the start of the projection (BM, BF, CM, ...)
double MaleBirthPropn[3]; // Propn of births that are male, by race
double LE_West26[91][2]; // Life expectancy, West level 26 life table, at 4% discount rate
double ThembisaTot[18][2][81]; // Thembisa totals for population, by age, sex and year
double RaceWeights[18][2][41][3]; // Propn in each race group (final index) by age, sex & year

// ----------------------------------------------------------------------------------------
// STD parameters not defined in the classes below
// ----------------------------------------------------------------------------------------

double OIincidence[7]; // Annual rate of OI incidence, by HIV stage
double HSVsheddingIncrease[6]; // % increase in HSV-2 shedding by HIV stage
double HSVrecurrenceIncrease[6]; // % increase in recurrence rate by HIV stage
double HSVsymptomInfecIncrease; // Multiple by which HSV-2 infectiousness increased
								// when symptomatic
double InfecIncreaseSyndrome[3][2]; // % by which HIV infectiousness increases when
									// experiencing syndrome
double SuscepIncreaseSyndrome[3][2]; // % by which HIV susceptibility increases when
									// experiencing syndrome
double MaxHIVprev; // Maximum HIV prevalence in any cohort in current STD cycle
double RelHIVfertility[6]; // Factor by which fertility rate is multiplied in each HIV stage
double PropnInfectedAtBirth; // Propn of children born to HIV+ mothers who are infected at
							 // or before birth
double PropnInfectedAfterBirth; // Propn of children born to HIV+ mothers who are infected 
								// after birth (through breastfeeding)
double MaleRxRate; // Rate at which adults males seek STD treatment
double MaleTeenRxRate; // Rate at which males aged <20 seek STD treatment
double FemRxRate; // Rate at which adults females seek STD treatment
double FemTeenRxRate; // Rate at which females aged <20 seek STD treatment
double FSWRxRate; // Rate at which female sex workers seek STD treatment
double InitMaleRxRate;
double InitMaleTeenRxRate;
double InitFemRxRate;
double InitFemTeenRxRate;
double InitRecurrenceRateM;
double InitRecurrenceRateF;
double PropnTreatedPublicM; // Propn of male STD cases treated in public health sector
double PropnTreatedPublicF; // Propn of female STD cases treated in public health sector
double PropnTreatedPrivateM; // Propn of male STD cases treated in private health sector
double PropnTreatedPrivateF; // Propn of female STD cases treated in private health sector
double PropnPublicUsingSM[81]; // Propn of providers in public sector using syndromic mngt
double PropnPrivateUsingSM[81]; // Propn of providers in private sector using syndromic mngt
double DrugShortage[81]; // % redn in public sector treatment effectiveness due to drug
						 // shortages, by year
double HCT1stTime[81][2]; // Annual rates of 1st-time HIV testing at age 25, by sex
double OIsDiagnosed[81]; // Fraction of OI patients tested for HIV
double ProbOIprecedesAIDSmort; // Fraction of HIV-related deaths where OI immediately precedes
double HAARTaccess[81]; // % of children progressing to AIDS who start HAART, by year
double ARTeligiblePreg[81][2]; // % of pregnant women eligible to start ART, CD4 200-349 & 350+
double ARTeligibleOI[81][4]; // % of OI patients eligible to start ART, by CD4 stage
double ARTeligibleGen[81][4]; // % of HIV+ adults eligible to start ART, by CD4 stage
double ARTuptakePreg[81]; // % of eligible pregnant women who start ART soon after diagnosis
double ARTuptakeOI[81]; // % of ART-eligible OI patients who start ART soon after diagnosis
double ARTinitiationU200[81][2]; // Ann ART initiation in diagnosed adults with CD4 <200 (M, F)
double ARTinterruption; // Annual rate of ART interruption (in women aged 35+)
double ARTresumption; // Annual rate of resuming ART after an interruption
double RRinterruptionM; // RR of ART interruption in men (compared to women)
double ARTinterruption20mult; // RR of ART interruption at age 20 relative to adult minimum
double AgeAdjInterrupt[81][2]; // RR of ART interruption by age and sex
double PMTCTaccess[81]; // % of pregnant women who are offered PMTCT for HIV, by year
double PropnCiproResistant[81]; // Propn of NG cases that are ciprofloxacin-resistant, by year
double PropnCiproTreated[81]; // Propn of treated NG cases that are given ciprofloxacin
double RxPhaseIn[81]; // Unscaled rates of phase-in for future STD interentions, by year
double MMCtoM_ratio[81]; // Ratio of # MMC operations to men aged 15-49
double AnnPrEPuptake[81]; // Annual rate of PrEP uptake in individuals who are eligible
double PrEPuptakeMSM[81]; // Annual rate of PrEP uptake in MSM
double AnnVaccineUptake[81]; // Annual rate of vaccination in eligible individuals
double HIVdiagNoBF[81]; // % of HIV-diagnosed mothers who choose not to breastfeed
double InitDrugEffNG; // Initial effectiveness of drugs for treating NG
double InitCorrectRxHSV;
double InitCorrectRxTVM;
double VCTageMean[2]; // Average age at which adults are tested for HIV, by sex
double VCTageSD[2]; // Std deviation of ages at which adults are tested for HIV, by sex
double VCTageEffect[81][2]; // Relative rates of HIV testing by age (relative to age 25)
double RetestAdj; // Relative rate of HIV testing if previously tested
double VCTadjEdu; // RR of HIV testing per year of increase in edu attainment
double RR_ARTstartPer100CD4; // Relative rates of ART initiation per 100 increase in CD4
double AcceptScreening; // % of pregnant women offered HIV screening who accept
double AcceptNVP; // % of women testing positive who agree to receive nevirapine
double RednNVP; // % reduction in HIV transmission at/before birth if woman receives NVP
double RednFF; // % reduction in HIV transmission after birth if woman receives formula
			   // feeding OR exclusive breasfeeding
double SecondaryRxMult; // Factor by which the rate of treatment seeking is multiplied
						// when experiencing symptoms of secondary syphilis
double SecondaryCureMult; // Factor by which the probability of cure is multiplied
						  // when treated for secondary syphilis
double FSWasympRxRate; // Rate at which female sex workers seek STD treatment when 
					   // asymptomatically infected with an STD
double FSWasympCure; // Prob that treatment for symptomatic STD in FSW cures other 
					 // asymptomatic STDs (expressed as multiple of the prob of cure if the
					 // STD was symptomatic)
double InitFSWasympRxRate;
double InitFSWasympCure; 
double InitANCpropnScreened;
double InitANCpropnTreated;
double InitHIVprevHigh; // % of high risk group initially HIV-positive (assumed to be asymp)
double InitHIVprevRace[3]; // Adjustment factors to initial HIV prevalence by race
double HlabisaRatio[7][2]; // Ratio of age- and sex-specific initial prev to 15-49 fem prev
double InitHIVtransm[3][2]; // HIV transm probs per sex act at start of epidemic, by
							// nature of rel (1st index) and sex (2nd index)
double RatioUltToInitHIVtransm; // Ratio of ultimate HIV transmission prob to initial HIV
								// transm prob, in the 'no STD cofactor' scenario
double SDsuscepHIVadj; // Coefficient of variation in HIV susceptibility (after controlling 
					   // for other factors in the model)
double CurrPerinatal; // Propn of HIV-exposed kids infected perinatally in current year
double CurrPostnatal; // Propn of HIV-exposed kids infected postnatally in current year
double VLeffectTransm[2]; // Parameters determining effect of HIV VL on transmission prob
double HIVsusceptInjectable; // Factor by which M-to-F transm risk increases if woman uses
							 // injectable contraception
double HIVinfectHormonal; // Factor by which F-to-M transm risk increases if woman uses
						  // hormonal contraception
double PregEffectHIVtransm[2]; // Effect of pregnancy on M->F and F->M transmission of HIV
double HSVprevRace[3]; // RR of HSV-2 by race, in 1985
double PrEPefficacy; // PrEP efficacy if adherence is perfect (heterosexual couples)
double PrEPefficacyMSM; // PrEP efficacy if adherence is perfect (same-sex couples)
double PrEPdiscontinue; // Annual rate of discontinuing PrEP (not due to eligibility changes)
double VaccineWane; // Rate of waning of HIV vaccine protection (per month)
double HIVvaccProt[2]; // Reduction in HIV risk after prime and after prime+boost
double VaccAdherence; // Probability of return for each follow-up vaccine dose after first
double MMCuptake[16]; // RR of MMC uptake among uncircumcised men (relative to 10-14 age group)
double RR_VMMClogIncome; // RR of MMC uptake per log increase in household income
double MCprevBaseline[91][3]; // Fraction of men who are circumcised in 1985, by age & race
double PublicANCuse[3]; // % of pregnant women who use public antenatal clinics, by race
double StillbirthRate[2]; // % of pregnancies resulting in stillbirth, by HIV status

// Additional parameters for HIV testing Modelling Consortium project
double HBCTuptake[2]; // % of household members at home and agreeing to testing through HBCT
double OR_STuptake; // OR for uptake of ST compared to uptake of regular testing
double MobileTestUptake[2]; // Annual rate of uptake of mobile HCT if living in a community
							// where mobile HCT is offered (without/with community mobilization)
double SchoolTestUptake[2]; // Prob of HIV testing if school offers HIV testing (M and F)
double SchoolTestVirginRR; // RR of HIV testing in schools if not sexually experienced
double ANCpartnerTested[2][2]; // % of male partners of pregnant women tested (by female HIV 
							   // status (1st index) and marital status (2nd index)
double WorkTestUptake[2]; // % accepting offer of HIV testing through workplace (M and F)
double EmployedReachable; // % of employed population reachable through workforce testing
double FPCtestUptake; // annual testing rate in FPCs that offer HIV testing
double RetestAdjDiagnosed[15]; // RR of HIV testing in previously-diagnosed, by testing modality
							// 0=general, 1=OI, 2=ANC, 3=HBCT, 4=mobile, 5=STI, 6=FPC, 7=FSW,
							// 8=MSM, 9=school, 10=workplace, 11=partner, 12=prison
double RetestAdjART[15]; // RR of HIV testing in patients on ART (same coding as for Diagnosed)
double RR_ARTstartCommunity; // RR of ART initiation if diagnosed through community-based testing
double RR_ARTstartST; // RR of ART initiation if diagnosed through self-testing
double PropnSTreferral[81]; // Propn pregnant women given ST kits to give to partners
double PropnAssistedNotif[81]; // Propn newly-diagnosed offered assisted partner notification
double HBCTfreq[81][2]; // # times per annum home-based HCT is offered (rural then urban)
double PropnST_HBCT[81]; // Propn offered ST kits as part of home-based HCT
double MobileTestCoverage[81][2]; // % of people living in communities where mobile HCT is offered
double CommunityMobilization[81]; // % of mobile HCT services with community mobilization
double FSWtestUptake[81]; // Annual rate of uptake of HCT through FSW-specific programmes
double MSMtestUptake[81]; // Annual rate of uptake of HCT through MSM-specific programmes
double SchoolTestFreq[81]; // % of high schools in which HCT is offered
double ANCpartnersInvited[81][2]; // % of partners of pregnant women (HIV-/HIV+) invited for HCT
double TestingPrisons[81]; // Annual rate of testing in prisons
double TestingHIVinSTIs[81]; // % of STI patients tested for HIV
double TestingWorkplace[81]; // % of employed individuals offered HIV testing on an annual basis
double TestingFPC[81]; // % of family planning clinics offering women annual HIV testing
double ANCretestFreq[81]; // Rate of retesting in pregnant women in late pregnancy
double BaseModalityCost[8][2]; // Cost per HIV test by baseline modality and test outcome (+ or -)
double NewModalityCost[12][2]; // Cost per HIV test by new modality and test outcome (+ or -)
double AnnARTcost; // Annual cost of ART provision per patient on ART ($)
double AnnCTXcost; // Annual cost of cotrimoxazole provision per patient ($)
double CondomCost; // Cost per condom ($)
double MMCcost; // Cost per MMC operation ($)
double PMTCTcost[2]; // Cost of short-course PMTCT ($), depending on diagnosis in labour/before
double AnnPrEPcost; // Annual cost of PrEP provision ($) per person
double InfantTestCost[2]; // Cost of PCR test in infants at birth or at 6 weeks ($)
double PalliativeCareCost; // Annual cost of palliative care per patient ($)

double NewHIVbySex[2]; // Expected new HIV cases in current year, by sex
double NewCTbySex[2]; // New chlamydia cases in current year, by sex
double NewNGbySex[2]; // New gonorrhoea cases in current year, by sex
double NewTVbySex[2]; // New trichomoniasis cases in current year, by sex
double SummOut[1000][56]; // Summary outputs
int SummOutRow; // Row ID for Summary outputs

// --------------------------------------------------------------------------------------
// Parameters related to conception and contraception
// --------------------------------------------------------------------------------------

double MeanGestation; // Average time from conception to birth (in years)
double SDgestation; // Std deviation of times from conception to birth (in years)
double SDfecundability; // Std deviation of fecundability adjustment factor (mean = 1)
double RateFecund[8]; // Annual rates of becoming fecund, ages 12-19
double RateInfecund[15]; // Annual rates of becoming infecund, ages 35-49

// Note that 'contraception' here refers to injectable/pill/female sterilization.

// Contraception prevalence assumptions: 1985
double InitContr[3][3]; // Prevalence of contraception in 1985 by age (15-24, 25-34, 35-49) & race
						// among married women with previous birth, not currently pregnant
double InitInjectable[3][3]; // Injectable use among contracepters in 1985, by age & race
double ORcontrUnmarriedNeverPreg[2][2][3]; // OR for contraceptive use if unmarried (1st index),  
										   // never pregnant (2nd index), by race (3rd index)
double ORcontrEdu; // OR for contraceptive use per additional grade completed
double ORcontrAbstinent; // OR for contraceptive use if not sexually active
double ContrVirgin; // Prevalence of contraceptive use among virgins aged 15-49, in 1985
double ORinjectableEdu;  // OR for injectable use per additional grade completed
double ORinjectableNeverPreg;  // OR for injectable use if never pregnant
double SterilizationNotPill[7]; // % of women sterilized by age (as % of sterilized + pill)

// Contraception uptake and discontinuation
double SterilizationRates[7]; // Annual rate of sterilization in women, by age group
double RaceEffectSteril[3]; // Race-specific adjustments to rates of sterilization
double StartContrNewRel[3]; // Prob of starting contraception at time of new partner, by race
double StartContrPostnatal[3]; // Prob of starting contraception postnatally, by race
double StartContrOther; // Annual rate of starting hormonal contraception for other reasons
double PrevBirthEffectContr[3]; // OR of starting contraception if no previous birth, by race
double AgeEffectContr[7][2]; // OR of starting hormonal contraception by age group
double EduEffectContr; // OR for hormonal contraceptive use per additional grade completed
double HIVeffectContr; // OR for hormonal contraceptive use if HIV+
double CondomEffectContr; // OR for hormonal contraceptive use if consistently using condoms
double NoPrevUseContr; // OR for hormonal contraceptive use if never previously used
double NoPrevUseContr1; // As before, for period 1985-96
double NoPrevUseContr2; // As before, for period 2000 & after
double InjectablePref[3]; // % of women starting hormonal contr who choose injections, by race
double ORinjectableAge; // OR for injectable use per 5-year increase in age
double ORinjectableInit; // OR for injectable use comparing 1985 to 1998
double NewContrMethodWeight; // Prob of considering dif method from previous hormonal method
double StopContr[2]; // Monthly rate of stopping contraception (not sexually active, active)
double StopContrAbsentPartner; // Monthly rate of stopping contraception if partner away
double ORhormonalByYear[81][3]; // OR of hormonal use by year and age (15-24, 25-34, 35-49)
double RRsterilizationByYear[81]; // Relative rate of sterilization by year

// Contraceptive efficacy
double CondomEffPreg; // Condom effectiveness in preventing pregnancy
double ContrEffPreg[3]; // Effectiveness of injections, pills, sterilization in preventing preg 

// --------------------------------------------------------------------------------------
// Parameters related to breastfeeding 
// --------------------------------------------------------------------------------------

double EverFeed; // propn of HIV-neg/undiagnosed mothers who ever breastfeed
double MedianFeed[2]; // median duration of feeding (in months) by maternal HIV diagnosis
					  // (undiagnosed/HIV-neg, diagnosed positive)
double ShapeFeed[2]; // shape parameter for duration of feeding by maternal HIV diagnosis

// --------------------------------------------------------------------------------------
// Parameters related to urbanization and migration
// --------------------------------------------------------------------------------------

double InitUrbanPropn[18][8][2]; // % of 1985 pop urbanized, by age, edu/race, sex
double UrbanToRural[18][2]; // Annual rate of urban-to-rural migration, by age and sex
double RuralToUrban[18][2]; // Annual rate of rural-to-urban migration, by age and sex
double UrbanAdjRace[3]; // Racial adjustments to urban-to-rural migration
double RuralAdjRace[3]; // Racial adjustments to rural-to-urban migration
double UrbanAdjEdu[6]; // Adjustments to urban-to-rural migration by highest grade passed
double RuralAdjEdu[6]; // Adjustments to rural-to-urban migration by highest grade passed
double VisitFreq[2]; // Gamma dbn parameters for desired visit frequency
double LocationMixing; // RR of choosing indiv as partner if not in same location
double MaxVisitLength; // Maximum ave duration of visit to partner (in days)
double MaxCoitalFreqWeight; // Determines relative feq of sex during partner visits
double AbsentPartnerConcurAdj; // Determines adj to desired # partners if primary sexual
							   // partner is absent
double MarriedMigAdj[3]; // RR of migration for married cohabiting people, by race
double SeparatedMigAdj[3]; // RR of migration for married non-cohabiting people, by race

// --------------------------------------------------------------------------------------
// Parameters related to HIV disclosure
// --------------------------------------------------------------------------------------

double DiscloseProb; // Base prob of disclosure to partner immediately after HIV diagnosis
double DiscloseRate; // Base rate of disclosure to partner if individual did NOT
					 // disclose their HIV status immediately after HIV diagnosis
double DiscloseMale; // Increase to base prob/rate if the diagnosed individual is male
double DiscloseART; // Increase to base prob/rate if the diagnosed individual is on ART
double DiscloseMarried; // Increase to base prob/rate if the couple is married/cohabiting
double DiscloseEffectRelDur; // Effect of disclosure on relationship duration, if F+ & M-
double DiscloseEffectCondom; // Increase in condom use in serodiscordant relationships
							 // after HIV-positive partner discloses HIV status
double ARTdiscordance; // Relative prob that an ART patient selects an indiv as a partner
					   // if that indiv is not on ART
double ProbReferral[2]; // Prob partner seeks testing after disclosure (passive, assisted)

// --------------------------------------------------------------------------------------
// Parameters related to incarceration
// --------------------------------------------------------------------------------------

double PrisonEntryRate[16]; // Rates of entering prison by age (never prev. in prison)
double PrisonEntryRace[3]; // Race-specific adjustments to rate of entering prison
double PrisonEntryEdu; // Adjustment to prison entry rate if high school completed
double PrisonReentryAdj; // Adjustment to prison entry rate if previously in prison
double MeanUnsentencedDur; // Average time to sentencing
double UnsentencedRelease; // % of prisoners released without being convicted/sentenced
double SentenceLength[2]; // Gamma dbn parameters for newly-sentenced prisoners
double SentenceLength85[2]; // Gamma dbn parameters for prisoners in 1985
double ProbParole; // Prob that a newly-sentenced prisoner is ultimately granted parole
double PrisonEntryRisk[2]; // RR prisone entry for high-risk and low-risk men (incl virgins)

// --------------------------------------------------------------------------------------
// Parameters related to employment
// --------------------------------------------------------------------------------------

// Parameters that determine the fraction employed in 1985

double BaseEmployed; // Odds of employment in baseline covariate category
double EduEffectEmployed[14]; // Effect of edu on odds of employment
double RaceEffectEmployed[3]; // Effect of race on odds of employment
double UrbanEffectEmployed[2]; // Effect of urban/rural location on odds of employment
double SexEffectEmployed[2]; // Effect of sex on odds of employment
double AgeEffectEmployed[11]; // Effect of age group on odds of employment

// Parameters that determine changes in rates of employment thereafter

double BaseEmployed2[81]; // Odds of employment in baseline covariate category
double EduEffectEmployed2[5]; // Effect of edu on odds of employment
double RaceEffectEmployed2[3]; // Effect of race on odds of employment up to 2008
double RaceEffectEmployed3[3]; // Effect of race on odds of employment after 2008
double UrbanEffectEmployed2[2]; // Effect of urban/rural location on odds of employment
double SexEffectEmployed2[3]; // Effect of female sex on odds of employment, by race
double AgeEffectEmployed2[10]; // Effect of age group on odds of employment
double HIVeffectEmployed[6]; // Effect of HIV stage on odds of employment, grade <12
double ChildEffectEmployed; // Effect of being a child carer in women
double PrevEmployedEffect; // Effect of being employed 1 year previously on odds employed
double ORunemployedTraining; // OR for effect of training intervention on unemployment

// --------------------------------------------------------------------------------------
// Parameters related to income and grants
// --------------------------------------------------------------------------------------

double BaseLogIncome[5]; // Base log income if employed in 2002, 2006, 2010, 2014, 2018
double FemaleEffectLogIncome[5]; // Effect of female sex on log income, in each of 5 yrs
double EduEffectLogIncome[5][5]; // Effect of edu (5 levels) on log income, in each of 5 yrs
double RaceEffectLogIncome[3][5]; // Effect of race on log income, in each of 5 yrs
double AgeEffectLogIncome[9][5]; // Effect of age (9 levels) on log income, in each of 5 yrs
double ConscientiousEffectLogIncome; // Change in log income per SD change conscientiousness
double StdDevLogIncome; // Standard deviation of log incomes, controlling for other factors
double AnnRealEarningsGrowth; // Annual growth in earnings, adjusted for inflation
double CurrAveIncome[10][2][3][5]; // Current ave monthly log income by age, sex, race, edu

double BaseOddsPension; // Base odds of having a private pension
double EffectEduPension[5]; // Effect of education on odds of having a private pension
double EffectRacePension[3]; // Effect of race on odds of having a private pension
double PrivatePension2019; // Monthly amount (ZAR) of average private pension in 2019
double AnnRealPensionGrowth; // Annual growth in private pensions, adjusted for inflation
double CurrPrivatePension; // Monthly amount of private pension in current year

double ConsPriceIndex[81]; // Consumer price index at middle of each year
double CSGamount[81]; // Amount of child support grant, by year
double CSGageLimit[81]; // Age limit for child support grant eligibility, by year
double CSGincomeLimit[81]; // Monthly amount (ZAR) of child support grant, by year
double OAPamount[81]; // Monthly amount (ZAR) of state old age pension, by year

double ChildWeightEquivScale; // Child weight (relative to adult) in equivalence scale
double HHsizeAdjEquivScale; // Exponent in equivalence scale
double MeanLogIncome; // Based on equivalence scale
double MeanLogIncomeStore[81];
double MedianIncomeRace[3]; // Median household per capita income, per month, by race
double UnsortedIncomes[MaxPop][4]; // 2nd index: total, African, coloured, white
double SortedIncomes[MaxPop];
int TotalIncomes[4];

// --------------------------------------------------------------------------------------
// Parameters related to men who have sex with men (MSM)
// --------------------------------------------------------------------------------------

double MSMfraction; // Fraction of men who EVER have a propensity for same-sex behaviour
double BiFraction; // Fraction of MSMfraction who EVER have propensity for sex with women
double InitMalePrefBeta[2]; // Beta parameters for assigning desired % of partners who are 
							// male to bisexual men aged <= 20
double AnnChangeMalePref[2]; // Mean and SD of annual change in % of partners who are male
double InsertRecep[2]; // % of MSM who are exclusively insertive and exclusively receptive
double SameSexMarried; // Relative rate of marriage to partner of same sex
double AgePrefMSM[16][16]; // Proportion of MSM partners in each age group
double CasualEntry; // Annual rate of entry to casual sex, high risk MSM aged 20 (excl. bi)
double CasualAgeAdj; // Factor by which entry to casual sex reduces per yr of age
double CasualLowAdj; // Factor by which entry to casual sex is reduced in low risk MSM
double CasualHighAdj; // Factor by which entry to casual sex is reduced in high risk MSM
					  // who have a regular partner
double CasualExit; // Annual exit rate from casual sex in M who only have sex with other M
double CasualSexFreq; // Monthly frequency of casual sex for men engaging in casual sex
double FreqSexST_MSM; // Freq of sex per month in MSM in ST relationships
double ORcondomCasual; // Increase in odds of condom use in casual rels (relative to ST rels)
double InitHIVtransmMSM[3][2]; // Prob of HIV transm per act of anal sex, by rel type (1st
							   // index) & type of sex (receptive=0, insertive=1)
double RolePrefBi[3]; // Propn of bisexual men with RolePrefs of 0, 0.5 and 1
double RolePrefHomo[3]; // Propn of exclusively homosexual men with RolePrefs of 0, 0.5 and 1
double MeanDurST_MSM; // Average duration of ST relationships
double ActualPropnLTH_MSM[2]; // % of LT partners who are high risk
double ActualPropnSTH_MSM[2]; // % of ST partners who are high risk
double CumAgePrefMSM[16][16];
double DesiredPartnerRiskMSM[2][2]; // Desired proportion of partners in high and low risk
									// groups, by indiv risk group (1st index)
double AdjSTrateMSM[2][2]; // Adjustment to rate at which males form partnerships with males
int TotCurrCasual; // Total # MSM currently engaging in casual sex
double CondomUseCasual[16]; // Rates of condom use in casual sex, by age
double RatioInitPrevMSM; // Ratio of initial HIV prevalence in MSM to that in heterosexual men
double IntraEventVersatility; // Prob of receptive & insertive AI in same episode if both versatile

//----------------------------------------------------------------------------------------
// Education parameters
//----------------------------------------------------------------------------------------

double Grade1entry[5][3]; // Annual prob of entry to grade 1, ages 5-9, by race
double GradeRepetition[13][3][2]; // Prob of repeating current grade, by grade, race & sex
double ConscientiousEffectGradeRep; // Increase in odds grade rep per SD increase conscientiousness
double ParentEduEffectGradeRep; // Increase in odds grade rep per year increase in parent education
double Dropout[13][3][2]; // Annual prob of dropout, by grade, race and sex
double ParentEduEffectDropout; // Increase in odds dropout per year increase in parent education
double PregDropout[3]; // Prob of permanent dropout due to pregnancy, by race
double TertiaryEnrol[3]; // Prob of tertiary enrolment if matric completed, by race
double ConscientiousEffectTertiaryEd; // Increase in odds tertiary enrolment per SD conscientiousness
double CumGrade1985[91][13][3]; // Cumulative % with highest grade <=x, by age, x & race
double InSchool1985[31][13][3]; // % currently in school, by age, grade & race
double DropoutRedn[3]; // Factor by which dropout rate reduces per annum, by race
double EduAssort; // Assortativeness of mixing by edu attainment (0 => random mixing)
double EduMixing[14][14]; // Relative prob that indiv with edu i selects partner with edu j
double RRdropoutSupport; // RR of school dropout if receiving educational support
double RRdropoutIncSupport; // RR of school dropout for R800 income support (2005 ZAR)
double ProbReturnSupport; // Prob of return to school if receiving educational support

//----------------------------------------------------------------------------------------
// Parameters related to inequitable gender norms
//----------------------------------------------------------------------------------------

double AgeEffectIneqGender; // OR for effect of age (years) on inequitable gender norms
double EduEffectIneqGender[3]; // OR for effect of education on inequitable gender norms
double RuralEffectIneqGender; // OR for effect of rural location on inequitable gender norms
double RaceEffectIneqGender[3]; // OR for effect of race on inequitable gender norms
double HighRiskEffectIneqGender; // OR for effect of being high risk on inequitable gender norms
double EffectIneqGenderConcurrency; // Increase in concurrency rate if endorse ineq gender norms
double EffectIneqGenderCasual; // RR entry into casual sex per 0.1 decrease in ineq gender norms
double ORcondomGenderIneq; // OR for condom use if always endorse inequitable gender norms (vs never)
double RRgenderIneqIndiv; // RR inequitable norms after indiv-level gender-transformative intervention 
double RRgenderIneqComm; // RR inequitable norms after community-level intervention
double ProbGenderIneqReversion; // Annual prob of reverting to baseline gender norms

//----------------------------------------------------------------------------------------
// Parameters related to alcohol consumption
//----------------------------------------------------------------------------------------

double BaseDrinkProb[2]; // # drinking days per week (log scale), by sex
double UrbanEffectDrinkProb[2]; // Increase in # drinking days/week in urban dwellers, log scale
double AgeEffectDrinkProb[5][2]; // Age effect on # drinking days/week, log scale
double MarriedEffectDrinkProb[2]; // Increase in # drinking days/week if married, log scale
double EduEffectDrinkProb[4][2]; // Educational effect on # drinking days/week, log scale
double RaceEffectDrinkProb[3][2]; // Race effect on # drinking days/week, log scale
double StdDevDrinkProb[2]; // Standard deviation in daily drinking prob, controlling for covariates

double BaseDrinksPerDD[2]; // # drinks per drinking day (log scale), by sex
double UrbanEffectDrinksPerDD[2]; // Increase in # drinks/drinking day in urban dwellers, log scale
double AgeEffectDrinksPerDD[5][2]; // Age effect on # drinks/drinking day, log scale
double EmployedEffectDrinksPerDD[2]; // Increase in # drinks/drinking day if employed, log scale
double EduEffectDrinksPerDD[4][2]; // Educational effect on # drinks/drinking day, log scale
double RaceEffectDrinksPerDD[3][2]; // Race effect on # drinks/drinking day, log scale
double ConscientiousEffectDrinksPerDD; // Effect of conscientiousness score
double GenderIneqEffectDrinksPerDD; // Effect of endorsing inequitable gender norms
double StdDevDrinksPerDD[2]; // Standard deviation in # drinks/DD, controlling for covariates
double ConfoundingAlcSex; // Extent of confounding in assn between alcohol and high-risk group (0-1)

double AlcPropnReported; // Proportion of alcohol consumption that is reported
double RatioMaxDrinksToAve; // Ratio of maximum drinks on a DD to average drinks per DD
double ProbMaxDrinksPerDD; // Prob a drinking day is one on which max drinks are consumed
double RRalcoholSingle; // RR alcohol consumption after single-session counselling
double RRalcoholMultiple; // RR alcohol consumption after multi-session counselling
double ProbAlcoholReversion; // Annual prob of reverting to baseline alcohol consumption

//----------------------------------------------------------------------------------------
// Parameters related to household (HH) formation
//----------------------------------------------------------------------------------------

// Parameters for assigning headship at start of simulation

double BaseOddsHead[2]; // Odds of being HH head in baseline category, for M and F
double MarriedEffectHead[2]; // Effect of being married on odds of being head, in M and F
double AgeEffectHead[11][2]; // Effect of adult age group on odds of being head, in M and F
double RaceEffectHead[3][2]; // Effect of race on odds of being head, in M and F
double EmployedEffectHead[2]; // Effect of employment on odds of being head, in M and F
double InSchoolEffectHead[2]; // Effect of being in school on odds of being head, in M and F
double EduEffectHead[5][2]; // Effect of edu attainment on odds of being head, in M and F
double BirthEffectHead; // Effect of having given birth on odds of being head, in F

// Parameters for assessing if non-head adults are living with parents (ALWP) at start

double BaseOddsALWP[2]; // Odds of a non-head adult living with parents in baseline category
double MarriedEffectALWP[2]; // Effect of being married on odds of living with parent
double AgeEffectALWP[8][2]; // Effect of age group on odds of living with parent
double RaceEffectALWP[3][2]; // Effect of race on odds of living with parent
double EmployedEffectALWP[2]; // Effect of being employed on odds of living with parent
double InSchoolEffectALWP[2]; // Effect of being in school on odds of living with parent
double EduEffectALWP[5][2]; // Effect of edu attainment on odds of living with parents

// Parameters for assessing if children are living with parents (CLWP) at start

double BaseOddsCLWP; // Base odds of child living with parents
double AgeEffectCLWP[3]; // Effect of age group on odds of child living with parents
double RaceEffectCLWP[3]; // Effect of race on odds of child living with parents
double InSchoolEffectCLWP; // Effect of being in school on odds of child living with parents

// Parameters for assessing if children are living with surviving parent (CLSP): 1st index 
// is for mother, 2nd index is for father conditional on child not living with mother

double BaseOddsCLSP[2]; // Base odds of child living with parent, if parent is alive
double AgeEffectCLSP[2]; // Increase in odds of living with parent per year of age
double Age2EffectCLSP[2]; // Effect of age squared
double RaceEffectCLSP[3][2]; // Effect of race on odds of living with surviving parent
double SchoolingEffectCLSP[2]; // Increase in odds of living with parent if in school

// Parameters for assessing if a migrant from urban to rural (or rural to urban)
// establishes a new household instead of joining an existing household

double BaseOddsNewHH[2]; // Odds of starting new HH in baseline category, for M and F
double MarriedEffectNewHH[2]; // Effect of being married on odds of starting new HH, in M and F
double AgeEffectNewHH[7][2]; // Effect of adult age group on odds of starting new HH, in M and F
double RaceEffectNewHH[3][2]; // Effect of race on odds of starting new HH, in M and F
double EmployedEffectNewHH[2]; // Effect of employment on odds of starting new HH, in M and F
double InSchoolEffectNewHH[2]; // Effect of being in school on odds of starting new HH, in M and F
double EduEffectNewHH[5][2]; // Effect of edu attainment on odds of starting new HH, in M and F
double UrbanEffectNewHH[2]; // Effect of urban location on odds of starting new HH, in F

// Parameters for assessing if divorced woman establishes a new household after leaving partner

double BaseOddsNHADS; // Base odds of new household after divorce/separation (NHADS)
double AgeEffectNHADS; // Increase in odds of NHADS per year of age
double Age2EffectNHADS; // Increase in odds of NHADS per unit increase in age^2
double RaceEffectNHADS[3]; // Effect of race on odds of NHADS
double EmployedEffectNHADS; // Effect of employment on odds of NHADS
double ParentAliveNHADS; // Effect of having a surviving parent on odds of NHADS

// Parameters for leaving the nest

double BaseOddsLeaveNest; // Base odds of living in a house where head is not parent/grandparent
double ORleaveNestAge; // Increase in odds of having left the nest per year of age
double ORleaveNestAge2; // Increase in odds of having left the nest per unit increase in age^2
double RRleaveNestFemale; // RR of leaving nest per year, in females
double RRleaveNestRace[3]; // RR of leaving nest per year, by race
double ORleaveNestEmployed; // OR of having left the nest if employed 
double RRleaveNestInSchool; // RR of leaving nest per year, if currently studying
double RRleaveNestHH10plus; // RR of leaving nest per year, if in house of 10 or more people
double BaseOddsHHALN; // Base odds of being household head after leaving the nest
double OR_HHALNfemale; // OR for being household head after leaving nest, in females
double OR_HHALNage; // OR for being household head after leaving nest, per year of age
double OR_HHALNage2; // OR for being household head after leaving nest, per unit increase in age^2
double OR_HHALNrace[3]; // OR for being household head after leaving nest, by race
double OR_HHALNemployed; // OR for being household head after leaving nest, if employed
double OR_HHALNinSchool; // OR for being household head after leaving nest, if currently studying
double AnnProbLeaveNest[35]; // Annual prob of leaving the nest, by age (15-49)

// Parameters for homelessness

double ProbHomeless; // Prob that indiv with no local relative becomes homeless
double BaseOddsNotHomeless; // Annual odds of ceasing to be homeless (unemployed men)
double ORnotHomelessFem; // OR for ceasing to be homeless in women (relative to men)
double ORnotHomelessEmployed; // OR for ceasing to be homeless if employed

// ---------------------------------------------------------------------------------------
// Arrays for random numbers and doing multiple runs, likelihood calculations
// ---------------------------------------------------------------------------------------

double r[InitPop], rprisk[InitPop], rpID[InitPop], rSTI[MaxPop][8];
double r2[MaxPop], revent[MaxPop], rpAge[MaxPop], rpID2[MaxPop], rpref[MaxPop];
double rbirth[MaxPop], rconcep[MaxPop];
const int ParamCombs = 50; // number of input parameter combinations
const int IterationsPerPC = 1; // number of iterations per parameter combination
const int samplesize = 50; // number of simulations (must = ParamCombs * IterationsPerPC)
int SeedRecord[ParamCombs][2]; // seeds used when FixedUncertainty = 1
int GetSDfromData = 1; // Set to 1 if you want to calculate the standard deviation in the
					   // likelihood function for the national prevalence data based on the
					   // model fit to the data rather than the 95% CIs around the survey
					   // estimates (similar to the approach in Sevcikova et al (2006)).
int VaryParameters = 1; // 1 if parameters vary across simulations
int FixedUncertainty = 1; // 0 if parameters vary stochastically, 1 if previously-
						  // generated parameters are read in from input files
int CurrSim; // Counter for the current simulation
int FixedPeriod = 0; // # years after start of projection in which we use the same seeds
					 // (across all simulations), to limit variability in HIV simulations
double TotalLogL; // The log of the likelihood for the current simulation (based on HIV and
				  // STD prevalence data)
double StructLogL; // The log likelihood in respect of HIV structural driver data
double ANCbiasVar = 0.0; // Variance of the ANC bias term (previously set to 0.0225)
double HSRCbiasVar = 0.0; // Variance of the HSRC bias term (previously set to 0.0225)
double IneqGenderBingeSvy[4][2], IneqGenderMultSvy[4][2], BingeMultSvyM[4][2], BingeMultSvyF[4][2];
double BingeCondomSvyM[4][2], BingeCondomSvyF[4][2], EmployedTransSvy[4][2], EmployedMultSvyM[4][2];
double EmployedMultSvyF[4][2], EmployedHIVsvyM[4][2], EmployedHIVsvyF[4][2], EduCondomSvyM[4][2];
double EduCondomSvyF[4][2], EduHIVsvyM[4][2], EduHIVsvyF[4][2], SchoolMarriageSvy[4][2];

// ---------------------------------------------------------------------------------------
// Classes for STD prevalence data and likelihoods
// ---------------------------------------------------------------------------------------

class PrevalenceData
{
public:
	PrevalenceData();

	double LogL;
	int Observations;
	int SexInd; // 0 = males, 1 = females
};

class NationalData : public PrevalenceData
{
public:
	NationalData();

	// Note that in defining the arrays below we are assuming that there would not be more
	// than 100 data points to calibrate to. If for some reason there are (i.e. Observations
	// >100), then one should change the dimensions of the arrays below.

	int StudyYear[100];
	double StudyPrev[100];
	double PrevSE[100];
	int AgeStart[100];
	double ExpSe[100];
	double ExpSp[100];
	double ModelPrev[100];
	double BiasMult[19];
	double ModelVarEst;

	void CalcModelVar();
	void CalcLogL();
};

class AntenatalN : public NationalData
{
public:
	AntenatalN();
};

class HouseholdN : public NationalData
{
public:
	HouseholdN();
};

class SentinelData : public PrevalenceData
{
public:
	SentinelData();

	// Note that in defining the arrays below we are assuming that there would not be more
	// than 25 data points to calibrate to. If for some reason there are (i.e. Observations
	// >25), then one should change the dimensions of the arrays below.

	int StudyYear[25];
	int StudyN[25];
	double StudyPrev[25];
	double StudyPos[25];
	double ExpSe[25];
	double ExpSp[25];
	double VarSe[25];
	double VarSp[25];
	int HIVprevInd[25]; // 1 if HIV prevalence was measured in the study, 0 otherwise
	double HIVprev[25];
	double ModelPrev[25];
	double CurrentModelPrev;
	double BiasMult;
	double VarStudyEffect; // Sigma(b) squared

	void CalcLogL();
};

class Household : public SentinelData
{
public:
	Household();

	int AgeStart[25];
	int AgeEnd[25];
	int ExclVirgins[25];
};

class NonHousehold : public SentinelData
{
public:
	NonHousehold();
};

class ANC : public NonHousehold
{
public:
	ANC();
};

class FPC : public NonHousehold
{
public:
	FPC();
};

class GUD : public NonHousehold
{
public:
	GUD();
};

class CSW : public NonHousehold
{
public:
	CSW();
};

class RCTdata
{
public:
	RCTdata(int n, int code);

	int DataPoints;
	int InterventionCode;
	double StudyDetails[100][10]; // 100 upper limit is arbitrary - may need to expand this
								  // if it's less than DataPoints
	double ModelEstimates[100][2]; // Model estimates of outcome in control and intervention arms
	double ModelEffects[100];
	double AltLogL; // Only for the purpose of sensitivity analysis
	int MaxTerm;

	void CalcModelEst(int ii);
	void GetCurrModelEsts();
	void CalcModelEffects();
	void GetMaxTerm();
	double GetLikelihood();
};

// ---------------------------------------------------------------------------------------
// Classes for outputs
// ---------------------------------------------------------------------------------------

class PostOutputArray
{
	// We use this to record outputs for different years, from multiple simulations. 

public:
	PostOutputArray(int n);

	int columns;
	double out[samplesize][64]; // None of the arrays require > 64 columns. 
	double Means[64]; // Averages across samplesize
	double StdErrors[64]; // Standard errors of means
	double PeriodTot; // Sum of output across multiple years (mean)
	double PeriodTotSE; // SE of sum of output across multiple years 
	double PeriodIncrease; // Increase over a period, relative to baseline scenario
	double PeriodIncreaseSE; // SE of saving relative to baseline

	void RecordSample(const char* filout);
	void ReadBaseline(const char* input);
	void CalcMeans(); // Also calculates StdError
	void CalcPeriodTot(int Y1, int Y2);
	void OutputRatio(PostOutputArray* array1, PostOutputArray* array2);
	void CalcPeriodRatio(int Y1, int Y2, PostOutputArray* array1, PostOutputArray* array2);
	void CalcPeriodIncrease(int Y1, int Y2, PostOutputArray* array1);
	void GetSummGenOut(int col); 
	void GetSummHCTout1(int row); // Specific to Modelling Consortium HCT modelling project
	void GetSummHCTout2(int row, PostOutputArray* array1, PostOutputArray* array2); 
	void CalcSvyAssnLogL(double SvyMat[4][2], int year);
};

class PostOutputArray2
{
	// Same as PostOutputArray, but with dimension ParamCombs instead of samplesize. 

public:
	PostOutputArray2(int n);

	int columns;
	double out[ParamCombs][64]; // None of the arrays require > 64 columns. 

	void RecordSample(const char* filout);
};

class PostOutputArrayB
{
	// Similar to PostOutputArray, but for calibration to behavioural data (MSM). 

public:
	PostOutputArrayB(int n, int n2, int n3);

	int studies; // Adjust dimension of CalibData if n2>10
	int columns;
	int AgeStd; // 1 if calibration should be age-standardized, 0 otherwise
	double out[samplesize][81]; // None of the arrays require > 81 columns.
	double out2[samplesize][81][2]; // Outputs stratified by age (3rd index): <25 or 25+
	double CalibData[15][4]; // Calib data (columns for mean, SE, % aged 25+, year)
	double RandomEffectVar;
	double ModelAve[15];
	double LogL;

	void CalcRandomEffectVar();
	void CalcLogL();
	void RecordSample(const char* filout);
};

class Database
{
	// For recording outputs at an individual level

public:
	Database(int n);

	int TotIndivs; // total individuals in database
	int columns;
	double out[MaxPop][10]; // Increase to >10 columns if necessary

	void RecordSample(const char* filout);
};

// ---------------------------------------------------------------------------------------
// Classes for STD parameters and probabilities of transition between STD disease states
// ---------------------------------------------------------------------------------------

class STDtransition
{
public:
	STDtransition();

	// NB: The CondomEff and SuscepIncrease members actually belong to the
	// opposite sex. This is potentially confusing, but we have done things this way 
	// to be consistent with the deterministic model.

	int nStates;
	int SexInd; // 0 = male, 1 = female
	double AveDuration[6]; // Average number of weeks spent in state if untreated
	double CondomEff; // Probability that condom prevents transmission
	double CircEff; // Probability that male circumcision prevents transmission
	double SuscepIncrease[16]; // Multiple by which susceptibility to STD increases, by age
	double HIVinfecIncrease[6]; // % by which HIV infectiousness increases, by HIV/STD stage
	double RelTransmCSW; // Ratio of M->F transm prob per sex act in CSW-client relationships 
						 // to that in non-spousal relationships
	double RelTransmLT; // Ratio of transm prob in spousal relationships to that in
						// non-spousal relationships

	// Objects used for calibration purposes
	ANC ANClogL;
	FPC FPClogL;
	GUD GUDlogL;
	CSW CSWlogL;
	Household HouseholdLogL;
	AntenatalN AntenatalNlogL;
	HouseholdN HouseholdNlogL;
	double CSWprevUnsmoothed[41];
	double ANCprevUnsmoothed[81];

	void ReadPrevData(const char* input);
	void GetCSWprev();
	void GetANCprevSmooth();
	void SetVarStudyEffect(double Variance);
};

class HIVtransition: public STDtransition
{
public:
	HIVtransition(int Sex, int ObsANC, int ObsFPC, int ObsGUD, int ObsCSW,
		int ObsHH, int ObsANCN, int ObsHHN);

	double TransmProb[7]; // Transmission probability per sex act (by relationship type)
	double TransmProbMSM[7][2]; // Only relevant to MSM (2nd index: 0=receptive)
	double From1to2; // Prob of progressing from acute HIV to asymp per STD cycle
	double From2to3; // Prob of progressing from asymp to pre-AIDS symptoms, per STD cycle
	double From3to4; // Prob of progressing from pre-AIDS to untreated AIDS per STD cycle
	double From3to5; // Prob of progressing to AIDS & starting HAART per STD cycle
	double From5to6; // Prob of interrupting ART per STD cycle
	double From6to5; // Prob of resuming ART after an interruption
	double From4toDead; // Prob of dying from AIDS if not receiving HAART
	double From5toDead; // Prob of dying from AIDS while receiving HAART
	double From6toDead; // Prob of dying from AIDS while interrupting HAART

	// New variables for modelling CD4 and VL changes
	double MeanInitCD4; // Mean parameter for lognormal dbn of initial CD4 counts
	double SDinitCD4; // Std deviation parameter for lognormal dbn of initial CD4 counts
	double MeanInitVL; // Mean initial VL (after acute HIV) on log10 scale
	double SDinitVL; // Std deviation of initial VLs (after acute HIV) on log10 scale
	double AgeEffectInitVL; // Increase in log of init VL per year of increase in age
	double ChangeCD4; // Annual average change in CD4 on ln scale
	double VLeffectCD4; // Effect of logVL on change in annual CD4 (ln scale)
	double AgeEffectCD4; // Effect of age on change in annual CD4 (ln scale)
	double ChangeVL; // Annual average change in VL on log10 scale
	double VLeffectVL; // Effect of logVL on change in logVL
	double CD4effectVL; // Effect of 100-unit CD4 change on change in logVL
	double AgeEffectVL; // Effect of age on change in logVL
	double MortZeroCD4; // Annual AIDS mortality if CD4 count is zero
	double RednMortCD4; // Reduction in AIDS mortality per unit increase in CD4
	double ARTmortRedn; // RR HIV-related mortality per year increase in ART duration
	double ARTmortConstant; // Offset to annual reduction factor

	void CalcTransitionProbs();
	void GetPrev();
	double GetTransmProb(int ID);
	double GetTransmProbMSM(int ID);
	double GetCofactorType2(int ID, int ID2);
	double GetInfectivityMult(int ID);
	void GetNewStage(int ID, double p);
};

class NonHIVtransition: public STDtransition
{
public:
	NonHIVtransition();

	double CorrectRxPreSM; // % of cases correctly treated before syndromic management 
	double CorrectRxWithSM; // % of cases correctly treated under syndromic management
	double DrugEff; // Prob of complete cure if effective drugs are prescribed
	double TradnalEff; // Prob of cure if treated by a traditional healer
	double ProbCompleteCure; // Prob that individual seeking treatment is cured
	double HIVsuscepIncrease[6]; // % by which HIV susceptibility increases, by STD stage
	double PropnByStage[320][7]; // 20 behav states (max) x 16 age groups = 320 rows
								 // Max of 7 disease states. This matrix determines initial
								 // % of individuals in different STD stages.

	void CalcProbCure();
	int GetSTDstage(int offset, double r);
};

class TradnalSTDtransition: public NonHIVtransition
{
public:
	TradnalSTDtransition();

	int ImmuneState; // The index of the state representing individuals immune to reinfection
					 // (value is 0 if there is no immune state)
	double TransmProb; // Probability of transmission per act of sex
	double TransmProbSW; // Probability of transmission per act of sex (note that this is
						 // only used for client-to-sex worker transmission)
};

class SyphilisTransition: public TradnalSTDtransition
{
public:
	SyphilisTransition(int Sex, int ObsANC, int ObsFPC, int ObsGUD, int ObsCSW, int ObsHH,
		int ObsANCN, int ObsHHN);

	double ANCpropnScreened; // % of women attending ANCs who are screened for syphilis
	double ANCpropnTreated; // % of women testing positive who receive treatment
	double ANCpropnCured; // % of women attending ANCs who get cured
	double PropnSuscepAfterRx; // Proportion of individuals who are susceptible to reinfection
							   // after cure of primary syphilis

	double From1to2; // Prob that incubating syphilis becomes primary syphilis per STD cycle
	double From2to0; // Prob that primary syphilis is cured and RPR- per STD cycle
	double From2to0T; // Prob that primary syphilis is cured and RPR- per STD cycle, age <20
	double From2to0C; // Prob that primary syphilis is cured and RPR- per STD cycle, in CSWs
	double From2to3; // Prob that primary syphilis becomes secondary syphilis per STD cycle
	double From2to3T; // Prob that primary syphilis becomes secondary syphilis, age <20
	double From2to3C; // Prob that primary syphilis becomes secondary syphilis, in CSWs
	double From2to5; // Prob that primary syphilis is cured but still RPR+ per STD cycle
	double From2to5T; // Prob that primary syphilis is cured but still RPR+, age <20
	double From2to5C; // Prob that primary syphilis is cured but still RPR+, in CSWs
	double From3to4; // Prob that secondary syphilis becomes latent per STD cycle
	double From3to4T; // Prob that secondary syphilis becomes latent per STD cycle, age <20
	double From3to4C; // Prob that secondary syphilis becomes latent per STD cycle, in CSWs
	double From3to5; // Prob that secondary syphilis is cured per STD cycle
	double From3to5T; // Prob that secondary syphilis is cured per STD cycle, age <20
	double From3to5C; // Prob that secondary syphilis is cured per STD cycle, in CSWs
	double From4to6; // Prob that latent syphilis resolves per STD cycle
	double From4to6C; // Prob that latent syphilis resolves per STD cycle, in CSWs
	double From5to0; // Prob of seroreversion if cured in early syphilis, per STD cycle
	double From6to0; // Prob of seroreversion if resolved in latent syphilis, per STD cycle

	void CalcTransitionProbs();
	double GetTransmProb(int ID);
	void GetNewStage(int ID, double p);
};

class HerpesTransition: public TradnalSTDtransition
{
public:
	HerpesTransition(int Sex, int ObsANC, int ObsFPC, int ObsGUD, int ObsCSW, int ObsHH,
		int ObsANCN, int ObsHHN);

	double RecurrenceRate; // Rate at which symptomatic recurrences occur
	double SymptomaticPropn; // % of people who develop primary ulcer
	double From1to2; // Prob that primary ulcer resolves per STD cycle
	double From1to2T; // Prob that primary ulcer resolves per STD cycle, age <20
	double From1to2C; // Prob that primary ulcer resolves per STD cycle, in CSWs
	double From2to3[7]; // Prob of symptomatic recurrence per STD cycle, by HIV stage
	double From3to2; // Prob of resolution of recurrent ulcer per STD cycle
	double From3to2T; // Prob of resolution of recurrent ulcer per STD cycle, age <20
	double From3to2C; // Prob of resolution of recurrent ulcer per STD cycle, in CSWs
	double From2to4; // Prob that transiently asymp indiv becomes permanently asymp

	void CalcTransitionProbs();
	double GetTransmProb(int ID);
	void GetNewStage(int ID, double p);
};

class OtherSTDtransition: public TradnalSTDtransition
{
public:
	OtherSTDtransition(int Sex, int ObsANC, int ObsFPC, int ObsGUD, int ObsCSW, int ObsHH,
		int ObsANCN, int ObsHHN);

	double SymptomaticPropn; // % of individuals who develop symptoms
	double PropnImmuneAfterRx; // % of individuals who are immune to reinfection after 
							   // successful treatment
	double PropnImmuneAfterSR; // % of individuals who are immune to reinfection after 
							   // spontaneous resolution

	double From1to0; // Prob of symptomatic indiv reverting to susceptible per STD cycle
	double From1to0T; // Prob of symptomatic indiv reverting to susceptible per STD cycle, 
					  // age <20
	double From1to0C; // Prob of symptomatic indiv reverting to susceptible per STD cycle, 
					  // in CSWs
	double From1to3; // Prob of symptomatic indiv becoming immune per STD cycle
	double From1to3T; 
	double From1to3C; 
	double From2to0; // Prob of asymptomatic indiv reverting to susceptible per STD cycle
	double From2to0C; // Prob of asymptomatic indiv reverting to susceptible per STD cycle, 
					  // in CSWs
	double From2to3; // Prob of asymptomatic indiv becoming immune per STD cycle
	double From2to3C; 
	double From3to0; // Prob of immune individual reverting to susceptible per STD cycle

	void CalcTransitionProbs();
	double GetTransmProb(int ID, int STD); // STD 1 for NG, 2 for CT, 3 for TV and 4 for HD 
	void GetNewStage(int ID, double p, int STD); // Same STD codes as before
};

class NonSTDtransition: public NonHIVtransition
{
public:
	NonSTDtransition();

	double DrugPartialEff; // Prob of partial cure if effective drugs are prescribed
	double ProbPartialCure; // Prob that individual seeking treatment is partially cured

	void CalcProbPartialCure();
};

class BVtransition: public NonSTDtransition
{
public:
	BVtransition(int Sex, int ObsANC, int ObsFPC, int ObsGUD, int ObsCSW, int ObsHH, int ObsANCN,
		int ObsHHN);

	double SymptomaticPropn; // Propn of women developing BV who experience symptoms
	double Incidence1; // Weekly incidence of BV in women with intermediate flora, with
					   // one current partner
	double IncidenceMultTwoPartners; // Multiple by which incidence is increased in women
									 // with more than one partner
	double IncidenceMultNoPartners; // Multiple by which incidence is decreased in women
									// with no partners
	double CtsTransition[4][4]; // Continuous rates of transition between states
	double From1to2; // Prob of going from normal vaginal flora to intermediate per STD cycle
	double From2to1ind; // Independent prob of reverting to normal flora per STD cycle
	double From3to1; // Prob of returning to normal flora if BV is currently symptomatic
	double From3to1T; // Prob of returning to normal flora, age <20
	double From3to1C; // Prob of returning to normal flora, in CSWs
	double From3to2; // Prob of returning to intermediate flora if BV is currently symptomatic
	double From3to2T; // Prob of returning to intermediate flora, age <20
	double From3to2C; // Prob of returning to intermediate flora, in CSWs
	double From4to1; // Prob of returning to normal flora if BV is currently asymp
	double From4to1C; // Prob of returning to normal flora if BV is currently asymp, in CSWs
	double From4to2; // Prob of returning to intermediate flora if BV is currently asymp
	double From4to2C; // Prob of returning to intermediate flora if BV is asymp, in CSWs
	double From2to3ind[3]; // Independent prob of developing symptomatic BV per STD cycle, if 
						   // currently intermediate, by # current partners
	double From2to4ind[3]; // Independent prob of developing asymp BV per STD cycle, if 
						   // currently intermediate, by # current partners
	double From2to1dep[3]; // Dependent prob of reverting to normal flora per STD cycle, if 
						   // currently intermediate, by # current partners
	double From2to3dep[3]; // Dependent prob of developing symptomatic BV per STD cycle, if 
						   // currently intermediate, by # current partners
	double From2to4dep[3]; // Dependent prob of developing asymp BV per STD cycle, if 
						   // currently intermediate, by # current partners

	void CalcTransitionProbs();
	void GetNewStage(int ID, double p);
};

class VCtransition: public NonSTDtransition
{
public:
	VCtransition(int Sex, int ObsANC, int ObsFPC, int ObsGUD, int ObsCSW, int ObsHH, int ObsANCN,
		int ObsHHN);

	double RecurrenceRate; // Rate at which symptoms develop in women with yeast colonization
	double Incidence; // Rate at which asymptomatic yeast colonization occurs
	double IncidenceIncrease[6]; // % increase in incidence of yeasts by HIV stage
	double From1to2; // Prob that asymp colonization becomes symptomatic per STD cycle
	double From1to2C; // Prob that asymp colonization becomes symptomatic, in CSWs
	double From1to0; // Prob that asymp colonization resolves per STD cycle
	double From1to0C; // Prob that asymp colonization resolves per STD cycle, in CSWs
	double From2to1; // Prob that symptomatic infection becomes asymp per STD cycle
	double From2to1T; // Prob that symptomatic infection becomes asymp per STD cycle, age <20
	double From2to1C; // Prob that symptomatic infection becomes asymp per STD cycle, in CSWs
	double From2to0; // Prob that symptomatic infection is cured per STD cycle
	double From2to0T; // Prob that symptomatic infection is cured per STD cycle, age <20
	double From2to0C; // Prob that symptomatic infection is cured per STD cycle, in CSWs
	double From0to1[7][7]; // Prob that uninfected woman gets asymptomatically colonized per
						   // STD cycle (by age and HIV stage)
	
	void CalcTransitionProbs();
	void GetNewStage(int ID, double p);
};

class PaedHIV
{
public:
	PaedHIV();

	double PreAIDSmedian;
	double PreAIDSshape;
	double MeanAIDSsurvival;
	double AveYrsOnART; // Average survival after ART initiation (ignoring non-HIV mortality)
	double NonAIDSmort[15];
	double MortProb1st6m; // Independent prob of non-AIDS mortality in 1st 6 months of life

	void GetNewStage(int ID, double p);
};

class Child
{
public:
	Child(int Sex);

	int SexInd; // 0 = male, 1 = female
	double PropnBirths;

	PaedHIV Perinatal;
	PaedHIV Breastmilk;
};

// --------------------------------------------------------------------------------------
// Classes for individuals
// --------------------------------------------------------------------------------------

class Indiv
{
public:
	Indiv();

	int AliveInd; // 0 = dead, 1 = alive
	int SexInd; // 0 = male, 1 = female (same as in TSHISA)
	int RiskGroup; // 1 = high risk, 2 = low risk (same as in TSHISA)
	double PopWeight; // Scaling factor to get SA totals, in current year
	double DOB; // Date of birth (calendar year, with decimal indicating timing of birth)
	double DOD; // Date of death
	double DOLC; // Date of last conception
	double DOLB; // Date of last birth
	double DOLW; // Date of last weaning
	double DOAI; // Date of ART initiation
	double DOUD; // Date of union dissolution (due to divorce/widowhood)
	double Fecundability; // Fertility adjustment factor
	double Conscientiousness; // z-score for personality
	int AgeGroup; // 5-year age group (0 corresponds to 0-4)
	int CurrAge; // Age at start of year
	int Race; // 0 = black African, 1 = coloured, 2 = other (white/Asian)
	int HouseholdID;
	double DateHomeless; // Date of becoming homeless, if HouseholdID = 0
	int ParentID[2]; // IDs of father and mother
	int ChildIDs[20]; // IDs of children
	int FatherBaby; // ID of father (only relevant to women currently pregnant)
	int CircInd; // 0 = uncircumcised, 1 = circumcised (males only)
	int HighestGrade; // Highest grade completed (13 = tertiary qualification)
	int InSchool; // 1 => currently enrolled at an educational institution (excl. PG)
	int CurrContr; // Current contraception (0=none, 1=injectable, 2=pill, 3=sterilization)
	int PrevContr; // Previous contraception (0=never, 1=injectable, 2=pill, 3=sterilization)
	int CurrUrban; // 0 = rural, 1 = urban
	int Imprisoned; // 0 = not in prison, 1 = awaiting trial in prison, 2 = sentenced in prison
	double ReleaseDate; // Date of release from prison
	int Employed; // 0 = unemployed, 1 = employed
	double LogIncomeDif; // Difference between log income and average for age/sex/race/edu
	int PrivatePension; // 1 if person is receiving a private pension, 0 otherwise
	int MarriedInd; // 0 = unmarried, 1 = married (same as in TSHISA)
	int FSWind; // 1 => Female sex worker (same as in TSHISA)
	int CasualInd; // 1 = currently available for casual sex (only relevant to MSM)
	int HetCasualInd; // 1 = available for casual sex with people of the opposite sex
	int VirginInd; // 0 = sexually experienced, 1 = virgin (same as in TSHISA)
	int IDprimary; // ID of primary partner (0 => no primary partner)
	int ID2ndary; // ID of secondary partner (0 => no secondary partner)
	// Note that the individual's ID is always 1 + their row number in the array (so that the
	// first individual simulated has ID number 1, not 0).
	int CurrPartners; // Number of current sexual partners (max 2)
	int LifetimePartners; // Cumulative number of sexual partners (ever)
	int AnnPartners; // Number of sexual partners in the last year
	double IneqGender; // Prob of endorsing inequitable gender norms (applies to males only)
	double NonHIVmortProb; // Prob of dying from causes unrelated to HIV, per behaviour cycle
	double NonHIVfertRate; // Rate of giving birth in HIV-negative women, per behaviour cycle
	double CumSelectionProb; // Prob that this indiv (or someone with lower ID) is selected
	double MalePref; // Desired fraction of partners who are men (only relevant to men)
	double ChangeMalePref; // Annual change in fraction of partners who are male, for bisexuals
	double InsertivePref; // 0 = exclusively receptive, 1 = exclusively insertive, 0.5 = no pref
	double CD4; // Current CD4 count (only relevant for HIV-positive individuals)
	double BaselineCD4; // CD4 at first ART initiation
	double logVL; // Current log10 viral load (only relevant for HIV-positive individuals)
	double SPVL; // Set point viral load (log10 scale)
	int HIVstage; // 0 for HIV-neg, 1=acute, 2=CD4 350+, 3=CD4 200-349, 4=CD4<200, 5=ART, 6=interrupting
	int HIVstageE; // HIV stage at end of STD cycle
	int CTstage; // 0 for susceptible, 1 = symptomatic, 2 = asymptomatic infected, 3 = immune
	int CTstageE; // CT stage at end of STD cycle
	int HDstage; // 0 for susceptible, 1 = symptomatic, 2 = asymptomatic infected, 3 = immune
	int HDstageE; // HD stage at end of STD cycle
	int NGstage; // 0 for susceptible, 1 = symptomatic, 2 = asymptomatic infected, 3 = immune
	int NGstageE; // NG stage at end of STD cycle
	int TVstage; // 0 for susceptible, 1 = symptomatic, 2 = asymptomatic infected, 3 = immune
	int TVstageE; // TV stage at end of STD cycle
	int TPstage; // 0 for susceptible, 1 = incubation, 2 = primary, 3 = 2ndary, 4 = latent, 5/6 = immune
	int TPstageE; // TP stage at end of STD cycle
	int HSVstage; // 0 for susceptible, 1 = primary ulcer, 2 = transiently asymptomatic,
				  // 3 = recurrent ulcer, 4 = long-term asymptomatic
	int HSVstageE; // HSV stage at end of STD cycle
	int BVstage; // 0 = normal, 1 = intermediate, 2 = symptomatic BV, 3 = asymptomatic BV
	int BVstageE; // BV stage at end of cycle
	int VCstage; // 0 = uninfected, 1 = asymptomatic VC, 2 = symptomatic VC
	int VCstageE; // VC stage at end of cycle
	double SuscepHIVadj; // Adjustment to individual's HIV susceptibility
	double DateInfect; // Date at which HIV was acquired
	double PartnerRateAdj; // adjustment to individual's rate of partnership formation and
						   // dissolution (allowing for variation between individuals)
	double DesiredNewPartners; // rate at which new partners are acquired
	double CondomPref; // adjustment to individual's odds of condom use based on edu, race
	double TempCondomAdj; // adjustment to reflect current level of binge drinking, income
						  // and (in men) inequitable gender norms
	double CondomPrimary; // fraction of sex acts with primary partner that are protected
	double Condom2ndary; // fraction of sex acts with 2ndary partner that are protected
	int DisclosedPrimary; // 1 if disclosed HIV+ status to primary partner, 0 otherwise
	int Disclosed2ndary; // 1 if disclosed HIV+ status to 2ndary partner, 0 otherwise
	int NewStatus; // 0 if indiv has not yet had status updated in current behav cycle
	int UVIprimary; // # unprotected acts of vaginal/anal intercourse with primary partner
	int PVIprimary; // # protected acts of vaginal/anal intercourse with primary partner
	int UVI2ndary; // # unprotected acts of vaginal/anal intercourse with 2ndary partner
	int PVI2ndary; // # protected acts of vaginal/anal intercourse with 2ndary partner
	int UVICSW; // # unprotected acts of vaginal intercourse with CSWs
	int PVICSW; // # protected acts of vaginal intercourse with CSWs
	int UAIcasual; // # unprotected acts of vaginal/anal intercoure with casual partners
	int PAIcasual; // # protected acts of vaginal/anal intercoure with casual partners 
	int IDofCSW; // ID of CSW (if man visits CSW)
	int IDofCasual; // ID of casual sex partner 
	int VCThistory; // 0 = never HIV tested, 1 = last tested neg, 2 = diagnosed positive
	int Visiting; // 1 if currently away from usual place of residence to visit partner
	double VisitFreq; // Annual frequency of partner visits if geographically separated
	double OnPrEP; // Level of PrEP use (0 = none, 1 = perfect adherence)
	double VaccineEff; // Reduction in HIV transmission risk due to vaccination
	double Date1stVacc; // Date of first HIV vaccination
	double DailyDrinkProb; // Daily prob of consuming any alcohol
	double DrinksPerDD; // Average drinks per drinking day
	double DrinkProbConstant; // DailyDrinkProb component that is constant over time
	double DrinksPerDDconstant; // DrinksPerDD component that is constant over time

	// MSM output variables
	int EverMSM; // 1 for men who have ever had sex with other men, 0 otherwise
	int EverBi; // 1 for bisexual men who have ever had sex with women
	int RecentMSM; // # MSM partners in last 6 months
	int RecentBi; // # female partners in last 6 months (if bisexual)
	// Other output variables (only needed for calibration)
	int EverInjectable; // temporary, only necessary for calibration
	int EverPill; // temporary, only necessary for calibration
	int PrevImprisoned; // 0 = never previously in prison, 1 otherwise
	double BirthInterval; // temporary, only necessary for calibration
	double DateMig; // Date of last migration
	int EverCasual; // 1 if ever had once-off casual sex
	int CasualLast3mo; // # once-off casual sex partners in last 3 months
	int PartnersLast3mo; // # ST and LT partners in last 3 months (excl. casual, commercial sex)

	// Structural RCT calibration variables
	int BaselineAge;
	int BaselinePoor;
	int BaselineSchool;
	int BaselineBinge;
	int BaselineUnemployed;
	int BaselineVirgin;
	int CumCasual;
	int CumSTIs;
	int CumHSV2;
	int CumTeenPreg;

	int SelectEvent(int ID, double rnd, int sexpref); // Output is 0 if no new event, 1 if acquire S1,  
								 // 2 if acquire S2, 3 if marry primary, 4 if marry 2ndary, 5 if    
								 // divorce, 6 if end primary S, 7 if end 2ndary S
	double ConvertToDependent1(double rate1);
	double ConvertToDependent2(double rate1, double rate2, int ord);
	double ConvertToDependent3(double rate1, double rate2, double rate3, int ord);
	double ConvertToDependent4(double rate1, double rate2, double rate3, double rate4, int ord);
	void AssignCondom(int RelType, int ID, int PID);
	void SimulateSexActs(int ID);
	int RandomSexGenerator(double p, double lambda);
	void GetInitCD4andVL();
	void GetNewHIVstate(int ID, double p);
	void GetNewHSVstate(int ID, double p);
	void GetNewTPstate(int ID, double p);
	void GetNewHDstate(int ID, double p);
	void GetNewNGstate(int ID, double p);
	void GetNewCTstate(int ID, double p);
	void GetNewTVstate(int ID, double p);
	void TestStartContr(int EventType, int condom);
	void GetDiagnosed(int ID, int Setting);
	void GetRediagnosed(int ID, int Setting); // For individuals who were previously diagnosed positive
	void TestDisclose(int ID, int Partner, int Timing);
	void TestCondomChange(int ID);
	void UpdateChildren(int ChildID);
	void FindMother(int ID);
	void FindFather(int ID);
	void AssignChildToHH(int ID);
	void ChangeHHafterMigration(int ID);
	void ChangeHHleavingNest(int ID);
	void CheckReturnToNest(int ID);
	int FindRelativeToLiveWith(int ID, int Sibling); // Returns 0 if no relative could be found, 1 otherwise
	int TestHomeless(int ID); // Returns 1 if individual becomes homeless, 0 otherwise
};

class HouseholdGroup
{
public:
	HouseholdGroup(int head);

	int IDhead; // ID of the household head
	int Members[MaxHHsize]; // IDs of all household members
	int Size; // Number of household members
	int Kids; // Number of children under age 15
	int Active; // 0 if the household no longer exists, 1 otherwise
	int Urban; // 0 for a rural household, 1 for an urban household
	double PerCapitaIncome; // Monthly household income, divided by # household members
	double PerCapitaIncomeAdj; // Using # HH members adjusted by equivalence scale

	void AddMember(int ID);
	void RemoveMember(int ID); // Always call AddMember/NewHousehold BEFORE RemoveMember.
	void GetPerCapitaIncome(int hhid);
	void GetUnsortedIncomes();
};

class Pop
{
public:
	Pop(int xxx);

	int PopPyramid[18][2];
	int BehavPyramid[18][6]; // virgin M, experienced unmarried M, married M (then F)
	int PrefMatrix[16][16]; // 1st index is age of male, 2nd index is age of female
	double NumberPartners[16][10]; // unmarried 0, unmarried 1, unmarried 2, married 1,
								// married 2 (males then females)
	double FSWcontacts2016[16]; // # male contacts with FSWs in 2016, by age
	int TotSex[16][8]; // ST UVI, ST PVI, LT UVI, LT PVI (males then females)
	double AdultHIVstageTrend[81][14]; // 2nd index is HIV stage in M (0-6) and F (7-13)
	int CSWregister[MaxCSWs][3]; // Records IDs of women who are currently sex workers
	int CasualRegister[MaxCasual]; // Records IDs of MSM currently engaging in casual sex
	int CasualRegisterM[MaxCasual][3]; // IDs of men engaging in casual sex with women, by race
	int CasualRegisterF[MaxCasual][3]; // IDs of women engaging in casual sex with men, by race
	double EligibleByAgeM[16][2], EligibleByAgeF[16][2]; // Sum of DesiredNewPartners
		// in individuals eligible to acquire new partners, by age and risk group

	// Functions to assign initial profile
	void AssignAgeSex();
	void AssignPersonality();
	void GetAllPartnerRates();
	void AssignBehav();
	double GetMarriedPropnAtStart(int Sex, int Age, int Race);
	void AssignVisitFreq();
	void AssignEdu();
	void AssignBaseIncome();
	void AssignUrban();
	void AssignParents();
	void AssignBirth();
	void AssignIneqGender();
	void AssignAlcohol();
	void AssignInitHousehold();
	void AssignHeadship();
	void AssignAdultsHH();
	void AssignChildrenHH();
	void AssignContraception();
	void AssignPrison();
	void SetEmployment();
	void ChooseLTpartner(int ID, double rnd1, double rnd2);
	void ChooseSTpartner(int ID, double rnd1, double rnd2);
	void AssignHIV();
	void AssignHIV1990(int race);
	void AssignSTIs();
	void ResetToBaselineStart();

	// Output functions
	void GetPopPyramid();
	void GetBehavPyramid();
	void GetSexuallyExp();
	void GetAgeDisparate();
	void GetPrefMatrix(int type); // type = 1 for married, 2 for unmarried
	void GetNumberPartners();
	void GetTotSex();
	void GetCondomUse(); // For effects of edu, race, diagnosis
	void GetCondomUse2(); // For trends in condom use
	void GetPrEP(); // For trends in condom use
	void GetNumbersByHIVstage();
	void GetInitHIVprevH();
	void GetHHprevByAge();
	void GetANCprevByAge();
	void GetHIVprevByEdu();
	void GetCD4profile();
	void GetInitSTIprev();
	void GetCurrSTIprev();
	double GetQstatistic(int Cmatrix[3][3], int MatDim);
	void GetSTIconcordance();
	void SavePopPyramid(const char* filout);
	void SaveBehavPyramid(const char* filout);
	void SavePrefMatrix(const char* filout);
	void SaveNumberPartners(const char* filout);
	void SaveFSWcontacts2016(const char* filout);
	void SaveLifetimePartners(const char* filout);
	void SaveTotSex(const char* filout);
	void SaveMarriageRates(const char* filout);
	void SaveAdultHIVstage(const char* filout);
	void GetHIVsurvivalTimes();
	void GetVLs();
	void GetUnprotected();
	void GetMigHIV();
	void GetPrisons();
	void GetEduProfile();
	void GetAssortativeness();
	void GetCohabitTemp();
	void GetConcurrency();
	void GetMarriage();
	void GetMigrationOutput();
	void GetAlcoholOutput();
	void GetGenderNormOutput();
	void GetHouseholdOutput();
	void GetDiscordance();
	void GetDiscordance2();
	void GetCasualCalib();

	// Functions updated every year
	void OneYear();
	void ResetFlow();
	void ResetAnnPartners();
	void UpdateAgeGroup();
	void UpdateMaleCirc();
	void UpdateUrban();
	void ExternalSTIintroductions();
	void UpdateEduProfile();
	void UpdateIneqGender();
	void UpdateAlcohol();
	void UpdateEmployment();
	void UpdateHHincome();
	void UpdateHouseholds();
	void UpdatePopWeights();
	void CalcNonHIVmort();
	void UpdateFecundability();
	void CalcNonHIVfert();
	void CalcContrOutput();
	void UpdateProbCure();
	void UpdateSTDtransitionProbs();
	void UpdateCondomUse();
	void UpdateCondomEffects();
	void CalcInitFSWageDbn(int iy);
	void CalcCurrMarriageRates();
	void GetANCprev(STDtransition* a, int STDind); // 1 = HSV, 2 = TP, 3 = HD, 4 = NG, 5 = CT, 6 = TV, 7 = BV, 8 = VC
	void GetHHprev(STDtransition* a, int STDind);
	void GetFPCprev(STDtransition* a, int STDind);
	void GetCSWprev(STDtransition* a, int STDind);
	void CalcHIVandSTIincidence();
	void CalcStockCosts();
	void CalcFlowCosts();
	void ResetStructOutcomes();

	// Functions updated on a monthly basis
	void UpdateConception(int month);
	void UpdateContraception(int month);
	void UpdateDisclosure();
	void UpdateIncarceration(int month);
	void UpdatePrEP(int month);
	void UpdateVaccination(int month);
	void StructInterventionWaning(int month);

	// Functions updated every sexual behaviour cycle
	void OneBehavCycle();
	void UpdateBirths();
	void UpdateNonHIVmort();
	void UpdateFSWprofile();
	void UpdateMarriageIncidence();
	void UpdateVisitStatus();
	void BalanceSexualPartners();
	void CalcPartnerTransitions();
	void UpdateMSMcasual();
	void UpdateHetCasual();
	void GetNewPartner(int ID, double rnd, int rsk);
	void GetNewMSMpartner(int ID, double rnd, int rsk);
	void SetNewStatusTo1(int ID);
	void SetToDead(int ID, int Cause); // Cause is 1 for HIV, 0 otherwise
	void NewBirth(int ID);
	void NewHousehold(int ID);
	void ChangeHHafterMarriage(int ID1, int ID2);
	void ChangeHHafterDivorce(int ID1, int ID2);
	void GetStructRCTests();

	// Functions updated every STD cycle
	void OneSTDcycle();
	void GetSexActs();
	void CalcHIVtesting();
	void GetHIVtransitions();
	void GetSTDtransitions();
};

class Partner
{
public:
	Partner();

	int ID;
	double DesiredNewPartners; // rate at which new partners are acquired
};

class PartnerCohort
{
	// This class defines a collection of individuals of the same age group, risk group 
	// and sex, who are eligible to form new partnerships.

public:
	PartnerCohort();
	
	vector<Partner> Pool;
	double TotalDesire;

	void Reset(); // Note that this function should be redundant in the new model.
	void AddMember(int PID, int Pref); // Pref = 0 for males, 1 for females
	int SamplePool(double rand1); // Returns the ID of a randomly chosen partner
	int SamplePool2(int ID1); // Returns the ID of a randomly chosen partner
};

// --------------------------------------------------------------------------------------
// General functions
// --------------------------------------------------------------------------------------

// Functions for reading input parameters

void ReadSexAssumps(const char* input);
void ReadSTDepi(const char* input);
void ReadConception();
void ReadBreastfeeding();
void ReadMigration();
void ReadDisclosure();
void ReadGenderNormAssumps();
void ReadAlcoholAssumps();
void ReadHouseholdAssumps();
void ReadPrisonAssumps();
void ReadEmployment();
void ReadIncome();
void ReadMSM();
void ReadMSMbehavData();
void ReadMSM_HIVdata();
void ReadRatesByYear();
void ReadMortTables();
void ReadFertTables();
void ReadPopWeights();
void ReadOneStartProfileM(ifstream* file, int group);
void ReadOneStartProfileF(ifstream* file, int group);
void ReadStartProfile(const char* input);
void ReadInitHIV();
void ReadStartPop();
void ReadSTDprev();
void ReadSTDparameters();
void ReadMCprev();
void ReadEduAssumps();
void ReadStartEdu();
void ReadCostAssumps();
void ReadHCTparameters();
void ReadStructuralDriverData();
void ReadStructuralRCTdata();
void ReadAllInputFiles();

// Mathematical functions

//double GetQstatistic(int Cmatrix[3][3], int MatDim);

// Calibration/likelihood functions

void SetCalibParameters();
void CalcStructuralLogL();
void CalcTotalLogL();

// Functions for running simulations and uncertainty analysis

void OneSimulation();
void RunSimulations(int sims); // sims = # completed simulations
void StoreOutputs();
void StoreTempOutputs(); // Backing up, only relevant when FixedUncertainty=0
void RestartFromStore(int sims); // Restoring after a shutdown
void StoreHCToutputs(const char* filout); // For Modelling Consortium HCT project
void StoreGenOutputs(const char* filout);
void StoreStructuralOutputs(const char* filout);
void AggregateSims();
void SimulateParameters();
void SimulateHIVparams();
void SimulateTPparameters();
void SimulateHSVparameters();
void SimulateNGparameters();
void SimulateCTparameters();
void SimulateTVparameters();
void SimulateBVparameters();
void SimulateVCparameters();
void SimulateHCTparameters();
void SimulateStructParams();
void SaveBaselinePop();
void RestartFromBaseline();

// New MSM functions

void UpdateMalePref();
void ResetRecentMSM();
void GetMSMcalib();
void GetMSMcalib2();
void SimulateMSMparameters();
void SimulateMSM_HIV();

// New functions for structural drivers/interventions

void CalcStructuralAssns();
void RunStructRCTs();

// --------------------------------------------------------------------------------------
// Objects created from classes defined above
// --------------------------------------------------------------------------------------

// STI prevalence and incidence
PostOutputArray BVprev15to49F(56);
PostOutputArray CTprev15to49F(56);
PostOutputArray CTprev15to49M(56);
PostOutputArray CTprevCSW(56);
PostOutputArray HDprev15to49F(21);
PostOutputArray HDprev15to49M(21);
PostOutputArray HIVprev15to49F(56);
PostOutputArray HIVprev15to49M(56);
PostOutputArray HIVprev15to49(56);
PostOutputArray HIVprevCSW(56);
PostOutputArray HSVprev15to49F(56);
PostOutputArray HSVprev15to49M(56);
PostOutputArray HSVprevCSW(56);
PostOutputArray HSVprevANC(56);
PostOutputArray NGprev15to49F(56);
PostOutputArray NGprev15to49M(56);
PostOutputArray NGprevCSW(56);
PostOutputArray TPprev15to49F(56);
PostOutputArray TPprev15to49M(56);
PostOutputArray TPprevCSW(56);
PostOutputArray TPprevANC(56);
PostOutputArray TVprev15to49F(56);
PostOutputArray TVprev15to49M(56);
PostOutputArray TVprevCSW(56);
PostOutputArray VCprev15to49F(56);
PostOutputArray NewCT(56);
PostOutputArray NewHIV(56);
PostOutputArray NewHIV_U(56);
PostOutputArray NewHIV_R(56);
PostOutputArray NewHIVexp(56);
PostOutputArray NewHIVexpF(56);
PostOutputArray NewHIVexpM(56);
PostOutputArray HIVinc15to49(56);
PostOutputArray HIVinc15to49_U(56);
PostOutputArray HIVinc15to49_R(56);
PostOutputArray NewHSV(56);
PostOutputArray NewNG(56);
PostOutputArray NewTP(56);
PostOutputArray NewTV(56);
PostOutputArray NewCTNGTV(56);
PostOutputArray NewCTNGTV_F(56);
PostOutputArray NewCTNGTV_M(56);

PostOutputArray PrevPreg15(28);
PostOutputArray PrevPreg20(28);
PostOutputArray PrevPreg25(28);
PostOutputArray PrevPreg30(28);
PostOutputArray PrevPreg35(28);
PostOutputArray PrevPreg40(28);
PostOutputArray PrevPregTotal(28);
PostOutputArray PrevFertile15(30);
PostOutputArray PrevFertile20(30);
PostOutputArray PrevFertile25(30);
PostOutputArray PrevFertile30(30);
PostOutputArray PrevFertile35(30);
PostOutputArray PrevPregNoEdu(26);
PostOutputArray PrevPregPrimary(26);
PostOutputArray PrevPreg2ndary(26);
PostOutputArray PrevPregTertiary(26);
PostOutputArray PrevPregAllEdu(26);
PostOutputArray PrevPregUrban(28);
PostOutputArray PrevPregRural(28);
PostOutputArray PrevPregB(28);
PostOutputArray PrevPregC(28);
PostOutputArray PrevPregW(28);
PostOutputArray PrevHH2005(18);
PostOutputArray PrevHH2008(18);
PostOutputArray PrevHH2012(18);
PostOutputArray PrevByEdu2002(5);
PostOutputArray DiscordantPropnM(31);
PostOutputArray DiscordantPropnF(31);
PostOutputArray DiscordantPropn(31);
PostOutputArray DiscordantFtoM(31);
PostOutputArray DiscordantMposLT(36);
PostOutputArray DiscordantMposST(36);
PostOutputArray DiscordantFposLT(36);
PostOutputArray DiscordantFposST(36);
PostOutputArray ConcordantPosLT(36);
PostOutputArray ConcordantPosST(36);
PostOutputArray ConcordantNegLT(36);
PostOutputArray ConcordantNegST(36);

// STI prevalence by race, urban/rural
PostOutputArray HIVprev15to49B(51);
PostOutputArray HIVprev15to49C(51);
PostOutputArray HIVprev15to49W(51);
PostOutputArray HSVprev15to49B(51);
PostOutputArray HSVprev15to49C(51);
PostOutputArray HSVprev15to49W(51);
PostOutputArray UrbanHIV(51);
PostOutputArray RuralHIV(51);

/*PostOutputArray CTconcordance(21);
PostOutputArray HDconcordance(21);
PostOutputArray HIVconcordance(21);
PostOutputArray HSVconcordance(21);
PostOutputArray NGconcordance(21);
PostOutputArray TPconcordance(21);
PostOutputArray TVconcordance(21);*/

// MSM outputs/outputs by sexual orientation
PostOutputArrayB MSMmarried(40, 8, 1);
PostOutputArrayB MSMrecentBi(40, 6, 1);
PostOutputArrayB MSMeverBi(40, 12, 1);
PostOutputArrayB MSMaged25plus(40, 12, 0);
PostOutputArrayB MSMpositive(40, 13, 1);
PostOutputArrayB MSMcurrentReg(40, 3, 1);
PostOutputArrayB MSMmultPartners(40, 5, 1);
PostOutputArray EverMSMpropn(30);
PostOutputArrayB EverMSMcurrentReg(40, 1, 1);
PostOutputArray EverMSMeverBi(30);
PostOutputArrayB EverMSMpositive(40, 1, 1);
PostOutputArray NeverMSMpositive(30);
PostOutputArray HIVageProfileMSM(36);
PostOutputArray CasualAgeProfile(9);
PostOutputArray BiAgeProfile(9);
PostOutputArray HighRiskMSM(2);
PostOutputArray NewHIV_MSM(30);
PostOutputArray MSMpositiveBi(30);
PostOutputArray MSMpositiveGay(30);
PostOutputArray MSMposInsertive(30);
PostOutputArray MSMposReceptive(30);
PostOutputArray MSMposVersatile(30);
PostOutputArray HomoPrefPos(30);
PostOutputArray BiPrefPos(30);
PostOutputArray HeteroPrefPos(30);
PostOutputArray MSMprev(56); // Same as MSMpositive

// HIV disease progression outputs
PostOutputArray CD4500plus(30);
PostOutputArray CD4350to499(30);
PostOutputArray CD4200to349(30);
PostOutputArray CD4under200(30);
PostOutputArray UndiagnosedAdults(56);
PostOutputArray DiagnosedUntreated(56);
PostOutputArray TreatedAdults(56);
PostOutputArray AdultARTcoverage(56);
PostOutputArray AveHIVdurPreg(56);
PostOutputArray AveHIVdurPreg15(56);
PostOutputArray AveHIVdurPreg20(56);
PostOutputArray AveHIVdurPreg25(56);
PostOutputArray AveHIVdurPreg30(56);
PostOutputArray AveHIVdurPreg35(56);

// Education outputs
PostOutputArray CompletedMatric(56);
PostOutputArray CompletedMatricAM(56);
PostOutputArray CompletedMatricAF(56);
PostOutputArray CompletedMatricCM(56);
PostOutputArray CompletedMatricCF(56);
PostOutputArray CompletedMatricWM(56);
PostOutputArray CompletedMatricWF(56);
PostOutputArray TertiaryEnrolRatioA(56);
PostOutputArray TertiaryEnrolRatioC(56);
PostOutputArray TertiaryEnrolRatioW(56);
/*PostOutputArray AttainmentByAge2007(16);
PostOutputArray AttainmentByGrade2007(12);
PostOutputArray Enrolled2001A(14);
PostOutputArray Enrolled2001C(14);
PostOutputArray Enrolled2001W(14);*/
PostOutputArray EduByParentEdu(5);
PostOutputArray DropoutDuePregnancy(56);

// Employment outputs
PostOutputArray EmployedPropn(56);
PostOutputArray EmployedPropnAM(56);
PostOutputArray EmployedPropnAF(56);
PostOutputArray EmployedPropnCM(56);
PostOutputArray EmployedPropnCF(56);
PostOutputArray EmployedPropnWM(56);
PostOutputArray EmployedPropnWF(56);

// Income outputs
PostOutputArray MedianIncome(56);
PostOutputArray ZeroIncomeHH(56);
PostOutputArray GiniIncome(56);
PostOutputArray GiniIncomeB(56);
PostOutputArray GiniIncomeC(56);
PostOutputArray GiniIncomeW(56);
PostOutputArray PalmaRatio(56);
PostOutputArray MeanHHincomePC(56);

// Sexual behaviour outputs
PostOutputArray SexuallyExperienced(22);
PostOutputArray SexuallyExpAdol(56);
PostOutputArray EduConcordance(25);
PostOutputArray RaceConcordance(9);
PostOutputArray ConcurrencyPrev(33);
PostOutputArray ConcurrencyTot(56);
PostOutputArray MarriedTrendM(56);
PostOutputArray MarriedTrendF(56);
PostOutputArray MarriedYouth(56);
PostOutputArray Married2001M(30);
PostOutputArray Married2001F(30);
PostOutputArray Married2011M(30);
PostOutputArray Married2011F(30);
PostOutputArray MarriedEdu2011(10);
PostOutputArray CondomUseByEdu(30);
PostOutputArray CondomUseByRace(3);
PostOutputArray CondomUse15to24(56);
PostOutputArray CondomUse25to49(56);
PostOutputArray CondomUse15to49(56);
PostOutputArray SWsexActs(56);
PostOutputArray SWsexActsProt(56);
PostOutputArray CohabitTemp(6);
PostOutputArray CohabitTempTrend(51);
PostOutputArray PartnerAgeDifM(22);
PostOutputArray PartnerAgeDifF(22);
PostOutputArray AgeDisparate2000(2);
PostOutputArray CasualCalib2005(9);
PostOutputArray CasualPrevM(56);
PostOutputArray CasualPrevF(56);
PostOutputArray MultPartnersM(51);
PostOutputArray MultPartnersF(51);

// Intervention outputs
PostOutputArray TakingPrEP(56);
PostOutputArray VaccineDoses(56);

// Demographic, migration, contraception outputs
PostOutputArray TotBirths(56);
PostOutputArray TeenBirths(56);
PostOutputArray TotPop(56);
PostOutputArray TotPop15to49_U(56);
PostOutputArray TotPop15to49_R(56);
PostOutputArray UrbanTrend(30);
PostOutputArray UrbanAge(36);
PostOutputArray UrbanAge1996(36);
PostOutputArray UrbanRace(3);
PostOutputArray UrbanRace1996(3);
PostOutputArray UrbanEdu(6);
PostOutputArray ContrPrev(36);
PostOutputArray EverUseContr98(8);
PostOutputArray EverUseContr03(8);
PostOutputArray CurrUseContr98(21);
PostOutputArray CurrUseContr03(21);
PostOutputArray CurrUseContr16(21);
PostOutputArray CurrContrRace98(9);
PostOutputArray CurrContrRace03(9);
PostOutputArray CurrContrEdu98(18);
PostOutputArray CurrContrEdu03(18);
PostOutputArray CurrInjectable(56);
PostOutputArray CurrPill(56);
PostOutputArray YouthHormonal(56);
PostOutputArray BirthIntervals(15);
PostOutputArray SexuallyExpFertB(35);
PostOutputArray SexuallyExpFertC(35);
PostOutputArray SexuallyExpFertW(35);
PostOutputArray FertByMarriage(28);

// Prison outputs
PostOutputArray FractionInPrison(51);
PostOutputArray PrisonYouth(51);
PostOutputArray PrisonCompletedSchool(51);
PostOutputArray PrisonHIV(51);
PostOutputArray PrisonPrevious(51);
PostOutputArray PrisonSentenced(51);
PostOutputArray PrisonRace(3);

// Alcohol outputs
PostOutputArray BingeDrinking15to24F(51);
PostOutputArray BingeDrinking15to24M(51);
PostOutputArray BingeDrinking15to49F(51);
PostOutputArray BingeDrinking15to49M(51);
PostOutputArray BingeDrinkingHSchoolF(51);
PostOutputArray BingeDrinkingHSchoolM(51);
PostOutputArray AlcoholPast12moF(51);
PostOutputArray AlcoholPast12moM(51);
PostOutputArray AveDrinksPerDay(51);

// Household outputs
PostOutputArray IndivsInHHsized1or2(51);
PostOutputArray IndivsInHHsized3or4(51);
PostOutputArray IndivsInHHsized5or6(51);
PostOutputArray IndivsInHHsized7to9(51);
PostOutputArray IndivsInHHsized10plus(51);
PostOutputArray HeadOfHH(51);
PostOutputArray PartnerOfHead(51);
PostOutputArray ChildOfHead(51);
PostOutputArray GrandchildOfHead(51);
PostOutputArray Homeless(51);
PostOutputArray HomelessAdult(51);
PostOutputArray HomelessChild(51);
PostOutputArray HomelessChildPropn(51);
PostOutputArray HomelessFem(51);
PostOutputArray HomelessEmployed(51);
PostOutputArray HomelessAveDur(51);
PostOutputArray HomelessDrinkGT2perWeek(51);
PostOutputArray HHmembersByAge(64);

// Gender norm outputs
PostOutputArray IneqGenderAllM(51);
PostOutputArray IneqGender15to24M(51);
PostOutputArray IneqGender25to34M(51);
PostOutputArray IneqGender35plusM(51);

// Additional HCT outputs for Modelling Consortium project
PostOutputArray TotalTestsOI(56);
PostOutputArray PosTestsOI(56);
PostOutputArray NewDiagOI(56);
PostOutputArray YieldOI(56);
PostOutputArray TotalTestsPrEP(56);
PostOutputArray PosTestsPrEP(56);
PostOutputArray YieldPrEP(56);
PostOutputArray TotalTestsANC(56);
PostOutputArray PosTestsANC(56);
PostOutputArray NewDiagANC(56);
PostOutputArray YieldANC(56);
PostOutputArray TotalTestsGen(56);
PostOutputArray PosTestsGen(56);
PostOutputArray NewDiagGen(56);
PostOutputArray YieldGen(56);
PostOutputArray TotalTestsMMC(56);
PostOutputArray PosTestsMMC(56);
PostOutputArray YieldMMC(56);
PostOutputArray TotalTestsPartner(56);
PostOutputArray PosTestsPartner(56);
PostOutputArray NewDiagPartner(56);
PostOutputArray YieldPartner(56);
PostOutputArray TotalTestsHH_U(56);
PostOutputArray PosTestsHH_U(56);
PostOutputArray NewDiagHH_U(56);
PostOutputArray YieldHH_U(56);
PostOutputArray TotalTestsHH_R(56);
PostOutputArray PosTestsHH_R(56);
PostOutputArray NewDiagHH_R(56);
PostOutputArray YieldHH_R(56);
PostOutputArray TotalTestsMobile_U(56);
PostOutputArray PosTestsMobile_U(56);
PostOutputArray NewDiagMobile_U(56);
PostOutputArray YieldMobile_U(56);
PostOutputArray TotalTestsMobile_R(56);
PostOutputArray PosTestsMobile_R(56);
PostOutputArray NewDiagMobile_R(56);
PostOutputArray YieldMobile_R(56);
PostOutputArray TotalTestsFSW(56);
PostOutputArray PosTestsFSW(56);
PostOutputArray NewDiagFSW(56);
PostOutputArray YieldFSW(56);
PostOutputArray TotalTestsMSM(56);
PostOutputArray PosTestsMSM(56);
PostOutputArray NewDiagMSM(56);
PostOutputArray YieldMSM(56);
PostOutputArray TotalTestsSchool(56);
PostOutputArray PosTestsSchool(56);
PostOutputArray NewDiagSchool(56);
PostOutputArray YieldSchool(56);
PostOutputArray TotalTestsANCpartner0(56);
PostOutputArray PosTestsANCpartner0(56);
PostOutputArray NewDiagANCpartner0(56);
PostOutputArray YieldANCpartner0(56);
PostOutputArray TotalTestsANCpartner1(56);
PostOutputArray PosTestsANCpartner1(56);
PostOutputArray NewDiagANCpartner1(56);
PostOutputArray YieldANCpartner1(56);
PostOutputArray TotalTestsPrison(56);
PostOutputArray PosTestsPrison(56);
PostOutputArray NewDiagPrison(56);
PostOutputArray YieldPrison(56);
PostOutputArray TotalTestsSTI(56);
PostOutputArray PosTestsSTI(56);
PostOutputArray NewDiagSTI(56);
PostOutputArray YieldSTI(56);
PostOutputArray TotalTestsWork(56);
PostOutputArray PosTestsWork(56);
PostOutputArray NewDiagWork(56);
PostOutputArray YieldWork(56);
PostOutputArray TotalTestsFPC(56);
PostOutputArray PosTestsFPC(56);
PostOutputArray NewDiagFPC(56);
PostOutputArray YieldFPC(56);
PostOutputArray TotalTests(56);
PostOutputArray TotPosTests(56);
PostOutputArray AcuteTests(56);
PostOutputArray NewDiagnoses(56);
PostOutputArray MaleDiagnosed(56);
PostOutputArray FemDiagnosed(56);
PostOutputArray Diagnosed15to24(56);
PostOutputArray Diagnosed25to49(56);
PostOutputArray Diagnosed50plus(56);
PostOutputArray FSWdiagnosed(56);
PostOutputArray MSMdiagnosed(56);
PostOutputArray LYsLost(56);
PostOutputArray LYsLostExp(56);
PostOutputArray ARTdeathsExp(56);
PostOutputArray NewART200(56);
PostOutputArray NewART350(56);
PostOutputArray NewART500(56);
PostOutputArray NewART500plus(56);
PostOutputArray NewARTexp(56);
PostOutputArray MMCoperations(56);
PostOutputArray DiagNoART200(56);
PostOutputArray DiagNoART350(56);
PostOutputArray DiagNoART500(56);
PostOutputArray DiagNoART500plus(56);
PostOutputArray ProtSexActs(56);
PostOutputArray TreatedAdults200(56);
PostOutputArray TreatedAdults350(56);
PostOutputArray TreatedAdults500(56);
PostOutputArray TreatedAdults500plus(56);
PostOutputArray TotBirthsHIV(56);
PostOutputArray TotBirthsART(56);
PostOutputArray AIDSdeaths(56);
PostOutputArray NonAIDSdeaths(56);
PostOutputArray TotalCosts(56);
PostOutputArray TotalTestCosts(56);
PostOutputArray NewHIV_B(56);
PostOutputArray LYsLostB(56);
PostOutputArray TotalTestsB(56);
PostOutputArray NewDiagnosesB(56);
PostOutputArray TotalCostsB(56);

// Structural driver outputs
PostOutputArray IneqGenderBingeAssn(56);
PostOutputArray IneqGenderMultAssn(56);
PostOutputArray BingeMultAssnM(56);
PostOutputArray BingeMultAssnF(56);
PostOutputArray BingeCondomAssnM(56);
PostOutputArray BingeCondomAssnF(56);
PostOutputArray EmployedTransAssn(56);
PostOutputArray EmployedMultAssnM(56);
PostOutputArray EmployedMultAssnF(56);
PostOutputArray EmployedHIVassnM(56);
PostOutputArray EmployedHIVassnF(56);
PostOutputArray EduCondomAssnM(56);
PostOutputArray EduCondomAssnF(56);
PostOutputArray EduHIVassnM(56);
PostOutputArray EduHIVassnF(56);
PostOutputArray SchoolMarriageAssn(56);
PostOutputArray SingleSessionAlcohol(10);
PostOutputArray MultiSessionAlcohol(17);
PostOutputArray CashTransferOut(51);
PostOutputArray SchoolSupportOut(30);
PostOutputArray VocationalTrainingOut(17);
PostOutputArray GenderCommun(6);
PostOutputArray GenderIndiv(14);

// Parameter/log likelihood outputs
PostOutputArray2 HIVparamsLogL(10);
PostOutputArray2 TPparamsLogL(11);
PostOutputArray2 HSVparamsLogL(12);
PostOutputArray2 NGparamsLogL(11);
PostOutputArray2 CTparamsLogL(11);
PostOutputArray2 TVparamsLogL(13);
PostOutputArray2 BVparamsLogL(9);
PostOutputArray2 VCparamsLogL(7);
PostOutputArray2 MSMparamsLogL(12);
PostOutputArray2 MSMTparamsLogL(9);
PostOutputArray2 StructParamsLogL(30);
PostOutputArray2 RandomUniformHIV(9);
PostOutputArray2 RandomUniformTP(10);
PostOutputArray2 RandomUniformHSV(11);
PostOutputArray2 RandomUniformNG(10);
PostOutputArray2 RandomUniformCT(10);
PostOutputArray2 RandomUniformTV(12);
PostOutputArray2 RandomUniformBV(8);
PostOutputArray2 RandomUniformVC(6);
PostOutputArray2 RandomUniformMSM(11);
PostOutputArray2 RandomUniformMSMT(8);
PostOutputArray2 RandomUniformHCT(13);
PostOutputArray2 RandomUniformStruct(23);
PostOutputArray2 OutANCbias(1);
PostOutputArray2 OutModelVarANC(1);
PostOutputArray2 OutModelVarHH(2);

Database HIVsurvivalTimes(6);
Database ViralLoads(3);
Database UnprotectedRisk(7);
Database MigrationHIVassn(9);

HIVtransition HIVtransitionM(0, 0, 0, 0, 0, 0, 0, 27);
HIVtransition HIVtransitionF(1, 0, 0, 0, 11, 0, 95, 27);
SyphilisTransition TPtransitionM(0, 0, 0, 0, 0, 6, 0, 0);
SyphilisTransition TPtransitionF(1, 22, 6, 0, 7, 5, 0, 0);
HerpesTransition HSVtransitionM(0, 0, 0, 0, 0, 10, 0, 0);
HerpesTransition HSVtransitionF(1, 2, 1, 0, 1, 10, 0, 0);
OtherSTDtransition NGtransitionM(0, 0, 0, 0, 0, 9, 0, 0);
OtherSTDtransition NGtransitionF(1, 11, 7, 0, 7, 9, 0, 0);
OtherSTDtransition CTtransitionM(0, 0, 0, 0, 0, 10, 0, 0);
OtherSTDtransition CTtransitionF(1, 11, 7, 0, 7, 10, 0, 0);
OtherSTDtransition HDtransitionM(0, 0, 0, 11, 0, 0, 0, 0);
OtherSTDtransition HDtransitionF(1, 0, 0, 3, 0, 0, 0, 0);
OtherSTDtransition TVtransitionM(0, 0, 0, 0, 0, 2, 0, 0);
OtherSTDtransition TVtransitionF(1, 11, 7, 0, 2, 3, 0, 0);
BVtransition BVtransitionF(1, 5, 4, 0, 0, 0, 0, 0);
VCtransition VCtransitionF(1, 5, 3, 0, 2, 0, 0, 0);

RCTdata SingleSessionAlcoholCounselling(12, 1);
RCTdata MultiSessionAlcoholCounselling(17, 2);
RCTdata CashTransfers(51, 3);
RCTdata SchoolSupport(30, 4);
RCTdata VocationalTraining(17, 5);
RCTdata GenderTransformCommun(6, 6);
RCTdata GenderTransformIndiv(14, 7);

Child MaleChild(0);
Child FemChild(1);

vector<Indiv> Register(InitPop);
vector<HouseholdGroup> HHregister;
Pop RSApop(0);
PartnerCohort MalePartners[16][2]; // Cohorts of men eligible to enter relationships, by age and risk group
PartnerCohort FemPartners[16][2]; // Cohorts of women eligible to enter relationships, by age and risk group
PartnerCohort MSMpartners[16][2]; // Cohorts of MSM eligible to enter relationships

// Objects for simulating RCTs
Pop TempRSApop(0);
vector<Indiv> TempRegister(InitPop);
vector<HouseholdGroup> TempHHregister;
int TempTotCurrFSW[3];
int TempTotCurrCasualHet[3][2];
int TempTotCurrCasual;