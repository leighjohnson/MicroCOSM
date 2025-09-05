// This is the main project file for VC++ application project 
// generated using an Application Wizard.

//using std::cout;
//using std::vector;
//using std::ifstream;
using namespace std;

#include "Microsimulation.h"
#include "StatFunctions.h"
#include "randomc.h"
#include <time.h>
#include <cstring>
#include <sstream>
#include <string>
#include <stddef.h>

CRandomMersenne rg(111);
int process_num; 

int main(int argc, const char *argv[])
{
    clock_t start, finish;
	double elapsed_time;

	start = clock();

	if (argc > 1) { // There is at least one argument
		process_num = atoi(argv[1]); // Convert the first argument to an integer
		rg = CRandomMersenne(process_num); // Reinitialize the RNG
	}
	else {
		process_num = 0; // No arguments, so set the process_num to zero.
	}

	RunSimulations(0);

	/*for (int i = 0; i < Register.size(); ++i){
		memset(&Register[i], 0, sizeof(Register[0]));}
	ReadAllInputFiles();
	RSApop.AssignAgeSex();
	RSApop.GetAllPartnerRates();
	RSApop.AssignBehav();
	cout << "Propn of MSM aged <25 is " << PropnMSMagedU25() << endl;
	RSApop.AssignHIV();
	if(HSVind==1 || TPind==1 || HDind==1 || NGind==1 || CTind==1 || 
		TVind==1 || BVind==1 || VCind==1){
			RSApop.AssignSTIs();}
	CurrYear = StartYear;
	for(int ii=0; ii<ProjectionTerm; ii++){
		RSApop.OneYear();}
	RSApop.UpdateAgeGroup();*/
	//RSApop.GetPrefMatrix(1);
	//RSApop.SavePrefMatrix("PrefMatrixL.txt");
	//RSApop.GetPrefMatrix(2);
	//RSApop.SavePrefMatrix("PrefMatrixS.txt");
	//RSApop.GetPopPyramid();
	//RSApop.SavePopPyramid("PopPyramid.txt");
	//RSApop.GetBehavPyramid();
	//RSApop.SaveBehavPyramid("BehavPyramid.txt");
	//RSApop.GetNumberPartners();
	//RSApop.SaveNumberPartners("NumberPartners.txt");
	//RSApop.SaveFSWcontacts2016("FSWcontacts2016.txt");
	//RSApop.SaveLifetimePartners("LifetimePartners.txt");
	//RSApop.GetTotSex();
	//RSApop.SaveTotSex("TotSex.txt");
	//RSApop.SaveMarriageRates("MarriageRates.txt");
	//RSApop.SaveAdultHIVstage("AdultHIVstage.txt");
	//RSApop.GetInitSTIprev();

	finish = clock();
	elapsed_time = (finish - start);
	cout<<"Time taken: "<<elapsed_time<<endl;
	system("PAUSE");
	return 0;
}


PrevalenceData::PrevalenceData(){
	LogL = 0.0;
}

NationalData::NationalData()
{
	int ia;

	for (ia = 0; ia<19; ia++){
		BiasMult[ia] = 1.0;
	}
	ModelVarEst = 0.0;
}

void NationalData::CalcModelVar()
{
	// This function should ONLY get called if GetSDfromData = 1

	int ia, ic;
	double AdjPrev, StochasticVar;
	double VarLogitPrev[100]; // Variance of the logit-transformed prevalence estimates 
	// (sampling variation only). Change dimension if n>100.
	double SpSum, BiasSum, BiasLevel;

	BiasSum = 0.0;
	for (ia = 0; ia<19; ia++){
		BiasSum += log(BiasMult[ia]);
	}
	SpSum = 0.0;
	for (ic = 0; ic<Observations; ic++){
		SpSum += 1.0 - ExpSp[ic];
	}

	BiasLevel = 0.0;
	/*if (BiasSum > 0.0 || SpSum > 0.0){ // ANC data
		for (ic = 0; ic<Observations; ic++){
			BiasLevel += log(StudyPrev[ic] / (1.0 - StudyPrev[ic])) - log(ModelPrev[ic] / (1.0 -
				ModelPrev[ic]));
		}
		BiasLevel = BiasLevel / Observations;
		if (BiasLevel < 0.0){ BiasLevel = 0.0; }
	}*/

	ModelVarEst = 0.0;
	for (ic = 0; ic<Observations; ic++){
		AdjPrev = 1.0 / (1.0 + (1.0 / ModelPrev[ic] - 1.0) * exp(-BiasLevel));
		StochasticVar = 0.0;
		//StochasticVar = exp(2.0 * (111.271 - 0.0561 * StudyYear[ic] - 4.299 * ModelPrev[ic] +
		//	3.342 * pow(ModelPrev[ic], 2.0))) / IterationsPerPC;
		VarLogitPrev[ic] = pow(PrevSE[ic] / (StudyPrev[ic] * (1.0 - StudyPrev[ic])), 2.0) + StochasticVar;
		if (BiasSum > 0.0 || SpSum > 0.0){ VarLogitPrev[ic] += ANCbiasVar; }
		else{ VarLogitPrev[ic] += HSRCbiasVar; }
		ModelVarEst += pow(log(StudyPrev[ic] / (1.0 - StudyPrev[ic])) - log(AdjPrev / (1.0 -
			AdjPrev)), 2.0) - VarLogitPrev[ic];
	}
	ModelVarEst = ModelVarEst / Observations;
}

void NationalData::CalcLogL()
{
	int ia, ic, SimCount2;
	double AdjPrev, StochasticVar;
	double VarLogitPrev[100]; // Variance of the logit-transformed prevalence estimates 
	// (sampling variation only). Change dimension if n>100.
	double SpSum, BiasSum, BiasLevel;

	SimCount2 = (CurrSim - 1) / IterationsPerPC;
	BiasSum = 0.0;
	for (ia = 0; ia<19; ia++){
		BiasSum += log(BiasMult[ia]);
	}
	SpSum = 0.0;
	for (ic = 0; ic<Observations; ic++){
		SpSum += 1.0 - ExpSp[ic];
	}

	BiasLevel = 0.0;
	/*if (BiasSum > 0.0 || SpSum > 0.0){ // ANC data
		for (ic = 0; ic<Observations; ic++){
			BiasLevel += log(StudyPrev[ic] / (1.0 - StudyPrev[ic])) - log(ModelPrev[ic] / (1.0 -
				ModelPrev[ic]));
		}
		BiasLevel = BiasLevel / Observations;
		if (BiasLevel < 0.0){ BiasLevel = 0.0; }
		if (FixedUncertainty == 1){ OutANCbias.out[SimCount2][0] = exp(BiasLevel); }
	}*/

	LogL = 0.0;
	for (ic = 0; ic<Observations; ic++){
		AdjPrev = 1.0 / (1.0 + (1.0 / ModelPrev[ic] - 1.0) * exp(-BiasLevel));
		if (GetSDfromData == 0){
			LogL += -0.5 * log(2.0 * 3.141592654) - log(PrevSE[ic])
				- 0.5 * pow((StudyPrev[ic] - AdjPrev) / PrevSE[ic], 2.0);
		}
		else{
			StochasticVar = 0.0;
			//StochasticVar = exp(2.0 * (111.271 - 0.0561 * StudyYear[ic] - 4.299 * ModelPrev[ic] +
			//	3.342 * pow(ModelPrev[ic], 2.0))) / IterationsPerPC;
			VarLogitPrev[ic] = pow(PrevSE[ic] / (StudyPrev[ic] * (1.0 - StudyPrev[ic])), 2.0) + StochasticVar;
			if (BiasSum > 0.0 || SpSum > 0.0){ VarLogitPrev[ic] += ANCbiasVar; }
			else{ VarLogitPrev[ic] += HSRCbiasVar; }
		}
		ModelPrev[ic] = AdjPrev;
	}

	if (GetSDfromData == 1){
		if (FixedUncertainty == 1){
			if (BiasSum > 0.0 || SpSum > 0.0){
				OutModelVarANC.out[SimCount2][0] = ModelVarEst;
			}
			else if (SexInd == 1){
				OutModelVarHH.out[SimCount2][0] = ModelVarEst;
			}
			else{
				OutModelVarHH.out[SimCount2][1] = ModelVarEst;
			}
		}
		for (ic = 0; ic<Observations; ic++){
			LogL += -0.5 * (log(2.0 * 3.141592654 * (VarLogitPrev[ic] + ModelVarEst)) +
				pow(log(StudyPrev[ic] / (1.0 - StudyPrev[ic])) - log(ModelPrev[ic] / (1.0 -
				ModelPrev[ic])), 2.0) / (VarLogitPrev[ic] + ModelVarEst));
		}
	}
}

AntenatalN::AntenatalN(){}

HouseholdN::HouseholdN(){}

SentinelData::SentinelData()
{
	BiasMult = 1.0;
}

void SentinelData::CalcLogL()
{
	int ic;
	double Mean, Var; // The mean and variance of theta(i), the modelled prevalence
	// of the STD in study i
	double MeanFb, VarFb; // The mean and variance of f(bi), where bi is the 'random
	// effect' for study i
	double alpha, beta; // The alpha and beta parameters for the beta prior on theta(i)
	double a, b, c, d; // Arguments for the gamma functions
	double LogLikelihood;

	LogL = 0.0;
	for (ic = 0; ic<Observations; ic++){
		if (ModelPrev[ic]<0.0001){ ModelPrev[ic] = 0.0001; } // To be consistent with GetPrev function for HIV
		if (ModelPrev[ic]>0.9999){ ModelPrev[ic] = 0.9999; }
		ModelPrev[ic] *= BiasMult;
		MeanFb = ModelPrev[ic] + VarStudyEffect * ModelPrev[ic] * (1.0 - ModelPrev[ic]) *
			(0.5 - ModelPrev[ic]);
		VarFb = pow(ModelPrev[ic] * (1.0 - ModelPrev[ic]), 2.0) * (VarStudyEffect +
			(1.5 - 8.0 * ModelPrev[ic] + 8.0 * pow(ModelPrev[ic], 2.0)) * pow(VarStudyEffect,
			2.0) + pow(pow(ModelPrev[ic], 2.0) - ModelPrev[ic] + 1.0 / 6.0, 2.0) * 15.0 *
			pow(VarStudyEffect, 3.0));
		Mean = 1.0 - ExpSp[ic] + MeanFb * (ExpSe[ic] + ExpSp[ic] - 1.0);
		Var = VarFb * (pow(ExpSe[ic] + ExpSp[ic] - 1.0, 2.0) + VarSp[ic] + VarSe[ic]) +
			VarSe[ic] * pow(MeanFb, 2.0) + VarSp[ic] * pow(1.0 - MeanFb, 2.0);
		StudyPos[ic] = StudyN[ic] * StudyPrev[ic];

		if (Var>0){
			alpha = Mean * (Mean * (1.0 - Mean) / Var - 1.0);
			beta = (1.0 - Mean) * (Mean * (1.0 - Mean) / Var - 1.0);

			a = alpha + beta;
			b = alpha + StudyPos[ic];
			c = beta + StudyN[ic] - StudyPos[ic];
			d = alpha + beta + StudyN[ic];

			LogLikelihood = gamma_log(&a) + gamma_log(&b) + gamma_log(&c) - gamma_log(&d) -
				gamma_log(&alpha) - gamma_log(&beta);
		}
		else{
			// In this case, the mean is fixed, so the likelihood is just the binomial.
			a = StudyN[ic] + 1.0;
			b = StudyPos[ic] + 1.0;
			c = StudyN[ic] - StudyPos[ic] + 1.0;
			LogLikelihood = gamma_log(&a) - gamma_log(&b) - gamma_log(&c) + StudyPos[ic] *
				log(Mean) + (StudyN[ic] - StudyPos[ic]) * log(1.0 - Mean);
		}
		LogL += LogLikelihood;
	}
}

Household::Household(){}

NonHousehold::NonHousehold(){}

ANC::ANC(){}

FPC::FPC(){}

GUD::GUD(){}

CSW::CSW(){}

RCTdata::RCTdata(int n, int code)
{
	DataPoints = n;
	InterventionCode = code;
}

void RCTdata::CalcModelEst(int ii)
{
	int ic, ia, ia1, ia2, ig, ig2, ie, ir, ih, PID, AgeDif, eligible;
	double temp1, temp2, AUDITC, HHID;
	
	ia1 = StudyDetails[ii][2];
	ia2 = StudyDetails[ii][3];
	ig = StudyDetails[ii][1];
	ie = StudyDetails[ii][4];

	temp1 = 0.0;
	temp2 = 0.0;

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].BaselineAge >= ia1 && Register[ic].BaselineAge <= ia2 &&
			(Register[ic].SexInd == ig || ig == 2) && (Register[ic].BaselineSchool == ie || ie == 3 ||
			(ie == 2 && Register[ic].BaselineSchool == 0 && Register[ic].BaselineUnemployed == 1))){
		//if (Register[ic].AliveInd == 1 && Register[ic].CurrAge >= ia1 && Register[ic].CurrAge <= ia2 &&
		//	(Register[ic].SexInd == ig || ig == 2) && (Register[ic].InSchool == ie || ie == 3 ||
		//	(ie == 2 && Register[ic].InSchool == 0 && Register[ic].Employed == 0))){
			eligible = 1;
			if (InterventionCode == 1 || InterventionCode == 2){
				// Alcohol counselling interventions limit to regular binge drinkers
				if (Register[ic].BaselineBinge == 0){ eligible = 0; }
			}
			if (InterventionCode == 3 || InterventionCode == 4 || InterventionCode == 5){
				// Economic interventions limited to people in households with log income below average
				if (Register[ic].BaselinePoor == 0){ eligible = 0; }
				// School support interventions limited to youth in school
				if (Register[ic].BaselineSchool == 0 && InterventionCode == 4 && ie==1){ eligible = 0; }
				if (Register[ic].BaselineSchool == 1 && InterventionCode == 4 && (ie==0 || ie==2)){ eligible = 0; }
				// Vocational training limited to people out of school
				if (Register[ic].BaselineSchool == 1 && InterventionCode == 5){ eligible = 0; }
			}
			if (eligible == 1){
				temp2 += Register[ic].PopWeight;
				if (StudyDetails[ii][0] == 1){ // Proportion of drinking days
					temp1 += Register[ic].DailyDrinkProb * Register[ic].PopWeight;
				}
				if (StudyDetails[ii][0] == 2){ // Drinks per drinking day
					temp1 += Register[ic].DrinksPerDD * Register[ic].PopWeight;
				}
				if (StudyDetails[ii][0] == 3){ // Binge drinking in last month
					if (Register[ic].DrinksPerDD >= 5.0 && Register[ic].DailyDrinkProb > 1.0 / 30.0){
						temp1 += Register[ic].PopWeight;
					}
				}
				if (StudyDetails[ii][0] == 4 || StudyDetails[ii][0] == 13){ // AUDIT or AUDIT-C 
					AUDITC = 0.0;
					if (Register[ic].DailyDrinkProb >= 0.5){ AUDITC += 4.0; }
					else if (Register[ic].DailyDrinkProb >= 0.21){ AUDITC += 3.0; }
					else if (Register[ic].DailyDrinkProb >= 0.05){ AUDITC += 2.0; }
					else if (Register[ic].DailyDrinkProb >= 0.01){ AUDITC += 1.0; }
					if (Register[ic].DrinksPerDD >= 9.5){ AUDITC += 4.0; }
					else if (Register[ic].DrinksPerDD >= 6.5){ AUDITC += 3.0; }
					else if (Register[ic].DrinksPerDD >= 4.5){ AUDITC += 2.0; }
					else if (Register[ic].DrinksPerDD >= 2.5){ AUDITC += 1.0; }
					if (Register[ic].DrinksPerDD >= 5.5){
						if (Register[ic].DailyDrinkProb >= 0.5){ AUDITC += 4.0; }
						else if (Register[ic].DailyDrinkProb >= 0.09){ AUDITC += 3.0; }
						else if (Register[ic].DailyDrinkProb >= 0.03){ AUDITC += 2.0; }
						else if (Register[ic].DailyDrinkProb >= 0.01){ AUDITC += 1.0; }
					}
					if (StudyDetails[ii][0] == 4){
						// Adjustment to take into account that change in AUDIT is greater than 
						// change in AUDIT-C (based on Wandera et al, 2017).
						if (StudyDetails[ii][7] >= 12){ AUDITC *= 1.7; }
						else{ AUDITC *= (1.0 + 0.7 * (StudyDetails[ii][7] / 12.0)); }
					}
					temp1 += AUDITC * Register[ic].PopWeight;
				}
				if (StudyDetails[ii][0] == 5){ // No drinking in last month
					temp1 += exp(-30.0 * Register[ic].DailyDrinkProb) * Register[ic].PopWeight;
				}
				if (StudyDetails[ii][0] == 6){ // Employment
					if (Register[ic].Employed == 1){ temp1 += Register[ic].PopWeight; }
				}
				if (StudyDetails[ii][0] == 7){ // Monthly earnings on log(x+1) scale
					if (Register[ic].Employed == 1){
						ia = Register[ic].AgeGroup - 3;
						if (ia > 9){ ia = 9; }
						ig2 = Register[ic].SexInd;
						ir = Register[ic].Race;
						if (Register[ic].HighestGrade == 0){ ih = 0; }
						else if (Register[ic].HighestGrade <= 7){ ih = 1; }
						else if (Register[ic].HighestGrade <= 11){ ih = 2; }
						else if (Register[ic].HighestGrade == 12){ ih = 3; }
						else { ih = 4; }
						temp1 += log(1.0 + exp(CurrAveIncome[ia][ig2][ir][ih] + Register[ic].LogIncomeDif)) *
							Register[ic].PopWeight;
					}
				}
				if (StudyDetails[ii][0] == 8){ // School dropout
					if (Register[ic].InSchool == 0 && Register[ic].HighestGrade < 12){
						temp1 += Register[ic].PopWeight;
					}
				}
				if (StudyDetails[ii][0] == 9){ // Prob of endorsing inequitable gender norms
					temp1 += Register[ic].IneqGender * Register[ic].PopWeight;
				}
				if (StudyDetails[ii][0] == 10){ // Concurrency
					// This output code used to be used for 'cumulative unprotected sex' but
					// this output is no longer defined.
					//if (Register[ic].CumUnprotected>0){ temp1 += Register[ic].PopWeight; }
					if (Register[ic].CurrPartners > 0){
						if (Register[ic].CurrPartners == 2){ temp1 += Register[ic].PopWeight; }
					}
					else{ temp2 = temp2 - Register[ic].PopWeight; }
				}
				if (StudyDetails[ii][0] == 11){ // Unprotected sex with recent partners
					if (Register[ic].CurrPartners > 0 && Register[ic].CondomPrimary == 0){
						temp1 += Register[ic].PopWeight;
					}
					else if (Register[ic].CurrPartners == 2 && Register[ic].Condom2ndary == 0){
						temp1 += Register[ic].PopWeight;
					}
					else if (Register[ic].UAIcasual > 0 || Register[ic].UVICSW > 0){
						temp1 += Register[ic].PopWeight;
					}
				}
				if (StudyDetails[ii][0] == 12){ // Condom use at last sex
					if (Register[ic].CurrPartners == 1){
						temp1 += Register[ic].CondomPrimary * Register[ic].PopWeight;
					}
					else if (Register[ic].CurrPartners == 2){
						temp1 += 0.5 * (Register[ic].CondomPrimary + Register[ic].Condom2ndary) * Register[ic].PopWeight;
					}
					else{ temp2 = temp2 - Register[ic].PopWeight; }
				}
				if (StudyDetails[ii][0] == 13){ // Current partners
					// I decided not to include this because only Kalichman (2007) reported this and
					// there wasn't enough data presented to calculate the SE properly.
					// Subsequently code 13 was re-assigned to the AUDIT-C outcome (see code for outcome 4).
				}
				if (StudyDetails[ii][0] == 14){ // Cumulative partners (over last 12 months)
					if (Register[ic].AnnPartners > 1){ temp1 += Register[ic].PopWeight; }
				}
				if (StudyDetails[ii][0] == 15){ // Cumulative casual sex
					if (Register[ic].CumCasual>0){ temp1 += Register[ic].PopWeight; }
				}
				if (StudyDetails[ii][0] == 16){ // Current casual sex
					if (Register[ic].CasualInd == 1 || Register[ic].HetCasualInd == 1){ 
						temp1 += Register[ic].PopWeight; }
				}
				if (StudyDetails[ii][0] == 17){ // Sexual debut
					if (Register[ic].VirginInd==0){ temp1 += Register[ic].PopWeight; }
				}
				if (StudyDetails[ii][0] == 18){ // Recent sexual activity
					if (Register[ic].CurrPartners > 0 || Register[ic].CasualInd == 1 || 
						Register[ic].HetCasualInd == 1 || Register[ic].FSWind == 1){
						temp1 += Register[ic].PopWeight;
					}
				}
				if (StudyDetails[ii][0] == 19){ // Partner age difference >5 years
					if (Register[ic].CurrPartners > 0 || Register[ic].HetCasualInd == 1 || Register[ic].CasualInd == 1){
						PID = 0; // default
						if (Register[ic].HetCasualInd == 1 || Register[ic].CasualInd == 1){
							PID = Register[ic].IDofCasual;
						}
						if(PID == 0){ PID = Register[ic].IDprimary; }
						if (PID > 0){
							AgeDif = Register[PID - 1].CurrAge - Register[ic].CurrAge;
							if (AgeDif >= 5 || AgeDif <= -5){ temp1 += Register[ic].PopWeight; }
						}
					}
				}
				if (StudyDetails[ii][0] == 20){ // Cumulative incidence of marriage
					if (Register[ic].MarriedInd == 1){ temp1 += Register[ic].PopWeight; }
				}
				if (StudyDetails[ii][0] == 21){ // Cumulative HIV incidence
					if (Register[ic].HIVstage > 0 && Register[ic].DateInfect > (BaselineStart + 0.5)){ 
						temp1 += Register[ic].PopWeight; }
				}
				if (StudyDetails[ii][0] == 22){ // Cumulative HSV-2 incidence
					if (Register[ic].HSVstage > 0 && Register[ic].CumHSV2 > 0){
						temp1 += Register[ic].PopWeight;
					}
				}
				if (StudyDetails[ii][0] == 23){ // Cumulative NG/CT/TV incidence
					if (Register[ic].CumSTIs > 0){ temp1 += Register[ic].PopWeight; }
				}
				if (StudyDetails[ii][0] == 24){ // Teenage pregnancy
					if (Register[ic].BaselineAge < 20 && Register[ic].CumTeenPreg == 1){ 
						temp1 += Register[ic].PopWeight; }
					else if (Register[ic].BaselineAge >= 20){
						temp2 = temp2 - Register[ic].PopWeight; }
				}
				if (StudyDetails[ii][0] == 25){ // HIV prevalence
					if (Register[ic].HIVstage > 0){ temp1 += Register[ic].PopWeight; }
				}
				if (StudyDetails[ii][0] == 26){ // HSV-2 prevalence
					if (Register[ic].HSVstage > 0){ temp1 += Register[ic].PopWeight; }
				}
				if (StudyDetails[ii][0] == 27){ // Male circumcision
					if (Register[ic].CircInd==1){ temp1 += Register[ic].PopWeight; }
				}
			}
		}
	}
	if (StructIntScenario == 0){ ModelEstimates[ii][0] = temp1 / temp2; }
	else{ ModelEstimates[ii][1] = temp1 / temp2; }
}

void RCTdata::GetCurrModelEsts()
{
	int ii;
	double duration;

	duration = (12.0 * BehavCycleCount / CycleS) + 12.0 * (CurrYear - BaselineStart);
	
	for (ii = 0; ii < DataPoints; ii++){
		if (StudyDetails[ii][7] == duration){ 
			CalcModelEst(ii); 
		}
	}
}

void RCTdata::CalcModelEffects()
{
	int ii;
	double ScalingFactor;

	for (ii = 0; ii < DataPoints; ii++){
		if (StudyDetails[ii][6] == 1){
			ModelEffects[ii] = log(ModelEstimates[ii][1] / (1.0 - ModelEstimates[ii][1])) -
				log(ModelEstimates[ii][0] / (1.0 - ModelEstimates[ii][0]));
		}
		else if (StudyDetails[ii][6] == 2){
			ModelEffects[ii] = ModelEstimates[ii][1] / ModelEstimates[ii][0];
		}
		else if (StudyDetails[ii][6] == 3){
			ModelEffects[ii] = log(ModelEstimates[ii][1] / ModelEstimates[ii][0]);
		}
		else if (StudyDetails[ii][6] == 4){
			ModelEffects[ii] = ModelEstimates[ii][1] - ModelEstimates[ii][0];
		}
		else if (StudyDetails[ii][6] == 5){
			ModelEffects[ii] = ModelEstimates[ii][1] - ModelEstimates[ii][0];
		}
		if (InterventionCode == 3 || InterventionCode == 4){
			ScalingFactor = pow(StudyDetails[ii][5] / 800.0, 0.5);
			if (InterventionCode == 4 && StudyDetails[ii][0] == 8){
				ScalingFactor = (log(RRdropoutSupport) + ScalingFactor * log(RRdropoutIncSupport)) / 
					(log(RRdropoutSupport) + log(RRdropoutIncSupport));
			}
			if (StudyDetails[ii][6] != 2){ ModelEffects[ii] *= ScalingFactor; }
			else{ ModelEffects[ii] = pow(ModelEffects[ii], ScalingFactor); }
		}
	}
}

void RCTdata::GetMaxTerm()
{
	int ii;
	double Term, intpart1, DoubleIntDif1;

	MaxTerm = StudyDetails[0][7];
	for (ii = 0; ii < DataPoints; ii++){
		if (StudyDetails[ii][7] > MaxTerm){ MaxTerm = StudyDetails[ii][7]; }
	}
	Term = 1.0 * MaxTerm / 12.0;
	DoubleIntDif1 = modf(Term, &intpart1);
	if (DoubleIntDif1 == 0){ MaxTerm = Term; }
	else{ MaxTerm = Term + 1; }
}

double RCTdata::GetLikelihood()
{
	int ii;
	double dif, LogL;

	LogL = 0.0;
	AltLogL = 0.0;
	//CalcModelEffects();
	for (ii = 0; ii < DataPoints; ii++){
		if (ModelEstimates[ii][0] > 0.0 || ModelEstimates[ii][1] > 0.0){
			dif = StudyDetails[ii][8] - ModelEffects[ii];
			LogL += -0.5 * (log(2.0 * 3.141592654 * pow(StudyDetails[ii][9], 2.0)) +
				pow(dif, 2.0) / pow(StudyDetails[ii][9], 2.0));
			// Sensitivity analysis to assess the effect of excluding the Kalichman RCTs
			if (InterventionCode == 1 && StudyDetails[ii][0] != 11 &&
				StudyDetails[ii][0] != 12 && StudyDetails[ii][0] != 16){
				AltLogL += -0.5 * (log(2.0 * 3.141592654 * pow(StudyDetails[ii][9], 2.0)) +
					pow(dif, 2.0) / pow(StudyDetails[ii][9], 2.0));
			}
		}
		/*if (ii == 0){
			cout << "InterventionCode: " << InterventionCode << endl;
		}
		cout << "ModelEffects for study " << ii + 1 << ": " << ModelEffects[ii] << endl;*/
	}
	//cout << "LogL: " <<LogL << endl;

	return LogL;
}

PostOutputArray::PostOutputArray(int n)
{
	columns = n;
}

void PostOutputArray::RecordSample(const char* filout)
{
	int i, c;
	ostringstream s; 

	if (process_num > 0){
		s << process_num << filout;
	}
	else{
		s << filout;
	}

	//ofstream file(filout);
	ofstream file(s.str().c_str()); // Converts s to a C string

	for(i=0; i<samplesize; i++){
		file<<setw(6)<<right<<i<<"	";
		for(c=0; c<columns; c++){
			file << "	" << setw(10) << right << out[i][c];
		}
		file<<endl;
	}
	file.close();
}

void PostOutputArray::ReadBaseline(const char* input)
{
	int ir, ic, dummy;
	ifstream file;

	file.open(input);
	if (file.fail()) {
		cerr << "Could not open input file.txt\n";
		exit(1);
	}

	for (ir = 0; ir < samplesize; ir++){
		file >> dummy;
		for (ic = 0; ic < columns; ic++){
			file >> out[ir][ic]; }
	}
	file.close();
}

void PostOutputArray::CalcMeans()
{
	int iy, ic, count;
	double sumx, sumx2, tempsum;

	for (iy = 0; iy <= CurrYear - StartYear; iy++){
		sumx = 0.0;
		sumx2 = 0.0;
		count = 0;
		for (ic = 0; ic < samplesize; ic++){
			if (out[ic][iy] < 0.0 || out[ic][iy] >= 0.0){
				count += 1;
				sumx += out[ic][iy];
				sumx2 += out[ic][iy] * out[ic][iy];
			}
		}
		Means[iy] = 1.0 * sumx / count;
		StdErrors[iy] = 1.0 * pow((count * sumx2 - sumx * sumx) / (count - 1), 0.5) / count;
	}
}

void PostOutputArray::CalcPeriodTot(int Y1, int Y2)
{
	// Calculate mean and standard deviation of aggregated total for years Y1-Y2

	int iy, ic;
	double sumx, sumx2, tempsum;

	sumx = 0.0;
	sumx2 = 0.0;
	for (ic = 0; ic < samplesize; ic++){
		tempsum = 0.0;
		for (iy = Y1 - StartYear; iy <= Y2 - StartYear; iy++){
			tempsum += out[ic][iy];
		}
		sumx += tempsum;
		sumx2 += tempsum * tempsum;
	}
	PeriodTot = 1.0 * sumx / samplesize;
	PeriodTotSE = 1.0 * pow((samplesize * sumx2 - sumx * sumx) / (samplesize - 1), 0.5) / samplesize;
}

void PostOutputArray::OutputRatio(PostOutputArray* array1, PostOutputArray* array2)
{
	int ic, iy;

	for (ic = 0; ic < samplesize; ic++){
		for (iy = 0; iy <= CurrYear - StartYear; iy++){
			if (array2->out[ic][iy] > 0.0){
				out[ic][iy] = 1.0 * array1->out[ic][iy] / array2->out[ic][iy];}
			else{ out[ic][iy] = 0.0; }
		}
	}
}

void PostOutputArray::CalcPeriodRatio(int Y1, int Y2, PostOutputArray* array1, PostOutputArray* array2)
{
	int iy, ic;
	double sumx, sumx2, tempsum, tempsum2;

	sumx = 0.0;
	sumx2 = 0.0;
	for (ic = 0; ic < samplesize; ic++){
		tempsum = 0.0;
		tempsum2 = 0.0;
		for (iy = Y1 - StartYear; iy <= Y2 - StartYear; iy++){
			tempsum += array1->out[ic][iy];
			tempsum2 += array2->out[ic][iy];
		}
		if (tempsum2 > 0.0){
			sumx += 1.0 * tempsum / tempsum2;
			sumx2 += 1.0 * pow(tempsum / tempsum2, 2.0);
		}
	}
	PeriodTot = 1.0 * sumx / samplesize;
	PeriodTotSE = 1.0 * pow((samplesize * sumx2 - sumx * sumx) / (samplesize - 1), 0.5) / samplesize;
}

void PostOutputArray::CalcPeriodIncrease(int Y1, int Y2, PostOutputArray* array1)
{
	int iy, ic;
	double sumx, sumy, sumxmy, sumxmy2;

	sumxmy = 0.0;
	sumxmy2 = 0.0;
	for (ic = 0; ic < samplesize; ic++){
		sumx = 0.0;
		sumy = 0.0;
		for (iy = Y1 - StartYear; iy <= Y2 - StartYear; iy++){
			sumx += out[ic][iy];
			sumy += array1->out[ic][iy];
		}
		sumxmy += sumx - sumy;
		sumxmy2 += pow(sumx - sumy, 2.0);
	}
	PeriodIncrease = 1.0 * sumxmy / samplesize;
	PeriodIncreaseSE = 1.0 * pow((samplesize * sumxmy2 - sumxmy * sumxmy) / (samplesize - 1), 0.5) / samplesize;
}

void PostOutputArray::GetSummGenOut(int col)
{
	int iy;

	CalcMeans();

	for (iy = col; iy < CurrYear - StartYear; iy++){
		SummOut[SummOutRow][iy] = Means[iy-col];
		SummOut[SummOutRow + 1][iy] = StdErrors[iy-col];
	}
	SummOutRow += 2;
}

void PostOutputArray::GetSummHCTout1(int row)
{
	int iy;

	CalcMeans();
	CalcPeriodTot(2019, 2038);

	for (iy = 0; iy < CurrYear - 2000; iy++){
		SummOut[row - 1][iy] = Means[iy + 15];}

	SummOut[row - 1][41] = PeriodTot;
	SummOut[row - 1][42] = PeriodTotSE;
}

void PostOutputArray::GetSummHCTout2(int row, PostOutputArray* array1, PostOutputArray* array2)
{
	int iy;

	OutputRatio(array1, array2);
	CalcMeans();
	CalcPeriodRatio(2019, 2038, array1, array2);

	for (iy = 0; iy < CurrYear - 2000; iy++){
		SummOut[row - 1][iy] = Means[iy + 15];
	}

	SummOut[row - 1][41] = PeriodTot;
	SummOut[row - 1][42] = PeriodTotSE;
}

void PostOutputArray::CalcSvyAssnLogL(double SvyMat[4][2], int year)
{
	double TempLogL, SvyVar, ModelLogOR;
	int iy, ii;

	if (year == 0){ iy = 21; } // 2 years prior to 2008
	if (year == 1){ iy = 25; } // 2 years prior to 2012
	if (year == 2){ iy = 29; } // 2 years prior to 2016
	if (year == 3){ iy = 30; } // 2 years prior to 2017

	if (SvyMat[year][1] > 0.0){
		SvyVar = pow(SvyMat[year][1] / 1.96, 2.0);
		ModelLogOR = 0.0;
		for (ii = iy; ii < iy + 5; ii++){
			if (out[CurrSim - 1][ii] > 0.0 || out[CurrSim - 1][ii] < 0.0){ ModelLogOR += out[CurrSim - 1][ii]; }
		}
		ModelLogOR = ModelLogOR / 5.0;
		TempLogL = -0.5 * (log(2.0 * 3.141592654 * SvyVar) + pow(SvyMat[year][0] -
			ModelLogOR, 2.0) / SvyVar);
		StructLogL += TempLogL;
	}
}

PostOutputArray2::PostOutputArray2(int n)
{
	columns = n;
}

void PostOutputArray2::RecordSample(const char* filout)
{
	int i, c;
	ostringstream s;

	if (process_num > 0){
		s << process_num << filout;
	}
	else{
		s << filout;
	}

	//ofstream file(filout);
	ofstream file(s.str().c_str()); // Converts s to a C string

	for (i = 0; i<CurrSim; i++){
		file << setw(6) << right << i << "	";
		for (c = 0; c<columns; c++){
			file << "	" << setw(10) << setprecision(10) << right << out[i][c];
		}
		file << endl;
	}
	file.close();
}

PostOutputArrayB::PostOutputArrayB(int n, int n2, int n3)
{
	columns = n;
	studies = n2;
	AgeStd = n3;
}

void PostOutputArrayB::CalcRandomEffectVar()
{
	int ir, iy, BeginYr;
	double dif, AveByAge[15][2];

	RandomEffectVar = 0.0;
	if (AgeStd == 0){
		for (ir = 0; ir < studies; ir++){
			ModelAve[ir] = 0.0;
			BeginYr = CalibData[ir][3] - StartYear - 2;
			for (iy = BeginYr - 2; iy < BeginYr + 5; iy++){
				ModelAve[ir] += out[CurrSim - 1][iy] * 0.2;}
			dif = log(ModelAve[ir] / (1.0 - ModelAve[ir])) - log(CalibData[ir][0] / (1.0 - CalibData[ir][0]));
			RandomEffectVar += pow(dif, 2.0) - pow(CalibData[ir][1], 2.0);
		}
	}
	else{
		for (ir = 0; ir < studies; ir++){
			AveByAge[ir][0] = 0.0;
			AveByAge[ir][1] = 0.0;
			BeginYr = CalibData[ir][3] - StartYear - 2;
			for (iy = BeginYr - 2; iy < BeginYr + 5; iy++){
				AveByAge[ir][0] += out2[CurrSim - 1][iy][0] * 0.2;
				AveByAge[ir][1] += out2[CurrSim - 1][iy][1] * 0.2;
			}
			ModelAve[ir] = (1.0 - CalibData[ir][2]) * AveByAge[ir][0] + CalibData[ir][2] * AveByAge[ir][1];
			dif = log(ModelAve[ir] / (1.0 - ModelAve[ir])) - log(CalibData[ir][0] / (1.0 - CalibData[ir][0]));
			RandomEffectVar += pow(dif, 2.0) - pow(CalibData[ir][1], 2.0);
		}
	}
	RandomEffectVar = RandomEffectVar / studies;
}

void PostOutputArrayB::CalcLogL()
{
	int ir;
	double dif;

	LogL = 0.0;
	if (RandomEffectVar < 0.0){ RandomEffectVar = 0.0; }
	for (ir = 0; ir < studies; ir++){
		dif = log(ModelAve[ir] / (1.0 - ModelAve[ir])) - log(CalibData[ir][0] / (1.0 - CalibData[ir][0]));
		LogL += -0.5 * (log(2.0 * 3.141592654 * (RandomEffectVar + pow(CalibData[ir][1], 2.0))) + 
			pow(dif, 2.0) / (RandomEffectVar + pow(CalibData[ir][1], 2.0)));
	}
}

void PostOutputArrayB::RecordSample(const char* filout)
{
	int i, c;
	ostringstream s;

	if (process_num > 0){
		s << process_num << filout;
	}
	else{
		s << filout;
	}

	//ofstream file(filout);
	ofstream file(s.str().c_str()); // Converts s to a C string

	for (i = 0; i<samplesize; i++){
		file << setw(6) << right << i << "	";
		for (c = 0; c<columns; c++){
			file << "	" << setw(10) << right << out[i][c];
		}
		file << endl;
	}
	file.close();
}

Database::Database(int n)
{
	columns = n;
}

void Database::RecordSample(const char* filout)
{
	int i, c;
	ostringstream s;

	if (process_num > 0){
		s << process_num << filout;
	}
	else{
		s << filout;
	}

	//ofstream file(filout);
	ofstream file(s.str().c_str()); // Converts s to a C string

	for (i = 0; i<TotIndivs; i++){
		file << setw(6) << right << i << "	";
		for (c = 0; c<columns; c++){
			file << "	" << setw(10) << right << out[i][c];
		}
		file << endl;
	}
	file.close();
}

STDtransition::STDtransition(){}

void STDtransition::ReadPrevData(const char* input)
{
	int ic;
	ifstream file;

	file.open(input);
	if (file.fail()) {
		cerr << "Could not open input file.txt\n";
		exit(1);
	}
	file.ignore(255, '\n');
	if (ANClogL.Observations>0){
		for (ic = 0; ic<ANClogL.Observations; ic++){
			file >> ANClogL.StudyYear[ic] >> ANClogL.StudyN[ic] >> ANClogL.StudyPrev[ic] >>
				ANClogL.ExpSe[ic] >> ANClogL.ExpSp[ic] >> ANClogL.VarSe[ic] >> ANClogL.VarSp[ic] >>
				ANClogL.HIVprevInd[ic] >> ANClogL.HIVprev[ic];
		}
		file.ignore(255, '\n');
	}
	file.ignore(255, '\n');
	if (FPClogL.Observations>0){
		for (ic = 0; ic<FPClogL.Observations; ic++){
			file >> FPClogL.StudyYear[ic] >> FPClogL.StudyN[ic] >> FPClogL.StudyPrev[ic] >>
				FPClogL.ExpSe[ic] >> FPClogL.ExpSp[ic] >> FPClogL.VarSe[ic] >> FPClogL.VarSp[ic] >>
				FPClogL.HIVprevInd[ic] >> FPClogL.HIVprev[ic];
		}
		file.ignore(255, '\n');
	}
	file.ignore(255, '\n');
	if (GUDlogL.Observations>0){
		for (ic = 0; ic<GUDlogL.Observations; ic++){
			file >> GUDlogL.StudyYear[ic] >> GUDlogL.StudyN[ic] >> GUDlogL.StudyPrev[ic] >>
				GUDlogL.ExpSe[ic] >> GUDlogL.ExpSp[ic] >> GUDlogL.VarSe[ic] >> GUDlogL.VarSp[ic] >>
				GUDlogL.HIVprevInd[ic] >> GUDlogL.HIVprev[ic];
		}
		file.ignore(255, '\n');
	}
	file.ignore(255, '\n');
	if (CSWlogL.Observations>0){
		for (ic = 0; ic<CSWlogL.Observations; ic++){
			file >> CSWlogL.StudyYear[ic] >> CSWlogL.StudyN[ic] >> CSWlogL.StudyPrev[ic] >>
				CSWlogL.ExpSe[ic] >> CSWlogL.ExpSp[ic] >> CSWlogL.VarSe[ic] >> CSWlogL.VarSp[ic] >>
				CSWlogL.HIVprevInd[ic] >> CSWlogL.HIVprev[ic];
		}
		file.ignore(255, '\n');
	}
	file.ignore(255, '\n');
	if (HouseholdLogL.Observations>0){
		for (ic = 0; ic<HouseholdLogL.Observations; ic++){
			file >> HouseholdLogL.StudyYear[ic] >> HouseholdLogL.StudyN[ic] >>
				HouseholdLogL.StudyPrev[ic] >> HouseholdLogL.ExpSe[ic] >>
				HouseholdLogL.ExpSp[ic] >> HouseholdLogL.VarSe[ic] >> HouseholdLogL.VarSp[ic] >>
				HouseholdLogL.HIVprevInd[ic] >> HouseholdLogL.HIVprev[ic] >>
				HouseholdLogL.AgeStart[ic] >> HouseholdLogL.AgeEnd[ic] >>
				HouseholdLogL.ExclVirgins[ic];
		}
		file.ignore(255, '\n');
	}
	file.ignore(255, '\n');
	if (AntenatalNlogL.Observations>0){
		for (ic = 0; ic<AntenatalNlogL.Observations; ic++){
			file >> AntenatalNlogL.StudyYear[ic] >> AntenatalNlogL.AgeStart[ic] >>
				AntenatalNlogL.StudyPrev[ic] >> AntenatalNlogL.PrevSE[ic] >>
				AntenatalNlogL.ExpSe[ic] >> AntenatalNlogL.ExpSp[ic];
		}
		file.ignore(255, '\n');
	}
	file.ignore(255, '\n');
	if (HouseholdNlogL.Observations>0){
		for (ic = 0; ic<HouseholdNlogL.Observations; ic++){
			file >> HouseholdNlogL.StudyYear[ic] >> HouseholdNlogL.AgeStart[ic] >>
				HouseholdNlogL.StudyPrev[ic] >> HouseholdNlogL.PrevSE[ic] >>
				HouseholdNlogL.ExpSe[ic] >> HouseholdNlogL.ExpSp[ic];
		}
	}
	file.close();
}

void STDtransition::GetCSWprev()
{
	int ic, iy;

	if (CSWlogL.Observations>0){
		for (ic = 0; ic<CSWlogL.Observations; ic++){
			iy = CSWlogL.StudyYear[ic] - StartYear;
			CSWlogL.ModelPrev[ic] = 0.05 * CSWprevUnsmoothed[iy - 3];
			CSWlogL.ModelPrev[ic] += 0.12 * CSWprevUnsmoothed[iy - 2];
			CSWlogL.ModelPrev[ic] += 0.20 * CSWprevUnsmoothed[iy - 1];
			CSWlogL.ModelPrev[ic] += 0.26 * CSWprevUnsmoothed[iy];
			CSWlogL.ModelPrev[ic] += 0.20 * CSWprevUnsmoothed[iy + 1];
			CSWlogL.ModelPrev[ic] += 0.12 * CSWprevUnsmoothed[iy + 2];
			CSWlogL.ModelPrev[ic] += 0.05 * CSWprevUnsmoothed[iy + 3];
		}
	}
}

void STDtransition::GetANCprevSmooth()
{
	int ic, iy;

	if (ANClogL.Observations>0){
		for (ic = 0; ic<ANClogL.Observations; ic++){
			iy = ANClogL.StudyYear[ic] - StartYear;
			ANClogL.ModelPrev[ic] = 0.05 * ANCprevUnsmoothed[iy - 3];
			ANClogL.ModelPrev[ic] += 0.12 * ANCprevUnsmoothed[iy - 2];
			ANClogL.ModelPrev[ic] += 0.20 * ANCprevUnsmoothed[iy - 1];
			ANClogL.ModelPrev[ic] += 0.26 * ANCprevUnsmoothed[iy];
			ANClogL.ModelPrev[ic] += 0.20 * ANCprevUnsmoothed[iy + 1];
			ANClogL.ModelPrev[ic] += 0.12 * ANCprevUnsmoothed[iy + 2];
			ANClogL.ModelPrev[ic] += 0.05 * ANCprevUnsmoothed[iy + 3];
		}
	}
}

void STDtransition::SetVarStudyEffect(double Variance)
{
	ANClogL.VarStudyEffect = Variance;
	FPClogL.VarStudyEffect = Variance;
	CSWlogL.VarStudyEffect = Variance;
	GUDlogL.VarStudyEffect = Variance;
	HouseholdLogL.VarStudyEffect = Variance;
}

NonHIVtransition::NonHIVtransition(){}

HIVtransition::HIVtransition(int Sex, int ObsANC, int ObsFPC, int ObsGUD, int ObsCSW,
	int ObsHH, int ObsANCN, int ObsHHN)
{
	SexInd = Sex;
	nStates = 6;

	ANClogL.Observations = ObsANC;
	FPClogL.Observations = ObsFPC;
	GUDlogL.Observations = ObsGUD;
	CSWlogL.Observations = ObsCSW;
	HouseholdLogL.Observations = ObsHH;
	AntenatalNlogL.Observations = ObsANCN;
	HouseholdNlogL.Observations = ObsHHN;
	HouseholdNlogL.SexInd = Sex;
}

void HIVtransition::CalcTransitionProbs()
{
	int iy;

	iy = CurrYear - StartYear;
	From1to2 = 1.0 - exp(-(1.0/AveDuration[0]) * 52.0/CycleD);
	From2to3 = 1.0 - exp(-(1.0/AveDuration[1]) * 52.0/CycleD);
	From3to4 = (1.0 - exp(-(1.0/AveDuration[2]) * 52.0/CycleD)) * (1.0 - HAARTaccess[iy]);
	if (HAARTaccess[iy] > 1.0){
		From3to4 = 0.0;}
	From3to5 = (1.0 - exp(-(1.0/AveDuration[2]) * 52.0/CycleD)) * HAARTaccess[iy];
	From5to6 = 1.0 - exp(-ARTinterruption / CycleD);
	From6to5 = 1.0 - exp(-ARTresumption / CycleD);
	From4toDead = (1.0 - exp(-(1.0/AveDuration[3]) * 52.0/CycleD));
	From5toDead = (1.0 - exp(-(1.0/AveDuration[4]) * 52.0/CycleD));
	From6toDead = (1.0 - exp(-(1.0/AveDuration[4]) * 52.0/CycleD));
}

void HIVtransition::GetPrev()
{
	// Aggregates the results for the relevant parameter combination

	int ii, ia, iy;

	if (SexInd == 1){
		for (ii = 0; ii < AntenatalNlogL.Observations; ii++){
			AntenatalNlogL.ModelPrev[ii] = 0.0;
		}
	}
	for (ii = 0; ii < HouseholdNlogL.Observations; ii++){
		HouseholdNlogL.ModelPrev[ii] = 0.0;
	}

	for (ii = 0; ii < IterationsPerPC; ii++){
		// Antenatal results
		if (SexInd == 1){
			for (iy = 0; iy < 16; iy++){
				AntenatalNlogL.ModelPrev[iy * 5] += 0.05 * PrevPreg15.out[CurrSim - 1 - ii][iy + 4] +
					0.12 * PrevPreg15.out[CurrSim - 1 - ii][iy + 5] + 0.20 * PrevPreg15.out[CurrSim - 1 - ii][iy + 6] +
					0.26 * PrevPreg15.out[CurrSim - 1 - ii][iy + 7] + 0.20 * PrevPreg15.out[CurrSim - 1 - ii][iy + 8] +
					0.12 * PrevPreg15.out[CurrSim - 1 - ii][iy + 9] + 0.05 * PrevPreg15.out[CurrSim - 1 - ii][iy + 10];
				AntenatalNlogL.ModelPrev[iy * 5 + 1] += 0.05 * PrevPreg20.out[CurrSim - 1 - ii][iy + 4] +
					0.12 * PrevPreg20.out[CurrSim - 1 - ii][iy + 5] + 0.20 * PrevPreg20.out[CurrSim - 1 - ii][iy + 6] +
					0.26 * PrevPreg20.out[CurrSim - 1 - ii][iy + 7] + 0.20 * PrevPreg20.out[CurrSim - 1 - ii][iy + 8] +
					0.12 * PrevPreg20.out[CurrSim - 1 - ii][iy + 9] + 0.05 * PrevPreg20.out[CurrSim - 1 - ii][iy + 10];
				AntenatalNlogL.ModelPrev[iy * 5 + 2] += 0.05 * PrevPreg25.out[CurrSim - 1 - ii][iy + 4] +
					0.12 * PrevPreg25.out[CurrSim - 1 - ii][iy + 5] + 0.20 * PrevPreg25.out[CurrSim - 1 - ii][iy + 6] +
					0.26 * PrevPreg25.out[CurrSim - 1 - ii][iy + 7] + 0.20 * PrevPreg25.out[CurrSim - 1 - ii][iy + 8] +
					0.12 * PrevPreg25.out[CurrSim - 1 - ii][iy + 9] + 0.05 * PrevPreg25.out[CurrSim - 1 - ii][iy + 10];
				AntenatalNlogL.ModelPrev[iy * 5 + 3] += 0.05 * PrevPreg30.out[CurrSim - 1 - ii][iy + 4] +
					0.12 * PrevPreg30.out[CurrSim - 1 - ii][iy + 5] + 0.20 * PrevPreg30.out[CurrSim - 1 - ii][iy + 6] +
					0.26 * PrevPreg30.out[CurrSim - 1 - ii][iy + 7] + 0.20 * PrevPreg30.out[CurrSim - 1 - ii][iy + 8] +
					0.12 * PrevPreg30.out[CurrSim - 1 - ii][iy + 9] + 0.05 * PrevPreg30.out[CurrSim - 1 - ii][iy + 10];
				AntenatalNlogL.ModelPrev[iy * 5 + 4] += 0.05 * PrevPreg35.out[CurrSim - 1 - ii][iy + 4] +
					0.12 * PrevPreg35.out[CurrSim - 1 - ii][iy + 5] + 0.20 * PrevPreg35.out[CurrSim - 1 - ii][iy + 6] +
					0.26 * PrevPreg35.out[CurrSim - 1 - ii][iy + 7] + 0.20 * PrevPreg35.out[CurrSim - 1 - ii][iy + 8] +
					0.12 * PrevPreg35.out[CurrSim - 1 - ii][iy + 9] + 0.05 * PrevPreg35.out[CurrSim - 1 - ii][iy + 10];
			}
		}
		// Household survey results
		for (ia = 0; ia < 9; ia++){
			HouseholdNlogL.ModelPrev[ia] += PrevHH2005.out[CurrSim - 1 - ii][ia + 9 * SexInd];
			HouseholdNlogL.ModelPrev[ia+9] += PrevHH2008.out[CurrSim - 1 - ii][ia + 9 * SexInd];
			HouseholdNlogL.ModelPrev[ia+18] += PrevHH2012.out[CurrSim - 1 - ii][ia + 9 * SexInd];
		}
	}

	// Take averages
	if (SexInd == 1){
		for (ii = 0; ii < AntenatalNlogL.Observations; ii++){
			AntenatalNlogL.ModelPrev[ii] = AntenatalNlogL.ModelPrev[ii]/IterationsPerPC; 
			if (AntenatalNlogL.ModelPrev[ii] == 0.0){ AntenatalNlogL.ModelPrev[ii] = 0.0001; }
		}
	}
	for (ii = 0; ii < HouseholdNlogL.Observations; ii++){
		HouseholdNlogL.ModelPrev[ii] = HouseholdNlogL.ModelPrev[ii]/IterationsPerPC;
		if (HouseholdNlogL.ModelPrev[ii] == 0.0){ HouseholdNlogL.ModelPrev[ii] = 0.0001; }
	}
}

double HIVtransition::GetTransmProb(int ID)
{
	int PID1, PID2, CSWID, is, ig, ia, ic;
	int IRisk, PRisk, SameSex;
	double NoTransmProb, SingleActProb;

	IRisk = Register[ID-1].RiskGroup;
	SameSex = 0;
	NoTransmProb = 1.0;

	// Infection from primary partner
	if(Register[ID-1].CurrPartners>0){
		PID1 = Register[ID-1].IDprimary;
		if (Register[PID1 - 1].SexInd == Register[ID - 1].SexInd){ SameSex = 1; }
		if(Register[PID1-1].HIVstage>0 && SameSex==0){
			// Get base probability depending on partnership type, risk group
			PRisk = Register[PID1-1].RiskGroup;
			if(Register[ID-1].MarriedInd==1){
				if(PRisk==1 && IRisk==1){
					SingleActProb = TransmProb[4];}
				else if(PRisk==2 && IRisk==2){
					SingleActProb = TransmProb[6];}
				else{
					SingleActProb = TransmProb[5];}
				if(PRisk==2 && NoViralTransm12==1){
					SingleActProb = 0.0;}
				if(PRisk==2 && IRisk==2 && NoViralTransm22==1){
					SingleActProb = 0.0;}
			}
			else{
				if(PRisk==1 && IRisk==1){
					SingleActProb = TransmProb[1];}
				else if(PRisk==2 && IRisk==2){
					SingleActProb = TransmProb[3];}
				else{
					SingleActProb = TransmProb[2];}
			}
			// Make adjustment for HIV stage of partner
			SingleActProb *= GetInfectivityMult(PID1);
			// Make adjustment for age of individual
			ia = Register[ID-1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			// Make adjustment for contraception
			if (Register[PID1 - 1].SexInd == 1 && (Register[PID1 - 1].CurrContr == 1 ||
				Register[PID1 - 1].CurrContr == 2)){
				SingleActProb *= HIVinfectHormonal;
			}
			if (Register[ID - 1].SexInd == 1 && Register[ID - 1].CurrContr == 1){
				SingleActProb *= HIVsusceptInjectable;}
			// Make adjustment for pregnancy
			if (Register[ID - 1].SexInd == 1 && Register[ID - 1].DOLB > (1.0*CurrYear + 0.5 + 
				1.0*(BehavCycleCount - 1.0) / CycleS)){
				SingleActProb *= PregEffectHIVtransm[0];
			}
			if (Register[PID1 - 1].SexInd == 1 && Register[PID1 - 1].DOLB > (1.0*CurrYear + 0.5 +
				1.0*(BehavCycleCount - 1.0) / CycleS)){
				SingleActProb *= PregEffectHIVtransm[1];
			}
			// Make adjustment for STI cofactors
			if (CofactorType == 2){
				SingleActProb *= GetCofactorType2(ID, PID1); }
			// Make adjustment for male circumcision and PrEP
			if (Register[ID - 1].CircInd == 1){
				SingleActProb *= (1.0 - CircEff); }
			if (Register[ID - 1].OnPrEP > 0.0){
				SingleActProb *= 1.0 - PrEPefficacy * Register[ID - 1].OnPrEP;}
			// Make adjustment for HIV vaccines
			SingleActProb *= (1.0 - Register[ID - 1].VaccineEff);
			// Make adjustment for other sources of heterogeneity
			if(AllowHIVsuscepAdj==1){
				SingleActProb *= Register[ID-1].SuscepHIVadj;}
			if(SingleActProb>1.0){
				SingleActProb = 1.0;}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID-1].UVIprimary);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff), 
				Register[ID-1].PVIprimary);
		}
	}

	// Infection from 2ndary partner
	if(Register[ID-1].CurrPartners==2){
		PID2 = Register[ID-1].ID2ndary;
		if (Register[PID2 - 1].SexInd == Register[ID - 1].SexInd){ SameSex = 1; }
		else{ SameSex = 0; }
		if(Register[PID2-1].HIVstage>0 && SameSex==0){
			// Get base probability depending on risk group
			PRisk = Register[PID2-1].RiskGroup;
			if(PRisk==1 && IRisk==1){
				SingleActProb = TransmProb[1];}
			else if(PRisk==2 && IRisk==2){
				SingleActProb = TransmProb[3];}
			else{
				SingleActProb = TransmProb[2];}
			// Make adjustment for HIV stage of partner
			SingleActProb *= GetInfectivityMult(PID2);
			// Make adjustment for age of individual
			ia = Register[ID-1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			// Make adjustment for contraception
			if (Register[PID2 - 1].SexInd == 1 && (Register[PID2 - 1].CurrContr == 1 ||
				Register[PID2 - 1].CurrContr == 2)){
				SingleActProb *= HIVinfectHormonal;
			}
			if (Register[ID - 1].SexInd == 1 && Register[ID - 1].CurrContr == 1){
				SingleActProb *= HIVsusceptInjectable;}
			// Make adjustment for pregnancy
			if (Register[ID - 1].SexInd == 1 && Register[ID - 1].DOLB > (1.0*CurrYear + 0.5 +
				1.0*(BehavCycleCount - 1.0) / CycleS)){
				SingleActProb *= PregEffectHIVtransm[0];
			}
			if (Register[PID2 - 1].SexInd == 1 && Register[PID2 - 1].DOLB > (1.0*CurrYear + 0.5 +
				1.0*(BehavCycleCount - 1.0) / CycleS)){
				SingleActProb *= PregEffectHIVtransm[1];
			}
			// Make adjustment for STI cofactors
			if (CofactorType == 2){
				SingleActProb *= GetCofactorType2(ID, PID2); }
			// Make adjustment for male circumcision and PrEP
			if (Register[ID - 1].CircInd == 1){
				SingleActProb *= (1.0 - CircEff); }
			if (Register[ID - 1].OnPrEP > 0.0){
				SingleActProb *= 1.0 - PrEPefficacy * Register[ID - 1].OnPrEP;}
			// Make adjustment for HIV vaccines
			SingleActProb *= (1.0 - Register[ID - 1].VaccineEff);
			// Make adjustment for other sources of heterogeneity
			if(AllowHIVsuscepAdj==1){
				SingleActProb *= Register[ID-1].SuscepHIVadj;}
			if(SingleActProb>1.0){
				SingleActProb = 1.0;}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID-1].UVI2ndary);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff), 
				Register[ID-1].PVI2ndary);
		}
	}

	// Infection from CSW (relevant only to high-risk men)
	if(Register[ID-1].UVICSW + Register[ID-1].PVICSW > 0){
		PID2 = Register[ID-1].IDofCSW;
		if(Register[PID2-1].HIVstage>0){
			SingleActProb = TransmProb[0];
			SingleActProb *= GetInfectivityMult(PID2);
			// Make adjustment for age of individual
			ia = Register[ID-1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			// Make adjustment for contraception
			if (Register[PID2 - 1].CurrContr == 1 || Register[PID2 - 1].CurrContr == 2){
				SingleActProb *= HIVinfectHormonal; }
			// Make adjustment for pregnancy
			if (Register[PID2 - 1].DOLB > (1.0*CurrYear + 0.5 + 1.0*(BehavCycleCount - 1.0) / CycleS)){
				SingleActProb *= PregEffectHIVtransm[1];}
			// Make adjustment for STI cofactors
			if (CofactorType == 2){
				SingleActProb *= GetCofactorType2(ID, PID2); }
			// Make adjustment for male circumcision and PrEP
			if (Register[ID - 1].CircInd == 1){
				SingleActProb *= (1.0 - CircEff); }
			if (Register[ID - 1].OnPrEP > 0.0){
				SingleActProb *= 1.0 - PrEPefficacy * Register[ID - 1].OnPrEP;}
			// Make adjustment for HIV vaccines
			SingleActProb *= (1.0 - Register[ID - 1].VaccineEff);
			// Make adjustment for other sources of heterogeneity
			if(AllowHIVsuscepAdj==1){
				SingleActProb *= Register[ID-1].SuscepHIVadj;}
			if(SingleActProb>1.0){
				SingleActProb = 1.0;}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID-1].UVICSW);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff), 
				Register[ID-1].PVICSW);
		}
	}

	// Infection from client (relevant only to CSWs)
	if(Register[ID-1].FSWind==1){
		ia = Register[ID-1].AgeGroup - 2;
		for(ic=0; ic<Register.size(); ic++){
			if(Register[ic].IDofCSW==ID && (Register[ic].UVICSW + Register[ic].PVICSW) > 0){
				if(Register[ic].HIVstage>0){
					SingleActProb = TransmProb[0];
					SingleActProb *= GetInfectivityMult(ic+1);
					if(AllowHIVsuscepAdj==1){
						SingleActProb *= Register[ID-1].SuscepHIVadj;}
					// Make adjustment for contraception
					if (Register[ID - 1].CurrContr == 1){
						SingleActProb *= HIVsusceptInjectable; }
					// Make adjustment for pregnancy
					if (Register[ID - 1].DOLB > (1.0*CurrYear + 0.5 + 1.0*(BehavCycleCount - 1.0) / CycleS)){
						SingleActProb *= PregEffectHIVtransm[0];}
					// Make adjustment for STI cofactors
					if (CofactorType == 2){
						SingleActProb *= GetCofactorType2(ID, ic+1);}
					// Make adjustment for PrEP
					if (Register[ID - 1].OnPrEP > 0.0){
						SingleActProb *= 1.0 - PrEPefficacy * Register[ID - 1].OnPrEP;}
					// Make adjustment for HIV vaccines
					SingleActProb *= (1.0 - Register[ID - 1].VaccineEff);
					// Make adjustment for other sources of heterogeneity
					if (AllowHIVsuscepAdj == 1){
						SingleActProb *= Register[ID - 1].SuscepHIVadj;}
					if(SingleActProb>1.0){
						SingleActProb = 1.0;}
					NoTransmProb *= pow(1.0 - SingleActProb, Register[ic].UVICSW);
					NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff), 
						Register[ic].PVICSW);
				}
			}
		}
	}

	// Infection from casual partner
	if (Register[ID - 1].HetCasualInd == 1 && (Register[ID - 1].UAIcasual + Register[ID - 1].PAIcasual) > 0){
		PID2 = Register[ID - 1].IDofCasual;
		if (Register[PID2 - 1].HIVstage>0){
			// Get base probability depending on risk group (same as for ST partners)
			PRisk = Register[PID2 - 1].RiskGroup;
			if (PRisk == 1 && IRisk == 1){
				SingleActProb = TransmProb[1]; }
			else if (PRisk == 2 && IRisk == 2){
				SingleActProb = TransmProb[3]; }
			else{
				SingleActProb = TransmProb[2]; }
			// Make adjustment for HIV stage of partner
			SingleActProb *= GetInfectivityMult(PID2);
			// Make adjustment for STI cofactors
			if (CofactorType == 2){
				SingleActProb *= GetCofactorType2(ID, PID2); }
			// Make adjustment for PrEP
			if (Register[ID - 1].OnPrEP > 0.0){
				SingleActProb *= 1.0 - PrEPefficacy * Register[ID - 1].OnPrEP; }
			// Make adjustment for HIV vaccines
			SingleActProb *= (1.0 - Register[ID - 1].VaccineEff);
			// Make adjustment for other sources of heterogeneity
			if (AllowHIVsuscepAdj == 1){
				SingleActProb *= Register[ID - 1].SuscepHIVadj; }
			if (SingleActProb>1.0){
				SingleActProb = 1.0; }
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID - 1].UAIcasual);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
				Register[ID - 1].PAIcasual);
		}
	}

	return 1.0 - NoTransmProb;
}

double HIVtransition::GetTransmProbMSM(int ID)
{
	int PID1, PID2, CSWID, is, ig, ia, ic;
	int PRisk, IRisk, SameSex;
	double IInsert, PInsert, NoTransmProb, SingleActProb;

	IRisk = Register[ID - 1].RiskGroup;
	IInsert = Register[ID - 1].InsertivePref;
	SameSex = 0;
	NoTransmProb = 1.0;

	// Infection from primary partner
	if (Register[ID - 1].CurrPartners>0){
		PID1 = Register[ID - 1].IDprimary;
		if (Register[PID1 - 1].SexInd == Register[ID - 1].SexInd){ SameSex = 1; }
		if (Register[PID1 - 1].HIVstage>0 && SameSex == 1){
			// Get base probability depending on partnership type, role pref
			PRisk = Register[PID1 - 1].RiskGroup;
			PInsert = Register[PID1 - 1].InsertivePref;
			if (Register[ID - 1].MarriedInd == 1){
				if (PInsert == 0.5 && IInsert == 0.5){ // Both versatile
					SingleActProb = 1.0 - pow((1.0 - TransmProbMSM[5][0]) * 
						(1.0 - TransmProbMSM[5][1]), (1.0+IntraEventVersatility)/2.0);
				}
				else if ((PInsert == 1.0 && IInsert<1.0) || (PInsert == 0.5 && IInsert == 0.0)){
					SingleActProb = TransmProbMSM[5][0]; // Receptive
				}
				else if ((PInsert == 0.0 && IInsert>0.0) || (PInsert == 0.5 && IInsert == 1.0)){
					SingleActProb = TransmProbMSM[5][1]; // Insertive
				}
				else{ SingleActProb = 0.0; } // No anal sex assumed if role prefs incompatible.
				if (PRisk == 2 && NoViralTransm12 == 1){
					SingleActProb = 0.0;
				}
				if (PRisk == 2 && IRisk == 2 && NoViralTransm22 == 1){
					SingleActProb = 0.0;
				}
			}
			else{
				if (PInsert == 0.5 && IInsert == 0.5){ // Both versatile
					SingleActProb = 1.0 - pow((1.0 - TransmProbMSM[2][0]) * 
						(1.0 - TransmProbMSM[2][1]), (1.0 + IntraEventVersatility) / 2.0);
				}
				else if ((PInsert == 1.0 && IInsert<1.0) || (PInsert == 0.5 && IInsert == 0.0)){
					SingleActProb = TransmProbMSM[2][0]; // Receptive
				}
				else if ((PInsert == 0.0 && IInsert>0.0) || (PInsert == 0.5 && IInsert == 1.0)){
					SingleActProb = TransmProbMSM[2][1]; // Insertive
				}
				else{ SingleActProb = 0.0; } // No anal sex assumed if role prefs incompatible.
			}
			// Make adjustment for HIV stage of partner
			SingleActProb *= GetInfectivityMult(PID1);
			// Make adjustment for STI cofactors
			if (CofactorType==2){
				SingleActProb *= GetCofactorType2(ID, PID1); }
			// Make adjustment for PrEP
			if (Register[ID - 1].OnPrEP > 0.0){
				SingleActProb *= 1.0 - PrEPefficacyMSM * Register[ID - 1].OnPrEP;}
			// Make adjustment for HIV vaccines
			SingleActProb *= (1.0 - Register[ID - 1].VaccineEff);
			// Make adjustment for other sources of heterogeneity
			if (AllowHIVsuscepAdj == 1){
				SingleActProb *= Register[ID - 1].SuscepHIVadj;
			}
			if (SingleActProb>1.0){
				SingleActProb = 1.0;
			}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID - 1].UVIprimary);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
				Register[ID - 1].PVIprimary);
		}
	}

	// Infection from 2ndary partner
	if (Register[ID - 1].CurrPartners == 2){
		PID2 = Register[ID - 1].ID2ndary;
		if (Register[PID2 - 1].SexInd == Register[ID - 1].SexInd){ SameSex = 1; }
		else{ SameSex = 0; }
		if (Register[PID2 - 1].HIVstage>0 && SameSex == 1){
			// Get base probability depending on risk group
			PRisk = Register[PID2 - 1].RiskGroup;
			PInsert = Register[PID2 - 1].InsertivePref;
			if (PInsert == 0.5 && IInsert == 0.5){ // Both versatile
				SingleActProb = 1.0 - pow((1.0 - TransmProbMSM[2][0]) *
					(1.0 - TransmProbMSM[2][1]), (1.0 + IntraEventVersatility) / 2.0);
			}
			else if ((PInsert == 1.0 && IInsert<1.0) || (PInsert == 0.5 && IInsert == 0.0)){
				SingleActProb = TransmProbMSM[2][0]; // Receptive
			}
			else if ((PInsert == 0.0 && IInsert>0.0) || (PInsert == 0.5 && IInsert == 1.0)){
				SingleActProb = TransmProbMSM[2][1]; // Insertive
			}
			else{ SingleActProb = 0.0; } // No anal sex assumed if role prefs incompatible.
			// Make adjustment for HIV stage of partner
			SingleActProb *= GetInfectivityMult(PID2);
			// Make adjustment for STI cofactors
			if (CofactorType == 2){
				SingleActProb *= GetCofactorType2(ID, PID2); }
			// Make adjustment for PrEP
			if (Register[ID - 1].OnPrEP > 0.0){
				SingleActProb *= 1.0 - PrEPefficacyMSM * Register[ID - 1].OnPrEP;}
			// Make adjustment for HIV vaccines
			SingleActProb *= (1.0 - Register[ID - 1].VaccineEff);
			// Make adjustment for other sources of heterogeneity
			if (AllowHIVsuscepAdj == 1){
				SingleActProb *= Register[ID - 1].SuscepHIVadj;
			}
			if (SingleActProb>1.0){
				SingleActProb = 1.0;
			}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID - 1].UVI2ndary);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
				Register[ID - 1].PVI2ndary);
		}
	}

	// Infection from casual partner
	if (Register[ID - 1].CasualInd == 1 && (Register[ID - 1].UAIcasual + Register[ID - 1].PAIcasual) > 0){
		PID2 = Register[ID - 1].IDofCasual;
		if (Register[PID2 - 1].HIVstage>0){
			// Get base probability depending on risk group
			PRisk = Register[PID2 - 1].RiskGroup;
			PInsert = Register[PID2 - 1].InsertivePref;
			if (PInsert == 0.5 && IInsert == 0.5){ // Both versatile
				SingleActProb = 1.0 - pow((1.0 - TransmProbMSM[0][0]) *
					(1.0 - TransmProbMSM[0][1]), (1.0 + IntraEventVersatility) / 2.0);
			}
			else if ((PInsert == 1.0 && IInsert<1.0) || (PInsert == 0.5 && IInsert == 0.0)){
				SingleActProb = TransmProbMSM[0][0]; // Receptive
			}
			else if ((PInsert == 0.0 && IInsert>0.0) || (PInsert == 0.5 && IInsert == 1.0)){
				SingleActProb = TransmProbMSM[0][1]; // Insertive
			}
			else{ SingleActProb = 0.0; } // No anal sex assumed if role prefs incompatible.
			// Make adjustment for HIV stage of partner
			SingleActProb *= GetInfectivityMult(PID2);
			// Make adjustment for STI cofactors
			if (CofactorType == 2){
				SingleActProb *= GetCofactorType2(ID, PID2); }
			// Make adjustment for PrEP
			if (Register[ID - 1].OnPrEP > 0.0){
				SingleActProb *= 1.0 - PrEPefficacyMSM * Register[ID - 1].OnPrEP;}
			// Make adjustment for HIV vaccines
			SingleActProb *= (1.0 - Register[ID - 1].VaccineEff);
			// Make adjustment for other sources of heterogeneity
			if (AllowHIVsuscepAdj == 1){
				SingleActProb *= Register[ID - 1].SuscepHIVadj;
			}
			if (SingleActProb>1.0){
				SingleActProb = 1.0;
			}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID - 1].UAIcasual);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
				Register[ID - 1].PAIcasual);
		}
	}

	return 1.0 - NoTransmProb;
}

double HIVtransition::GetCofactorType2(int ID, int ID2)
{
	int ig, ig2;
	double IndivSTI, PartnerSTI;

	IndivSTI = 1.0;
	PartnerSTI = 1.0;
	ig = Register[ID-1].SexInd;
	ig2 = Register[ID2-1].SexInd;

	// Calculate effect of HIV-neg individual's STIs
	if ((HSVind == 1 && (Register[ID - 1].HSVstage == 1 || Register[ID - 1].HSVstage == 3)) ||
		(TPind == 1 && Register[ID - 1].TPstage == 2) || (HDind == 1 && Register[ID - 1].HDstage == 1)){
		IndivSTI *= (1.0 + SuscepIncreaseSyndrome[0][ig]);
	}
	else if ((NGind == 1 && Register[ID - 1].NGstage == 1) || (CTind == 1 && Register[ID - 1].CTstage == 1) ||
		(TVind == 1 && Register[ID - 1].TVstage == 1) || (VCind == 1 && Register[ID - 1].VCstage == 2) ||
		(BVind == 1 && Register[ID - 1].BVstage == 2)){
		IndivSTI *= (1.0 + SuscepIncreaseSyndrome[1][ig]);
	}
	else if ((HSVind == 1 && (Register[ID - 1].HSVstage == 2 || Register[ID - 1].HSVstage == 4)) ||
		(TPind == 1 && (Register[ID - 1].TPstage == 1 || Register[ID - 1].TPstage == 3)) || 
		(HDind == 1 && Register[ID - 1].HDstage == 2) || (NGind == 1 && Register[ID - 1].NGstage == 2) ||
		(CTind == 1 && Register[ID - 1].CTstage == 2) || (TVind == 1 && Register[ID - 1].TVstage == 2) ||
		(VCind == 1 && Register[ID - 1].VCstage == 1) || (BVind == 1 && (Register[ID - 1].BVstage == 1 ||
		Register[ID - 1].BVstage == 3))){
		IndivSTI *= (1.0 + SuscepIncreaseSyndrome[2][ig]);
	}

	// Calculate effect of HIV-pos individual's STIs
	if ((HSVind == 1 && (Register[ID2 - 1].HSVstage == 1 || Register[ID2 - 1].HSVstage == 3)) ||
		(TPind == 1 && Register[ID2 - 1].TPstage == 2) || (HDind == 1 && Register[ID2 - 1].HDstage == 1)){
		PartnerSTI *= (1.0 + InfecIncreaseSyndrome[0][ig2]);
	}
	else if ((NGind == 1 && Register[ID2 - 1].NGstage == 1) || (CTind == 1 && Register[ID2 - 1].CTstage == 1) ||
		(TVind == 1 && Register[ID2 - 1].TVstage == 1) || (VCind == 1 && Register[ID2 - 1].VCstage == 2) ||
		(BVind == 1 && Register[ID2 - 1].BVstage == 2)){
		PartnerSTI *= (1.0 + InfecIncreaseSyndrome[1][ig2]);
	}
	else if ((HSVind == 1 && (Register[ID2 - 1].HSVstage == 2 || Register[ID2 - 1].HSVstage == 4)) ||
		(TPind == 1 && (Register[ID2 - 1].TPstage == 1 || Register[ID2 - 1].TPstage == 3)) ||
		(HDind == 1 && Register[ID2 - 1].HDstage == 2) || (NGind == 1 && Register[ID2 - 1].NGstage == 2) ||
		(CTind == 1 && Register[ID2 - 1].CTstage == 2) || (TVind == 1 && Register[ID2 - 1].TVstage == 2) ||
		(VCind == 1 && Register[ID2 - 1].VCstage == 1) || (BVind == 1 && (Register[ID2 - 1].BVstage == 1 ||
		Register[ID2 - 1].BVstage == 3))){
		PartnerSTI *= (1.0 + InfecIncreaseSyndrome[2][ig2]);
	}

	return IndivSTI * PartnerSTI;
}

double HIVtransition::GetInfectivityMult(int ID)
{
	int is;
	double Multiplier, VLoffset;

	is = Register[ID - 1].HIVstage;

	if (is == 1){
		VLoffset = 6.0 - Register[ID - 1].SPVL;
		if (VLoffset < 0.0){ VLoffset = 0.0; }
		Multiplier = (1.0 + HIVinfecIncrease[0]) * exp(-VLeffectTransm[0] * 
			(pow(VLoffset, VLeffectTransm[1]) - pow(1.5, VLeffectTransm[1])));
	}
	else{
		VLoffset = 6.0 - Register[ID - 1].logVL;
		if (VLoffset < 0.0){ VLoffset = 0.0; }
		Multiplier = exp(-VLeffectTransm[0] *
			(pow(VLoffset, VLeffectTransm[1]) - pow(1.5, VLeffectTransm[1])));
		if (is == 5){ Multiplier *= (1.0 + HIVinfecIncrease[4]); }
	}
	if (is == 4 || Register[ID - 1].CD4<200){ Multiplier *= (1.0 + HIVinfecIncrease[3]); }

	return Multiplier;
}

void HIVtransition::GetNewStage(int ID, double p)
{
	int is, iy, ia;
	double ExactAge, temp1, temp2, MortProb, duration;

	// Calculate CD4 & VL changes and mortality
	is = Register[ID-1].HIVstage;
	iy = CurrYear - StartYear;
	temp2 = 0.0;
	if(is==1){
		if(p<From1to2){Register[ID-1].HIVstageE = 2;}
		else{Register[ID-1].HIVstageE = 1;}
	}
	if (is > 1){
		MortProb = MortZeroCD4 * pow(RednMortCD4, Register[ID - 1].CD4);
		if (is == 6){ MortProb = MortZeroCD4 * pow(RednMortCD4, Register[ID - 1].BaselineCD4); }
		if (is >= 5){
			duration = 1.0 * CurrYear + 0.5 + 1.0 * (BehavCycleCount - 1)/CycleS + 
				1.0 * (STDcycleCount - 1) / CycleD - Register[ID - 1].DOAI;
			MortProb *= (ARTmortConstant * pow(Register[ID - 1].BaselineCD4, 0.5) + (1.0 -
				ARTmortConstant * pow(Register[ID - 1].BaselineCD4, 0.5)) * pow(ARTmortRedn, duration));
		}
		MortProb = 1.0 - exp(-MortProb / CycleD);
		if (is < 5){
			temp2 = MortProb * ProbOIprecedesAIDSmort * ARTuptakeOI[iy];
			if (Register[ID - 1].VCThistory < 2){ temp2 *= OIsDiagnosed[iy]; }
			MortProb = MortProb - temp2;
		}
		if (FixedUncertainty == 1){
			ia = Register[ID - 1].CurrAge;
			if (ia >= 15){ 
				LYsLostExp.out[CurrSim - 1][iy] += MortProb * LE_West26[ia][SexInd] * Register[ID - 1].PopWeight;
				if (is >= 5){ ARTdeathsExp.out[CurrSim - 1][iy] += MortProb * Register[ID-1].PopWeight; }
			}
		}
	}
	if (is == 2 || is == 3 || is == 4 || is == 6){
		if (p<MortProb){ 
			RSApop.SetToDead(ID, 1); 
		}
		ExactAge = 1.0 * (CurrYear + 0.5 + (BehavCycleCount - 1.0) / CycleS - Register[ID - 1].DOB);
		temp1 = exp(-(ChangeCD4 + VLeffectCD4 * (Register[ID - 1].logVL - 4.0) +
			AgeEffectCD4 * (ExactAge - 35.0)) / CycleD);
		if (temp1 < 1.0){ Register[ID - 1].CD4 *= temp1; }
		temp1 = (ChangeVL + VLeffectVL * (Register[ID - 1].logVL - 4.0) + CD4effectVL *
			(Register[ID - 1].CD4 - 500.0) / 100.0 + AgeEffectVL * (ExactAge - 35.0)) / CycleD;
		if (temp1>0.0){ Register[ID - 1].logVL += temp1; }
	}
	if(is==2){
		if (Register[ID - 1].CD4<350){ Register[ID - 1].HIVstageE = 3; }
		else{Register[ID-1].HIVstageE = 2;}
	}
	if(is==3){
		if (Register[ID - 1].CD4<200){ Register[ID - 1].HIVstageE = 4; }
		else{Register[ID-1].HIVstageE = 3;}
	}
	if(is==4){Register[ID-1].HIVstageE = 4;}
	if(is==5){
		ia = Register[ID - 1].CurrAge - 10;
		if(p<MortProb){RSApop.SetToDead(ID, 1);}
		else if (p<MortProb + From5to6*AgeAdjInterrupt[ia][SexInd]){ Register[ID - 1].HIVstageE = 6; }
		else{Register[ID-1].HIVstageE = 5;}
	}
	if (is == 6){
		if (p<MortProb + From6to5){ Register[ID - 1].HIVstageE = 5; }
		else{ Register[ID - 1].HIVstageE = 6; }
	}

	// Determine whether ART is initiated
	if (Register[ID - 1].VCThistory == 2 && is < 5 && Register[ID - 1].AliveInd == 1){
		temp1 = ARTinitiationU200[iy][SexInd] / CycleD;
		if (is > 1){ temp1 *= pow(RR_ARTstartPer100CD4, (Register[ID - 1].CD4 - 100.0) / 100.0); }
		if (is == 4){ temp1 *= ARTeligibleGen[iy][0]; }
		if (is == 3){ temp1 *= ARTeligibleGen[iy][1]; }
		if (is == 2){ 
			if (Register[ID - 1].CD4 < 500){ temp1 *= ARTeligibleGen[iy][2]; }
			else{ temp1 *= ARTeligibleGen[iy][3]; }
		}
		if (is == 1){ temp1 = 0.0; }
		if (temp1 + temp2 > 0.0){
			temp1 = 1.0 - exp(-temp1);
			if (temp2 > temp1){ temp1 = temp2; }
			if (FixedUncertainty == 1){ NewARTexp.out[CurrSim - 1][iy] += temp1; }
			p = (p - MortProb) / (1.0 - MortProb);
			if (p < temp1){ 
				Register[ID - 1].HIVstageE = 5; 
				if (Register[ID - 1].VCThistory < 2){ Register[ID - 1].GetDiagnosed(ID, 1); }
				Register[ID - 1].DOAI = 1.0 * CurrYear + 0.5 + 1.0 * (BehavCycleCount - 1) / CycleS +
					1.0 * (STDcycleCount - 1) / CycleD;
				Register[ID - 1].BaselineCD4 = Register[ID - 1].CD4;
				if (FixedUncertainty == 1){
					if (is == 4){ NewART200.out[CurrSim - 1][CurrYear - StartYear] += Register[ID - 1].PopWeight; }
					if (is == 3){ NewART350.out[CurrSim - 1][CurrYear - StartYear] += Register[ID - 1].PopWeight; }
					if (is == 2){
						if (Register[ID - 1].CD4 < 500){ NewART500.out[CurrSim - 1][CurrYear - StartYear] += Register[ID - 1].PopWeight; }
						else{ NewART500plus.out[CurrSim - 1][CurrYear - StartYear] += Register[ID - 1].PopWeight; }
					}
				}
			}
		}
	}
}

void NonHIVtransition::CalcProbCure()
{
	double PropnTreatedPublic, PropnTreatedPrivate;
	int iy;
	
	if(SexInd==0){
		PropnTreatedPublic = PropnTreatedPublicM;
		PropnTreatedPrivate = PropnTreatedPrivateM;}
	else{
		PropnTreatedPublic = PropnTreatedPublicF;
		PropnTreatedPrivate = PropnTreatedPrivateF;}
	iy = CurrYear - StartYear;
	
	ProbCompleteCure = (PropnTreatedPublic * (PropnPublicUsingSM[iy] * CorrectRxWithSM +
		(1.0 - PropnPublicUsingSM[iy]) * CorrectRxPreSM) * (1.0 - DrugShortage[iy]) +
		PropnTreatedPrivate * (PropnPrivateUsingSM[iy] * CorrectRxWithSM +
		(1.0 - PropnPrivateUsingSM[iy]) * CorrectRxPreSM)) * DrugEff + TradnalEff *
		(1.0 - PropnTreatedPublic - PropnTreatedPrivate);
}

int NonHIVtransition::GetSTDstage(int offset, double r)
{
	int is, stage;
	double CumProb;

	CumProb = 0.0;
	stage = 0;
	for(is=0; is<nStates; is++){
		CumProb += PropnByStage[offset][is];
		if(r<CumProb){
			stage = is;
			break;
		}
	}

	return stage;
}

TradnalSTDtransition::TradnalSTDtransition(){}

SyphilisTransition::SyphilisTransition(int Sex, int ObsANC, int ObsFPC, int ObsGUD,
	int ObsCSW, int ObsHH, int ObsANCN, int ObsHHN)
{
	SexInd = Sex;
	nStates = 7;
	ImmuneState = 1;

	ANClogL.Observations = ObsANC;
	FPClogL.Observations = ObsFPC;
	GUDlogL.Observations = ObsGUD;
	CSWlogL.Observations = ObsCSW;
	HouseholdLogL.Observations = ObsHH;
	AntenatalNlogL.Observations = ObsANCN;
	HouseholdNlogL.Observations = ObsHHN;
}

void SyphilisTransition::CalcTransitionProbs()
{
	int ia;
	double RxRate, TeenRxRate, RednAsympDurFSW;
	double Adj2ndary; // Adjustment to prob of cure (for primary syphilis) in individuals 
					  // with 2ndary syphilis

	if(SexInd==0){
		RxRate = MaleRxRate;
		TeenRxRate = MaleTeenRxRate;
		RednAsympDurFSW = 0.0;
	}
	else{
		RxRate = FemRxRate;
		TeenRxRate = FemTeenRxRate;
		RednAsympDurFSW = 1.0/(1.0 + 1.0/(FSWasympRxRate * ProbCompleteCure * FSWasympCure *
			AveDuration[3]));
	}

	ANCpropnCured = ANCpropnScreened * ANCpropnTreated * DrugEff;

	From1to2 = 1.0 - exp(-(1.0/AveDuration[0]) * 52.0/CycleD);
	From2to3 = (1.0 - exp(-(1.0/AveDuration[1]) * 52.0/CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-RxRate * ProbCompleteCure * 52.0/CycleD)));
	From2to3T = (1.0 - exp(-(1.0/AveDuration[1]) * 52.0/CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-TeenRxRate * ProbCompleteCure * 52.0/CycleD)));
	From2to3C = (1.0 - exp(-(1.0/AveDuration[1]) * 52.0/CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-FSWRxRate * ProbCompleteCure * 52.0/CycleD)));
	From2to5 = (1.0 - exp(-RxRate * ProbCompleteCure * 52.0/CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-(1.0/AveDuration[1]) * 52.0/CycleD)));
	From2to0 = From2to5 * PropnSuscepAfterRx;
	From2to5 *= (1.0 - PropnSuscepAfterRx);
	From2to5T = (1.0 - exp(-TeenRxRate * ProbCompleteCure * 52.0/CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-(1.0/AveDuration[1]) * 52.0/CycleD)));
	From2to0T = From2to5T * PropnSuscepAfterRx;
	From2to5T *= (1.0 - PropnSuscepAfterRx);
	From2to5C = (1.0 - exp(-FSWRxRate * ProbCompleteCure * 52.0/CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-(1.0/AveDuration[1]) * 52.0/CycleD)));
	From2to0C = From2to5C * PropnSuscepAfterRx;
	From2to5C *= (1.0 - PropnSuscepAfterRx);
	Adj2ndary = (1.0 - SecondaryRxMult) * (1.0 - SecondaryCureMult);
	From3to4 = (1.0 - exp(-(1.0/AveDuration[2]) * 52.0/CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-RxRate * ProbCompleteCure * Adj2ndary * 52.0/CycleD)));
	From3to4T = (1.0 - exp(-(1.0/AveDuration[2]) * 52.0/CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-TeenRxRate * ProbCompleteCure * Adj2ndary * 52.0/CycleD)));
	From3to4C = (1.0 - exp(-(1.0/AveDuration[2]) * 52.0/CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-FSWRxRate * ProbCompleteCure * Adj2ndary * 52.0/CycleD)));
	From3to5 = (1.0 - exp(-RxRate * ProbCompleteCure * Adj2ndary * 52.0/CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-(1.0/AveDuration[2]) * 52.0/CycleD)));
	From3to5T = (1.0 - exp(-TeenRxRate * ProbCompleteCure * Adj2ndary * 52.0/CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-(1.0/AveDuration[2]) * 52.0/CycleD)));
	From3to5C = (1.0 - exp(-FSWRxRate * ProbCompleteCure * Adj2ndary * 52.0/CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-(1.0/AveDuration[2]) * 52.0/CycleD)));
	From4to6 = 1.0 - exp(-(1.0/AveDuration[3]) * 52.0/CycleD);
	From4to6C = 1.0 - exp(-(1.0/(AveDuration[3] * (1 - RednAsympDurFSW))) * 52.0/CycleD);
	From5to0 = 1.0 - exp(-(1.0/AveDuration[4]) * 52.0/CycleD);
	From6to0 = 1.0 - exp(-(1.0/AveDuration[5]) * 52.0/CycleD);
}

double SyphilisTransition::GetTransmProb(int ID)
{
	// ID is the ID of the susceptible partner, but the SyphilisTransition object is for
	// the opposite sex (the sex of the infected partner).

	int PID1, PID2, CSWID, is, ig, ia, ic;
	int IRisk, PRisk;
	double NoTransmProb, SingleActProb;

	IRisk = Register[ID-1].RiskGroup;
	NoTransmProb = 1.0;

	// Infection from primary partner
	if(Register[ID-1].CurrPartners>0){
		PID1 = Register[ID-1].IDprimary;
		if(Register[PID1-1].TPstage==2 || Register[PID1-1].TPstage==3){
			// Get base probability depending on partnership type
			SingleActProb = TransmProb;
			PRisk = Register[PID1-1].RiskGroup;
			if(Register[ID-1].MarriedInd==1){
				SingleActProb *= RelTransmLT;
				if(PRisk==2 && NoBacterialTransm12==1){
					SingleActProb = 0.0;}
				if(PRisk==2 && IRisk==2 && NoBacterialTransm22==1){
					SingleActProb = 0.0;}
			}
			// Make adjustment for age of individual
			ia = Register[ID-1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			if(SingleActProb>1.0){
				SingleActProb = 1.0;}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID-1].UVIprimary);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff), 
				Register[ID-1].PVIprimary);
		}
	}

	// Infection from 2ndary partner
	if(Register[ID-1].CurrPartners==2){
		PID2 = Register[ID-1].ID2ndary;
		if(Register[PID2-1].TPstage==2 || Register[PID2-1].TPstage==3){
			SingleActProb = TransmProb;
			// Make adjustment for age of individual
			ia = Register[ID-1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			if(SingleActProb>1.0){
				SingleActProb = 1.0;}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID-1].UVI2ndary);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff), 
				Register[ID-1].PVI2ndary);
		}
	}

	// Infection from CSW (relevant only to high-risk men)
	if(Register[ID-1].UVICSW + Register[ID-1].PVICSW > 0){
		PID2 = Register[ID-1].IDofCSW;
		if(Register[PID2-1].TPstage==2 || Register[PID2-1].TPstage==3){
			SingleActProb = TransmProb;
			// F-to-M transmission prob is assumed to be the same in commercial sex
			// Make adjustment for age of individual
			ia = Register[ID-1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			if(SingleActProb>1.0){
				SingleActProb = 1.0;}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID-1].UVICSW);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff), 
				Register[ID-1].PVICSW);
		}
	}

	// Infection from client (relevant only to CSWs)
	if(Register[ID-1].FSWind==1){
		ia = Register[ID-1].AgeGroup - 2;
		for(ic=0; ic<Register.size(); ic++){
			if(Register[ic].IDofCSW==ID && (Register[ic].UVICSW + Register[ic].PVICSW) > 0){
				if(Register[ic].TPstage==2 || Register[ic].TPstage==3){
					SingleActProb = TransmProbSW;
					SingleActProb *= SuscepIncrease[ia];
					if(SingleActProb>1.0){
						SingleActProb = 1.0;}
					NoTransmProb *= pow(1.0 - SingleActProb, Register[ic].UVICSW);
					NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff), 
						Register[ic].PVICSW);
				}
			}
		}
	}

	// Infection from casual partner 
	if ((Register[ID - 1].CasualInd == 1 || Register[ID - 1].HetCasualInd == 1) && 
		(Register[ID - 1].UAIcasual + Register[ID - 1].PAIcasual) > 0){
		PID2 = Register[ID - 1].IDofCasual;
		if (Register[PID2 - 1].TPstage == 2 || Register[PID2 - 1].TPstage == 3){
			SingleActProb = TransmProb;
			ia = Register[ID - 1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			if (SingleActProb>1.0){SingleActProb = 1.0;}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID - 1].UAIcasual);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
				Register[ID - 1].PAIcasual);
		}
	}

	return 1.0 - NoTransmProb;
}

void SyphilisTransition::GetNewStage(int ID, double p)
{
	int is;

	is = Register[ID-1].TPstage;
	if(is==1){
		if(p<From1to2){Register[ID-1].TPstageE = 2;}
		else{Register[ID-1].TPstageE = 1;}
	}
	if(is==2){
		if (Register[ID - 1].FSWind == 1){
			if(p<From2to0C){Register[ID-1].TPstageE = 0;}
			else if(p<From2to0C + From2to3C){Register[ID-1].TPstageE = 3;}
			else if(p<From2to0C + From2to3C + From2to5C){
				Register[ID-1].TPstageE = 5;}
			else{Register[ID-1].TPstageE = 2;}
		}
		else if (Register[ID - 1].AgeGroup<4){
			if(p<From2to0T){Register[ID-1].TPstageE = 0;}
			else if(p<From2to0T + From2to3T){Register[ID-1].TPstageE = 3;}
			else if(p<From2to0T + From2to3T + From2to5T){
				Register[ID-1].TPstageE = 5;}
			else{Register[ID-1].TPstageE = 2;}
		}
		else{
			if(p<From2to0){Register[ID-1].TPstageE = 0;}
			else if(p<From2to0 + From2to3){Register[ID-1].TPstageE = 3;}
			else if(p<From2to0 + From2to3 + From2to5){
				Register[ID-1].TPstageE = 5;}
			else{Register[ID-1].TPstageE = 2;}
		}
	}
	if(is==3){
		if (Register[ID - 1].FSWind == 1){
			if(p<From3to4C){Register[ID-1].TPstageE = 4;}
			else if(p<From3to4C + From3to5C){Register[ID-1].TPstageE = 5;}
			else{Register[ID-1].TPstageE = 3;}
		}
		else if (Register[ID - 1].AgeGroup<4){
			if(p<From3to4T){Register[ID-1].TPstageE = 4;}
			else if(p<From3to4T + From3to5T){Register[ID-1].TPstageE = 5;}
			else{Register[ID-1].TPstageE = 3;}
		}
		else{
			if(p<From3to4){Register[ID-1].TPstageE = 4;}
			else if(p<From3to4 + From3to5){Register[ID-1].TPstageE = 5;}
			else{Register[ID-1].TPstageE = 3;}
		}
	}
	if(is==4){
		if(Register[ID-1].FSWind==1){
			if(p<From4to6C){Register[ID-1].TPstageE = 6;}
			else{Register[ID-1].TPstageE = 4;}
		}
		else{
			if(p<From4to6){Register[ID-1].TPstageE = 6;}
			else{Register[ID-1].TPstageE = 4;}
		}
	}
	if(is==5){
		if(p<From5to0){Register[ID-1].TPstageE = 0;}
		else{Register[ID-1].TPstageE = 5;}
	}
	if(is==6){
		if(p<From6to0){Register[ID-1].TPstageE = 0;}
		else{Register[ID-1].TPstageE = 6;}
	}
}

HerpesTransition::HerpesTransition(int Sex, int ObsANC, int ObsFPC, int ObsGUD, int ObsCSW,
	int ObsHH, int ObsANCN, int ObsHHN)
{
	SexInd = Sex;
	nStates = 5;

	ANClogL.Observations = ObsANC;
	FPClogL.Observations = ObsFPC;
	GUDlogL.Observations = ObsGUD;
	CSWlogL.Observations = ObsCSW;
	HouseholdLogL.Observations = ObsHH;
	AntenatalNlogL.Observations = ObsANCN;
	HouseholdNlogL.Observations = ObsHHN;
}

void HerpesTransition::CalcTransitionProbs()
{
	int is;
	double RxRate, TeenRxRate;

	if(SexInd==0){
		RxRate = MaleRxRate;
		TeenRxRate = MaleTeenRxRate;}
	else{
		RxRate = FemRxRate;
		TeenRxRate = FemTeenRxRate;}

	From1to2 = 1.0 - exp(-((1.0/AveDuration[0]) + RxRate * ProbCompleteCure)*52.0/CycleD);
	From1to2T = 1.0 - exp(-((1.0/AveDuration[0]) + TeenRxRate * ProbCompleteCure)*
		52.0/CycleD);
	From2to3[0] = (1.0 - exp(-RecurrenceRate*52.0/CycleD)) *
		(1.0 - 0.5*(1.0 - exp(-(1.0/AveDuration[1])*52.0/CycleD)));
	for(is=1; is<7; is++){
		From2to3[is] = 1.0 - exp(-RecurrenceRate * (1.0 + HSVrecurrenceIncrease[is-1])*
			52.0/CycleD);}
	From3to2 = 1.0 - exp(-((1.0/AveDuration[2]) + RxRate * ProbCompleteCure)*52.0/CycleD);
	From3to2T = 1.0 - exp(-((1.0/AveDuration[2]) + TeenRxRate * ProbCompleteCure)*
		52.0/CycleD);
	From2to4 = (1 - exp(-(1.0/AveDuration[1])*52.0/CycleD)) *
		(1.0 - 0.5*(1.0 - exp(-RecurrenceRate*52.0/CycleD)));
	if(SexInd==1){
		From1to2C = 1.0 - exp(-((1.0/AveDuration[0]) + FSWRxRate * ProbCompleteCure)*
			52.0/CycleD);
		From3to2C = 1.0 - exp(-((1.0/AveDuration[2]) + FSWRxRate * ProbCompleteCure)*
			52.0/CycleD);}
}

double HerpesTransition::GetTransmProb(int ID)
{
	// ID is the ID of the suscpetible partner, but the HerpesTransition object is for
	// the opposite sex (the sex of the infected partner).

	int PID1, PID2, CSWID, is, ig, ia, ic;
	int IRisk, PRisk;
	double NoTransmProb, SingleActProb;

	IRisk = Register[ID-1].RiskGroup;
	NoTransmProb = 1.0;

	// Infection from primary partner
	if(Register[ID-1].CurrPartners>0){
		PID1 = Register[ID-1].IDprimary;
		if(Register[PID1-1].HSVstage>0){
			// Get base probability depending on partnership type
			SingleActProb = TransmProb;
			PRisk = Register[PID1-1].RiskGroup;
			if(Register[ID-1].MarriedInd==1){
				SingleActProb *= RelTransmLT;
				if(PRisk==2 && NoViralTransm12==1){
					SingleActProb = 0.0;}
				if(PRisk==2 && IRisk==2 && NoViralTransm22==1){
					SingleActProb = 0.0;}
			}
			// Make adjustment for higher HSV-2 shedding when symptomatic
			if(Register[PID1-1].HSVstage==1 || Register[PID1-1].HSVstage==3){
				SingleActProb *= HSVsymptomInfecIncrease;}
			// Make adjustment for higher HSV-2 shedding if HIV-positive
			if(Register[PID1-1].HIVstage>0){
				is = Register[PID1-1].HIVstage;
				SingleActProb *= (1.0 + HSVsheddingIncrease[is-1]);
			}
			// Make adjustment for age of individual
			ia = Register[ID-1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			if(SingleActProb>1.0){
				SingleActProb = 1.0;}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID-1].UVIprimary);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff), 
				Register[ID-1].PVIprimary);
		}
	}

	// Infection from 2ndary partner
	if(Register[ID-1].CurrPartners==2){
		PID2 = Register[ID-1].ID2ndary;
		if(Register[PID2-1].HSVstage>0){
			SingleActProb = TransmProb;
			// Make adjustment for higher HSV-2 shedding when symptomatic
			if(Register[PID2-1].HSVstage==1 || Register[PID2-1].HSVstage==3){
				SingleActProb *= HSVsymptomInfecIncrease;}
			// Make adjustment for higher HSV-2 shedding if HIV-positive
			if(Register[PID2-1].HIVstage>0){
				is = Register[PID2-1].HIVstage;
				SingleActProb *= (1.0 + HSVsheddingIncrease[is-1]);
			}
			// Make adjustment for age of individual
			ia = Register[ID-1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			if(SingleActProb>1.0){
				SingleActProb = 1.0;}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID-1].UVI2ndary);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff), 
				Register[ID-1].PVI2ndary);
		}
	}

	// Infection from CSW (relevant only to high-risk men)
	if(Register[ID-1].UVICSW + Register[ID-1].PVICSW > 0){
		PID2 = Register[ID-1].IDofCSW;
		if(Register[PID2-1].HSVstage>0){
			SingleActProb = TransmProb;
			// F-to-M transmission prob is assumed to be the same in commercial sex
			// Make adjustment for higher HSV-2 shedding when symptomatic
			if(Register[PID2-1].HSVstage==1 || Register[PID2-1].HSVstage==3){
				SingleActProb *= HSVsymptomInfecIncrease;}
			// Make adjustment for higher HSV-2 shedding if HIV-positive
			if(Register[PID2-1].HIVstage>0){
				is = Register[PID2-1].HIVstage;
				SingleActProb *= (1.0 + HSVsheddingIncrease[is-1]);
			}
			// Make adjustment for age of individual
			ia = Register[ID-1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			if(SingleActProb>1.0){
				SingleActProb = 1.0;}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID-1].UVICSW);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff), 
				Register[ID-1].PVICSW);
		}
	}

	// Infection from client (relevant only to CSWs)
	if(Register[ID-1].FSWind==1){
		ia = Register[ID-1].AgeGroup - 2;
		for(ic=0; ic<Register.size(); ic++){
			if(Register[ic].IDofCSW==ID && (Register[ic].UVICSW + Register[ic].PVICSW) > 0){
				if(Register[ic].HSVstage>0){
					SingleActProb = TransmProbSW;
					// Make adjustment for higher HSV-2 shedding when symptomatic
					if(Register[ic].HSVstage==1 || Register[ic].HSVstage==3){
						SingleActProb *= HSVsymptomInfecIncrease;}
					// Make adjustment for higher HSV-2 shedding if HIV-positive
					if(Register[ic].HIVstage>0){
						is = Register[ic].HIVstage;
						SingleActProb *= (1.0 + HSVsheddingIncrease[is-1]);
					}
					// Make adjustment for age of individual
					SingleActProb *= SuscepIncrease[ia];
					if(SingleActProb>1.0){
						SingleActProb = 1.0;}
					NoTransmProb *= pow(1.0 - SingleActProb, Register[ic].UVICSW);
					NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff), 
						Register[ic].PVICSW);
				}
			}
		}
	}

	// Infection from casual partner 
	if ((Register[ID - 1].CasualInd == 1 || Register[ID - 1].HetCasualInd == 1) && 
		(Register[ID - 1].UAIcasual + Register[ID - 1].PAIcasual) > 0){
		PID2 = Register[ID - 1].IDofCasual;
		if (Register[PID2 - 1].HSVstage > 0){
			SingleActProb = TransmProb;
			// Make adjustment for higher HSV-2 shedding when symptomatic
			if (Register[PID2 - 1].HSVstage == 1 || Register[PID2 - 1].HSVstage == 3){
				SingleActProb *= HSVsymptomInfecIncrease;}
			// Make adjustment for higher HSV-2 shedding if HIV-positive
			if (Register[PID2 - 1].HIVstage>0){
				is = Register[PID2 - 1].HIVstage;
				SingleActProb *= (1.0 + HSVsheddingIncrease[is - 1]);
			}
			ia = Register[ID - 1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			if (SingleActProb>1.0){ SingleActProb = 1.0; }
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID - 1].UAIcasual);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
				Register[ID - 1].PAIcasual);
		}
	}

	return 1.0 - NoTransmProb;
}

void HerpesTransition::GetNewStage(int ID, double p)
{
	int is, is2;

	is = Register[ID-1].HSVstage;
	if(is==1){
		if (Register[ID - 1].FSWind == 1){
			if(p<From1to2C){Register[ID-1].HSVstageE = 2;}
			else{Register[ID-1].HSVstageE = 1;}
		}
		else if (Register[ID - 1].AgeGroup<4){
			if(p<From1to2T){Register[ID-1].HSVstageE = 2;}
			else{Register[ID-1].HSVstageE = 1;}
		}
		else{
			if(p<From1to2){Register[ID-1].HSVstageE = 2;}
			else{Register[ID-1].HSVstageE = 1;}
		}
	}
	if(is==2){
		is2 = Register[ID-1].HIVstage;
		if(p<From2to3[is2]){Register[ID-1].HSVstageE = 3;}
		else if(is2==0 && (p<From2to3[0]+From2to4)){
			Register[ID-1].HSVstageE = 4;}
		else{
			Register[ID-1].HSVstageE = 2;}
	}
	if(is==3){
		if (Register[ID - 1].FSWind == 1){
			if(p<From3to2C){Register[ID-1].HSVstageE = 2;}
			else{Register[ID-1].HSVstageE = 3;}
		}
		else if (Register[ID - 1].AgeGroup<4){
			if(p<From3to2T){Register[ID-1].HSVstageE = 2;}
			else{Register[ID-1].HSVstageE = 3;}
		}
		else{
			if(p<From3to2){Register[ID-1].HSVstageE = 2;}
			else{Register[ID-1].HSVstageE = 3;}
		}
	}
	if(is==4){
		Register[ID-1].HSVstageE = 4;}
}

OtherSTDtransition::OtherSTDtransition(int Sex, int ObsANC, int ObsFPC, int ObsGUD,
	int ObsCSW, int ObsHH, int ObsANCN, int ObsHHN)
{
	SexInd = Sex;
	nStates = 4;
	ImmuneState = 3;

	ANClogL.Observations = ObsANC;
	FPClogL.Observations = ObsFPC;
	GUDlogL.Observations = ObsGUD;
	CSWlogL.Observations = ObsCSW;
	HouseholdLogL.Observations = ObsHH;
	AntenatalNlogL.Observations = ObsANCN;
	HouseholdNlogL.Observations = ObsHHN;
}

void OtherSTDtransition::CalcTransitionProbs()
{
	double RxRate, TeenRxRate, RednAsympDurFSW;

	if(SexInd==0){
		RxRate = MaleRxRate;
		TeenRxRate = MaleTeenRxRate;}
	else{
		RxRate = FemRxRate;
		TeenRxRate = FemTeenRxRate;
		//RednAsympDurFSW = 1.0/(1.0 + 1.0/(FSWasympRxRate * ProbCompleteCure * FSWasympCure *
		//	AveDuration[1]));
	}

	From1to0 = 1.0 - exp(-(((1.0 - PropnImmuneAfterSR)/AveDuration[0]) + RxRate * 
		ProbCompleteCure * (1.0 - PropnImmuneAfterRx))*52.0/CycleD);
	From1to0T = 1.0 - exp(-(((1.0 - PropnImmuneAfterSR)/AveDuration[0]) + TeenRxRate * 
		ProbCompleteCure * (1.0 - PropnImmuneAfterRx))*52.0/CycleD);
	From2to0 = 1.0 - exp(-((1.0 - PropnImmuneAfterSR)/AveDuration[1])*52.0/CycleD);
	From1to3 = 1.0 - exp(-((PropnImmuneAfterSR/AveDuration[0]) + RxRate * 
		ProbCompleteCure * PropnImmuneAfterRx)*52.0/CycleD);
	From1to3T = 1.0 - exp(-((PropnImmuneAfterSR/AveDuration[0]) + TeenRxRate * 
		ProbCompleteCure * PropnImmuneAfterRx)*52.0/CycleD);
	From2to3 = 1.0 - exp(-(PropnImmuneAfterSR/AveDuration[1])*52.0/CycleD);
	From3to0 = 1.0 - exp(-(1.0/AveDuration[2])*52.0/CycleD);
	// Transition probabilities for sex workers
	if(SexInd==1){
		From1to0C = 1.0 - exp(-(((1.0 - PropnImmuneAfterSR)/AveDuration[0]) + FSWRxRate * 
			ProbCompleteCure * (1.0 - PropnImmuneAfterRx))*52.0/CycleD);
		From2to0C = 1.0 - exp(-(((1.0 - PropnImmuneAfterSR)/AveDuration[1]) + FSWasympRxRate *
			ProbCompleteCure * FSWasympCure * (1.0 - PropnImmuneAfterRx))*52.0/CycleD);
		From1to3C = 1.0 - exp(-((PropnImmuneAfterSR/AveDuration[0]) + FSWRxRate * 
			ProbCompleteCure * PropnImmuneAfterRx)*52.0/CycleD);
		From2to3C = 1.0 - exp(-((PropnImmuneAfterSR/AveDuration[1]) + FSWasympRxRate *
			ProbCompleteCure * FSWasympCure * PropnImmuneAfterRx)*52.0/CycleD);
	}
}

double OtherSTDtransition::GetTransmProb(int ID, int STD)
{
	// ID is the ID of the suscpetible partner, but the OtherSTDtransition object is for
	// the opposite sex (the sex of the infected partner).
	// STD indicator is 1 for NG, 2 for CT, 3 for TV and 4 for HD

	int PID1, PID2, CSWID, is, ig, ia, ic;
	int IRisk, PRisk;
	double NoTransmProb, SingleActProb;

	IRisk = Register[ID-1].RiskGroup;
	NoTransmProb = 1.0;

	// Infection from primary partner
	if(Register[ID-1].CurrPartners>0){
		PID1 = Register[ID-1].IDprimary;
		if((STD==1 && Register[PID1-1].NGstage>0 && Register[PID1-1].NGstage<3) ||
			(STD==2 && Register[PID1-1].CTstage>0 && Register[PID1-1].CTstage<3) ||
			(STD==3 && Register[PID1-1].TVstage>0 && Register[PID1-1].TVstage<3) ||
			(STD==4 && Register[PID1-1].HDstage>0 && Register[PID1-1].HDstage<3)){
			// Get base probability depending on partnership type
			SingleActProb = TransmProb;
			PRisk = Register[PID1-1].RiskGroup;
			if(Register[ID-1].MarriedInd==1){
				SingleActProb *= RelTransmLT;
				if(PRisk==2 && NoBacterialTransm12==1){
					SingleActProb = 0.0;}
				if(PRisk==2 && IRisk==2 && NoBacterialTransm22==1){
					SingleActProb = 0.0;}
			}
			// Make adjustment for age of individual
			ia = Register[ID-1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			if(SingleActProb>1.0){
				SingleActProb = 1.0;}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID-1].UVIprimary);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff), 
				Register[ID-1].PVIprimary);
		}
	}

	// Infection from 2ndary partner
	if(Register[ID-1].CurrPartners==2){
		PID2 = Register[ID-1].ID2ndary;
		if((STD==1 && Register[PID2-1].NGstage>0 && Register[PID2-1].NGstage<3) ||
			(STD==2 && Register[PID2-1].CTstage>0 && Register[PID2-1].CTstage<3) ||
			(STD==3 && Register[PID2-1].TVstage>0 && Register[PID2-1].TVstage<3) ||
			(STD==4 && Register[PID2-1].HDstage>0 && Register[PID2-1].HDstage<3)){
			SingleActProb = TransmProb;
			// Make adjustment for age of individual
			ia = Register[ID-1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			if(SingleActProb>1.0){
				SingleActProb = 1.0;}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID-1].UVI2ndary);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff), 
				Register[ID-1].PVI2ndary);
		}
	}

	// Infection from CSW (relevant only to high-risk men)
	if(Register[ID-1].UVICSW + Register[ID-1].PVICSW > 0){
		PID2 = Register[ID-1].IDofCSW;
		if((STD==1 && Register[PID2-1].NGstage>0 && Register[PID2-1].NGstage<3) ||
			(STD==2 && Register[PID2-1].CTstage>0 && Register[PID2-1].CTstage<3) ||
			(STD==3 && Register[PID2-1].TVstage>0 && Register[PID2-1].TVstage<3) ||
			(STD==4 && Register[PID2-1].HDstage>0 && Register[PID2-1].HDstage<3)){
			SingleActProb = TransmProb;
			// F-to-M transmission prob is assumed to be the same in commercial sex
			// Make adjustment for age of individual
			ia = Register[ID-1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			if(SingleActProb>1.0){
				SingleActProb = 1.0;}
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID-1].UVICSW);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff), 
				Register[ID-1].PVICSW);
		}
	}

	// Infection from client (relevant only to CSWs)
	if(Register[ID-1].FSWind==1){
		ia = Register[ID-1].AgeGroup - 2;
		for(ic=0; ic<Register.size(); ic++){
			if(Register[ic].IDofCSW==ID && (Register[ic].UVICSW + Register[ic].PVICSW) > 0){
				if((STD==1 && Register[ic].NGstage>0 && Register[ic].NGstage<3) ||
					(STD==2 && Register[ic].CTstage>0 && Register[ic].CTstage<3) ||
					(STD==3 && Register[ic].TVstage>0 && Register[ic].TVstage<3) ||
					(STD==4 && Register[ic].HDstage>0 && Register[ic].HDstage<3)){
					SingleActProb = TransmProbSW;
					SingleActProb *= SuscepIncrease[ia];
					if(SingleActProb>1.0){
						SingleActProb = 1.0;}
					NoTransmProb *= pow(1.0 - SingleActProb, Register[ic].UVICSW);
					NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff), 
						Register[ic].PVICSW);
				}
			}
		}
	}

	// Infection from casual partner 
	if ((Register[ID - 1].CasualInd == 1 || Register[ID - 1].HetCasualInd == 1) && 
		(Register[ID - 1].UAIcasual + Register[ID - 1].PAIcasual) > 0){
		PID2 = Register[ID - 1].IDofCasual;
		if ((STD == 1 && Register[PID2 - 1].NGstage>0 && Register[PID2 - 1].NGstage<3) ||
			(STD == 2 && Register[PID2 - 1].CTstage>0 && Register[PID2 - 1].CTstage<3) ||
			(STD == 3 && Register[PID2 - 1].TVstage>0 && Register[PID2 - 1].TVstage<3) ||
			(STD == 4 && Register[PID2 - 1].HDstage>0 && Register[PID2 - 1].HDstage<3)){
			SingleActProb = TransmProb;
			ia = Register[ID - 1].AgeGroup - 2;
			SingleActProb *= SuscepIncrease[ia];
			if (SingleActProb>1.0){ SingleActProb = 1.0; }
			NoTransmProb *= pow(1.0 - SingleActProb, Register[ID - 1].UAIcasual);
			NoTransmProb *= pow(1.0 - SingleActProb * (1.0 - CondomEff),
				Register[ID - 1].PAIcasual);
		}
	}

	return 1.0 - NoTransmProb;
}

void OtherSTDtransition::GetNewStage(int ID, double p, int STD)
{
	int StartStage, EndStage;

	if(STD==1){StartStage = Register[ID-1].NGstage;}
	if(STD==2){StartStage = Register[ID-1].CTstage;}
	if(STD==3){StartStage = Register[ID-1].TVstage;}
	if(STD==4){StartStage = Register[ID-1].HDstage;}

	if(StartStage==1){
		if (Register[ID - 1].FSWind == 1){
			if(p<From1to0C){EndStage = 0;}
			else if(p<From1to0C+From1to3C){EndStage = 3;}
			else{EndStage = 1;}
		}
		else if (Register[ID - 1].AgeGroup<4){
			if(p<From1to0T){EndStage = 0;}
			else if(p<From1to0T+From1to3T){EndStage = 3;}
			else{EndStage = 1;}
		}
		else{
			if(p<From1to0){EndStage = 0;}
			else if(p<From1to0+From1to3){EndStage = 3;}
			else{EndStage = 1;}
		}
	}
	if(StartStage==2){
		if(Register[ID-1].FSWind==1){
			if(p<From2to0C){EndStage = 0;}
			else if(p<From2to0C+From2to3C){EndStage = 3;}
			else{EndStage = 2;}
		}
		else{
			if(p<From2to0){EndStage = 0;}
			else if(p<From2to0+From2to3){EndStage = 3;}
			else{EndStage = 2;}
		}
	}
	if(StartStage==3){
		if(p<From3to0){EndStage = 0;}
		else{EndStage = 3;}
	}

	if(STD==1){Register[ID-1].NGstageE = EndStage;}
	if(STD==2){Register[ID-1].CTstageE = EndStage;}
	if(STD==3){Register[ID-1].TVstageE = EndStage;}
	if(STD==4){Register[ID-1].HDstageE = EndStage;}
}

NonSTDtransition::NonSTDtransition(){}

void NonSTDtransition::CalcProbPartialCure()
{
	int iy;
	
	iy = CurrYear - StartYear;
	
	ProbPartialCure = (PropnTreatedPublicF * (PropnPublicUsingSM[iy] * CorrectRxWithSM +
		(1.0 - PropnPublicUsingSM[iy]) * CorrectRxPreSM) * (1.0 - DrugShortage[iy]) +
		PropnTreatedPrivateF * (PropnPrivateUsingSM[iy] * CorrectRxWithSM +
		(1.0 - PropnPrivateUsingSM[iy]) * CorrectRxPreSM)) * DrugPartialEff + TradnalEff *
		(1.0 - PropnTreatedPublicF - PropnTreatedPrivateF);
}

BVtransition::BVtransition(int Sex, int ObsANC, int ObsFPC, int ObsGUD, int ObsCSW, int ObsHH,
	int ObsANCN, int ObsHHN)
{
	SexInd = Sex;
	nStates = 4;

	ANClogL.Observations = ObsANC;
	FPClogL.Observations = ObsFPC;
	GUDlogL.Observations = ObsGUD;
	CSWlogL.Observations = ObsCSW;
	HouseholdLogL.Observations = ObsHH;
	AntenatalNlogL.Observations = ObsANCN;
	HouseholdNlogL.Observations = ObsHHN;
}

void BVtransition::CalcTransitionProbs()
{
	int ij;
	double RednAsympDurFSW;

	RednAsympDurFSW = 1.0/(1.0 + 1.0/(FSWasympRxRate * (ProbCompleteCure + ProbPartialCure) * 
		FSWasympCure * AveDuration[3]));
	From1to2 = 1.0 - exp(-CtsTransition[0][1] * 52.0/CycleD);
	From2to1ind = 1.0 - exp(-CtsTransition[1][0] * 52.0/CycleD);
	From3to1 = (1.0 - exp(-(CtsTransition[2][0] + FemRxRate * ProbCompleteCure) * 52.0/
		CycleD)) * (1.0 - 0.5*(1.0 - exp(-(CtsTransition[2][1] + FemRxRate * 
		ProbPartialCure) * 52.0/CycleD)));
	From3to1T = (1.0 - exp(-(CtsTransition[2][0] + FemTeenRxRate * ProbCompleteCure) * 52.0/
		CycleD)) * (1.0 - 0.5*(1.0 - exp(-(CtsTransition[2][1] + FemTeenRxRate * 
		ProbPartialCure) * 52.0/CycleD)));
	From3to1C = (1.0 - exp(-(CtsTransition[2][0] + FSWRxRate * ProbCompleteCure) * 52.0/
		CycleD)) * (1.0 - 0.5*(1.0 - exp(-(CtsTransition[2][1] + FSWRxRate * 
		ProbPartialCure) * 52.0/CycleD)));
	From3to2 = (1.0 - exp(-(CtsTransition[2][1] + FemRxRate * ProbPartialCure) * 52.0/
		CycleD)) * (1.0 - 0.5*(1.0 - exp(-(CtsTransition[2][0] + FemRxRate * 
		ProbCompleteCure) * 52.0/CycleD)));
	From3to2 = (1.0 - exp(-(CtsTransition[2][1] + FemRxRate * ProbPartialCure) * 52.0/
		CycleD)) * (1.0 - 0.5*(1.0 - exp(-(CtsTransition[2][0] + FemRxRate * 
		ProbCompleteCure) * 52.0/CycleD)));
	From3to2T = (1.0 - exp(-(CtsTransition[2][1] + FemTeenRxRate * ProbPartialCure) * 52.0/
		CycleD)) * (1.0 - 0.5*(1.0 - exp(-(CtsTransition[2][0] + FemTeenRxRate * 
		ProbCompleteCure) * 52.0/CycleD)));
	From3to2C = (1.0 - exp(-(CtsTransition[2][1] + FSWRxRate * ProbPartialCure) * 52.0/
		CycleD)) * (1.0 - 0.5*(1.0 - exp(-(CtsTransition[2][0] + FSWRxRate * 
		ProbCompleteCure) * 52.0/CycleD)));
	From4to1 = (1.0 - exp(-CtsTransition[3][0] * 52.0/CycleD)) *
		(1.0 - 0.5 * (1.0 - exp(-CtsTransition[3][1] * 52.0/CycleD)));
	From4to1C = (1.0 - exp(-(CtsTransition[3][0]/(1.0 - RednAsympDurFSW))*52.0/CycleD)) *
		(1.0 - 0.5 * (1.0 - exp(-(CtsTransition[3][1]/(1.0 - RednAsympDurFSW))*52.0/CycleD)));
	From4to2 = (1.0 - exp(-CtsTransition[3][1] * 52.0/CycleD)) *
		(1.0 - 0.5 * (1.0 - exp(-CtsTransition[3][0] * 52.0/CycleD)));
	From4to2C = (1.0 - exp(-(CtsTransition[3][1]/(1.0 - RednAsympDurFSW))*52.0/CycleD)) *
		(1.0 - 0.5 * (1.0 - exp(-(CtsTransition[3][0]/(1.0 - RednAsympDurFSW))*52.0/CycleD)));

	From2to3ind[0] = 1.0 - exp(-CtsTransition[1][2] * (1.0 - IncidenceMultNoPartners)*52.0/
		CycleD);
	From2to3ind[1] = 1.0 - exp(-CtsTransition[1][2] * 52.0/CycleD);
	From2to3ind[2] = 1.0 - exp(-CtsTransition[1][2] * (1.0 + IncidenceMultTwoPartners)*52.0/
		CycleD);
	From2to4ind[0] = 1.0 - exp(-CtsTransition[1][3] * (1.0 - IncidenceMultNoPartners)*52.0/
		CycleD);
	From2to4ind[1] = 1.0 - exp(-CtsTransition[1][3] * 52.0/CycleD);
	From2to4ind[2] = 1.0 - exp(-CtsTransition[1][3] * (1.0 + IncidenceMultTwoPartners)*52.0/
		CycleD);
	for(ij=0; ij<3; ij++){
		From2to1dep[ij] = From2to1ind * (1.0 - 0.5 * (From2to3ind[ij] + From2to4ind[ij]) +
			From2to3ind[ij] * From2to4ind[ij]/3.0);
		From2to3dep[ij] = From2to3ind[ij] * (1.0 - 0.5 * (From2to1ind + From2to4ind[ij]) +
			From2to1ind * From2to4ind[ij]/3.0);
		From2to4dep[ij] = From2to4ind[ij] * (1.0 - 0.5 * (From2to1ind + From2to3ind[ij]) +
			From2to1ind * From2to3ind[ij]/3.0);
	}
}

void BVtransition::GetNewStage(int ID, double p)
{
	int partners;

	if(Register[ID-1].BVstage==0){
		if(p<From1to2){Register[ID-1].BVstageE = 1;}
		else{Register[ID-1].BVstageE = 0;}
	}
	if(Register[ID-1].BVstage==1){
		partners = Register[ID-1].CurrPartners;
		if(p<From2to1dep[partners]){Register[ID-1].BVstageE = 0;}
		else if(p<From2to1dep[partners] + From2to3dep[partners]){
			Register[ID-1].BVstageE = 2;}
		else if(p<From2to1dep[partners] + From2to3dep[partners] + From2to4dep[partners]){
			Register[ID-1].BVstageE = 3;}
		else{Register[ID-1].BVstageE = 1;}
	}
	if(Register[ID-1].BVstage==2){
		if (Register[ID - 1].FSWind == 1){
			if(p<From3to1C){Register[ID-1].BVstageE = 0;}
			else if(p<From3to1C+From3to2C){Register[ID-1].BVstageE = 1;}
			else{Register[ID-1].BVstageE = 2;}
		}
		else if (Register[ID - 1].AgeGroup<4){
			if(p<From3to1T){Register[ID-1].BVstageE = 0;}
			else if(p<From3to1T+From3to2T){Register[ID-1].BVstageE = 1;}
			else{Register[ID-1].BVstageE = 2;}
		}
		else{
			if(p<From3to1){Register[ID-1].BVstageE = 0;}
			else if(p<From3to1+From3to2){Register[ID-1].BVstageE = 1;}
			else{Register[ID-1].BVstageE = 2;}
		}
	}
	if(Register[ID-1].BVstage==3){
		if(Register[ID-1].FSWind==1){
			if(p<From4to1C){Register[ID-1].BVstageE = 0;}
			else if(p<From4to1C+From4to2C){Register[ID-1].BVstageE = 1;}
			else{Register[ID-1].BVstageE = 3;}
		}
		else{
			if(p<From4to1){Register[ID-1].BVstageE = 0;}
			else if(p<From4to1+From4to2){Register[ID-1].BVstageE = 1;}
			else{Register[ID-1].BVstageE = 3;}
		}
	}
}

VCtransition::VCtransition(int Sex, int ObsANC, int ObsFPC, int ObsGUD, int ObsCSW, int ObsHH,
	int ObsANCN, int ObsHHN)
{
	SexInd = Sex;
	nStates = 3;

	ANClogL.Observations = ObsANC;
	FPClogL.Observations = ObsFPC;
	GUDlogL.Observations = ObsGUD;
	CSWlogL.Observations = ObsCSW;
	HouseholdLogL.Observations = ObsHH;
	AntenatalNlogL.Observations = ObsANCN;
	HouseholdNlogL.Observations = ObsHHN;
}

void VCtransition::CalcTransitionProbs()
{
	int ia, is;
	double RednAsympDurFSW; // Although it's not strictly necessary to include this (since
							// we set it to 0), we may want to change this in future, so I
							// have included it.

	RednAsympDurFSW = 0.0;
	From1to2 = (1.0 - exp(-RecurrenceRate * 52.0/CycleD)) * 
		(1.0 - 0.5 * (1.0 - exp(-(1.0/AveDuration[0]) * 52.0/CycleD)));
	From1to2C = (1.0 - exp(-RecurrenceRate * 52.0/CycleD)) * (1.0 - 0.5 *
		(1.0 - exp(-(1.0/(AveDuration[0] * (1.0 - RednAsympDurFSW))) * 52.0/CycleD)));
	From1to0 = (1.0 - exp(-(1.0/AveDuration[0]) * 52.0/CycleD)) *
		(1.0 - 0.5 * (1.0 - exp(-RecurrenceRate * 52.0/CycleD)));
	From1to0C = (1.0 - exp(-(1.0/(AveDuration[0] * (1.0 - RednAsympDurFSW))) * 52.0/CycleD))*
		(1.0 - 0.5 * (1.0 - exp(-RecurrenceRate * 52.0/CycleD)));
	From2to1 = (1.0 - exp(-((1.0/AveDuration[1]) + FemRxRate * ProbPartialCure)*52.0/
		CycleD))*(1.0 - 0.5 * (1.0 - exp(-FemRxRate * ProbCompleteCure * 52.0/CycleD)));
	From2to1T = (1.0 - exp(-((1.0/AveDuration[1]) + FemTeenRxRate * ProbPartialCure)*52.0/
		CycleD))*(1.0 - 0.5 * (1.0 - exp(-FemTeenRxRate * ProbCompleteCure * 52.0/CycleD)));
	From2to1C = (1.0 - exp(-((1.0/AveDuration[1]) + FSWRxRate * ProbPartialCure)*52.0/
		CycleD))*(1.0 - 0.5 * (1.0 - exp(-FSWRxRate * ProbCompleteCure * 52.0/CycleD)));
	From2to0 = (1.0 - exp(-FemRxRate * ProbCompleteCure * 52.0/CycleD)) * (1.0 - 0.5 * 
		(1.0 - exp(-((1.0/AveDuration[1]) + FemRxRate * ProbPartialCure)*52.0/CycleD)));
	From2to0T = (1.0 - exp(-FemTeenRxRate * ProbCompleteCure * 52.0/CycleD)) * (1.0 - 0.5 * 
		(1.0 - exp(-((1.0/AveDuration[1]) + FemTeenRxRate * ProbPartialCure)*52.0/CycleD)));
	From2to0C = (1.0 - exp(-FSWRxRate * ProbCompleteCure * 52.0/CycleD)) * (1.0 - 0.5 * 
		(1.0 - exp(-((1.0/AveDuration[1]) + FSWRxRate * ProbPartialCure)*52.0/CycleD)));

	for(ia=0; ia<7; ia++){
		From0to1[ia][0] = 1.0 - exp(-Incidence * (FertilityTable[ia*5+2][0][0]/
			FertilityTable[2][0][0])*52.0 / CycleD);
	}
	if(HIVind==1){
		for(ia=0; ia<7; ia++){
			for(is=1; is<7; is++){
				From0to1[ia][is] = 1.0 - exp(-Incidence * (FertilityTable[ia*5+2][0][0] /
					FertilityTable[2][0][0]) * (1.0 + IncidenceIncrease[is-1]) * 52.0 / CycleD);
			}
		}
	}
}

void VCtransition::GetNewStage(int ID, double p)
{
	int ia, is;
	
	if(Register[ID-1].VCstage==0){
		ia = Register[ID-1].AgeGroup - 3;
		is = Register[ID-1].HIVstage;
		if(ia>=0 && ia<7){
			if(p<From0to1[ia][is]){Register[ID-1].VCstageE = 1;}
			else{Register[ID-1].VCstageE = 0;}
		}
		else{Register[ID-1].VCstageE = 0;}
	}
	if(Register[ID-1].VCstage==1){
		if(Register[ID-1].FSWind==1){
			if(p<From1to0C){Register[ID-1].VCstageE = 0;}
			else if(p<From1to0C + From1to2C){Register[ID-1].VCstageE = 2;}
			else{Register[ID-1].VCstageE = 1;}
		}
		else{
			if(p<From1to0){Register[ID-1].VCstageE = 0;}
			else if(p<From1to0 + From1to2){Register[ID-1].VCstageE = 2;}
			else{Register[ID-1].VCstageE = 1;}
		}
	}
	if(Register[ID-1].VCstage==2){
		if (Register[ID - 1].FSWind == 1){
			if(p<From2to0C){Register[ID-1].VCstageE = 0;}
			else if(p<From2to0C + From2to1C){Register[ID-1].VCstageE = 1;}
			else{Register[ID-1].VCstageE = 2;}
		}
		else if (Register[ID - 1].AgeGroup<4){
			if(p<From2to0T){Register[ID-1].VCstageE = 0;}
			else if(p<From2to0T + From2to1T){Register[ID-1].VCstageE = 1;}
			else{Register[ID-1].VCstageE = 2;}
		}
		else{
			if(p<From2to0){Register[ID-1].VCstageE = 0;}
			else if(p<From2to0 + From2to1){Register[ID-1].VCstageE = 1;}
			else{Register[ID-1].VCstageE = 2;}
		}
	}
}

PaedHIV::PaedHIV(){}

void PaedHIV::GetNewStage(int ID, double p)
{
	double t1, t2;
	double AIDSprog, From2to4, From4toDead, From5toDead;

	if(Register[ID-1].HIVstage==2){
		t2 = CurrYear + 0.5 + 1.0 * (STDcycleCount/CycleD) + 
			1.0 * (BehavCycleCount - 1.0)/CycleS - Register[ID-1].DOB;
		t1 = t2 - 1.0/CycleD;
		if(t1<0.0){t1 = 0.0;}
		AIDSprog = 1.0 - pow(0.5, pow(t2/PreAIDSmedian, PreAIDSshape) - 
			pow(t1/PreAIDSmedian, PreAIDSshape));
		From2to4 = AIDSprog * (1.0 - HAARTaccess[CurrYear - StartYear]);
		if(p < From2to4){
			Register[ID-1].HIVstageE = 4;}
		else if(p < AIDSprog){
			Register[ID-1].HIVstageE = 5;
			Register[ID - 1].DOAI = 1.0 * CurrYear + 0.5 + 1.0 * (BehavCycleCount - 1) / CycleS +
				1.0 * (STDcycleCount - 1) / CycleD;
			Register[ID - 1].BaselineCD4 = Register[ID - 1].CD4;
			Register[ID - 1].VCThistory = 2;
			Register[ID - 1].CondomPref *= DiscloseEffectCondom;
		}
		else{
			Register[ID-1].HIVstageE = 2;}
	}
	if(Register[ID-1].HIVstage==4){
		From4toDead = 1.0 - exp(-1.0/(MeanAIDSsurvival * CycleD));
		if(p < From4toDead){RSApop.SetToDead(ID, 1);}
		else{Register[ID-1].HIVstageE = 4;}
	}
	if(Register[ID-1].HIVstage==5){
		From5toDead = 1.0 - exp(-1.0/(AveYrsOnART * CycleD));
		if(p < From5toDead){RSApop.SetToDead(ID, 1);}
		else{Register[ID-1].HIVstageE = 5;}
	}
}

Child::Child(int Sex)
{
	SexInd = Sex;
}

Indiv::Indiv()
{
	int ii;

	HighestGrade = 0;
	InSchool = 0;
	CurrContr = 0;
	PrevContr = 0;
	EverInjectable = 0;
	EverPill = 0;
	BirthInterval = 0.0;
	DOLC = 0.0;
	DOLB = 0.0;
	DOLW = 0.0;
	DOUD = 0.0;
	Fecundability = 0.0;
	HouseholdID = 0;
	DateHomeless = 0.0;
	ParentID[0] = 0;
	ParentID[1] = 0;
	for (ii = 0; ii < 20; ii++){ ChildIDs[ii] = 0; }
	FatherBaby = 0;
	HIVstage = 0;
	CTstage = 0;
	HDstage = 0;
	NGstage = 0;
	TVstage = 0;
	TPstage = 0;
	HSVstage = 0;
	BVstage = 0;
	VCstage = 0;
	DateInfect = 0.0;
	DesiredNewPartners = 0.0;
	CondomPrimary = 0.0;
	Condom2ndary = 0.0;
	TempCondomAdj = 1.0;
	DisclosedPrimary = 0;
	Disclosed2ndary = 0;
	NewStatus = 0;
	UVIprimary = 0;
	PVIprimary = 0;
	UVI2ndary = 0;
	PVI2ndary = 0;
	UVICSW = 0;
	PVICSW = 0;
	IDofCSW = 0;
	VCThistory = 0;
	SuscepHIVadj = 1.0;
	EverMSM = 0;
	EverBi = 0;
	RecentMSM = 0;
	RecentBi = 0;
	CircInd = 0;
	Visiting = 0;
	DateMig = 0.0;
	Imprisoned = 0;
	ReleaseDate = 0.0;
	PrevImprisoned = 0;
	Employed = 0;
	LogIncomeDif = 0.0;
	PrivatePension = 0;
	OnPrEP = 0.0;
	VaccineEff = 0.0;
	Date1stVacc = 0.0;
	DailyDrinkProb = 0.0;
	DrinksPerDD = 1.0;
	EverCasual = 0;
	CasualLast3mo = 0;
	BaselineAge = 0;
	BaselinePoor = 0;
	BaselineSchool = 0;
	BaselineBinge = 0;
	BaselineUnemployed = 1;
	BaselineVirgin = 1;
	//CumUnprotected = 0;
	CumCasual = 0;
	CumSTIs = 0;
	CumHSV2 = 0;
	CumTeenPreg = 0;
}

int Indiv::SelectEvent(int ID, double rnd, int sexpref)
{
	int ia, ii, EventType;
	double AcquireS1, AcquireS2, Marry1, Marry2, ORadj, Divorce, LoseS1, LoseS2;
	// Dependent versions of above rates:
	double DAcquireS1, DAcquireS2, DMarry1, DMarry2, DDivorce, DLoseS1, DLoseS2;
	int prisk, temp;
	double CumEventProb[8];

	AcquireS1 = 0.0;
	AcquireS2 = 0.0;
	Marry1 = 0.0;
	Marry2 = 0.0;
	Divorce = 0.0;
	LoseS1 = 0.0;
	LoseS2 = 0.0;
	DAcquireS1 = 0.0;
	DAcquireS2 = 0.0;
	DMarry1 = 0.0;
	DMarry2 = 0.0;
	DDivorce = 0.0;
	DLoseS1 = 0.0;
	DLoseS2 = 0.0;

	// Calculate rates at which new ST partners are acquired
	if(DesiredNewPartners>0.0 && VirginInd==0){
		if(SexInd==0 && sexpref==0){ // Male chooses female
			AcquireS1 = DesiredNewPartners * DesiredPartnerRiskM[RiskGroup-1][0] * 
				AdjSTrateM[RiskGroup-1][0]/CycleS;
			AcquireS2 = DesiredNewPartners * DesiredPartnerRiskM[RiskGroup-1][1] * 
				AdjSTrateM[RiskGroup-1][1]/CycleS;
		}
		else if (SexInd == 0 && sexpref == 1){ // Male chooses male
			AcquireS1 = DesiredNewPartners * DesiredPartnerRiskMSM[RiskGroup - 1][0] *
				AdjSTrateMSM[RiskGroup - 1][0] / CycleS;
			AcquireS2 = DesiredNewPartners * DesiredPartnerRiskMSM[RiskGroup - 1][1] *
				AdjSTrateMSM[RiskGroup - 1][1] / CycleS;
		}
		else{ // Female chooses male
			AcquireS1 = DesiredNewPartners * DesiredPartnerRiskF[RiskGroup-1][0] * 
				AdjSTrateF[RiskGroup-1][0]/CycleS;
			AcquireS2 = DesiredNewPartners * DesiredPartnerRiskF[RiskGroup-1][1] * 
				AdjSTrateF[RiskGroup-1][1]/CycleS;
		}
	}
	else if(DesiredNewPartners>0.0 && VirginInd==1){ // No adjustment factor applied
		if(SexInd==0 && sexpref==0){ // Male chooses female
			AcquireS1 = DesiredNewPartners * DesiredPartnerRiskM[RiskGroup-1][0]/CycleS;
			AcquireS2 = DesiredNewPartners * DesiredPartnerRiskM[RiskGroup-1][1]/CycleS;
		}
		else if (SexInd == 0 && sexpref == 1){ // Male chooses male
			AcquireS1 = DesiredNewPartners * DesiredPartnerRiskMSM[RiskGroup - 1][0] / CycleS;
			AcquireS2 = DesiredNewPartners * DesiredPartnerRiskMSM[RiskGroup - 1][1] / CycleS;
		}
		else{ // Female chooses male
			AcquireS1 = DesiredNewPartners * DesiredPartnerRiskF[RiskGroup-1][0]/CycleS;
			AcquireS2 = DesiredNewPartners * DesiredPartnerRiskF[RiskGroup-1][1]/CycleS;
		}
	}

	// Calculate rates at which marriage occurs
	if(CurrPartners>0.0 && MarriedInd==0){
		ia = CurrAge - 10;
		Marry1 = MarriageIncidence[ia][SexInd][Race] / CycleS;
		ORadj = 1.0;
		if (HighestGrade == 12){ ORadj *= ORmarriage2ndaryEdu[SexInd]; }
		if (HighestGrade == 13){ ORadj *= ORmarriageTertiary[SexInd]; }
		if (InSchool == 1){ ORadj *= ORmarriageInSchool; }
		if (Race == 1 && CurrAge >= 45 && SexInd == 0){ ORadj *= 0.5; }
		if (Race == 2){
			if (CurrAge >= 25 && CurrAge<35){ ORadj *= 2.0; }
			else if(CurrAge < 25){ ORadj *= 0.5; }
		}
		if (DOUD + 1.0 > CurrYear + 0.5 + (BehavCycleCount - 1.0) / CycleS){ 
			ORadj *= ORremarriage[SexInd]; }
		Marry1 = 1.0 / (1.0 + (1.0 - Marry1) / (Marry1 * ORadj));
		if (CurrPartners == 2){
			if (Register[IDprimary - 1].MarriedInd == 0 && Register[IDprimary - 1].NewStatus == 0 &&
				Register[ID2ndary - 1].MarriedInd == 0 && Register[ID2ndary - 1].NewStatus == 0){
				Marry1 *= 0.5;
				Marry2 = Marry1;
			}
			else if (Register[IDprimary - 1].MarriedInd == 0 && Register[IDprimary - 1].NewStatus == 0){
				Marry2 = 0.0;
			}
			else if (Register[ID2ndary - 1].MarriedInd == 0 && Register[ID2ndary - 1].NewStatus == 0){
				Marry2 = Marry1;
				Marry1 = 0.0;
			}
			else{
				Marry1 = 0.0;
				Marry2 = 0.0;
			}
			if (SexInd == 0 && Register[ID2ndary - 1].SexInd == 0){ Marry2 *= SameSexMarried; }
		}
		if (CurrPartners == 1){
			if (Register[IDprimary - 1].MarriedInd == 1 || Register[IDprimary - 1].NewStatus == 1){
				Marry1 = 0.0;
			}
		}
		if (SexInd == 0 && Register[IDprimary - 1].SexInd == 0){ Marry1 *= SameSexMarried; }
		if (Register[IDprimary - 1].CurrAge < 15){ Marry1 = 0.0; }
		if (CurrPartners == 2){ if (Register[ID2ndary - 1].CurrAge < 15){ Marry2 = 0.0; } }
		if (Marry1 > 1.0){ Marry1 = 1.0; }
		if (Marry2 > 1.0){ Marry2 = 1.0; }
	}

	// Calculate rate at which divorce occurs
	if(MarriedInd==1){
		if(Register[IDprimary-1].NewStatus==0){ // Messy
			Divorce = LTseparation[AgeGroup-2][SexInd]/CycleS;
			Divorce *= DivorceAdj * pow(DivorceTrend, CurrYear - StartYear);
			// Adjust for effect of HIV disclosure
			if ((SexInd == 1 && DisclosedPrimary == 1 && Register[IDprimary - 1].VCThistory<2) || 
				(VCThistory<2 && Register[IDprimary - 1].SexInd == 1 && 
				Register[IDprimary - 1].DisclosedPrimary == 1)){
				Divorce *= 1.0 / DiscloseEffectRelDur;
			}
		}
	}

	// Calculate rate at which primary partnership is terminated if it is ST
	// (Note that the code differs from the corresponding function in TSHISA, where we weren't
	// worried about the order in which partnerships occur).
	if(CurrPartners>0 && MarriedInd==0){
		prisk = Register[IDprimary-1].RiskGroup - 1; // Messy
		if(Register[IDprimary-1].NewStatus==1){
			LoseS1 = 0.0;}
		else if(SexInd==0){
			if (Register[IDprimary - 1].SexInd == 1){ 
				LoseS1 = PartnerRateAdj / MeanDurSTrel[RiskGroup - 1][prisk]; }
			else{ LoseS1 = PartnerRateAdj / MeanDurST_MSM; }
		}
		else{
			LoseS1 = PartnerRateAdj/MeanDurSTrel[prisk][RiskGroup-1];}
		if (DrinksPerDD >= 5.0 && DailyDrinkProb > 1.0 / 30.0){ LoseS1 *= RRbreakupBinge; }
		if (Register[IDprimary - 1].DrinksPerDD >= 5.0 && Register[IDprimary - 1].DailyDrinkProb > 1.0 / 30.0){ 
			LoseS1 *= RRbreakupBinge; }
		//LoseS1 *= pow(0.90, HighestGrade + Register[IDprimary - 1].HighestGrade - 20.0);
		LoseS1 = LoseS1/CycleS;
		// Adjust for effect of HIV disclosure
		if ((SexInd == 1 && DisclosedPrimary == 1 && Register[IDprimary - 1].VCThistory<2) ||
			(VCThistory<2 && Register[IDprimary - 1].SexInd == 1 && ((Register[IDprimary - 1].IDprimary == ID &&
			Register[IDprimary - 1].DisclosedPrimary == 1) || (Register[IDprimary - 1].ID2ndary == ID &&
			Register[IDprimary - 1].Disclosed2ndary == 1)))){
			LoseS1 *= 1.0 / DiscloseEffectRelDur;
		}
	}

	// Calculate rate at which 2ndary partnership is terminated
	// (Note that the code differs from the corresponding function in TSHISA, where we weren't
	// worried about the order in which partnerships occur).
	if(CurrPartners==2){
		prisk = Register[ID2ndary-1].RiskGroup - 1; // Messy
		if(Register[ID2ndary-1].NewStatus==1){
			LoseS2 = 0.0;}
		else if(SexInd==0){
			if (Register[IDprimary - 1].SexInd == 1){ 
				LoseS2 = PartnerRateAdj / MeanDurSTrel[RiskGroup - 1][prisk]; }
			else{ LoseS2 = PartnerRateAdj / MeanDurST_MSM; }
		}
		else{
			LoseS2 = PartnerRateAdj/MeanDurSTrel[prisk][RiskGroup-1];}
		if (DrinksPerDD >= 5.0 && DailyDrinkProb > 1.0 / 30.0){ LoseS2 *= RRbreakupBinge; }
		if (Register[ID2ndary - 1].DrinksPerDD >= 5.0 && Register[ID2ndary - 1].DailyDrinkProb > 1.0 / 30.0){ 
			LoseS2 *= RRbreakupBinge; }
		//LoseS2 *= pow(0.90, HighestGrade + Register[ID2ndary - 1].HighestGrade - 20.0);
		LoseS2 = LoseS2/CycleS;
		// Adjust for effect of HIV disclosure
		if ((SexInd == 1 && Disclosed2ndary == 1 && Register[ID2ndary - 1].VCThistory<2) ||
			(VCThistory<2 && Register[ID2ndary - 1].SexInd == 1 && ((Register[ID2ndary - 1].IDprimary == ID &&
			Register[ID2ndary - 1].DisclosedPrimary == 1) || (Register[ID2ndary - 1].ID2ndary == ID &&
			Register[ID2ndary - 1].Disclosed2ndary == 1)))){
			LoseS2 *= 1.0 / DiscloseEffectRelDur;
		}
	}

	// Calculate dependent probabilities of different events
	if(CurrPartners==0){
		DAcquireS1 = ConvertToDependent2(AcquireS1, AcquireS2, 1);
		DAcquireS2 = ConvertToDependent2(AcquireS1, AcquireS2, 2);
	}
	if(CurrPartners==1 && MarriedInd==0){
		if(RiskGroup==1){
			DAcquireS1 = ConvertToDependent4(AcquireS1, AcquireS2, LoseS1, Marry1, 1);
			DAcquireS2 = ConvertToDependent4(AcquireS1, AcquireS2, LoseS1, Marry1, 2);
			DLoseS1 = ConvertToDependent4(AcquireS1, AcquireS2, LoseS1, Marry1, 3);
			DMarry1 = ConvertToDependent4(AcquireS1, AcquireS2, LoseS1, Marry1, 4);
		}
		else{
			DLoseS1 = ConvertToDependent2(LoseS1, Marry1, 1);
			DMarry1 = ConvertToDependent2(LoseS1, Marry1, 2);
		}
	}
	if(CurrPartners==1 && MarriedInd==1){
		if(RiskGroup==1){
			DAcquireS1 = ConvertToDependent3(AcquireS1, AcquireS2, Divorce, 1);
			DAcquireS2 = ConvertToDependent3(AcquireS1, AcquireS2, Divorce, 2);
			DDivorce = ConvertToDependent3(AcquireS1, AcquireS2, Divorce, 3);
		}
		else{
			DDivorce = ConvertToDependent1(Divorce);}
	}
	if(CurrPartners==2 && MarriedInd==0){
		DLoseS1 = ConvertToDependent4(LoseS1, LoseS2, Marry1, Marry2, 1);
		DLoseS2 = ConvertToDependent4(LoseS1, LoseS2, Marry1, Marry2, 2);
		DMarry1 = ConvertToDependent4(LoseS1, LoseS2, Marry1, Marry2, 3);
		DMarry2 = ConvertToDependent4(LoseS1, LoseS2, Marry1, Marry2, 4);
	}
	if(CurrPartners==2 && MarriedInd==1){
		DDivorce = ConvertToDependent2(Divorce, LoseS2, 1);
		DLoseS2 = ConvertToDependent2(Divorce, LoseS2, 2);
	}

	// Determine which event occurs
	CumEventProb[7] = 1.0;
	CumEventProb[6] = CumEventProb[7] - DLoseS2;
	CumEventProb[5] = CumEventProb[6] - DLoseS1;
	CumEventProb[4] = CumEventProb[5] - DDivorce;
	CumEventProb[3] = CumEventProb[4] - DMarry2;
	CumEventProb[2] = CumEventProb[3] - DMarry1;
	CumEventProb[1] = CumEventProb[2] - DAcquireS2;
	CumEventProb[0] = CumEventProb[1] - DAcquireS1;

	for(ii=0; ii<10; ii++){
		if(rnd<CumEventProb[ii]){
			EventType = ii;
			break;
		}
	}

	return EventType;
}

double Indiv::ConvertToDependent1(double rate1)
{
	double prob;
	prob = 1.0 - exp(-rate1);
	return prob;
}

double Indiv::ConvertToDependent2(double rate1, double rate2, int ord)
{
	double prob1, prob2, prob;

	if(TransitionCalc==0){
		prob1 = (1.0 - exp(-rate1)) * (1.0 - 0.5 * (1.0 - exp(-rate2)));
		prob2 = (1.0 - exp(-rate2)) * (1.0 - 0.5 * (1.0 - exp(-rate1)));
	}
	else{
		prob1 = (1.0 - exp(-rate1 - rate2)) * rate1/(rate1 + rate2);
		prob2 = (1.0 - exp(-rate1 - rate2)) * rate2/(rate1 + rate2);
	}

	if(ord==1){
		prob = prob1;}
	if(ord==2){
		prob = prob2;}

	return prob;
}

double Indiv::ConvertToDependent3(double rate1, double rate2, double rate3, int ord)
{
	double prob, totalrate;

	if(TransitionCalc==0){
		if(ord==1){
			prob = (1.0 - exp(-rate1)) * (1.0 - 0.5 * (2.0 - exp(-rate2) - exp(-rate3)) +
				(1.0 - exp(-rate2)) * (1.0 - exp(-rate3))/3.0);}
		if(ord==2){
			prob = (1.0 - exp(-rate2)) * (1.0 - 0.5 * (2.0 - exp(-rate1) - exp(-rate3)) +
				(1.0 - exp(-rate1)) * (1.0 - exp(-rate3))/3.0);}
		if(ord==3){
			prob = (1.0 - exp(-rate3)) * (1.0 - 0.5 * (2.0 - exp(-rate1) - exp(-rate2)) +
				(1.0 - exp(-rate1)) * (1.0 - exp(-rate2))/3.0);}
	}
	else{
		totalrate = rate1 + rate2 + rate3;
		if(ord==1){
			prob = (1.0 - exp(-totalrate)) * rate1/totalrate;}
		if(ord==2){
			prob = (1.0 - exp(-totalrate)) * rate2/totalrate;}
		if(ord==3){
			prob = (1.0 - exp(-totalrate)) * rate3/totalrate;}
	}

	return prob;
}

double Indiv::ConvertToDependent4(double rate1, double rate2, double rate3, double rate4, int ord)
{
	double prob, totalrate;
	double iprob1, iprob2, iprob3, iprob4;

	if(TransitionCalc==0){
		iprob1 = 1.0 - exp(-rate1);
		iprob2 = 1.0 - exp(-rate2);
		iprob3 = 1.0 - exp(-rate3);
		iprob4 = 1.0 - exp(-rate4);
		if(ord==1){
			prob = iprob1 * (1.0 - 0.5 * (iprob2 + iprob3 + iprob4) + (iprob2 * iprob3 +
				iprob2 * iprob4 + iprob3 * iprob4)/3.0 - 0.25 * iprob2 * iprob3 * iprob4);}
		if(ord==2){
			prob = iprob2 * (1.0 - 0.5 * (iprob1 + iprob3 + iprob4) + (iprob1 * iprob3 +
				iprob1 * iprob4 + iprob3 * iprob4)/3.0 - 0.25 * iprob1 * iprob3 * iprob4);}
		if(ord==3){
			prob = iprob3 * (1.0 - 0.5 * (iprob1 + iprob2 + iprob4) + (iprob1 * iprob2 +
				iprob1 * iprob4 + iprob2 * iprob4)/3.0 - 0.25 * iprob1 * iprob2 * iprob4);}
		if(ord==4){
			prob = iprob4 * (1.0 - 0.5 * (iprob1 + iprob2 + iprob3) + (iprob1 * iprob2 +
				iprob1 * iprob3 + iprob2 * iprob3)/3.0 - 0.25 * iprob1 * iprob2 * iprob3);}
	}
	else{
		totalrate = rate1 + rate2 + rate3 + rate4;
		if(ord==1){
			prob = (1.0 - exp(-totalrate)) * rate1/totalrate;}
		if(ord==2){
			prob = (1.0 - exp(-totalrate)) * rate2/totalrate;}
		if(ord==3){
			prob = (1.0 - exp(-totalrate)) * rate3/totalrate;}
		if(ord==4){
			prob = (1.0 - exp(-totalrate)) * rate4/totalrate;}
	}

	return prob;
}

void Indiv::AssignCondom(int RelType, int ID, int PID)
{
	int ia, ig;
	double randc, ICondom, PCondom, AveCondom;

	randc = rg.Random();

	if (RelType == 0){ // ST partnership (unmarried)
		ICondom = 1.0/(1.0 + (1.0/CondomUseST[AgeGroup - 2][SexInd] - 1.0)/(CondomPref * TempCondomAdj));
		ia = Register[PID - 1].AgeGroup - 2;
		ig = Register[PID - 1].SexInd;
		PCondom = 1.0 / (1.0 + (1.0 / CondomUseST[ia][ig] - 1.0) / (Register[PID-1].CondomPref *
			Register[PID - 1].TempCondomAdj));
		if (SexInd == 0 && ig == 0){ AveCondom = 0.5*(ICondom + PCondom); } // MSM
		else if (ig == 0){ AveCondom = GenderEquality * ICondom + (1.0 - GenderEquality) * PCondom; }
		else { AveCondom = GenderEquality * PCondom + (1.0 - GenderEquality) * ICondom; }
		if (randc < AveCondom){
			if (IDprimary == PID){ CondomPrimary = 1.0; }
			else{ Condom2ndary = 1.0; }
			if (Register[PID - 1].IDprimary == ID){ Register[PID - 1].CondomPrimary = 1.0; }
			else{ Register[PID - 1].Condom2ndary = 1.0; }
			if (Register[ID - 1].SexInd == 1 && Register[ID - 1].CurrContr == 0){
				Register[ID - 1].TestStartContr(1, 1); }
			if (Register[PID - 1].SexInd == 1 && Register[PID - 1].CurrContr == 0){ 
				Register[PID - 1].TestStartContr(1, 1); }
		}
		else{
			if (IDprimary == PID){ CondomPrimary = 0.0; }
			else{ Condom2ndary = 0.0; }
			if (Register[PID - 1].IDprimary == ID){ Register[PID - 1].CondomPrimary = 0.0; }
			else{ Register[PID - 1].Condom2ndary = 0.0; }
			if (Register[ID - 1].SexInd == 1 && Register[ID - 1].CurrContr == 0){ 
				Register[ID - 1].TestStartContr(1, 0); }
			if (Register[PID - 1].SexInd == 1 && Register[PID - 1].CurrContr == 0){ 
				Register[PID - 1].TestStartContr(1, 0); }
		}
	}
	else{ // LT partnership (married)
		if (VCThistory<2 && ((Register[PID - 1].IDprimary == ID && Register[PID - 1].DisclosedPrimary == 1) ||
			(Register[PID - 1].ID2ndary == ID && Register[PID - 1].Disclosed2ndary == 1))){
			randc = 0.0; // No discontinuation of condoms in serodiscordant relationships
		}
		if (IDprimary == PID){
			if (DisclosedPrimary == 1 && Register[PID - 1].VCThistory<2){ randc = 0.0; }
			if (CondomPrimary > 0.0){
				if (randc > RelEffectCondom[1]){
					CondomPrimary = 0.0;
					if (Register[PID - 1].IDprimary == ID){ Register[PID - 1].CondomPrimary = 0.0; }
					else{ Register[PID - 1].Condom2ndary = 0.0; }
					// Check if discontinuation of condom use leads to adoption of other contraception
					if (Register[ID - 1].SexInd == 1 && Register[ID - 1].CurrContr == 0){
						Register[ID - 1].TestStartContr(2, 0);}
					if (Register[PID - 1].SexInd == 1 && Register[PID - 1].CurrContr == 0){
						Register[PID - 1].TestStartContr(2, 0);}
				}
			}
		}
		else{
			if (Disclosed2ndary == 1 && Register[PID - 1].VCThistory<2){ randc = 0.0; }
			if (Condom2ndary > 0.0){
				if (randc > RelEffectCondom[1]){
					Condom2ndary = 0.0;
					if (Register[PID - 1].IDprimary == ID){ Register[PID - 1].CondomPrimary = 0.0; }
					else{ Register[PID - 1].Condom2ndary = 0.0; }
					// Check if discontinuation of condom use leads to adoption of other contraception
					if (Register[ID - 1].SexInd == 1 && Register[ID - 1].CurrContr == 0){
						Register[ID - 1].TestStartContr(2, 0);}
					if (Register[PID - 1].SexInd == 1 && Register[PID - 1].CurrContr == 0){
						Register[PID - 1].TestStartContr(2, 0);}
				}
			}
		}
	}
}

void Indiv::SimulateSexActs(int ID)
{
	// Note that this function only gets called if the individual is alive and is male.
	// Also note that in this function we're assuming that the specified frequencies of sex in
	// the SexAssumps sheet are monthly.
	// From Raikov's theorem we know that if the sum of 2 independent RVs is Poisson-distributed
	// then so is each of the RVs. So we can simulate the # unprotected and the # protected sex
	// acts as separate Poisson processes.

	double r[20]; // random numbers
	double ExpActs, ExpCondomUse, MaxFreq, AveLength;
	int ii, IndAge, PartnerAge, x, PartnerSex, Tries, TempID, ir;
	double xlam;

	// Note that I've arbitrarily fixed the seed for now.
	//int seed = 1272 + ID * 8 + (2025 - CurrYear) * 27 + BehavCycleCount * 31 + STDcycleCount * 95; 
	//if(CurrYear>=StartYear+FixedPeriod){
	//	seed += CurrSim;}
	//CRandomMersenne rg(seed);
	for(ii=0; ii<20; ii++){
		r[ii] = rg.Random();}

	IndAge = AgeGroup - 2;

	// Sex acts with primary partner
	if (CurrPartners>0 && Imprisoned==0){
		PartnerAge = Register[IDprimary-1].AgeGroup - 2;
		PartnerSex = Register[IDprimary - 1].SexInd;
		ExpCondomUse = CondomPrimary;
		if(MarriedInd==1){
			if (PartnerSex == 1){
				ExpActs = (GenderEquality * FreqSexLT[PartnerAge][1] + (1.0 - GenderEquality) *
					FreqSexLT[IndAge][0]) * 12.0 / CycleD;
			}
			else{ExpActs = 0.5 * (FreqSexLT[PartnerAge][0] + FreqSexLT[IndAge][0]) * 12.0 / CycleD;}
		}
		else{
			if (PartnerSex == 1){
				ExpActs = (GenderEquality * FreqSexST[PartnerAge][1] + (1.0 - GenderEquality) *
					FreqSexST[IndAge][0]) * 12.0 / CycleD;
			}
			else{ExpActs = FreqSexST_MSM * 12.0 / CycleD;}
		}
		// Adjust for partner absences
		// Note that we refer to CycleS not CycleD when making these adjustments because visits
		// to sexual partners are updated once every sexual behaviour cycle.
		if (CurrUrban != Register[IDprimary - 1].CurrUrban){
			if (Visiting == 1 && Register[IDprimary - 1].Visiting == 0){
				MaxFreq = 1.0 / (MaxVisitLength * pow(VisitFreq, 0.5));
				ExpActs *= pow(365.0, MaxCoitalFreqWeight) * pow(MaxFreq, MaxCoitalFreqWeight);
				AveLength = MaxVisitLength / pow(VisitFreq, 0.5);
				if (AveLength < 365.0 / CycleS){ ExpActs *= AveLength / (365.0 / CycleS); }
			}
			else if (Visiting == 0 && Register[IDprimary - 1].Visiting == 1){
				MaxFreq = 1.0 / (MaxVisitLength * pow(Register[IDprimary - 1].VisitFreq, 0.5));
				ExpActs *= pow(365.0, MaxCoitalFreqWeight) * pow(MaxFreq, MaxCoitalFreqWeight);
				AveLength = MaxVisitLength / pow(Register[IDprimary - 1].VisitFreq, 0.5);
				if (AveLength < 365.0 / CycleS){ ExpActs *= AveLength / (365.0 / CycleS); }
			}
			else{ ExpActs = 0.0; }
		}
		if (CurrUrban == Register[IDprimary - 1].CurrUrban && Visiting == 1){
			AveLength = MaxVisitLength / pow(VisitFreq, 0.5);
			if (AveLength < 365.0 / CycleS){
				ExpActs *= (1.0 - CycleS * AveLength / 365.0);}
			else{ ExpActs = 0.0; }
		}
		if (CurrUrban == Register[IDprimary - 1].CurrUrban && Register[IDprimary - 1].Visiting==1){
			AveLength = MaxVisitLength / pow(Register[IDprimary - 1].VisitFreq, 0.5);
			if (AveLength < 365.0 / CycleS){ 
				ExpActs *= (1.0 - CycleS * AveLength / 365.0); }
			else{ ExpActs = 0.0; }
		}

		// Simulate # unprotected sex acts
		xlam = ExpActs * (1.0 - ExpCondomUse);
		x = RandomSexGenerator(r[0], xlam);
		UVIprimary = x;
		if(Register[IDprimary-1].IDprimary==ID){
			Register[IDprimary-1].UVIprimary = x;}
		else{
			Register[IDprimary-1].UVI2ndary = x;}

		// Simulate # protected sex acts
		xlam = ExpActs * ExpCondomUse;
		x = RandomSexGenerator(r[1], xlam);
		PVIprimary = x;
		if(Register[IDprimary-1].IDprimary==ID){
			Register[IDprimary-1].PVIprimary = x;}
		else{
			Register[IDprimary-1].PVI2ndary = x;}
	}
	else{
		UVIprimary = 0;
		PVIprimary = 0;
	}
	/*if (UVIprimary > 0 && StructuralRCTcalib == 1 && CurrYear >= BaselineStart){
		CumUnprotected = 1;
		Register[IDprimary - 1].CumUnprotected = 1;
	}*/

	// Sex acts with 2ndary partner
	if (CurrPartners == 2 && Imprisoned == 0){
		PartnerAge = Register[ID2ndary-1].AgeGroup - 2;
		ExpCondomUse = Condom2ndary;
		if (PartnerSex == 1){
			ExpActs = (GenderEquality * FreqSexST[PartnerAge][1] + (1.0 - GenderEquality) *
				FreqSexST[IndAge][0]) * 12.0 / CycleD;
		}
		else{ExpActs = FreqSexST_MSM * 12.0 / CycleD;}
		// Adjust for partner absences
		if (CurrUrban != Register[ID2ndary - 1].CurrUrban){
			if (Visiting == 1 && Register[ID2ndary - 1].Visiting == 0){
				MaxFreq = 1.0 / (MaxVisitLength * pow(VisitFreq, 0.5));
				ExpActs *= pow(365.0, MaxCoitalFreqWeight) * pow(MaxFreq, MaxCoitalFreqWeight);
				AveLength = MaxVisitLength / pow(VisitFreq, 0.5);
				if (AveLength < 365.0 / CycleS){ ExpActs *= AveLength / (365.0 / CycleS); }
			}
			else if (Visiting == 0 && Register[ID2ndary - 1].Visiting == 1){
				MaxFreq = 1.0 / (MaxVisitLength * pow(Register[ID2ndary - 1].VisitFreq, 0.5));
				ExpActs *= pow(365.0, MaxCoitalFreqWeight) * pow(MaxFreq, MaxCoitalFreqWeight);
				AveLength = MaxVisitLength / pow(Register[ID2ndary - 1].VisitFreq, 0.5);
				if (AveLength < 365.0 / CycleS){ ExpActs *= AveLength / (365.0 / CycleS); }
			}
			else{ ExpActs = 0.0; }
		}
		if (CurrUrban == Register[ID2ndary - 1].CurrUrban && Visiting == 1){
			AveLength = MaxVisitLength / pow(VisitFreq, 0.5);
			if (AveLength < 365.0 / CycleS){
				ExpActs *= (1.0 - CycleS * AveLength / 365.0);}
			else{ ExpActs = 0.0; }
		}
		if (CurrUrban == Register[ID2ndary - 1].CurrUrban && Register[ID2ndary - 1].Visiting == 1){
			AveLength = MaxVisitLength / pow(Register[ID2ndary - 1].VisitFreq, 0.5);
			if (AveLength < 365.0 / CycleS){
				ExpActs *= (1.0 - CycleS * AveLength / 365.0);}
			else{ ExpActs = 0.0; }
		}

		// Simulate # unprotected sex acts
		xlam = ExpActs * (1.0 - ExpCondomUse);
		x = RandomSexGenerator(r[2], xlam);
		UVI2ndary = x;
		if(Register[ID2ndary-1].IDprimary==ID){
			Register[ID2ndary-1].UVIprimary = x;}
		else{
			Register[ID2ndary-1].UVI2ndary = x;}

		// Simulate # protected sex acts
		xlam = ExpActs * ExpCondomUse;
		x = RandomSexGenerator(r[3], xlam);
		PVI2ndary = x;
		if(Register[ID2ndary-1].IDprimary==ID){
			Register[ID2ndary-1].PVIprimary = x;}
		else{
			Register[ID2ndary-1].PVI2ndary = x;}
	}
	else{
		UVI2ndary = 0;
		PVI2ndary = 0;
	}
	/*if (UVI2ndary > 0 && StructuralRCTcalib == 1 && CurrYear >= BaselineStart){
		CumUnprotected = 1;
		Register[ID2ndary - 1].CumUnprotected = 1;
	}*/

	// Sex acts with sex workers
	if (RiskGroup == 1 && TotCurrFSW[Race]>0 && Imprisoned == 0){
		ExpActs = FSWcontactConstant * AgeEffectFSWcontact[CurrAge - 10]/CycleD;
		ExpActs *= (1.0 - MalePref);
		if (HIVstage > 0){ ExpActs *= HIVeffectPartners[HIVstage - 1]; }
		if (Employed == 1){ ExpActs *= RR_FSWcontactEmployedM; }
		//if (Employed == 1 || CurrYear>=2000){ ExpActs *= RR_FSWcontactEmployedM; }
		if(CurrPartners==0){
			ExpActs *= PartnerEffectFSWcontact[0];}
		if(CurrPartners==1 && MarriedInd==0){
			ExpActs *= PartnerEffectFSWcontact[1];}
		if(CurrPartners==1 && MarriedInd==1){
			ExpActs *= PartnerEffectFSWcontact[2];}
		if(CurrPartners==2 && MarriedInd==0){
			ExpActs *= PartnerEffectFSWcontact[3];}
		if(CurrPartners==2 && MarriedInd==1){
			ExpActs *= PartnerEffectFSWcontact[4];}
		ExpActs *= SWurban[1-CurrUrban]; 
		if (CurrPartners > 0){ ExpActs *= RaceEffectNew[Race]; }

		// Determine which sex worker the man has sex with.
		// (For simplicity we're assuming no man has contacts
		// with >1 CSW in a single STD cycle.)
		IDofCSW = r[6] * TotCurrFSW[Race];
		IDofCSW = RSApop.CSWregister[IDofCSW][Race];
		ExpCondomUse = 0.5 / (1.0 + (1.0 / CondomUseFSW - 1.0) / (Register[ID - 1].CondomPref * 
			Register[ID - 1].TempCondomAdj)) + 0.5 / (1.0 + (1.0 / CondomUseFSW - 1.0) /
			(Register[IDofCSW - 1].CondomPref * Register[IDofCSW - 1].TempCondomAdj));

		// Simulate # unprotected sex acts
		xlam = ExpActs * (1.0 - ExpCondomUse);
		x = RandomSexGenerator(r[4], xlam);
		UVICSW = x;
		
		// Simulate # protected sex acts
		xlam = ExpActs * ExpCondomUse;
		x = RandomSexGenerator(r[5], xlam);
		PVICSW = x;

		if ((MSMcalib == 1 || MSMcalibHIV == 1) && MalePref > 0.0){
			if (UVICSW > 0 || PVICSW > 0){
				EverBi = 1;
				RecentBi = 1;
			}
		}
	}
	else{
		UVICSW = 0;
		PVICSW = 0;
	}
	/*if (UVICSW > 0 && StructuralRCTcalib == 1 && CurrYear >= BaselineStart){
		CumUnprotected = 1;
		Register[IDofCSW - 1].CumUnprotected = 1;
	}*/

	// Sex acts with casual partners
	// First consider MSM relationships
	if (MalePref > 0.0 && CasualInd == 1 && IDofCasual == 0 && Imprisoned == 0){
		// Select ID of casual partner
		Tries = 0;
		while (Tries < 5 && IDofCasual==0){
			TempID = r[7 + Tries] * TotCurrCasual;
			TempID = RSApop.CasualRegister[TempID];
			if (Register[TempID - 1].IDofCasual == 0 && TempID != ID){
				IDofCasual = TempID;
				Register[TempID - 1].IDofCasual = ID;
			}
			Tries += 1; // Maximum of 5 attempts to sample casual partner
		}
		if (IDofCasual>0){
			ExpActs = CasualSexFreq * 12.0 / CycleD;
			PartnerAge = Register[IDofCasual - 1].AgeGroup - 2;
			ExpCondomUse = 0.5 / (1.0 + (1.0 / CondomUseCasual[PartnerAge] - 1.0) /
				(Register[IDofCasual - 1].CondomPref * Register[IDofCasual - 1].TempCondomAdj)) + 
				0.5 / (1.0 + (1.0 / CondomUseCasual[IndAge] - 1.0) /
				(Register[ID - 1].CondomPref * Register[ID - 1].TempCondomAdj));

			// Simulate # unprotected sex acts
			xlam = ExpActs * (1.0 - ExpCondomUse);
			x = RandomSexGenerator(r[12], xlam);
			UAIcasual = x;
			Register[IDofCasual - 1].UAIcasual = x;

			// Simulate # protected sex acts
			xlam = ExpActs * ExpCondomUse;
			x = RandomSexGenerator(r[13], xlam);
			PAIcasual = x;
			Register[IDofCasual - 1].PAIcasual = x;
		}
		else{
			UAIcasual = 0;
			PAIcasual = 0;
		}
	}
	// Then consider heterosexual casual relationships
	if (HetCasualInd == 1 && IDofCasual == 0 && Imprisoned == 0){
		// Select ID of casual partner
		Tries = 0;
		if(TotCurrCasualHet[Race][1]>0){
			while (Tries < 5 && IDofCasual == 0){
				TempID = r[7 + Tries] * TotCurrCasualHet[Race][1];
				TempID = RSApop.CasualRegisterF[TempID][Race];
				if (Register[TempID - 1].IDofCasual == 0 && TempID != ID){
					IDofCasual = TempID;
					Register[TempID - 1].IDofCasual = ID;
				}
				Tries += 1; // Maximum of 5 attempts to sample casual partner
			}
		}
		if (IDofCasual>0){
			ExpActs = CasualSexFreqHet * 12.0 / CycleD;
			PartnerAge = Register[IDofCasual - 1].AgeGroup - 2;
			ExpCondomUse = 0.5 / (1.0 + (1.0 / CondomUseCasualHet[PartnerAge][1] - 1.0) /
				(Register[IDofCasual - 1].CondomPref * Register[IDofCasual - 1].TempCondomAdj)) +
				0.5 / (1.0 + (1.0 / CondomUseCasualHet[IndAge][0] - 1.0) /
				(Register[ID - 1].CondomPref * Register[ID - 1].TempCondomAdj));

			// Simulate # unprotected sex acts
			xlam = ExpActs * (1.0 - ExpCondomUse);
			x = RandomSexGenerator(r[14], xlam);
			UAIcasual = x;
			Register[IDofCasual - 1].UAIcasual = x;

			// Simulate # protected sex acts
			xlam = ExpActs * ExpCondomUse;
			x = RandomSexGenerator(r[15], xlam);
			PAIcasual = x;
			Register[IDofCasual - 1].PAIcasual = x;
			if (UAIcasual + PAIcasual > 0 && CurrYear == 2004 && BehavCycleCount > CycleS*0.75){
				CasualLast3mo += 1;
			}
		}
		else{
			UAIcasual = 0;
			PAIcasual = 0;
		}
	}
	if (CasualInd == 0 && HetCasualInd == 0){
		UAIcasual = 0;
		PAIcasual = 0;
	}
	if (StructuralRCTcalib == 1 && CurrYear >= BaselineStart){
		/*if (UAIcasual > 0){
			CumUnprotected = 1;
			Register[IDofCasual - 1].CumUnprotected = 1;
		}*/
		if (UAIcasual + PAIcasual > 0){
			CumCasual = 1;
			Register[IDofCasual - 1].CumCasual = 1;
		}
	}

	if (FixedUncertainty == 1){ 
		ProtSexActs.out[CurrSim-1][CurrYear-StartYear] += (PVIprimary + PVI2ndary + PVICSW + PAIcasual) * PopWeight; 
		SWsexActs.out[CurrSim - 1][CurrYear - StartYear] += (UVICSW + PVICSW) * PopWeight;
		SWsexActsProt.out[CurrSim - 1][CurrYear - StartYear] += PVICSW * PopWeight;
	}
}

int Indiv::RandomSexGenerator(double p, double lambda)
{
	// Due to rounding there is a very small risk that the algorithm for returning
	// the Poisson-distributed RV might get into an infinite loop. To avoid this, I've
	// limited the # sex acts that can occur to at most 720 per year (which is a lot).

	int ii, SexActs, MaxSex;
	double CumProb, NextProb;

	MaxSex = 720/CycleD;
	NextProb = exp(-lambda);
	if(p<NextProb){
		SexActs = 0;}
	else{
		CumProb = NextProb;
		SexActs = MaxSex;
		for(ii=1; ii<MaxSex; ii++){
			NextProb *= lambda/ii;
			CumProb += NextProb;
			if(p<CumProb){
				SexActs = ii;
				break;
			}
		}
	}

	return SexActs;
}

void Indiv::GetInitCD4andVL()
{
	int ii, s, ind;
	double rand[2], a, b, p, q, x, bound, ExactAge;

	for (ii = 0; ii<2; ii++){rand[ii] = rg.Random();}
	bound = 0.0;
	ind = 2;
	s = 0;

	// Simulate initial CD4
	if (SexInd == 0){
		a = HIVtransitionM.MeanInitCD4;
		b = HIVtransitionM.SDinitCD4;
	}
	else{
		a = HIVtransitionF.MeanInitCD4;
		b = HIVtransitionF.SDinitCD4;
	}
	p = rand[0];
	q = 1.0 - rand[0];
	cdfnor(&ind, &p, &q, &x, &a, &b, &s, &bound);
	CD4 = exp(x);
	BaselineCD4 = 0.0; // Only relevant to those who initiate ART

	// Simulate initial VL (log scale)
	ExactAge = 1.0 * (CurrYear + 0.5 + (BehavCycleCount - 1.0) / CycleS - DOB);
	if (SexInd == 0){
		a = HIVtransitionM.MeanInitVL;
		a += HIVtransitionM.AgeEffectInitVL * (ExactAge - 28.0);
		b = HIVtransitionM.SDinitVL;
	}
	else{
		a = HIVtransitionF.MeanInitVL;
		a += HIVtransitionF.AgeEffectInitVL * (ExactAge - 28.0);
		b = HIVtransitionF.SDinitVL;
	}
	p = rand[1];
	q = 1.0 - rand[1];
	cdfnor(&ind, &p, &q, &x, &a, &b, &s, &bound);
	logVL = x;
	SPVL = x;
}

void Indiv::GetNewHIVstate(int ID, double p)
{
	double Prob1, Prob2;
	int MSMind;

	MSMind = 0;
	if (Register[ID - 1].SexInd == 0 && Register[ID - 1].MalePref > 0){ MSMind = 1; }

	Prob1 = 0.0;
	Prob2 = 0.0;
	if(HIVstage==0){
		// Note that we refer to opposite sex when getting transm prob
		if (SexInd == 0 && Register[ID - 1].MalePref < 1.0){
			Prob1 = HIVtransitionF.GetTransmProb(ID);}
		if (SexInd == 1){
			Prob1 = HIVtransitionM.GetTransmProb(ID);}
		if (SexInd == 0 && Register[ID - 1].MalePref > 0.0){
			Prob2 = HIVtransitionM.GetTransmProbMSM(ID);
			Prob1 = 1.0 - (1.0 - Prob1)*(1.0 - Prob2);
		}
		if (FixedUncertainty == 1){
			NewHIVexp.out[CurrSim - 1][CurrYear - StartYear] += Prob1; 
			NewHIVbySex[SexInd] += Prob1; 
		}
		if(p<Prob1){
			HIVstageE = 1;
			GetInitCD4andVL();
			DateInfect = 1.0 * (CurrYear + 0.5 + (STDcycleCount - 1.0)/CycleD +
				(BehavCycleCount - 1.0)/CycleS + (p/Prob1)/CycleD);
			if (FixedUncertainty == 1){ 
				NewHIV.out[CurrSim - 1][CurrYear - StartYear] += PopWeight; 
				if (AgeGroup>2 && AgeGroup<10 && CurrUrban == 0){ NewHIV_R.out[CurrSim - 1][CurrYear - StartYear] += PopWeight; }
				if (AgeGroup>2 && AgeGroup<10 && CurrUrban == 1){ NewHIV_U.out[CurrSim - 1][CurrYear - StartYear] += PopWeight; }
				if (SexInd == 0 && Register[ID - 1].MalePref > 0.0 && p<Prob2){
					NewHIV_MSM.out[CurrSim - 1][CurrYear - StartYear] += PopWeight;
				}
			}
		}
		else{
			HIVstageE = 0;}
	}
	else{
		if(SexInd==0){
			HIVtransitionM.GetNewStage(ID, p);}
		else{
			HIVtransitionF.GetNewStage(ID, p);}
	}
}

void Indiv::GetNewHSVstate(int ID, double p)
{
	double Prob1, Symptomatic;

	if(HSVstage==0){
		// Note that we refer to opposite sex when getting transm prob
		if(SexInd==0){
			Prob1 = HSVtransitionF.GetTransmProb(ID);}
		else{
			Prob1 = HSVtransitionM.GetTransmProb(ID);}
		if(p<Prob1){
			if(SexInd==0){Symptomatic = HSVtransitionM.SymptomaticPropn;}
			else{Symptomatic = HSVtransitionF.SymptomaticPropn;}
			if(p/Prob1 < Symptomatic){HSVstageE = 1;}
			else{HSVstageE = 4;}
			if (FixedUncertainty == 1){ NewHSV.out[CurrSim - 1][CurrYear - StartYear] += PopWeight; }
			if (StructuralRCTcalib == 1 && CurrYear >= BaselineStart){ CumHSV2 = 1; }
		}
		else{
			HSVstageE = 0;}
	}
	else{
		if(SexInd==0){
			HSVtransitionM.GetNewStage(ID, p);}
		else{
			HSVtransitionF.GetNewStage(ID, p);}
	}
}

void Indiv::GetNewTPstate(int ID, double p)
{
	double Prob1;

	if(TPstage==0){
		// Note that we refer to opposite sex when getting transm prob
		if(SexInd==0){
			Prob1 = TPtransitionF.GetTransmProb(ID);}
		else{
			Prob1 = TPtransitionM.GetTransmProb(ID);}
		if(p<Prob1){
			TPstageE = 1;
			if (FixedUncertainty == 1){ NewTP.out[CurrSim - 1][CurrYear - StartYear] += PopWeight; }
		}
		else{TPstageE = 0;}
	}
	else{
		if(SexInd==0){
			TPtransitionM.GetNewStage(ID, p);}
		else{
			TPtransitionF.GetNewStage(ID, p);}
	}
}

void Indiv::GetNewHDstate(int ID, double p)
{
	double Prob1, Symptomatic;

	if(HDstage==0){
		// Note that we refer to opposite sex when getting transm prob
		if(SexInd==0){
			Prob1 = HDtransitionF.GetTransmProb(ID, 4);}
		else{
			Prob1 = HDtransitionM.GetTransmProb(ID, 4);}
		if(p<Prob1){
			if(SexInd==0){Symptomatic = HDtransitionM.SymptomaticPropn;}
			else{Symptomatic = HDtransitionF.SymptomaticPropn;}
			if(p/Prob1 < Symptomatic){HDstageE = 1;}
			else{HDstageE = 2;}
		}
		else{
			HDstageE = 0;}
	}
	else{
		if(SexInd==0){
			HDtransitionM.GetNewStage(ID, p, 4);}
		else{
			HDtransitionF.GetNewStage(ID, p, 4);}
	}
}

void Indiv::GetNewNGstate(int ID, double p)
{
	double Prob1, Symptomatic;

	if(NGstage==0){
		// Note that we refer to opposite sex when getting transm prob
		if(SexInd==0){
			Prob1 = NGtransitionF.GetTransmProb(ID, 1);}
		else{
			Prob1 = NGtransitionM.GetTransmProb(ID, 1);}
		if(p<Prob1){
			if(SexInd==0){Symptomatic = NGtransitionM.SymptomaticPropn;}
			else{Symptomatic = NGtransitionF.SymptomaticPropn;}
			if(p/Prob1 < Symptomatic){NGstageE = 1;}
			else{NGstageE = 2;}
			if (FixedUncertainty == 1){
				NewNG.out[CurrSim - 1][CurrYear - StartYear] += PopWeight;
				NewNGbySex[SexInd] += PopWeight;
			}
			if (StructuralRCTcalib == 1 && CurrYear >= BaselineStart){ CumSTIs = 1; }
		}
		else{
			NGstageE = 0;}
	}
	else{
		if(SexInd==0){
			NGtransitionM.GetNewStage(ID, p, 1);}
		else{
			NGtransitionF.GetNewStage(ID, p, 1);}
	}
}

void Indiv::GetNewCTstate(int ID, double p)
{
	double Prob1, Symptomatic;

	if(CTstage==0){
		// Note that we refer to opposite sex when getting transm prob
		if(SexInd==0){
			Prob1 = CTtransitionF.GetTransmProb(ID, 2);}
		else{
			Prob1 = CTtransitionM.GetTransmProb(ID, 2);}
		if(p<Prob1){
			if(SexInd==0){Symptomatic = CTtransitionM.SymptomaticPropn;}
			else{Symptomatic = CTtransitionF.SymptomaticPropn;}
			if(p/Prob1 < Symptomatic){CTstageE = 1;}
			else{CTstageE = 2;}
			if (FixedUncertainty == 1){
				NewCT.out[CurrSim - 1][CurrYear - StartYear] += PopWeight;
				NewCTbySex[SexInd] += PopWeight;
			}
			if (StructuralRCTcalib == 1 && CurrYear >= BaselineStart){ CumSTIs = 1; }
		}
		else{
			CTstageE = 0;}
	}
	else{
		if(SexInd==0){
			CTtransitionM.GetNewStage(ID, p, 2);}
		else{
			CTtransitionF.GetNewStage(ID, p, 2);}
	}
}

void Indiv::GetNewTVstate(int ID, double p)
{
	double Prob1, Symptomatic;

	if(TVstage==0){
		// Note that we refer to opposite sex when getting transm prob
		if(SexInd==0){
			Prob1 = TVtransitionF.GetTransmProb(ID, 3);}
		else{
			Prob1 = TVtransitionM.GetTransmProb(ID, 3);}
		if(p<Prob1){
			if(SexInd==0){Symptomatic = TVtransitionM.SymptomaticPropn;}
			else{Symptomatic = TVtransitionF.SymptomaticPropn;}
			if(p/Prob1 < Symptomatic){TVstageE = 1;}
			else{TVstageE = 2;}
			if (FixedUncertainty == 1){
				NewTV.out[CurrSim - 1][CurrYear - StartYear] += PopWeight; 
				NewTVbySex[SexInd] += PopWeight;
			}
			if (StructuralRCTcalib == 1 && CurrYear >= BaselineStart){ CumSTIs = 1; }
		}
		else{
			TVstageE = 0;}
	}
	else{
		if(SexInd==0){
			TVtransitionM.GetNewStage(ID, p, 3);}
		else{
			TVtransitionF.GetNewStage(ID, p, 3);}
	}
}

void Indiv::TestStartContr(int EventType, int condom)
{
	// This function tests whether a woman starts hormonal contraception (and if so determines
	// the type of contraception). The function is only called if the woman is not currently
	// using hormonal contraception, and is not sterilized.
	// 1st argument is event type (1=new relationship, 2=stopping condom use, 3=postnatal).
	// 2nd argument is only relevant if event type is new relationship (0 implies no condom
	// use, 1 implies consistent condom use)

	double BaseProb, CumOR, AdjProb, AdjProb2, rc[2];

	rc[0] = rg.Random();
	rc[1] = rg.Random();

	// First test whether hormonal contraception is initiated
	CumOR = pow(EduEffectContr, HighestGrade - 10);
	//if (CurrYear >= 2000){ CumOR = pow(EduEffectContr, 3.0); }
	CumOR *= pow(Fecundability, 0.5);
	if (AgeGroup >= 10){ CumOR = 0.0; } // No hormonal contraception assumed in women aged 50+
	else if (AgeGroup < 5){ CumOR *= ORhormonalByYear[CurrYear - StartYear][0]; }
	else if (AgeGroup < 7){ CumOR *= ORhormonalByYear[CurrYear - StartYear][1]; }
	else { CumOR *= ORhormonalByYear[CurrYear - StartYear][2]; }
	if (HIVstage > 0){ CumOR *= HIVeffectContr; }
	//if (EventType==3 && VCThistory == 2){ CumOR *= HIVeffectContr; }
	if (PrevContr == 0){ CumOR *= NoPrevUseContr; }

	if (EventType == 1 || EventType == 2){ // New partner/cease condom use
		BaseProb = StartContrNewRel[Race];
		if (DOLB > 1.0 * CurrYear + 0.5 + 1.0 * (BehavCycleCount - 1.0)/ CycleS){ 
			BaseProb = 0.0; } // Still pregnant, so wouldn't start contraception
		if (AgeGroup > 3 && AgeGroup<10){ CumOR *= AgeEffectContr[AgeGroup - 3][0]; }
		CumOR *= PrevBirthEffectContr[Race];
		if (EventType == 1 && condom == 1){ CumOR *= CondomEffectContr; }
	}
	else{ // Postnatal initiation of contraception
		BaseProb = StartContrPostnatal[Race];
		if (CurrPartners == 0 && FSWind == 0){ BaseProb = 0.0; }
		if (AgeGroup > 3 && AgeGroup<10){ CumOR *= AgeEffectContr[AgeGroup - 3][1]; }
	}
	AdjProb = 1.0 / (1.0 - (1.0 - 1.0 / BaseProb) / CumOR);
	if (EventType == 2){
		AdjProb2 = 1.0 / (1.0 - (1.0 - 1.0 / AdjProb) / CondomEffectContr);
		AdjProb = (AdjProb - AdjProb2) / (1.0 - AdjProb2);
	}
	if (rc[0] < AdjProb){
		// Start hormonal contraception - select method
		BaseProb = InjectablePref[Race];
		CumOR = pow(ORinjectableEdu, HighestGrade - 10);
		//if (CurrYear >= 2000){ CumOR = pow(ORinjectableEdu, 3.0); }
		CumOR *= pow(ORinjectableAge, 1.0 * (CurrAge - 25.0) / 5.0);
		if (CurrYear<1995){
			CumOR *= ORinjectableInit + (1.0 - ORinjectableInit) *
				(CurrYear - StartYear) / 10.0;
		}
		if (DOLB>0){ CumOR *= 1.0 / ORinjectableNeverPreg; }
		AdjProb = 1.0 / (1.0 - (1.0 - 1.0 / BaseProb) / CumOR);
		if (PrevContr == 1){ AdjProb = 1.0 - (1.0 - AdjProb) * NewContrMethodWeight; }
		if (PrevContr == 2){ AdjProb = AdjProb * NewContrMethodWeight; }
		if (rc[1] < AdjProb){ CurrContr = 1; }
		else{ CurrContr = 2; }
	}
	else{CurrContr = 0;}
}

void Indiv::GetDiagnosed(int ID, int Setting)
{
	// Setting is 1 if diagnosed as result of an OI, 2 if diagnosed through PMTCT, 
	// 3 if diagnosed in community-based settings, 4 if diagnosed by ST, 0 otherwise

	int iy, is;
	double ARTprob, RandUni;

	iy = CurrYear - StartYear;
	RandUni = rg.Random();

	VCThistory = 2;
	CondomPref *= DiscloseEffectCondom;
	if (CurrPartners>0){
		TestDisclose(ID, 1, 0);}
	if (CurrPartners == 2){
		TestDisclose(ID, 2, 0);}

	// Calculate prob of ART initiation immediately after diagnosis
	if (HIVstage > 1 && HIVstage < 5 && HIVstageE < 5){
		if (Setting == 0 || Setting == 3 || Setting == 4){
			is = 4 - HIVstage;
			if (HIVstage == 2 && CD4 >= 500){ is = 3; }
			ARTprob = ARTeligibleGen[iy][is] * ARTinitiationU200[iy][SexInd] / 1.2;
			if (Setting == 3){ ARTprob *= RR_ARTstartCommunity; }
			if (Setting == 4){ ARTprob *= RR_ARTstartST; }
			ARTprob *= pow(RR_ARTstartPer100CD4, (CD4 - 100.0) / 100.0);
		}
		if (Setting == 1){
			is = 4 - HIVstage;
			if (HIVstage == 2 && CD4 >= 500){ is = 3; }
			ARTprob = ARTeligibleOI[iy][is] * ARTuptakeOI[iy];
		}
		if (Setting == 2){
			if (HIVstage == 4){ ARTprob = ARTeligibleGen[iy][0]; }
			else if (HIVstage == 3){ ARTprob = ARTeligiblePreg[iy][0]; }
			else { ARTprob = ARTeligiblePreg[iy][1]; }
			ARTprob = ARTprob * ARTuptakePreg[iy];
		}
		if (ARTprob > 1.0){ ARTprob = 1.0; }
		if (FixedUncertainty == 1 && AgeGroup > 2){ NewARTexp.out[CurrSim - 1][iy] += ARTprob * PopWeight; }
		if (RandUni < ARTprob){ 
			if (FixedUncertainty == 1){
				if (HIVstage == 4){ NewART200.out[CurrSim - 1][iy] += PopWeight; }
				if (HIVstage == 3){ NewART350.out[CurrSim - 1][iy] += PopWeight; }
				if (HIVstage == 2 && CD4 < 500){ NewART500.out[CurrSim - 1][iy] += PopWeight; }
				if (HIVstage == 2 && CD4 >= 500){ NewART500plus.out[CurrSim - 1][iy] += PopWeight; }
				if (Setting == 2){ TotBirthsART.out[CurrSim - 1][iy] += PopWeight; }
			}
			HIVstage = 5; 
			DOAI = 1.0 * CurrYear + 0.5 + 1.0 * (BehavCycleCount - 1) / CycleS +
				1.0 * (STDcycleCount - 1) / CycleD;
			BaselineCD4 = CD4;
		}
	}
}

void Indiv::GetRediagnosed(int ID, int Setting)
{
	// Setting is 1 if rediagnosed as result of an OI, 2 if rediagnosed through PMTCT, 
	// 3 if rediagnosed in community-based settings, 4 if rediagnosed by ST, 0 otherwise

	int iy, is;
	double ARTprob, RandUni;

	iy = CurrYear - StartYear;
	RandUni = rg.Random();

	// Calculate prob of ART initiation immediately after diagnosis
	if (HIVstage > 1 && HIVstage < 5){
		if (Setting == 0 || Setting == 3 || Setting == 4){
			is = 4 - HIVstage;
			if (HIVstage == 2 && CD4 >= 500){ is = 3; }
			ARTprob = ARTeligibleGen[iy][is] * ARTinitiationU200[iy][SexInd] / 1.2;
			if (Setting == 3){ ARTprob *= RR_ARTstartCommunity; }
			if (Setting == 4){ ARTprob *= RR_ARTstartST; }
			ARTprob *= pow(RR_ARTstartPer100CD4, (CD4 - 100.0) / 100.0);
		}
		if (Setting == 1){
			is = 4 - HIVstage;
			if (HIVstage == 2 && CD4 >= 500){ is = 3; }
			ARTprob = ARTeligibleOI[iy][is] * ARTuptakeOI[iy];
		}
		if (Setting == 2){
			if (HIVstage == 4){ ARTprob = ARTeligibleGen[iy][0]; }
			else if (HIVstage == 3){ ARTprob = ARTeligiblePreg[iy][0]; }
			else { ARTprob = ARTeligiblePreg[iy][1]; }
			ARTprob = ARTprob * ARTuptakePreg[iy];
		}
		if (ARTprob > 1.0){ ARTprob = 1.0; }
		if (FixedUncertainty == 1 && AgeGroup > 2){ NewARTexp.out[CurrSim - 1][iy] += ARTprob * PopWeight; }
		if (RandUni < ARTprob){
			if (FixedUncertainty == 1){
				if (HIVstage == 4){ NewART200.out[CurrSim - 1][iy] += PopWeight; }
				if (HIVstage == 3){ NewART350.out[CurrSim - 1][iy] += PopWeight; }
				if (HIVstage == 2 && CD4 < 500){ NewART500.out[CurrSim - 1][iy] += PopWeight; }
				if (HIVstage == 2 && CD4 >= 500){ NewART500plus.out[CurrSim - 1][iy] += PopWeight; }
				if (Setting == 2){ TotBirthsART.out[CurrSim - 1][iy] += PopWeight; }
			}
			HIVstage = 5;
			DOAI = 1.0 * CurrYear + 0.5 + 1.0 * (BehavCycleCount - 1) / CycleS +
				1.0 * (STDcycleCount - 1) / CycleD;
			BaselineCD4 = CD4;
		}
	}
}

void Indiv::TestDisclose(int ID, int Partner, int Timing)
{
	// This function tests whether or not an individual who has been diagnosed positive
	// discloses their HIV status to their primary (Partner=1) or secondary (Partner=2)
	// partner, either immediately after diagnosis (Timing=0) or thereafter (Timing=1).
	// If disclosure occurs, we also determine whether this leads to a change in condom
	// usage and whether the partner is tested.

	int ia, ig, iy;
	double BaseProb, CumOR, AdjProb, Condom1, Condom2, CondomP, rc[3];
	double OldProb, NewProb, TestProb, ProbReferral2;

	iy = CurrYear - StartYear;
	rc[0] = rg.Random();
	rc[1] = rg.Random();
	rc[2] = rg.Random();

	if (Timing == 0){ BaseProb = DiscloseProb; }
	else{ BaseProb = 1.0 - exp(-DiscloseRate/12.0); }
	if (Timing == 1){ ProbReferral2 = ProbReferral[0]; }
	else if (PropnAssistedNotif[iy] == 0.0){ ProbReferral2 = ProbReferral[0]; }
	else{ 
		ProbReferral2 = ProbReferral[0] * (1.0 - PropnAssistedNotif[iy]) + 
			ProbReferral[1] * PropnAssistedNotif[iy];
	}
	CumOR = 1.0;
	if (SexInd == 0){ CumOR *= DiscloseMale; }
	if (HIVstage == 5){ CumOR *= DiscloseART; }
	if (MarriedInd == 1 && Partner==1){ CumOR *= DiscloseMarried; }
	AdjProb = 1.0 / (1.0 - (1.0 - 1.0 / BaseProb) / CumOR);
	if (rc[0] < AdjProb){
		// Disclosure occurs - test if partner is tested and if there is a change in condom use.
		if (Partner == 1){ 
			DisclosedPrimary = 1; 
			if (rc[2] < ProbReferral2 && Register[IDprimary - 1].VCThistory<2){ // Partner gets tested
				TotalTestsPartner.out[CurrSim - 1][CurrYear - StartYear] += Register[IDprimary - 1].PopWeight;
				if (Register[IDprimary - 1].HIVstage > 1){
					Register[IDprimary - 1].GetDiagnosed(IDprimary, 0); 
					PosTestsPartner.out[CurrSim - 1][CurrYear - StartYear] += Register[IDprimary - 1].PopWeight;
					NewDiagPartner.out[CurrSim - 1][CurrYear - StartYear] += Register[IDprimary - 1].PopWeight;
				}
				else{ Register[IDprimary - 1].VCThistory = 1; }
				if (Register[IDprimary - 1].HIVstage == 1){ AcuteTests.out[CurrSim - 1][iy] += Register[IDprimary - 1].PopWeight; }
			}
			if (Register[IDprimary - 1].VCThistory < 2 && CondomPrimary==0.0){
				Condom1 = 1.0 / (1.0 + (1.0 / CondomUseST[AgeGroup - 2][SexInd] - 1.0) / (CondomPref *
					TempCondomAdj / DiscloseEffectCondom));
				Condom2 = 1.0 / (1.0 + (1.0 / CondomUseST[AgeGroup - 2][SexInd] - 1.0) / (CondomPref *
					TempCondomAdj));
				ia = Register[IDprimary - 1].AgeGroup - 2;
				ig = Register[IDprimary - 1].SexInd;
				CondomP = 1.0 / (1.0 + (1.0 / CondomUseST[ia][ig] - 1.0) / 
					(Register[IDprimary - 1].CondomPref * Register[IDprimary - 1].TempCondomAdj));
				if (MarriedInd == 1){
					Condom1 *= RelEffectCondom[1];
					Condom2 *= RelEffectCondom[1];
					CondomP *= RelEffectCondom[1];
				}
				if (ig == SexInd){ // MSM
					OldProb = 0.5 * (Condom1 + CondomP);
					NewProb = 0.5 * (Condom2 + CondomP);
				}
				else if (SexInd == 0){
					OldProb = (1.0 - GenderEquality) * Condom1 + GenderEquality * CondomP;
					NewProb = (1.0 - GenderEquality) * Condom2 + GenderEquality * CondomP;
				}
				else {
					OldProb = (1.0 - GenderEquality) * CondomP + GenderEquality * Condom1;
					NewProb = (1.0 - GenderEquality) * CondomP + GenderEquality * Condom2;
				}
				TestProb = (NewProb - OldProb) / (1.0 - OldProb);
				if (rc[1] < TestProb){ 
					CondomPrimary = 1.0; 
					if (Register[IDprimary - 1].IDprimary == ID){ 
						Register[IDprimary - 1].CondomPrimary = 1.0; }
					else{ Register[IDprimary - 1].Condom2ndary = 1.0; }
				}
			}
		}
		else{ 
			Disclosed2ndary = 1; 
			if (rc[2] < ProbReferral2 && Register[ID2ndary - 1].VCThistory<2){ // Partner gets tested
				TotalTestsPartner.out[CurrSim - 1][CurrYear - StartYear] += Register[ID2ndary - 1].PopWeight;
				if (Register[ID2ndary - 1].HIVstage > 1){
					Register[ID2ndary - 1].GetDiagnosed(ID2ndary, 0); 
					PosTestsPartner.out[CurrSim - 1][CurrYear - StartYear] += Register[ID2ndary - 1].PopWeight;
					NewDiagPartner.out[CurrSim - 1][CurrYear - StartYear] += Register[ID2ndary - 1].PopWeight;
				}
				else{ Register[ID2ndary - 1].VCThistory = 1; }
				if (Register[ID2ndary - 1].HIVstage == 1){ AcuteTests.out[CurrSim - 1][iy] += Register[ID2ndary - 1].PopWeight; }
			}
			if (Register[ID2ndary - 1].VCThistory < 2 && Condom2ndary == 0.0){
				Condom1 = 1.0 / (1.0 + (1.0 / CondomUseST[AgeGroup - 2][SexInd] - 1.0) / (CondomPref *
					TempCondomAdj / DiscloseEffectCondom));
				Condom2 = 1.0 / (1.0 + (1.0 / CondomUseST[AgeGroup - 2][SexInd] - 1.0) / (CondomPref *
					TempCondomAdj));
				ia = Register[ID2ndary - 1].AgeGroup - 2;
				ig = Register[ID2ndary - 1].SexInd;
				CondomP = 1.0 / (1.0 + (1.0 / CondomUseST[ia][ig] - 1.0) /
					(Register[ID2ndary - 1].CondomPref * Register[ID2ndary - 1].TempCondomAdj));
				if (ig == SexInd){ // MSM
					OldProb = 0.5 * (Condom1 + CondomP);
					NewProb = 0.5 * (Condom2 + CondomP);
				}
				else if (SexInd == 0){
					OldProb = (1.0 - GenderEquality) * Condom1 + GenderEquality * CondomP;
					NewProb = (1.0 - GenderEquality) * Condom2 + GenderEquality * CondomP;
				}
				else {
					OldProb = (1.0 - GenderEquality) * CondomP + GenderEquality * Condom1;
					NewProb = (1.0 - GenderEquality) * CondomP + GenderEquality * Condom2;
				}
				TestProb = (NewProb - OldProb) / (1.0 - OldProb);
				if (rc[1] < TestProb){
					Condom2ndary = 1.0;
					if (Register[ID2ndary - 1].IDprimary == ID){
						Register[ID2ndary - 1].CondomPrimary = 1.0;
					}
					else{ Register[ID2ndary - 1].Condom2ndary = 1.0; }
				}
			}
		}
	}
}

void Indiv::TestCondomChange(int ID)
{
	// This function tests whether an individual changes their condom use after a structural
	// intervention that changes either their level of binge drinking or their gender norms.
	// (I still need to extend the function to allow for the effect of economic
	// interventions that change a person's SES). The focus is on change in condom use with
	// CURRENT partners, not new partners

	int ia, ig, iy;
	double Condom1, Condom2, CondomP, ORcondomPostInt, NormAdj, AlcAdj;
	double OldProb, NewProb, TestProb;

	iy = CurrYear - StartYear;

	// First check if the individual received any of the interventions (code similar to that in the 
	// StructInterventionWaning function)
	ORcondomPostInt = 1.0;
	if ((CondomPrimary == 0 || Condom2ndary == 0) && AliveInd == 1){
		if (SexInd == 0 && ((CurrAge >= 18 && CurrAge < 50 && StructIntScenario == 6) ||
			(CurrAge >= 15 && CurrAge < 30 && StructIntScenario == 7))){
			// Inequitable norm adjustment
			NormAdj = 0.0;
			if (StructIntScenario == 6){ NormAdj = RRgenderIneqComm; }
			else if (StructIntScenario == 7){ NormAdj = RRgenderIneqIndiv; }
			ORcondomPostInt = pow(ORcondomGenderIneq, -IneqGender * (1.0 - NormAdj));
			// Not adjusting for the effect on binge-drinking, which is a second-order effect.
		}
		if (CurrAge >= 15 && DailyDrinkProb >= 1.0 / 30.0 && DrinksPerDD >= 5.0 &&
			(StructIntScenario == 1 || StructIntScenario == 2)){
			// Binge drinking adjustment
			if (StructIntScenario == 1){ AlcAdj = pow(RRalcoholSingle, 0.5); }
			else if (StructIntScenario == 2){ AlcAdj = pow(RRalcoholMultiple, 0.5); }
			if (DrinksPerDD >= 5.0){
				if (DrinksPerDD * AlcAdj < 5.0){
					ORcondomPostInt *= pow(ORcondomBingePW, -DailyDrinkProb * 7.0); }
				else{ 
					ORcondomPostInt *= pow(ORcondomBingePW, -DailyDrinkProb * (1.0 - pow(AlcAdj, 2.0)) * 7.0); }
			}
		}
	}
	// Then call additional code if their condom preference was affected by the intervention
	if (ORcondomPostInt > 1.0){
		if (CurrPartners == 1){
			if (CondomPrimary == 0.0){
				Condom1 = 1.0 / (1.0 + (1.0 / CondomUseST[AgeGroup - 2][SexInd] - 1.0) / (CondomPref *
					TempCondomAdj));
				Condom2 = 1.0 / (1.0 + (1.0 / CondomUseST[AgeGroup - 2][SexInd] - 1.0) / (CondomPref *
					TempCondomAdj * ORcondomPostInt));
				ia = Register[IDprimary - 1].AgeGroup - 2;
				ig = Register[IDprimary - 1].SexInd;
				CondomP = 1.0 / (1.0 + (1.0 / CondomUseST[ia][ig] - 1.0) /
					(Register[IDprimary - 1].CondomPref * Register[IDprimary - 1].TempCondomAdj));
				if (MarriedInd == 1){
					Condom1 *= RelEffectCondom[1];
					Condom2 *= RelEffectCondom[1];
					CondomP *= RelEffectCondom[1];
				}
				if (ig == SexInd){ // MSM
					OldProb = 0.5 * (Condom1 + CondomP);
					NewProb = 0.5 * (Condom2 + CondomP);
				}
				else if (SexInd == 0){
					OldProb = (1.0 - GenderEquality) * Condom1 + GenderEquality * CondomP;
					NewProb = (1.0 - GenderEquality) * Condom2 + GenderEquality * CondomP;
				}
				else {
					OldProb = (1.0 - GenderEquality) * CondomP + GenderEquality * Condom1;
					NewProb = (1.0 - GenderEquality) * CondomP + GenderEquality * Condom2;
				}
				TestProb = (NewProb - OldProb) / (1.0 - OldProb);
				if (r2[ID-1] < TestProb){
					CondomPrimary = 1.0;
					if (Register[IDprimary - 1].IDprimary == ID){
						Register[IDprimary - 1].CondomPrimary = 1.0;
					}
					else{ Register[IDprimary - 1].Condom2ndary = 1.0; }
				}
			}
		}
		if(CurrPartners == 2){
			if (Condom2ndary == 0.0){
				Condom1 = 1.0 / (1.0 + (1.0 / CondomUseST[AgeGroup - 2][SexInd] - 1.0) / (CondomPref *
					TempCondomAdj));
				Condom2 = 1.0 / (1.0 + (1.0 / CondomUseST[AgeGroup - 2][SexInd] - 1.0) / (CondomPref *
					TempCondomAdj * ORcondomPostInt));
				ia = Register[ID2ndary - 1].AgeGroup - 2;
				ig = Register[ID2ndary - 1].SexInd;
				CondomP = 1.0 / (1.0 + (1.0 / CondomUseST[ia][ig] - 1.0) /
					(Register[ID2ndary - 1].CondomPref * Register[ID2ndary - 1].TempCondomAdj));
				if (ig == SexInd){ // MSM
					OldProb = 0.5 * (Condom1 + CondomP);
					NewProb = 0.5 * (Condom2 + CondomP);
				}
				else if (SexInd == 0){
					OldProb = (1.0 - GenderEquality) * Condom1 + GenderEquality * CondomP;
					NewProb = (1.0 - GenderEquality) * Condom2 + GenderEquality * CondomP;
				}
				else {
					OldProb = (1.0 - GenderEquality) * CondomP + GenderEquality * Condom1;
					NewProb = (1.0 - GenderEquality) * CondomP + GenderEquality * Condom2;
				}
				TestProb = (NewProb - OldProb) / (1.0 - OldProb);
				if (r2[ID - 1] < TestProb){
					Condom2ndary = 1.0;
					if (Register[ID2ndary - 1].IDprimary == ID){
						Register[ID2ndary - 1].CondomPrimary = 1.0;
					}
					else{ Register[ID2ndary - 1].Condom2ndary = 1.0; }
				}
			}
		}
	}
}

void Indiv::UpdateChildren(int ChildID)
{
	int PrevKids;

	PrevKids = 0;
	while (ChildIDs[PrevKids] > 0){
		PrevKids += 1;
	}
	ChildIDs[PrevKids] = ChildID;
}

void Indiv::FindMother(int ID)
{
	// This function only gets called at the start of the simulation (1985).

	int ic, AgeDif;
	double MotherWeight, TotSelectionProb, CumProb;

	// Calculate sampling weights
	TotSelectionProb = 0.0;
	for (ic = 0; ic < InitPop; ic++){
		AgeDif = Register[ic].CurrAge - CurrAge;
		if (Register[ic].SexInd == 1 && Register[ic].AliveInd == 1 && AgeDif > 14 && AgeDif < 50 &&
			Register[ic].VirginInd == 0 && Register[ic].Race == Race && Register[ic].ChildIDs[19]==0 &&
			Register[ic].CurrUrban == CurrUrban){
			MotherWeight = FertilityTable[AgeDif - 15][0][Race];
			if (Register[ic].Fecundability > 0.0){ MotherWeight *= Register[ic].Fecundability; }
			if (Register[ic].Fecundability == 0.0 && Register[ic].CurrAge < 20){ MotherWeight = 0.0; }
			TotSelectionProb += MotherWeight;
		}
		Register[ic].CumSelectionProb = TotSelectionProb;
	}

	// Next select a mother
	//r[0] = rg.Random();
	//r[0] *= TotSelectionProb;
	CumProb = rpID2[ID - 1] * TotSelectionProb;
	if (TotSelectionProb > 0.0){
		for (ic = 0; ic < InitPop; ic++){
			//if (Register[ic].CumSelectionProb >= r[0]){
			if (Register[ic].CumSelectionProb >= CumProb){
				ParentID[1] = ic + 1;
				Register[ic].UpdateChildren(ID);
				break;
			}
		}
	}
	else{ ParentID[1] = 0; }
}

void Indiv::FindFather(int ID)
{
	// This function only gets called at the start of the simulation (1985).

	int ic, AgeDif;
	double FatherWeight, TotSelectionProb, CumProb;

	// Calculate sampling weights
	TotSelectionProb = 0.0;
	for (ic = 0; ic < InitPop; ic++){
		AgeDif = Register[ic].CurrAge - CurrAge;
		if (Register[ic].SexInd == 0 && Register[ic].AliveInd == 1 && AgeDif > 14 &&
			Register[ic].VirginInd == 0 && Register[ic].Race == Register[ID-1].Race &&
			Register[ic].ChildIDs[19] == 0 && Register[ic].CurrUrban == Register[ID-1].CurrUrban){
			FatherWeight = AgeEffectPartners[AgeDif - 10][0];
			FatherWeight *= (1.0 - Register[ic].MalePref);
			TotSelectionProb += FatherWeight;
		}
		Register[ic].CumSelectionProb = TotSelectionProb;
	}

	// Next select a father
	//r[0] = rg.Random();
	//r[0] *= TotSelectionProb;
	CumProb = TotSelectionProb * rpID2[ID - 1];
	if (TotSelectionProb > 0.0){
		for (ic = 0; ic < InitPop; ic++){
			//if (Register[ic].CumSelectionProb >= r[0]){
			if (Register[ic].CumSelectionProb >= CumProb){
				ParentID[0] = ic + 1;
				Register[ic].UpdateChildren(ID);
				break;
			}
		}
	}
	else{ ParentID[0] = 0; }
}

void Indiv::AssignChildToHH(int ID)
{
	int ic, MID, PID, MMID, MPID, PMID, PPID, NewHH;
	int SurvivingSiblings, AuntUncleIDs[20], TempID, ChosenID;
	double TempOdds, TempOdds2, RandHH[10];

	for (ic = 0; ic<10; ic++){
		RandHH[ic] = rg.Random(); }

	NewHH = 0;
	// First check if the mother is alive and calculate probability of staying with mother
	MID = ParentID[1];
	/*if (MID > 0 && CurrAge == 0 && Race == 0){
		if (Register[MID - 1].CurrUrban == 0){ NewHH = Register[MID - 1].HouseholdID; }
		if (Register[MID - 1].CurrUrban == 1){
			MMID = Register[MID - 1].ParentID[1];
			if (MMID > 0){
				if (Register[MMID - 1].AliveInd == 1 && Register[MMID - 1].CurrUrban == 0){
					NewHH = Register[MMID - 1].HouseholdID;
				}
			}
		}
		if (NewHH > 0){
			if (HHregister[NewHH - 1].Size < MaxHHsize){
				HHregister[NewHH - 1].AddMember(ID);
				HouseholdID = NewHH;
			}
			else{ NewHH = 0; }
		}
	}*/
	if (NewHH == 0 && MID > 0){
		if (Register[MID - 1].AliveInd == 1){
			TempOdds = BaseOddsCLSP[0];
			//TempOdds *= pow(AgeEffectCLSP[0], CurrAge);
			//TempOdds *= pow(Age2EffectCLSP[0], CurrAge * CurrAge);
			TempOdds *= RaceEffectCLSP[Race][0];
			if (InSchool == 1){ TempOdds *= SchoolingEffectCLSP[0]; }
			//if (Register[MID - 1].CurrUrban == 1 && Race == 0 && CurrAge == 0){ TempOdds *= 0.50; }
			//if (Register[MID - 1].CurrUrban == 0 && Race == 0 && CurrAge == 0){ TempOdds *= 2.00; }
			if (RandHH[0] < TempOdds / (1.0 + TempOdds)){
				NewHH = Register[MID - 1].HouseholdID;
				if (NewHH > 0){
					if (HHregister[NewHH - 1].Size < MaxHHsize && (HHregister[NewHH - 1].IDhead == 
						MID || HHregister[NewHH - 1].IDhead == Register[MID - 1].IDprimary)){
						HHregister[NewHH - 1].AddMember(ID);
						HouseholdID = NewHH;
					}
					else{ NewHH = 0; }
				}
			}
		}
	}

	// Failing that, check if father is alive and calculate probability of staying with father
	PID = ParentID[0];
	if (NewHH == 0 && PID > 0){
		if (Register[PID - 1].AliveInd == 1){
			TempOdds2 = BaseOddsCLSP[1];
			//TempOdds2 *= pow(AgeEffectCLSP[1], CurrAge);
			//TempOdds2 *= pow(Age2EffectCLSP[1], CurrAge * CurrAge);
			TempOdds2 *= RaceEffectCLSP[Race][1];
			if (InSchool == 1){ TempOdds2 *= SchoolingEffectCLSP[1]; }
			//if (Register[PID - 1].CurrUrban == 1 && Race == 0 && CurrAge == 0){ TempOdds *= 0.50; }
			//if (Register[PID - 1].CurrUrban == 0 && Race == 0 && CurrAge == 0){ TempOdds *= 2.00; }
			if (RandHH[1] < TempOdds2 / (1.0 + TempOdds2)){
				NewHH = Register[PID - 1].HouseholdID;
				if (NewHH > 0){
					if (HHregister[NewHH - 1].Size < MaxHHsize){
						HHregister[NewHH - 1].AddMember(ID);
						HouseholdID = NewHH;
					}
					else{ NewHH = 0; }
				}
			}
		}
	}

	// Failing that, check if maternal grandmother/grandfather is alive and assign them there
	// if they are  household head.
	if (MID > 0 && NewHH == 0){
		MMID = Register[MID - 1].ParentID[1];
		if (MMID > 0){
			if (Register[MMID - 1].AliveInd == 1){
				NewHH = Register[MMID - 1].HouseholdID;
				if (NewHH > 0){
					if (HHregister[NewHH - 1].Size < MaxHHsize && HHregister[NewHH - 1].IDhead == MMID){
						HHregister[NewHH - 1].AddMember(ID);
						HouseholdID = NewHH;
					}
					else{ NewHH = 0; }
				}
			}
		}
		MPID = Register[MID - 1].ParentID[0];
		if (MPID > 0 && NewHH == 0){
			if (Register[MPID - 1].AliveInd == 1){
				NewHH = Register[MPID - 1].HouseholdID;
				if (NewHH > 0){
					if (HHregister[NewHH - 1].Size < MaxHHsize && HHregister[NewHH - 1].IDhead == MPID){
						HHregister[NewHH - 1].AddMember(ID);
						HouseholdID = NewHH;
					}
					else{ NewHH = 0; }
				}
			}
		}
	}

	// Failing that, check if paternal grandmother/grandfather is alive and assign them there
	// if they are  household head.
	if (PID > 0 && NewHH == 0){
		PMID = Register[PID - 1].ParentID[1];
		if (PMID > 0){
			if (Register[PMID - 1].AliveInd == 1){
				NewHH = Register[PMID - 1].HouseholdID;
				if (NewHH > 0){
					if (HHregister[NewHH - 1].Size < MaxHHsize && HHregister[NewHH - 1].IDhead == PMID){
						HHregister[NewHH - 1].AddMember(ID);
						HouseholdID = NewHH;
					}
					else{ NewHH = 0; }
				}
			}
		}
		PPID = Register[PID - 1].ParentID[0];
		if (PPID > 0 && NewHH == 0){
			if (Register[PPID - 1].AliveInd == 1){
				NewHH = Register[PPID - 1].HouseholdID;
				if (NewHH > 0){
					if (HHregister[NewHH - 1].Size < MaxHHsize && HHregister[NewHH - 1].IDhead == PPID){
						HHregister[NewHH - 1].AddMember(ID);
						HouseholdID = NewHH;
					}
					else{ NewHH = 0; }
				}
			}
		}
	}

	// Failing that, assign them to a parent/grandparent who is NOT a household head
	if (MID > 0 && NewHH == 0){
		if (Register[MID - 1].AliveInd == 1){
			NewHH = Register[MID - 1].HouseholdID;
			if (NewHH > 0){
				if (HHregister[NewHH - 1].Size < MaxHHsize){
					HHregister[NewHH - 1].AddMember(ID);
					HouseholdID = NewHH;
				}
				else{ NewHH = 0; }
			}
		}
	}
	/*if (PID > 0 && NewHH == 0){
		if (Register[PID - 1].AliveInd == 1){
			NewHH = Register[PID - 1].HouseholdID;
			if (NewHH > 0){
				if (HHregister[NewHH - 1].Size < MaxHHsize){
					HHregister[NewHH - 1].AddMember(ID);
					HouseholdID = NewHH;
				}
				else{ NewHH = 0; }
			}
		}
	}*/
	if (MID > 0 && NewHH == 0){
		if (MMID > 0){
			if (Register[MMID - 1].AliveInd == 1){
				NewHH = Register[MMID - 1].HouseholdID;
				if (NewHH > 0){
					if (HHregister[NewHH - 1].Size < MaxHHsize){
						HHregister[NewHH - 1].AddMember(ID);
						HouseholdID = NewHH;
					}
					else{ NewHH = 0; }
				}
			}
		}
		if (MPID > 0 && NewHH == 0){
			if (Register[MPID - 1].AliveInd == 1){
				NewHH = Register[MPID - 1].HouseholdID;
				if (NewHH > 0){
					if (HHregister[NewHH - 1].Size < MaxHHsize){
						HHregister[NewHH - 1].AddMember(ID);
						HouseholdID = NewHH;
					}
					else{ NewHH = 0; }
				}
			}
		}
	}
	if (PID > 0 && NewHH == 0){
		if (PMID > 0){
			if (Register[PMID - 1].AliveInd == 1){
				NewHH = Register[PMID - 1].HouseholdID;
				if (NewHH > 0){
					if (HHregister[NewHH - 1].Size < MaxHHsize){
						HHregister[NewHH - 1].AddMember(ID);
						HouseholdID = NewHH;
					}
					else{ NewHH = 0; }
				}
			}
		}
		if (PPID > 0 && NewHH == 0){
			if (Register[PPID - 1].AliveInd == 1){
				NewHH = Register[PPID - 1].HouseholdID;
				if (NewHH > 0){
					if (HHregister[NewHH - 1].Size < MaxHHsize){
						HHregister[NewHH - 1].AddMember(ID);
						HouseholdID = NewHH;
					}
					else{ NewHH = 0; }
				}
			}
		}
	}

	// Failing that assign them to a maternal aunt/uncle
	if (MID > 0 && NewHH == 0){
		if (MMID > 0){
			SurvivingSiblings = 0;
			for (ic = 0; ic < 20; ic++){
				TempID = Register[MMID - 1].ChildIDs[ic];
				if (TempID > 0){
					if (Register[TempID - 1].AliveInd == 1 && TempID != MID &&
						Register[TempID - 1].CurrAge >= 15){
						AuntUncleIDs[SurvivingSiblings] = TempID;
						SurvivingSiblings += 1;
					}
				}
				else{ break; }
			}
			if (SurvivingSiblings > 0){
				ChosenID = SurvivingSiblings * RandHH[2];
				TempID = AuntUncleIDs[ChosenID];
				NewHH = Register[TempID - 1].HouseholdID;
				if (NewHH > 0){
					if (HHregister[NewHH - 1].Size < MaxHHsize){
						HHregister[NewHH - 1].AddMember(ID);
						HouseholdID = NewHH;
					}
					else{ NewHH = 0; }
				}
			}
		}
	}

	// Failing that assign them to a paternal aunt/uncle
	if (PID > 0 && NewHH == 0){
		if (PMID > 0){
			SurvivingSiblings = 0;
			for (ic = 0; ic < 20; ic++){
				TempID = Register[PMID - 1].ChildIDs[ic];
				if (TempID > 0){
					if (Register[TempID - 1].AliveInd == 1 && TempID != PID &&
						Register[TempID - 1].CurrAge >= 15){
						AuntUncleIDs[SurvivingSiblings] = TempID;
						SurvivingSiblings += 1;
					}
				}
				else{ break; }
			}
			if (SurvivingSiblings > 0){
				ChosenID = SurvivingSiblings * RandHH[3];
				TempID = AuntUncleIDs[ChosenID];
				NewHH = Register[TempID - 1].HouseholdID;
				if (NewHH > 0){
					if (HHregister[NewHH - 1].Size < MaxHHsize){
						HHregister[NewHH - 1].AddMember(ID);
						HouseholdID = NewHH;
					}
					else{ NewHH = 0; }
				}
			}
		}
	}

	// Failing that, assume the child is homeless/institutionalized.
	if (NewHH == 0){ 
		HouseholdID = 0; 
		DateHomeless = 1.0 * CurrYear + 0.5 + 1.0 * (BehavCycleCount - 1.0) / CycleS;
	}
}

void Indiv::ChangeHHafterMigration(int ID)
{
	int PID, HHID, NewHH, ia, ih;
	double OddsNewHH, ProbNewHH, RandHH;

	NewHH = 0; // Indicator of whether person has been assigned to a new HH yet.
	RandHH = rg.Random();

	if (MarriedInd == 1){
		PID = IDprimary;
		if (CurrUrban == Register[PID - 1].CurrUrban){
			// Migrant moves to same household as spouse if spouse is head
			HHID = Register[PID - 1].HouseholdID;
			if (HHID > 0){
				if (HHregister[HHID - 1].IDhead == PID && HHregister[HHID - 1].Size < MaxHHsize){
					NewHH = 1;
					HouseholdID = HHID;
					HHregister[HHID - 1].AddMember(ID);
				}
			}
		}
	}

	if (NewHH == 0){
		// Determine if individual forms a new household.
		ia = AgeGroup - 3;
		if (AgeGroup > 6){ AgeGroup = 6; }
		ih = 0;
		if (HighestGrade > 0 && HighestGrade <= 7){ ih = 1; }
		if (HighestGrade > 7 && HighestGrade <= 11){ ih = 2; }
		if (HighestGrade == 12){ ih = 3; }
		if (HighestGrade == 13){ ih = 4; }
		OddsNewHH = BaseOddsNewHH[SexInd] * AgeEffectNewHH[ia][SexInd] * 
			RaceEffectNewHH[Race][SexInd] * EduEffectNewHH[ih][SexInd];
		if (MarriedInd == 1){ OddsNewHH *= MarriedEffectNewHH[SexInd]; }
		if (Employed == 1){ OddsNewHH *= EmployedEffectNewHH[SexInd]; }
		if (InSchool == 1){ OddsNewHH *= InSchoolEffectNewHH[SexInd]; }
		if (CurrUrban == 1){ OddsNewHH *= UrbanEffectNewHH[SexInd]; }
		ProbNewHH = OddsNewHH / (1.0 + OddsNewHH);
		if (RandHH < ProbNewHH){
			// A new household is formed
			NewHH = 1;
			RSApop.NewHousehold(ID);
			if (MarriedInd == 1){
				if (CurrUrban == Register[PID - 1].CurrUrban){
					// Spouse moves to same household as migrant
					HHregister[HouseholdID - 1].AddMember(PID);
					Register[PID - 1].HouseholdID = HouseholdID;
					if (HHID > 0){ HHregister[HHID - 1].RemoveMember(PID); }
				}
			}
		}
	}

	// Failing that, assume they move in with non-head spouse if they're married
	if (NewHH == 0 && MarriedInd == 1){
		if (CurrUrban == Register[PID - 1].CurrUrban){
			if (HHID > 0){
				if (HHregister[HHID - 1].Size < MaxHHsize){
					NewHH = 1;
					HouseholdID = HHID;
					HHregister[HHID - 1].AddMember(ID);
				}
			}
		}
	}

	// Failing that, search for another relative to stay with, or test for homelessness
	if (NewHH == 0){ NewHH = FindRelativeToLiveWith(ID, 0); }
	if (NewHH == 0){ NewHH = TestHomeless(ID); }

	// Failing that, assume they form a new household
	if (NewHH == 0){
		NewHH = 1;
		RSApop.NewHousehold(ID);
		if (MarriedInd == 1){
			if (CurrUrban == Register[PID - 1].CurrUrban){
				// Spouse moves to same household as migrant
				HHregister[HouseholdID - 1].AddMember(PID);
				Register[PID - 1].HouseholdID = HouseholdID;
				if (HHID > 0){ HHregister[HHID - 1].RemoveMember(PID); }
			}
		}
	}
}

void Indiv::ChangeHHleavingNest(int ID)
{
	// This function only gets called for people aged 15 or older

	int LiveParent, HHID, SHID, MID, MMID, MPID, SMID, PID, PMID, PPID, SPID;
	int OldHH, LeaveInd, FoundRel;
	double LeaveProb, rand1, rand2, ProbHead;

	// First determine if the individual is living with their parent/grandparent

	LiveParent = 0;
	if (HouseholdID > 0){
		HHID = HHregister[HouseholdID - 1].IDhead;
		SHID = Register[HHID - 1].IDprimary;
		if (ParentID[0] > 0){
			if (ParentID[0] == HHID || ParentID[0] == SHID){ LiveParent = 1; }
		}
		if (LiveParent == 0 && ParentID[1] > 0){
			if (ParentID[1] == HHID || ParentID[1] == SHID){ LiveParent = 1; }
		}
		if (LiveParent == 0 && ParentID[1] > 0){
			MID = ParentID[1];
			MMID = Register[MID - 1].ParentID[1];
			MPID = Register[MID - 1].ParentID[0];
			SMID = Register[MID - 1].IDprimary;
			if (HHID == MMID || HHID == MPID || (SHID > 0 && (MMID == SHID || MPID == SHID))){
				LiveParent = 1;
			}
			else if (SMID > 0){
				if (Register[SMID - 1].ParentID[0] == HHID || Register[SMID - 1].ParentID[1] == HHID){
					LiveParent = 1;
				}
			}
		}
		if (LiveParent == 0 && ParentID[0] > 0){
			PID = ParentID[0];
			PMID = Register[PID - 1].ParentID[1];
			PPID = Register[PID - 1].ParentID[0];
			SPID = Register[PID - 1].IDprimary;
			if (HHID == PMID || HHID == PPID || (SHID > 0 && (PMID == SHID || PPID == SHID))){
				LiveParent = 1;
			}
			else if (SPID > 0){
				if (Register[SPID - 1].ParentID[0] == HHID || Register[SPID - 1].ParentID[1] == HHID){
					LiveParent = 1;
				}
			}
		}
	}

	// Next determine if they decide to leave the nest

	if (LiveParent == 1 && HouseholdID > 0){
		if (CurrAge < 50){ LeaveProb = AnnProbLeaveNest[CurrAge - 15]; }
		else{ LeaveProb = AnnProbLeaveNest[34]; }
		if (SexInd == 1){ LeaveProb *= RRleaveNestFemale; }
		if (Employed == 0 && InSchool == 0){ LeaveProb = 0.0; }
		else if (InSchool == 1 && HighestGrade >= 12){ 
			LeaveProb *= RRleaveNestInSchool / ORleaveNestEmployed; }
		if (HHregister[HouseholdID-1].Size >= 10){ LeaveProb *= RRleaveNestHH10plus; }
		LeaveProb *= RRleaveNestRace[Race];
		rand1 = rg.Random();
		if (rand1 < LeaveProb){ LeaveInd = 1; }
		else{ LeaveInd = 0; }
	}

	// Next determine which household they join if they leave the nest

	if (LiveParent == 1 && HouseholdID > 0 && LeaveInd == 1){
		OldHH = HouseholdID;
		ProbHead = BaseOddsHHALN * pow(OR_HHALNage, CurrAge) * 
			pow(OR_HHALNage2, CurrAge * CurrAge);
		if (SexInd == 1){ ProbHead *= OR_HHALNfemale; }
		if (Employed == 1){ ProbHead *= OR_HHALNemployed; }
		if (InSchool == 1){ ProbHead *= OR_HHALNinSchool; }
		ProbHead *= OR_HHALNrace[Race];
		ProbHead = ProbHead / (1.0 + ProbHead); // Converting OR to prob
		rand2 = rg.Random();
		if (rand2 < ProbHead){ RSApop.NewHousehold(ID); }
		else{
			FoundRel = FindRelativeToLiveWith(ID, 1);
			if (FoundRel == 0){ RSApop.NewHousehold(ID); }
		}
		if (OldHH > 0){ HHregister[OldHH - 1].RemoveMember(ID); }
		if (MarriedInd == 1){
			PID = IDprimary;
			if (CurrUrban == Register[PID - 1].CurrUrban &&
				OldHH == Register[PID-1].HouseholdID){
				// Spouse moves to same household 
				HHregister[HouseholdID - 1].AddMember(PID);
				Register[PID - 1].HouseholdID = HouseholdID;
				if (OldHH > 0){ HHregister[OldHH - 1].RemoveMember(PID); }
			}
		}
	}
}

void Indiv::CheckReturnToNest(int ID)
{
	int MID, PID, NewHH, TempHH, OldHH;

	NewHH = 0;
	if (MarriedInd == 0 && HouseholdID > 0 && CurrAge < 40){
		// Check if they move back to staying with their mother
		if (ParentID[1]>0){
			MID = ParentID[1];
			if (Register[MID - 1].AliveInd == 1 && Register[MID - 1].HouseholdID > 0){
				TempHH = Register[MID - 1].HouseholdID;
				if (HHregister[TempHH - 1].Size < MaxHHsize && (MID == HHregister[TempHH - 1].IDhead 
					|| Register[MID - 1].IDprimary == HHregister[TempHH - 1].IDhead) &&
					CurrUrban == Register[MID - 1].CurrUrban){
					NewHH = TempHH;
				}
			}
		}
		// Check if they move back to staying with their father
		if (ParentID[0]>0 && NewHH == 0){
			PID = ParentID[0];
			if (Register[PID - 1].AliveInd == 1 && Register[PID - 1].HouseholdID > 0){
				TempHH = Register[PID - 1].HouseholdID;
				if (HHregister[TempHH - 1].Size < MaxHHsize && (PID == HHregister[TempHH - 1].IDhead
					|| Register[PID - 1].IDprimary == HHregister[TempHH - 1].IDhead) &&
					CurrUrban == Register[PID - 1].CurrUrban){
					NewHH = TempHH;
				}
			}
		}
	}

	if (NewHH > 0){ // Move them to their parent's household
		OldHH = HouseholdID;
		HHregister[NewHH - 1].AddMember(ID);
		HouseholdID = NewHH;
		if (OldHH > 0){ HHregister[OldHH - 1].RemoveMember(ID); }
	}
}

int Indiv::FindRelativeToLiveWith(int ID, int Sibling)
{
	// This function follows a similar structure to AssignChildToHH, but 
	// with a different ordering of relatives

	// Sibling variable indicates if individual only seeks home of a sibling

	int ic, MID, PID, MMID, MPID, PMID, PPID, NewHH;
	int SurvivingSiblings, SiblingIDs[20], TempID, ChosenID;
	double RandHH[10];

	for (ic = 0; ic<10; ic++){
		RandHH[ic] = rg.Random();}
	NewHH = 0;

	// Check if individual can stay with mother
	MID = ParentID[1];
	if (MID > 0 && Sibling == 0){
		if (Register[MID - 1].AliveInd == 1 && CurrUrban == Register[MID - 1].CurrUrban &&
			Register[MID - 1].HouseholdID > 0 && Register[MID - 1].HouseholdID!=HouseholdID){
			NewHH = Register[MID - 1].HouseholdID;
			if (HHregister[NewHH - 1].Size < MaxHHsize){
				HHregister[NewHH - 1].AddMember(ID);
				HouseholdID = NewHH;
			}
			else{ NewHH = 0; }
		}
	}

	// Otherwise check if individual can stay with father
	PID = ParentID[0];
	if (PID > 0 && NewHH == 0 && Sibling == 0){
		if (Register[PID - 1].AliveInd == 1 && CurrUrban == Register[PID - 1].CurrUrban &&
			Register[PID - 1].HouseholdID > 0 && Register[PID - 1].HouseholdID != HouseholdID){
			NewHH = Register[PID - 1].HouseholdID;
			if (HHregister[NewHH - 1].Size < MaxHHsize){
				HHregister[NewHH - 1].AddMember(ID);
				HouseholdID = NewHH;
			}
			else{ NewHH = 0; }
		}
	}

	// Otherwise check if there is a maternal sibling the individual can stay with
	if (MID > 0 && NewHH == 0){
		SurvivingSiblings = 0;
		for (ic = 0; ic < 20; ic++){
			TempID = Register[MID - 1].ChildIDs[ic];
			if (TempID > 0){
				if (Register[TempID - 1].AliveInd == 1 && TempID != ID &&
					Register[TempID - 1].CurrAge >= 15 && CurrUrban == Register[TempID - 1].CurrUrban){
					SiblingIDs[SurvivingSiblings] = TempID;
					SurvivingSiblings += 1;
				}
			}
			else{ break; }
		}
		if (SurvivingSiblings > 0){
			ChosenID = SurvivingSiblings * RandHH[2];
			TempID = SiblingIDs[ChosenID];
			NewHH = Register[TempID - 1].HouseholdID;
			if (NewHH > 0){
				if (HHregister[NewHH - 1].Size < MaxHHsize){
					HHregister[NewHH - 1].AddMember(ID);
					HouseholdID = NewHH;
				}
				else{ NewHH = 0; }
			}
		}
	}

	// Otherwise check if there is a paternal sibling the individual can stay with
	if (PID > 0 && NewHH == 0){
		SurvivingSiblings = 0;
		for (ic = 0; ic < 20; ic++){
			TempID = Register[PID - 1].ChildIDs[ic];
			if (TempID > 0){
				if (Register[TempID - 1].AliveInd == 1 && TempID != ID &&
					Register[TempID - 1].CurrAge >= 15 && CurrUrban == Register[TempID - 1].CurrUrban){
					SiblingIDs[SurvivingSiblings] = TempID;
					SurvivingSiblings += 1;
				}
			}
			else{ break; }
		}
		if (SurvivingSiblings > 0){
			ChosenID = SurvivingSiblings * RandHH[3];
			TempID = SiblingIDs[ChosenID];
			NewHH = Register[TempID - 1].HouseholdID;
			if (NewHH > 0){
				if (HHregister[NewHH - 1].Size < MaxHHsize){
					HHregister[NewHH - 1].AddMember(ID);
					HouseholdID = NewHH;
				}
				else{ NewHH = 0; }
			}
		}
	}

	// Failing that, check if there is a child the individual can stay with
	if (NewHH == 0 && Sibling == 0){
		SurvivingSiblings = 0;
		for (ic = 0; ic < 20; ic++){
			TempID = ChildIDs[ic];
			if (TempID > 0){
				if (Register[TempID - 1].AliveInd == 1 && CurrUrban == Register[TempID - 1].CurrUrban){
					SiblingIDs[SurvivingSiblings] = TempID;
					SurvivingSiblings += 1;
				}
			}
			else{ break; }
		}
		if (SurvivingSiblings > 0){
			ChosenID = SurvivingSiblings * RandHH[4];
			TempID = SiblingIDs[ChosenID];
			NewHH = Register[TempID - 1].HouseholdID;
			if (NewHH > 0){
				if (HHregister[NewHH - 1].Size < MaxHHsize){
					HHregister[NewHH - 1].AddMember(ID);
					HouseholdID = NewHH;
				}
				else{ NewHH = 0; }
			}
		}
	}

	// Failing that check if there is a maternal grandparent the person can stay with
	if (MID > 0 && NewHH == 0 && Sibling == 0){
		MMID = Register[MID - 1].ParentID[1];
		if (MMID > 0){
			if (Register[MMID - 1].AliveInd == 1 && CurrUrban == Register[MMID - 1].CurrUrban){
				NewHH = Register[MMID - 1].HouseholdID;
				if (NewHH > 0){
					if (HHregister[NewHH - 1].Size < MaxHHsize){
						HHregister[NewHH - 1].AddMember(ID);
						HouseholdID = NewHH;
					}
					else{ NewHH = 0; }
				}
			}
		}
		MPID = Register[MID - 1].ParentID[0];
		if (MPID > 0 && NewHH == 0){
			if (Register[MPID - 1].AliveInd == 1 && CurrUrban == Register[MPID - 1].CurrUrban){
				NewHH = Register[MPID - 1].HouseholdID;
				if (NewHH > 0){
					if (HHregister[NewHH - 1].Size < MaxHHsize){
						HHregister[NewHH - 1].AddMember(ID);
						HouseholdID = NewHH;
					}
					else{ NewHH = 0; }
				}
			}
		}
	}

	// Failing that check if there is a paternal grandparent the person can stay with
	if (PID > 0 && NewHH == 0 && Sibling == 0){
		PMID = Register[PID - 1].ParentID[1];
		if (PMID > 0){
			if (Register[PMID - 1].AliveInd == 1 && CurrUrban == Register[PMID - 1].CurrUrban){
				NewHH = Register[PMID - 1].HouseholdID;
				if (NewHH > 0){
					if (HHregister[NewHH - 1].Size < MaxHHsize){
						HHregister[NewHH - 1].AddMember(ID);
						HouseholdID = NewHH;
					}
					else{ NewHH = 0; }
				}
			}
		}
		PPID = Register[PID - 1].ParentID[0];
		if (PPID > 0 && NewHH == 0){
			if (Register[PPID - 1].AliveInd == 1 && CurrUrban == Register[PPID - 1].CurrUrban){
				NewHH = Register[PPID - 1].HouseholdID;
				if (NewHH > 0){
					if (HHregister[NewHH - 1].Size < MaxHHsize){
						HHregister[NewHH - 1].AddMember(ID);
						HouseholdID = NewHH;
					}
					else{ NewHH = 0; }
				}
			}
		}
	}

	if (NewHH > 0){ NewHH = 1; }
	return NewHH;
}

int Indiv::TestHomeless(int ID)
{
	// This function determines whether an individual becomes homeless.

	int Homeless;
	double RandHH;

	RandHH = rg.Random();
	if ((RandHH < ProbHomeless && Employed == 0 && MarriedInd==0) || CurrAge < 15 || HouseholdID == 0){ 
		Homeless = 1; 
		if (HouseholdID > 0){
			// Don't call RemoveMember here because it already gets called in  
			// ChangeHHafterMigration, ChangeHHafterDivorce and UpdateIncarceration.
			HouseholdID = 0;
			DateHomeless = 1.0 * CurrYear + 0.5 + 1.0 * (BehavCycleCount - 1.0) / CycleS;
		}
	}
	else{ Homeless = 0; }

	return Homeless;
}

HouseholdGroup::HouseholdGroup(int head)
{
	IDhead = head;
	Size = 1;
	Members[0] = IDhead;
	Urban = Register[head - 1].CurrUrban;
	Active = 1;
	PerCapitaIncome = MeanLogIncome;
	PerCapitaIncomeAdj = MeanLogIncome;
}

void HouseholdGroup::AddMember(int ID)
{
	Members[Size] = ID;
	Size += 1;
	Register[ID - 1].CurrUrban = Urban;
}

void HouseholdGroup::RemoveMember(int ID)
{
	int ii, ij, ig, ir, ia, ih, NewHead, TempID, OldHH, NewHH;
	double OddsHead, MaxOdds;

	ij = 0;
	for (ii = 0; ii < Size; ii++){
		if (Members[ii] == ID){
			ij = ii;
			break;
		}
	}
	if (ij < Size - 1){
		for (ii = ij; ii < Size - 1; ii++){
			Members[ii] = Members[ii + 1];
		}
	}
	Members[Size - 1] = 0;
	Size = Size - 1;
	if (Size == 0){ Active = 0; }
	MaxOdds = 0.0;

	// Handle the case where removed member was head
	if (IDhead == ID && Size > 0){
		// Find the adult with the highest baseline odds of being the household head.
		for (ii = 0; ii < Size; ii++){
			TempID = Members[ii];
			if (Register[TempID - 1].CurrAge >= 15){
				ig = Register[TempID - 1].SexInd;
				ir = Register[TempID - 1].Race;
				ia = Register[TempID - 1].AgeGroup - 3;
				if (ia > 10){ ia = 10; }
				ih = 0;
				if (Register[TempID - 1].HighestGrade > 0 && Register[TempID - 1].HighestGrade <= 7){ ih = 1; }
				if (Register[TempID - 1].HighestGrade > 7 && Register[TempID - 1].HighestGrade < 12){ ih = 2; }
				if (Register[TempID - 1].HighestGrade == 12){ ih = 3; }
				if (Register[TempID - 1].HighestGrade == 13){ ih = 4; }
				OddsHead = BaseOddsHead[ig];
				if (Register[TempID - 1].MarriedInd == 1){ OddsHead *= MarriedEffectHead[ig]; }
				OddsHead *= AgeEffectHead[ia][ig];
				OddsHead *= RaceEffectHead[ir][ig];
				if (Register[TempID - 1].Employed == 1){ OddsHead *= EmployedEffectHead[ig]; }
				if (Register[TempID - 1].InSchool == 1){ OddsHead *= InSchoolEffectHead[ig]; }
				OddsHead *= EduEffectHead[ih][ig];
				if (ig == 1 && Register[TempID - 1].DOLB > 0.0){ OddsHead *= BirthEffectHead; }
				if (OddsHead > MaxOdds){
					IDhead = TempID;
					MaxOdds = OddsHead;
				}
			}
		}
		// Handle the case where there is no surviving adult
		if (MaxOdds == 0){
			for (ii = 0; ii < Size; ii++){
				TempID = Members[ii];
				Register[TempID - 1].AssignChildToHH(TempID);
				//Register[TempID - 1].HouseholdID = 0;
				Members[ii] = 0;
			}
			Size = 0;
			Active = 0;
		}
	}

	// Reassign children of the departed household member
	if (Size >= 1){
		for (ii = 0; ii < Size; ii++){
			TempID = Members[ii];
			if (Register[TempID - 1].ParentID[1] == ID && Register[TempID - 1].CurrAge < 15){
				OldHH = Register[TempID - 1].HouseholdID;
				Register[TempID - 1].AssignChildToHH(TempID);
				if (OldHH > 0){ HHregister[OldHH - 1].RemoveMember(TempID); }
			}
			/*else if (Register[TempID - 1].ParentID[1] == ID && Register[TempID - 1].MarriedInd==0
				&& Register[ID - 1].AliveInd == 1){
				OldHH = Register[TempID - 1].HouseholdID;
				NewHH = Register[ID - 1].HouseholdID;
				if (NewHH > 0 && HHregister[OldHH-1].IDhead!=TempID){
					HHregister[NewHH - 1].AddMember(TempID);
					Register[TempID - 1].HouseholdID = NewHH;
					HHregister[OldHH - 1].RemoveMember(TempID);
				}
			}*/
		}
	}
}

void HouseholdGroup::GetPerCapitaIncome(int hhid)
{
	int ic, ii, ih, ia, ir, ig, iy, ID, CID, CIDinHH, PID, tempID;
	double temp, r1, Income[MaxHHsize], PartnerIncome;

	PerCapitaIncome = 0.0;
	iy = CurrYear - StartYear;
	Kids = 0;

	// First calculate salary income
	for (ic = 0; ic < Size; ic++){
		ID = Members[ic];
		if (Register[ID - 1].AliveInd == 1 && Register[ID - 1].CurrAge < 15){ Kids += 1; }
		if (Register[ID - 1].AliveInd == 1 && Register[ID - 1].Employed == 1){
			ia = Register[ID - 1].AgeGroup - 3;
			if (ia > 9){ ia = 9; }
			ig = Register[ID - 1].SexInd;
			ir = Register[ID - 1].Race;
			if (Register[ID - 1].HighestGrade == 0){ ih = 0; }
			else if (Register[ID - 1].HighestGrade <= 7){ ih = 1; }
			else if (Register[ID - 1].HighestGrade <= 11){ ih = 2; }
			else if (Register[ID - 1].HighestGrade == 12){ ih = 3; }
			else { ih = 4; }
			Income[ic] = exp(CurrAveIncome[ia][ig][ir][ih] + Register[ID - 1].LogIncomeDif);
			PerCapitaIncome += Income[ic];
		}
		else{ Income[ic] = 0.0; }
	}

	// Next add private pensions
	for (ic = 0; ic < Size; ic++){
		ID = Members[ic];
		if (Register[ID - 1].AliveInd == 1 && Register[ID - 1].CurrAge >= 60){
			// First check if new 60-year olds receive private pension
			if (Register[ID - 1].CurrAge == 60){
				r1 = rg.Random();
				ir = Register[ID - 1].Race;
				if (Register[ID - 1].HighestGrade == 0){ ih = 0; }
				else if (Register[ID - 1].HighestGrade <= 7){ ih = 1; }
				else if (Register[ID - 1].HighestGrade <= 11){ ih = 2; }
				else if (Register[ID - 1].HighestGrade == 12){ ih = 3; }
				else { ih = 4; }
				temp = BaseOddsPension * EffectEduPension[ih] * EffectRacePension[ir];
				temp = temp / (1.0 + temp);
				if (r1 < temp || ir == 2){ Register[ID - 1].PrivatePension = 1; }
			}
			if (Register[ID - 1].PrivatePension == 1){ 
				PerCapitaIncome += CurrPrivatePension; 
				Income[ic] += CurrPrivatePension;
			}
		}
	}

	// Then add state pension
	for (ic = 0; ic < Size; ic++){
		ID = Members[ic];
		if (Register[ID - 1].AliveInd == 1 && Register[ID - 1].CurrAge >= 60 &&
			Register[ID - 1].PrivatePension == 0){
			if (Register[ID - 1].SexInd == 1){ 
				PerCapitaIncome += OAPamount[iy]; 
				Income[ic] += OAPamount[iy];
			}
			else if (Register[ID - 1].CurrAge >= 65 || CurrYear >= 2010){ 
				PerCapitaIncome += OAPamount[iy]; 
				Income[ic] += OAPamount[iy];
			}
		}
	}

	// Finally add the child support grant
	temp = 0.0;
	for (ic = 0; ic < Size; ic++){
		ID = Members[ic];
		if (Register[ID - 1].AliveInd == 1 && Register[ID - 1].CurrAge < CSGageLimit[iy]){
			if (CurrYear < 2008 && PerCapitaIncome < CSGincomeLimit[iy]){ temp += CSGamount[iy]; }
			if (CurrYear >= 2008){
				// Find the ID of the caregiver (CID)
				CID = 0;
				for (ii = 0; ii < Size; ii++){
					tempID = Members[ii];
					if (Register[tempID - 1].CurrAge >= 15){ 
						if (Register[ID - 1].ParentID[1] == tempID){
							CID = tempID;
							CIDinHH = ii;
							break;
						}
						else if (CID == 0){ 
							CID = tempID; 
							CIDinHH = ii;
						}
						else if (Register[tempID - 1].SexInd == 1){
							if (Register[CID - 1].SexInd == 0 || Register[CID - 1].CurrAge >
								Register[tempID - 1].CurrAge){
								CID = tempID;
								CIDinHH = ii;
							}
						}
						else{
							if (Register[CID - 1].SexInd == 0 && Register[CID - 1].CurrAge >
								Register[tempID - 1].CurrAge){
								CID = tempID;
								CIDinHH = ii;
							}
						}
					}
				}
				// Then check if the means test applies
				if (Register[CID - 1].MarriedInd == 0 && Income[CIDinHH] < CSGincomeLimit[iy]){
					PerCapitaIncome += CSGamount[iy];
				}
				else if (Register[CID - 1].MarriedInd == 1){
					PID = Register[CID - 1].IDprimary;
					PartnerIncome = 0.0;
					for (ii = 0; ii < Size; ii++){
						if (Members[ii] == PID){
							PartnerIncome = Income[ii];
							break;
						}
					}
					if (Income[CIDinHH] + PartnerIncome < CSGincomeLimit[iy] * 2.0){
						PerCapitaIncome += CSGamount[iy];
					}
				}
			}
		}
	}
	PerCapitaIncome += temp;
	
	// Finally calculate per capita income 
	if (Size > 0){ 
		temp = PerCapitaIncome;
		PerCapitaIncome = temp / Size; 
		PerCapitaIncomeAdj = temp / pow(Size - Kids * (1.0 - ChildWeightEquivScale),
			HHsizeAdjEquivScale);
		if (temp == 0){ 
			if (CurrYear > StartYear){ PerCapitaIncomeAdj = exp(MeanLogIncome - 3.0); }
			else{ PerCapitaIncomeAdj = 1.0; }
		}
	}
}

void HouseholdGroup::GetUnsortedIncomes()
{
	int ic, ii, ir, ID;

	if (Size > 0){
		for (ic = 0; ic < Size; ic++){
			ID = Members[ic];
			if (Register[ID - 1].AliveInd > 0){
				ir = Register[ID - 1].Race;
				ii = TotalIncomes[0];
				UnsortedIncomes[ii][0] = PerCapitaIncome;
				TotalIncomes[0] += 1;
				ii = TotalIncomes[ir + 1];
				UnsortedIncomes[ii][ir + 1] = PerCapitaIncome;
				TotalIncomes[ir + 1] += 1;
			}
		}
	}
}

Pop::Pop(int xxx){}

void Pop::AssignAgeSex()
{
	int ic, ir, ia, column;
	double TotalPop;
	double AgeExact;
	vector<Indiv>::const_iterator ic2;
	
	TotalPop = 0.0;
	for(ia=0; ia<91; ia++){
		for (ir = 0; ir < 6; ir++){
			TotalPop += StartPop[ia][ir];
		}
	}
	StartPop[0][0] = StartPop[0][0]/TotalPop;
	for(ia=1; ia<91; ia++){
		StartPop[ia][0] = StartPop[ia - 1][0] + StartPop[ia][0] / TotalPop;}
	for (ir = 1; ir < 6; ir++){
		StartPop[0][ir] = StartPop[90][ir-1] + StartPop[0][ir] / TotalPop;
		for (ia = 1; ia < 91; ia++){
			StartPop[ia][ir] = StartPop[ia - 1][ir] + StartPop[ia][ir] / TotalPop;}
	}

	//int seed = 5481; // Note that I've arbitrarily fixed the seed for now.
	//if(FixedPeriod==0){
	//	seed += CurrSim;}
	//CRandomMersenne rg(seed);
	for(ir=0; ir<InitPop; ir++){
		r[ir] = rg.Random();
	}

	int tpp = Register.size();

	// Simulate sex, race and date of birth for each individual, assign PopWeight
	for(ic=0; ic<tpp; ic++){
		Register[ic].AliveInd = 1;
		Register[ic].PopWeight = 1.0;
		for (ir = 0; ir < 6; ir++){
			if (r[ic] <= StartPop[90][ir]){
				column = ir;
				break;
			}
		}
		if (column == 0 || column == 2 || column == 4){
			Register[ic].SexInd = 0;}
		else{ Register[ic].SexInd = 1; }
		if (column < 2){ Register[ic].Race = 0; }
		else if (column < 4){ Register[ic].Race = 1; }
		else{ Register[ic].Race = 2; }
		if(r[ic] <= StartPop[0][column]){
			Register[ic].DOB = 0.5 + StartYear - r[ic]/StartPop[0][column];}
		else for(ir=1; ir<91; ir++){
			if (r[ic]>StartPop[ir - 1][column] && r[ic] <= StartPop[ir][column]){
				Register[ic].DOB = 0.5 + StartYear - ir - (r[ic] - 
					StartPop[ir - 1][column]) / (StartPop[ir][column] - 
					StartPop[ir - 1][column]);
				break;
			}
		}
	}

	// Calculate age group for each individual
	for(ic=0; ic<InitPop; ic++){
		AgeExact = StartYear + 0.5 - Register[ic].DOB;
		if(AgeExact < 90.0){
			Register[ic].CurrAge = AgeExact;
			Register[ic].AgeGroup = AgeExact/5;
		}
		else{
			Register[ic].CurrAge = 90;
			Register[ic].AgeGroup = 17;
		}
	}

	// Calculate male circumcision status
	for (ir = 0; ir<InitPop; ir++){
		r2[ir] = rg.Random();
	}
	for (ic = 0; ic < InitPop; ic++){
		if (Register[ic].SexInd == 0){
			ia = Register[ic].CurrAge;
			ir = Register[ic].Race;
			if (r2[ic] < MCprevBaseline[ia][ir]){Register[ic].CircInd = 1;}
			else{ Register[ic].CircInd = 0; }
		}
	}
}

void Pop::AssignPersonality()
{
	int ir, ic, ind, s;
	double x, a, b, p, q, bound;

	for (ir = 0; ir<InitPop; ir++){
		r[ir] = rg.Random();
	}

	a = 0.0;
	b = 1.0;
	ind = 2;
	s = 0;
	bound = 0.0;
	for (ic = 0; ic < InitPop; ic++){
		p = r[ic];
		q = 1.0 - r[ic];
		cdfnor(&ind, &p, &q, &x, &a, &b, &s, &bound);
		Register[ic].Conscientiousness = x;
	}
}

void Pop::GetAllPartnerRates()
{
	int ia, ig;
	double lambda, alpha;
	double GammaScaling[2];

	// Calculate the age adjustment factors and the high risk partner acquisition rates
	for(ig=0; ig<2; ig++){
		lambda = (GammaMeanST[ig] - 10.0)/pow(GammaStdDevST[ig], 2.0);
		alpha = (GammaMeanST[ig] - 10.0) * lambda;
		for(ia=0; ia<81; ia++){
			AgeEffectPartners[ia][ig] = pow(lambda, alpha) * pow(ia + 0.5, alpha - 1.0) *
				exp(-lambda * (ia + 0.5));
		}
		GammaScaling[ig] = BasePartnerAcqH[ig]/AgeEffectPartners[7][ig];
		for(ia=0; ia<81; ia++){
			AgeEffectPartners[ia][ig] *= GammaScaling[ig];}
		PartnershipFormation[0][ig] = 1.0;
	}
}

void Pop::AssignBehav()
{
	int ic, ia, ib, ig, ii, ij, ir, ih;
	int ExactAge, found, index1, s;
	double MarriageProb, High;
	int UnpartneredCount; // # married individuals not yet assigned a partner
	double TotalSelectionProb, partnerrisk, partnerID, temp1;
	double CasualEntryMSM, CasualExitMSM;
	int STpartners, ind, PID, TempTot;
	double x, y, a, b, p, q, bound;
	vector<int> indices(InitPop);

	// In the arrays defined below, the I- prefix indicates individual ages, to be
	// consistent with the TSHISA model. Age starts at 10 (not 0), since 10 is the
	// assumed youngest age of sexual activity.
	double IVirginPropnH[81][2][3], IVirginPropnL[81][2][3];
	double IVirginsH[81][2][3], IVirginsL[81][2][3];
	double ISexuallyExpH[81][2][3], ISexuallyExpL[81][2][3];
	double ISexuallyExpUnmarriedH[81][2][3], ISexuallyExpUnmarriedL[81][2][3];

	double cHD, cLD; // Rates of partnership formation x ave duration in high & low risk
	double GUnmarriedH0[16][2], GUnmarriedH1[16][2], GUnmarriedH2[16][2]; 
	double GUnmarriedL0[16][2], GUnmarriedL1[16][2];
	double GMarriedH1[16][2], GMarriedH2[16][2];

	double alpha, lambda; // Parameters for the gamma dbn used to determine relative
						  // frequencies of FSW contact at different ages
	double FSWcontactBase[81], EligibleForCSW[16];
	double BaseFSWdemand, MalePop15to49, FemPref;

	// First calculate IVirginPropn
	for (ig = 0; ig < 2; ig++){
		for (ia = 0; ia < 81; ia++){
			// Note that we use ia+1 in this equation, not ia, because ia is age last birthday at the 
			// START of the year and we want the average exact age over the course of the year.
			SexualDebut[ia][ig] = DebutShape[ig] / ((ia + 1.0) * (1.0 + pow((ia + 1.0)/(DebutMedian[ig] - 
				10.0), -DebutShape[ig])));
		}
	}
	for (ir = 0; ir < 3; ir++){
		IVirginPropnH[0][0][ir] = 1.0;
		IVirginPropnH[0][1][ir] = 1.0;
		IVirginPropnL[0][0][ir] = 1.0;
		IVirginPropnL[0][1][ir] = 1.0;
	}
	for(ia=1; ia<81; ia++){
		for (ir = 0; ir < 3; ir++){
			for (ig = 0; ig < 2; ig++){
				IVirginPropnH[ia][ig][ir] = IVirginPropnH[ia - 1][ig][ir] *
					exp(-SexualDebut[ia - 1][ig] * DebutAdjRace[ir]);
				IVirginPropnL[ia][ig][ir] = IVirginPropnL[ia - 1][ig][ir] *
					exp(-SexualDebut[ia - 1][ig] * DebutAdjLow[ig] * DebutAdjRace[ir]);
			}
		}
	}

	// Assign sexual orientation & role preference to each individual
	int tpp = Register.size();
	for (ic = 0; ic<tpp; ic++){ r[ic] = rg.Random();}
	for (ic = 0; ic<tpp; ic++){ r2[ic] = rg.Random(); }
	for (ic = 0; ic<tpp; ic++){ rprisk[ic] = rg.Random(); }
	for (ic = 0; ic < tpp; ic++){
		Register[ic].CasualInd = 0;
		Register[ic].HetCasualInd = 0;
		if (Register[ic].SexInd == 0){
			if (r[ic]>MSMfraction){ Register[ic].MalePref = 0.0; }
			else{
				if (r[ic]>(MSMfraction*BiFraction)){ 
					Register[ic].MalePref = 1.0; 
					if (r2[ic] < RolePrefHomo[0]){ Register[ic].InsertivePref = 0.0; }
					else if (r2[ic] < RolePrefHomo[0] + RolePrefHomo[1]){
						Register[ic].InsertivePref = 0.5;
					}
					else{ Register[ic].InsertivePref = 1.0; }
				}
				else{
					// Assign initial male preference
					a = InitMalePrefBeta[0];
					b = InitMalePrefBeta[1];
					p = r[ic] / (MSMfraction*BiFraction);
					q = 1 - p;
					bound = 0.0;
					ind = 2;
					cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
					Register[ic].MalePref = x; 
					// Assign annual change in male preference
					a = AnnChangeMalePref[0];
					b = AnnChangeMalePref[1];
					p = rprisk[ic];
					q = 1.0 - rprisk[ic];
					s = 0;
					cdfnor(&ind, &p, &q, &x, &a, &b, &s, &bound);
					Register[ic].ChangeMalePref = x;
					// Calculate CURRENT male preference
					ExactAge = StartYear + 0.5 - Register[ic].DOB;
					if (ExactAge > 20){
						Register[ic].MalePref += Register[ic].ChangeMalePref * (ExactAge - 20.0);
						if (Register[ic].MalePref < 0.001){ Register[ic].MalePref = 0.001; }
						if (Register[ic].MalePref > 0.999){ Register[ic].MalePref = 0.999; }
					}
					// Assign role preference
					if (r2[ic] < RolePrefBi[0]){ Register[ic].InsertivePref = 0.0; }
					else if (r2[ic] < RolePrefBi[0] + RolePrefBi[1]){
						Register[ic].InsertivePref = 0.5;
					}
					else{ Register[ic].InsertivePref = 1.0; }
				}
			}
		}
		else{Register[ic].MalePref = 1.0;} // Currently only modelling same-sex preference in men
	}

	// Now assign risk group to each individual
	
	for(ic=0; ic<tpp; ic++){
		r[ic] = rg.Random(); }
	for (ic = 0; ic<tpp; ic++){
		r2[ic] = rg.Random(); }

	for(ic=0; ic<tpp; ic++){
		if(Register[ic].SexInd==0){
			temp1 = 1.0 / (1.0 + ((1.0 - HighPropnM) / HighPropnM) / pow(ConscientiousEffectSex, 
				Register[ic].Conscientiousness));	
		}
		else{
			temp1 = 1.0 / (1.0 + ((1.0 - HighPropnF) / HighPropnF) / pow(ConscientiousEffectSex,
				Register[ic].Conscientiousness));
		}
		if (r[ic] < temp1){
			Register[ic].RiskGroup = 1;
		}
		else{
			Register[ic].RiskGroup = 2;
		}
		if (r2[ic] < ConfoundingAlcSex){ Register[ic].DrinksPerDDconstant = 1.0 - r[ic]; }
		if (r2[ic] < ConfoundingAlcSex){ Register[ic].DrinkProbConstant = 1.0 - r[ic]; }
	}

	// Next assign virgin status to each individual

	for(ic=0; ic<tpp; ic++){
		r[ic] = rg.Random();
	}

	for(ic=0; ic<tpp; ic++){
		if(Register[ic].AgeGroup < 2){
			Register[ic].VirginInd = 1;}
		else if(Register[ic].AgeGroup >= 6){
			Register[ic].VirginInd = 0;}
		else{
			ExactAge = StartYear + 0.5 - Register[ic].DOB - 10;
			ig = Register[ic].SexInd;
			ir = Register[ic].Race;
			if(Register[ic].RiskGroup==1){
				if(r[ic] < IVirginPropnH[ExactAge][ig][ir]){
					Register[ic].VirginInd = 1;}
				else{
					Register[ic].VirginInd = 0;}
			}
			else{
				if (r[ic] < IVirginPropnL[ExactAge][ig][ir]){
					Register[ic].VirginInd = 1;}
				else{
					Register[ic].VirginInd = 0;}
			}
		}
	}

	// Next assign marital status to each individual

	for(ic=0; ic<tpp; ic++){
		r[ic] = rg.Random();
	}

	for(ic=0; ic<tpp; ic++){
		if(Register[ic].VirginInd==1){
			Register[ic].MarriedInd = 0;}
		else{
			ExactAge = StartYear + 0.5 - Register[ic].DOB;
			ig = Register[ic].SexInd;
			ir = Register[ic].Race;
			if (Register[ic].HighestGrade == 13){ ih = 4; }
			else if (Register[ic].HighestGrade == 12){ ih = 3; }
			else if (Register[ic].HighestGrade >= 8){ ih = 2; }
			else if (Register[ic].HighestGrade >= 3){ ih = 1; }
			else{ ih = 0; }
			MarriageProb = BaseOddsMarriage[ig] * OR_MarriedEdu[ih][ig] * OR_MarriedRace[ir][ig];
			if (Register[ic].Employed == 1){
				MarriageProb *= OR_MarriedEmployed[ig];
				if (ExactAge < 65){ MarriageProb *= pow(OR_MarriedAgeEmployed[ig], 65.0 - ExactAge); }
			}
			MarriageProb *= pow(OR_MarriedAge[0][ig], ExactAge - 15.0) * pow(OR_MarriedAge[1][ig],
				(ExactAge - 15.0) * (ExactAge - 15.0)) * pow(OR_MarriedAge[2][ig], (ExactAge - 15.0) *
				(ExactAge - 15.0) * (ExactAge - 15.0));
			MarriageProb = MarriageProb / (1.0 + MarriageProb);
			if(ExactAge<30){
				if(Register[ic].RiskGroup==1){
					MarriageProb = MarriageProb / (1.0 - IVirginPropnH[ExactAge-10][ig][ir]);}
				else{
					MarriageProb = MarriageProb / (1.0 - IVirginPropnL[ExactAge-10][ig][ir]);}
			}
			if (Register[ic].SexInd == 0 && Register[ic].MalePref == 1.0){
				MarriageProb = 1.0 - pow(1.0 - MarriageProb, SameSexMarried);
			}
			if(r[ic] < MarriageProb){
				Register[ic].MarriedInd = 1;}
			else{
				Register[ic].MarriedInd = 0;}
		}
	}

	// Next calculate sexual mixing in married individuals (copied from code in TSHISA)

	DesiredMarriagesM[0][0] = ((1.0 - AssortativeM) + AssortativeM * HighPropnF) * HighPropnM;
	DesiredMarriagesM[0][1] = AssortativeM * (1.0 - HighPropnF) * HighPropnM;
	DesiredMarriagesM[1][1] = ((1.0 - AssortativeM) + AssortativeM * (1.0 - HighPropnF)) *
		(1.0 - HighPropnM);
	DesiredMarriagesM[1][0] = AssortativeM * HighPropnF * (1 - HighPropnM);
	DesiredMarriagesF[0][0] = ((1.0 - AssortativeF) + AssortativeF * HighPropnM) * HighPropnF;
	DesiredMarriagesF[0][1] = AssortativeF * (1.0 - HighPropnM) * HighPropnF;
	DesiredMarriagesF[1][1] = ((1.0 - AssortativeF) + AssortativeF * (1.0 - HighPropnM)) *
		(1.0 - HighPropnF);
	DesiredMarriagesF[1][0] = AssortativeF * HighPropnM * (1 - HighPropnF);

	for(ii=0; ii<2; ii++){
		for(ij=0; ij<2; ij++){
			if(AllowBalancing==1){
				AdjLTrateM[ii][ij] = (GenderEquality * DesiredMarriagesF[ij][ii] +
					(1.0 - GenderEquality) * DesiredMarriagesM[ii][ij]);
				AdjLTrateF[ii][ij] = (GenderEquality * DesiredMarriagesF[ii][ij] + 
					(1.0 - GenderEquality) * DesiredMarriagesM[ij][ii]);
			}
			else{
				AdjLTrateM[ii][ij] = DesiredMarriagesM[ii][ij];
				AdjLTrateF[ii][ij] = DesiredMarriagesF[ii][ij];
			}
		}
		ActualPropnLTH[ii][0] = AdjLTrateM[ii][0]/(AdjLTrateM[ii][0] + AdjLTrateM[ii][1]);
		ActualPropnLTH[ii][1] = AdjLTrateF[ii][0]/(AdjLTrateF[ii][0] + AdjLTrateF[ii][1]);
	}
	ActualPropnLTH_MSM[0] = (1.0 - AssortativeM) + AssortativeM * HighPropnM;
	ActualPropnLTH_MSM[1] = AssortativeM * HighPropnM;

	// Assign partner IDs to each married individual

	for(ic=0; ic<tpp; ic++){
		r[ic] = rg.Random();}
	for(ic=0; ic<tpp; ic++){
		rprisk[ic] = rg.Random();}
	for(ic=0; ic<tpp; ic++){
		rpID[ic] = rg.Random();}

	MaxNewPartnerInd = 0;
	//for(ii=0; ii<InitPop; ii++){
	//	indices[ii] = ii;}
	ii = 0;
	while(MaxNewPartnerInd==0 && ii<tpp){
		/*index1 = indices[r[ii] * indices.size()];
		ic = indices[index1];
		if(Register[ic].MarriedInd==1 && Register[ic].IDprimary==0){
			partnerrisk = rprisk[ii];
			partnerID = rpID[ii];
			ChooseLTpartner(ic+1, partnerrisk, partnerID);
		}
		indices[index1] = indices[indices.size()-1];
		indices.pop_back();
		ii += 1;*/
		UnpartneredCount = 0;
		for(ic=0; ic<tpp; ic++){
			if(Register[ic].MarriedInd==1 && Register[ic].IDprimary==0){
				UnpartneredCount += 1;}
		}
		TotalSelectionProb = 0.0;
		for(ic=0; ic<tpp; ic++){
			if(Register[ic].MarriedInd==1 && Register[ic].IDprimary==0){
				TotalSelectionProb += 1.0/UnpartneredCount;
				if(r[ii] <= TotalSelectionProb && r[ii] > (TotalSelectionProb -
					1.0/UnpartneredCount)){
						partnerrisk = rprisk[ii];
						partnerID = rpID[ii];
						ChooseLTpartner(ic+1, partnerrisk, partnerID);
					}
			}
		}
		ii += 1;
	}
	// For married people who haven't been assigned a partner ID, change status to unmarried.
	for(ic=0; ic<tpp; ic++){
		if(Register[ic].MarriedInd==1 && Register[ic].IDprimary==0){
			Register[ic].MarriedInd = 0;}
	}

	// Assign number of short-term partners to each individual
	for(ig=0; ig<2; ig++){ // Modified from GetStartProfile function in TSHISA
		for(ia=0; ia<16; ia++){
			cHD = PartnershipFormation[0][ig] * AgeEffectPartners[ia * 5 + 2][ig] * MeanDurSTrel[0][0];
			cLD = PartnershipFormation[1][ig] * AgeEffectPartners[ia * 5 + 2][ig] * MeanDurSTrel[0][0];
			GUnmarriedH0[ia][ig] = 1.0/(1.0 + cHD + 0.5 * PartnerEffectNew[0][ig] * pow(cHD, 2.0));
			GUnmarriedH1[ia][ig] = cHD/(1.0 + cHD + 0.5 * PartnerEffectNew[0][ig] * pow(cHD, 2.0));
			GUnmarriedH2[ia][ig] = 1.0 - GUnmarriedH0[ia][ig] - GUnmarriedH1[ia][ig];
			GUnmarriedL0[ia][ig] = 1.0/(1.0 + cLD);
			GUnmarriedL1[ia][ig] = cLD/(1.0 + cLD);
			GMarriedH1[ia][ig] = 1.0/(1.0 + cHD * PartnerEffectNew[1][ig]);
			GMarriedH2[ia][ig] = 1.0 - GMarriedH1[ia][ig];
		}
		// Calculate InSTrelationship
		if (ig == 0){ High = HighPropnM; }
		else{ High = HighPropnF; }
		for (ia = 0; ia < 81; ia++){
			cHD = PartnershipFormation[0][ig] * AgeEffectPartners[ia][ig] * MeanDurSTrel[0][0];
			cLD = PartnershipFormation[1][ig] * AgeEffectPartners[ia][ig] * MeanDurSTrel[0][0];
			for (ir = 0; ir < 3; ir++){
				InSTrelationship[ia][ig][ir] = 1.0 - High * (IVirginPropnH[ia][ig][ir] + (1.0 - 
					IVirginPropnH[ia][ig][ir]) / (1.0 + cHD + 0.5 * PartnerEffectNew[0][ig] * 
					RaceEffectNew[ir] * pow(cHD, 2.0))) - (1 - High) * (IVirginPropnL[ia][ig][ir] +
					(1.0 - IVirginPropnL[ia][ig][ir]) / (1.0 + cLD));
			}
		}
	}

	for(ic=0; ic<tpp; ic++){
		r[ic] = rg.Random();}

	for(ic=0; ic<tpp; ic++){
		if(Register[ic].VirginInd==1){
			Register[ic].CurrPartners = 0;}
		else{
			ia = Register[ic].AgeGroup - 2;
			ig = Register[ic].SexInd;
			if(Register[ic].RiskGroup==1){
				if(Register[ic].MarriedInd==1){
					if(r[ic]<GMarriedH1[ia][ig]){
						Register[ic].CurrPartners = 1;}
					else{
						Register[ic].CurrPartners = 2;}
				}
				else{
					if(r[ic]<GUnmarriedH0[ia][ig]){
						Register[ic].CurrPartners = 0;}
					else if(r[ic]<1.0-GUnmarriedH2[ia][ig]){
						Register[ic].CurrPartners = 1;}
					else{
						Register[ic].CurrPartners = 2;}
				}
			}
			else{
				if(Register[ic].MarriedInd==1){
					Register[ic].CurrPartners = 1;}
				else{
					if(r[ic]<GUnmarriedL0[ia][ig]){
						Register[ic].CurrPartners = 0;}
					else{
						Register[ic].CurrPartners = 1;}
				}
			}
		}
	}

	// Next calculate sexual mixing in ST partnerships (adapted from code in TSHISA)
	DesiredSTpartners[0][0] = 0;
	DesiredSTpartners[0][1] = 0;
	DesiredSTpartners[0][2] = 0;
	DesiredSTpartners[1][0] = 0;
	DesiredSTpartners[1][1] = 0;
	DesiredSTpartners[1][2] = 0;
	for(ic=0; ic<tpp; ic++){
		if(Register[ic].CurrPartners>0){
			STpartners = Register[ic].CurrPartners - Register[ic].MarriedInd;
			ib = Register[ic].RiskGroup - 1;
			ig = Register[ic].SexInd;
			DesiredSTpartners[ib][ig] += STpartners;
		}
	}

	DesiredPartnerRiskM[0][0] = (1.0 - AssortativeM) + AssortativeM * 
		DesiredSTpartners[0][1] / (DesiredSTpartners[0][1] + DesiredSTpartners[1][1]);
	DesiredPartnerRiskM[0][1] = 1.0 - DesiredPartnerRiskM[0][0];
	DesiredPartnerRiskM[1][1] = (1.0 - AssortativeM) + AssortativeM * 
		DesiredSTpartners[1][1] / (DesiredSTpartners[0][1] + DesiredSTpartners[1][1]);
	DesiredPartnerRiskM[1][0] = 1.0 - DesiredPartnerRiskM[1][1];
	DesiredPartnerRiskF[0][0] = (1.0 - AssortativeF) + AssortativeF * 
		DesiredSTpartners[0][0] / (DesiredSTpartners[0][0] + DesiredSTpartners[1][0]);
	DesiredPartnerRiskF[0][1] = 1.0 - DesiredPartnerRiskF[0][0];
	DesiredPartnerRiskF[1][1] = (1.0 - AssortativeF) + AssortativeF * 
		DesiredSTpartners[1][0] / (DesiredSTpartners[0][0] + DesiredSTpartners[1][0]);
	DesiredPartnerRiskF[1][0] = 1.0 - DesiredPartnerRiskF[1][1];

	for(ii=0; ii<2; ii++){
		for(ij=0; ij<2; ij++){
			if(AllowBalancing==1){
				AdjSTrateM[ii][ij] = (GenderEquality * DesiredSTpartners[ij][1] * 
					DesiredPartnerRiskF[ij][ii] + (1.0 - GenderEquality) * 
					DesiredSTpartners[ii][0] * DesiredPartnerRiskM[ii][ij]);
				AdjSTrateF[ii][ij] = (GenderEquality * DesiredSTpartners[ii][1] * 
					DesiredPartnerRiskF[ii][ij] + (1.0 - GenderEquality) * 
					DesiredSTpartners[ij][0] * DesiredPartnerRiskM[ij][ii]);
			}
			else{
				AdjSTrateM[ii][ij] = DesiredSTpartners[ii][0] * DesiredPartnerRiskM[ii][ij];
				AdjSTrateF[ii][ij] = DesiredSTpartners[ii][1] * DesiredPartnerRiskF[ii][ij];
			}
		}
		ActualPropnSTH[ii][0] = AdjSTrateM[ii][0]/(AdjSTrateM[ii][0] + AdjSTrateM[ii][1]);
		ActualPropnSTH[ii][1] = AdjSTrateF[ii][0]/(AdjSTrateF[ii][0] + AdjSTrateF[ii][1]);
	}
	ActualPropnSTH_MSM[0] = (1.0 - AssortativeM) + AssortativeM * HighPropnM;
	ActualPropnSTH_MSM[1] = AssortativeM * HighPropnM;

	// Assign partner ID(s) to each individual in a ST relationship

	for(ic=0; ic<tpp; ic++){ // It's theoretically possible that # ST partnerships may
								 // be > InitPop, but very unlikely.
		r[ic] = rg.Random();}
	for(ic=0; ic<tpp; ic++){
		rprisk[ic] = rg.Random();}
	for(ic=0; ic<tpp; ic++){
		rpID[ic] = rg.Random();}

	MaxNewPartnerInd = 0;
	ii = 0;
	while(MaxNewPartnerInd==0 && ii<tpp){
		UnpartneredCount = 0;
		for(ic=0; ic<tpp; ic++){
			if(Register[ic].CurrPartners==1 && Register[ic].IDprimary==0){
				UnpartneredCount += 1;}
			if(Register[ic].CurrPartners==2 && Register[ic].IDprimary==0){
				UnpartneredCount += 2;}
			if(Register[ic].CurrPartners==2 && Register[ic].IDprimary>0 
				&& Register[ic].ID2ndary==0){
					UnpartneredCount += 1;}
		}
		TotalSelectionProb = 0.0;
		found = 0;
		for(ic=0; ic<tpp && !found; ic++){
			if(Register[ic].CurrPartners==1 && Register[ic].IDprimary==0){
				ij = 1;}
			else if(Register[ic].CurrPartners==2 && Register[ic].IDprimary==0){
				ij = 2;}
			else if(Register[ic].CurrPartners==2 && Register[ic].ID2ndary==0){
				ij = 1;}
			else{
				ij = 0;}
			if(ij > 0){
				temp1 = TotalSelectionProb;
				TotalSelectionProb += 1.0 * ij/UnpartneredCount;
				if(r[ii]<=TotalSelectionProb && r[ii]>temp1){
						partnerrisk = rprisk[ii];
						partnerID = rpID[ii];
						ChooseSTpartner(ic+1, partnerrisk, partnerID);
						found = 1;
					}
			}
		}
		ii += 1;
	}
	// For people in ST relationships who haven't been assigned a partner ID, reduce # partners.
	for(ic=0; ic<tpp; ic++){
		if(Register[ic].CurrPartners>0 && Register[ic].IDprimary==0){
			Register[ic].CurrPartners = 0;}
		else if(Register[ic].CurrPartners==2 && Register[ic].ID2ndary==0){
			Register[ic].CurrPartners = 1;}
	}

	// Determine rate at which men visit CSWs (copied from code in TSHISA)
	lambda = GammaMeanFSW/pow(GammaStdDevFSW, 2.0);
	alpha = GammaMeanFSW * lambda;
	for(ia=0; ia<81; ia++){
		AgeEffectFSWcontact[ia] = pow(lambda, alpha) * pow(ia + 0.5, alpha - 1.0) *
			exp(-lambda * (ia + 0.5));
		FSWcontactBase[ia] = 0.0;
	}
	MalePop15to49 = 0.0;
	for(ic=0; ic<tpp; ic++){
		if(Register[ic].SexInd==0 && Register[ic].RiskGroup==1 && Register[ic].VirginInd==0){
			ia = Register[ic].CurrAge - 10;
			FemPref = 1.0 - Register[ic].MalePref;
			if (Register[ic].CurrUrban == 1){ FemPref *= SWurban[0]; }
			else{ FemPref *= SWurban[1]; }
			if (Register[ic].Employed == 1){ FemPref *= RR_FSWcontactEmployedM; }
			if (Register[ic].Imprisoned == 1 || Register[ic].Imprisoned == 2){ FemPref = 0.0; }
			if (Register[ic].CurrPartners>0){
				ir = Register[ic].Race;
				FemPref *= RaceEffectNew[ir]; 
			}
			if(Register[ic].CurrPartners==0){
				FSWcontactBase[ia] += AgeEffectFSWcontact[ia] * PartnerEffectFSWcontact[0] * FemPref;}
			if(Register[ic].CurrPartners==1 && Register[ic].MarriedInd==0){
				FSWcontactBase[ia] += AgeEffectFSWcontact[ia] * PartnerEffectFSWcontact[1] * FemPref;}
			if(Register[ic].CurrPartners==1 && Register[ic].MarriedInd==1){
				FSWcontactBase[ia] += AgeEffectFSWcontact[ia] * PartnerEffectFSWcontact[2] * FemPref;}
			if(Register[ic].CurrPartners==2 && Register[ic].MarriedInd==0){
				FSWcontactBase[ia] += AgeEffectFSWcontact[ia] * PartnerEffectFSWcontact[3] * FemPref;}
			if(Register[ic].CurrPartners==2 && Register[ic].MarriedInd==1){
				FSWcontactBase[ia] += AgeEffectFSWcontact[ia] * PartnerEffectFSWcontact[4] * FemPref;}
		}
		if(Register[ic].SexInd==0 && Register[ic].AgeGroup>=3 && Register[ic].AgeGroup<10){
			MalePop15to49 += 1.0;}
	}
	BaseFSWdemand = 0.0;
	/*for (ia = 0; ia<16; ia++){
		BaseFSWdemand += FSWcontactBase[ia];}
	FSWcontactConstant = MeanFSWcontacts / (pow(lambda, alpha) * pow(11.0, alpha - 1.0) * exp(-lambda * 11.0));*/
	for(ia=5; ia<40; ia++){
		BaseFSWdemand += FSWcontactBase[ia];}
	FSWcontactConstant = MeanFSWcontacts * MalePop15to49/BaseFSWdemand;
	BaseFSWdemand += FSWcontactBase[0];
	for(ia=40; ia<81; ia++){
		BaseFSWdemand += FSWcontactBase[ia];}
	DesiredFSWcontacts[0] = BaseFSWdemand * FSWcontactConstant;
	RequiredNewFSW = DesiredFSWcontacts[0]/AnnNumberClients;

	// Assign women to commercial sex work
	CalcInitFSWageDbn(0);
	for(ia=0; ia<16; ia++){
		EligibleForCSW[ia] = 0;}
	for(ic=0; ic<tpp; ic++){
		if(Register[ic].SexInd==1 && Register[ic].RiskGroup==1 && 
			Register[ic].VirginInd==0 && Register[ic].CurrPartners==0){
				ia = Register[ic].AgeGroup - 2;
				EligibleForCSW[ia] += 1;
		}
	}

	for(ic=0; ic<tpp; ic++){
		r[ic] = rg.Random();}
	for (ir = 0; ir < 3; ir++){ TotCurrFSW[ir] = 0; }
	for(ic=0; ic<tpp; ic++){
		if(Register[ic].SexInd==1 && Register[ic].RiskGroup==1 && 
			Register[ic].VirginInd==0 && Register[ic].CurrPartners==0){
				ia = Register[ic].AgeGroup - 2;
				ir = Register[ic].Race;
				if(r[ic] < RequiredNewFSW * InitFSWageDbn[ia] / EligibleForCSW[ia]){
					Register[ic].FSWind = 1;
					CSWregister[TotCurrFSW[ir]][ir] = ic + 1;
					TotCurrFSW[ir] += 1;
				}
		}
	}

	// I haven't copied the code from TSHISA for calculating the FSWentry array values,
	// because at baseline these would be very sensitive to the relative numbers of unpartnered
	// women in the high risk group at each age, and it would make more sense to calculate
	// the FSWentry values at the start of each sexual behaviour cycle (i.e. so that the results
	// are less sensitive to the baseline numbers of unpartnered women by age).

	// Assign MSM to casual sex
	TotCurrCasual = 0;
	if (InclMSM == 1){
		for (ic = 0; ic < tpp; ic++){
			r[ic] = rg.Random();}
		for (ic = 0; ic < tpp; ic++){
			if (Register[ic].SexInd == 0 && Register[ic].MalePref>0.0 &&
				(Register[ic].RiskGroup==1 || Register[ic].CurrPartners==0) &&
				Register[ic].VirginInd==0){
					ExactAge = StartYear + 0.5 - Register[ic].DOB;
					CasualEntryMSM = CasualEntry * pow(CasualAgeAdj, ExactAge - 20.0) *
						Register[ic].MalePref;
					if (Register[ic].RiskGroup == 1 && Register[ic].CurrPartners > 0){
						CasualEntryMSM *= CasualHighAdj;}
					if (Register[ic].RiskGroup == 2){ CasualEntryMSM *= CasualLowAdj; }
					CasualExitMSM = CasualExit / Register[ic].MalePref;
					if (r[ic] < (1.0 / (1.0 + CasualExitMSM / CasualEntryMSM))){
						Register[ic].CasualInd = 1;
						TotCurrCasual += 1;
						CasualRegister[TotCurrCasual - 1] = ic + 1;
					}
			}
		}
	}

	// And assign heterosexual casual sex indicators
	for (ir = 0; ir < 3; ir++){
		for (ig = 0; ig < 2; ig++){ TotCurrCasualHet[ir][ig] = 0; }
	}
	// Use the same set of random numbers as before since people can't have both MSM
	// and heterosexual casual sex contacts at the same time. (But we apply no MalePref
	// adjustment to avoid double-penalizing the bisexual men with high random numbers.)
	for (ic = 0; ic < tpp; ic++){
		if ((Register[ic].SexInd == 1 || Register[ic].MalePref < 1.0) &&
			(Register[ic].RiskGroup == 1 || Register[ic].CurrPartners == 0) &&
			Register[ic].VirginInd == 0){
				ig = Register[ic].SexInd;
				ir = Register[ic].Race;
				ExactAge = StartYear + 0.5 - Register[ic].DOB;
				if (ExactAge > 90){ ExactAge = 90; }
				//CasualEntryMSM = CasualEntryHet[ig] * pow(CasualAgeAdjHet[ig], ExactAge - 20.0);
				CasualEntryMSM = CasualEntryHet[ig] * AgeEffectPartners[ExactAge - 10][ig] /
					BasePartnerAcqH[ig];
				if (Register[ic].RiskGroup == 1 && Register[ic].CurrPartners > 0){
					CasualEntryMSM *= CasualHighAdjHet;
				}
				if (Register[ic].RiskGroup == 2){ CasualEntryMSM *= CasualLowAdjHet[ig]; }
				if (r[ic] < (1.0 / (1.0 + CasualExit / CasualEntryMSM))){
					Register[ic].HetCasualInd = 1;
					TotCurrCasualHet[ir][ig] += 1;
					TempTot = TotCurrCasualHet[ir][ig];
					if (ig == 0){ CasualRegisterM[TempTot - 1][ir] = ic + 1; }
					else{ CasualRegisterF[TempTot - 1][ir] = ic + 1; }
				}
		}
	}

	// MSM outputs
	for (ic = 0; ic < tpp; ic++){
		if (Register[ic].SexInd == 0 && Register[ic].MalePref>0.0 &&
			Register[ic].VirginInd == 0){
			if (Register[ic].CurrPartners > 0){ 
				PID = Register[ic].IDprimary;
				if (Register[PID - 1].SexInd == 0){ Register[ic].EverMSM = 1; }
				else{ Register[ic].EverBi = 1; }
			}
			if (Register[ic].CurrPartners == 2){
				PID = Register[ic].ID2ndary;
				if (Register[PID - 1].SexInd == 0){ Register[ic].EverMSM = 1; }
				else{ Register[ic].EverBi = 1; }
			}
			if (Register[ic].CasualInd == 1){ Register[ic].EverMSM = 1; }
		}
	}

	// Assign PartnerRateAdj
	
	if(AllowPartnerRateAdj==1){
		for(ic=0; ic<tpp; ic++){
			r[ic] = rg.Random();}
		ind = 2;
		// Note that the following formulas for a and b apply only when the gamma mean is 1.
		a = 1.0/pow(SDpartnerRateAdj, 2.0);
		b = a;
		for(ic=0; ic<tpp; ic++){
			p = r[ic];
			q = 1 - r[ic];
			cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
			Register[ic].PartnerRateAdj = x;
		}
	}
	else{
		for(ic=0; ic<tpp; ic++){
			Register[ic].PartnerRateAdj = 1.0;}
	}

	// Assign innate condom preference
	
	for (ic = 0; ic<tpp; ic++){
		r[ic] = rg.Random();}
	ind = 2;
	// Note that the following formulas for a and b apply only when the gamma mean is 1.
	a = 1.0 / pow(SDcondomPref, 2.0);
	b = a;
	for (ic = 0; ic<tpp; ic++){
		p = r[ic];
		q = 1 - r[ic];
		cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
		Register[ic].CondomPref = x;
		ir = Register[ic].Race;
		Register[ic].CondomPref *= CondomRace[ir];
	}
}

void Pop::AssignVisitFreq()
{
	int ic, ind, ir;
	double a, b, p, q, x, TempProb;

	for (ic = 0; ic<InitPop; ic++){
		r[ic] = rg.Random();}
	for (ic = 0; ic<InitPop; ic++){
		r2[ic] = rg.Random();}

	ind = 2;
	a = VisitFreq[0];
	b = VisitFreq[1];
	for (ic = 0; ic < InitPop; ic++){
		p = r[ic];
		q = 1 - r[ic];
		cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
		Register[ic].VisitFreq = x;
		ir = Register[ic].Race;
		TempProb = 8.0 * LocationMixing * MarriedMigAdj[ir]/SeparatedMigAdj[ir];
		// The constant adjustment of 8.0 was chosen to give a stable trend in
		// the fraction of married individuals not cohabiting (CohabitTempTrend)
		if (r2[ic] < TempProb){
			Register[ic].Visiting = 1; }
		else{ Register[ic].Visiting = 0; }
	}
}

void Pop::AssignEdu()
{
	int ii, ir, ia, ic;

	for (ic = 0; ic<InitPop; ic++){
		r[ic] = rg.Random();}
	for (ic = 0; ic<InitPop; ic++){
		r2[ic] = rg.Random();}

	// Simulate highest educational attainment for each individual
	for (ic = 0; ic<InitPop; ic++){
		ir = Register[ic].Race;
		ia = Register[ic].CurrAge;
		Register[ic].HighestGrade = 13;
		Register[ic].InSchool = 0;
		for (ii = 0; ii < 13; ii++){
			if (r[ic] < CumGrade1985[ia][ii][ir]){
				Register[ic].HighestGrade = ii;
				break;
			}
		}
		ii = Register[ic].HighestGrade;
		//Register[ic].CondomPref = CondomEdu[ii];
		if (ia <= 30 && Register[ic].HighestGrade < 13){
			ii = Register[ic].HighestGrade;
			if (r2[ic] < InSchool1985[ia][ii][ir]){
				Register[ic].InSchool = 1;}
		}
	}
}

void Pop::AssignBaseIncome()
{
	int ih, ir, ic, ind, s;
	double temp, x, a, b, p, q, bound;

	for (ic = 0; ic<InitPop; ic++){
		r[ic] = rg.Random();}
	for (ic = 0; ic<InitPop; ic++){
		r2[ic] = rg.Random();}
	
	// Assign log income difference
	a = 0.0;
	b = 1.0;
	ind = 2;
	s = 0;
	bound = 0.0;
	for (ic = 0; ic < InitPop; ic++){
		p = r[ic];
		q = 1.0 - r[ic];
		cdfnor(&ind, &p, &q, &x, &a, &b, &s, &bound);
		Register[ic].LogIncomeDif = x * StdDevLogIncome +
			Register[ic].Conscientiousness * ConscientiousEffectLogIncome;
	}

	// Assign private pension
	for (ic = 0; ic<InitPop; ic++){
		Register[ic].PrivatePension = 0;
		if (Register[ic].CurrAge >= 60){
			ir = Register[ic].Race;
			if (Register[ic].HighestGrade == 0){ ih = 0; }
			else if (Register[ic].HighestGrade <= 7){ ih = 1; }
			else if (Register[ic].HighestGrade <= 11){ ih = 2; }
			else if (Register[ic].HighestGrade == 12){ ih = 3; }
			else { ih = 4; }
			temp = BaseOddsPension * EffectRacePension[ir] * EffectEduPension[ih];
			temp = temp / (1.0 + temp);
			if (r2[ic] < temp || ir == 2){ Register[ic].PrivatePension = 1; }
		}
	}
}

void Pop::AssignUrban()
{
	int ii, ir, ia, ic, ig;

	for (ic = 0; ic<InitPop; ic++){
		r[ic] = rg.Random();}

	// Simulate urban-rural location for each individual
	for (ic = 0; ic<InitPop; ic++){
		ir = Register[ic].Race;
		if (ir == 0){
			if (Register[ic].HighestGrade == 0){ ii = 0; }
			else if (Register[ic].HighestGrade <=5){ ii = 1; }
			else if (Register[ic].HighestGrade <= 7){ ii = 2; }
			else if (Register[ic].HighestGrade <= 11){ ii = 3; }
			else if (Register[ic].HighestGrade == 12){ ii = 4; }
			else { ii = 5; }
		}
		else if (ir == 1){ ii = 6; }
		else{ ii = 7; }
		ia = Register[ic].AgeGroup;
		ig = Register[ic].SexInd;
		if (r[ic] < InitUrbanPropn[ia][ii][ig]){
			Register[ic].CurrUrban = 1;}
		else{ Register[ic].CurrUrban = 0; }
	}
}

void Pop::AssignParents()
{
	// This function only gets called at the start of the simulation (1985).

	int ic, ir, ia, ig, iu, MID, PID;
	double ParentSurvival[91][3][2]; // Prob parent is still alive, by child's
									 // current age & race, and by parent's sex

	// CalculateParentSurvival, assuming mother's age at birth was 27 and
	// father's age at birth was 30.
	for (ir = 0; ir < 3; ir++){
		for (ig = 0; ig < 2; ig++){ ParentSurvival[0][ir][ig] = 1.0; }
		for (ia = 1; ia < 61; ia++){
			ParentSurvival[ia][ir][0] = ParentSurvival[ia - 1][ir][0] *
				(1.0 - NonAIDSmortM[ia + 29][0][ir]);
		}
		for (ia = 61; ia < 91; ia++){ 
			ParentSurvival[ia][ir][0] = ParentSurvival[ia - 1][ir][0] *
				(1.0 - NonAIDSmortM[90][0][ir]);
		}
		for (ia = 1; ia < 64; ia++){
			ParentSurvival[ia][ir][1] = ParentSurvival[ia - 1][ir][1] *
				(1.0 - NonAIDSmortF[ia + 26][0][ir]);
		}
		for (ia = 64; ia < 91; ia++){ 
			ParentSurvival[ia][ir][1] = ParentSurvival[ia - 1][ir][1] *
				(1.0 - NonAIDSmortF[90][0][ir]);
		}
	}

	// Generate random numbers for assigning initial parent survival
	for (ic = 0; ic<InitPop; ic++){
		r[ic] = rg.Random();}
	for (ic = 0; ic<InitPop; ic++){
		r2[ic] = rg.Random();}
	for (ic = 0; ic<InitPop; ic++){
		rpID2[ic] = rg.Random();}

	// Assign parents/orphanhood status
	for (ic = 0; ic < Register.size(); ic++){
		ia = Register[ic].CurrAge;
		ir = Register[ic].Race;
		iu = Register[ic].CurrUrban;
		if (r[ic] < ParentSurvival[ia][ir][1]){
			Register[ic].FindMother(ic + 1);
		}
		else{ Register[ic].ParentID[1] = 0; }
		// For individuals aged <15 we require at least 1 parent alive initially.
		if (r2[ic] >= ParentSurvival[ia][ir][0] && (ia >= 15 || Register[ic].ParentID[1] > 0)){
			Register[ic].ParentID[0] = 0; 
		}
		else{
			if (Register[ic].ParentID[1] > 0){
				MID = Register[ic].ParentID[1];
				if (Register[MID - 1].MarriedInd == 1){ 
					PID = Register[MID - 1].IDprimary;
					if (Register[PID - 1].CurrAge - ia > 14 && Register[PID - 1].Race == ir &&
						Register[PID - 1].ChildIDs[19] == 0 && Register[PID - 1].CurrUrban == iu){
						Register[ic].ParentID[0] = PID;
						Register[PID - 1].UpdateChildren(ic + 1);
					}
					else{ Register[ic].FindFather(ic + 1); }
				}
				else{ Register[ic].FindFather(ic + 1); }
			}
			else{ Register[ic].FindFather(ic + 1); }
		}
	}
}

void Pop::AssignBirth()
{
	int ic, ia, ir, ind;
	double MostRecentBirth, Temp1, TempProb, a, b, p, q, x;

	for (ic = 0; ic<InitPop; ic++){
		revent[ic] = rg.Random();}
	for (ic = 0; ic<InitPop; ic++){
		r[ic] = rg.Random();}
	for (ic = 0; ic<InitPop; ic++){
		r2[ic] = rg.Random();}

	// Determine the initial fecundability parameters
	for (ic = 0; ic<InitPop; ic++){
		if (Register[ic].SexInd == 1 && Register[ic].CurrAge>=12 &&
			Register[ic].CurrAge < 50 && Register[ic].VirginInd == 0){
			ia = Register[ic].CurrAge;
			if (ia >= 12 && ia < 20){
				TempProb = 1.0 * (ia - 12) / 8.0;
				if (revent[ic] < TempProb){ Register[ic].Fecundability = 1.0; }
			}
			else if (ia <= 35){ 
				TempProb = 1.0;
				Register[ic].Fecundability = 1.0; 
			}
			else{
				TempProb = 1.0 - 1.0 * (ia - 35) / 15.0;
				if (revent[ic] < TempProb){ Register[ic].Fecundability = 1.0; }
			}
			if (Register[ic].Fecundability == 1){
				// Randomly assign initial fecundability
				ind = 2;
				a = 1.0 / pow(SDfecundability, 2.0);
				b = a;
				p = revent[ic] / TempProb;
				q = 1.0 - revent[ic] / TempProb;
				cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
				Register[ic].Fecundability = x;
			}
		}
	}

	// Determine which women are pregnant at baseline
	// Ideally we should adjust fertility rates to take into account virgins, but
	// this doesn't make a big difference.
	for (ic = 0; ic<InitPop; ic++){
		if (Register[ic].SexInd == 1 && Register[ic].AgeGroup>2 &&
			Register[ic].AgeGroup < 10 && Register[ic].VirginInd == 0 &&
			Register[ic].Fecundability > 0.0){
			ia = Register[ic].CurrAge - 15;
			ir = Register[ic].Race;
			if (r[ic] < FertilityTable[ia][0][ir] * MeanGestation){
				Register[ic].DOLB = 1985.5 + r[ic] / FertilityTable[ia][0][ir];
				if (Register[ic].CurrPartners > 0){ Register[ic].FatherBaby = Register[ic].IDprimary; }
				else{ Register[ic].FatherBaby = 0; }
			}
		}
	}
	
	// Determine which women have ever been pregnant
	// This is an improvement on version 9c, where we weren't linking mothers and children.
	for (ic = 0; ic<InitPop; ic++){
		if (Register[ic].SexInd == 1 && Register[ic].AgeGroup>2 &&
			Register[ic].VirginInd == 0 && Register[ic].DOLB==0 &&
			(Register[ic].Fecundability > 0.0 || Register[ic].AgeGroup > 3)){
			MostRecentBirth = 0.0;
			for (ia = 0; ia < 20; ia++){
				ir = Register[ic].ChildIDs[ia];
				if (ir > 0){ 
					Temp1 = Register[ir - 1].DOB; 
					if (Temp1 > MostRecentBirth){ MostRecentBirth = Temp1; }
				}
				else{ break; }
			}
			Register[ic].DOLB = MostRecentBirth;
		}
	}
}

void Pop::AssignIneqGender()
{
	int ic, ia, ir, ih;
	double TempOdds, AveRisk;

	AveRisk = HighPropnM + (1.0 - HighPropnM) * 2.0;
	for (ic = 0; ic<InitPop; ic++){
		r[ic] = rg.Random();}

	for (ic = 0; ic<InitPop; ic++){
		if (Register[ic].SexInd == 0){ // Only relevant to males
			TempOdds = exp(-r[ic]);
			ir = Register[ic].Race;
			TempOdds *= RaceEffectIneqGender[ir];
			TempOdds *= pow(HighRiskEffectIneqGender, AveRisk - Register[ic].RiskGroup);
			ia = Register[ic].CurrAge;
			if (ia >= 20){ 
				TempOdds *= pow(AgeEffectIneqGender, ia - 20); 
				if (Register[ic].CurrUrban == 0){ 
					TempOdds *= RuralEffectIneqGender; }
				ih = 0;
				if (Register[ic].HighestGrade >= 12){ ih = 2; }
				else if (Register[ic].HighestGrade >= 8){ ih = 1; }
				TempOdds *= EduEffectIneqGender[ih];
			}
			Register[ic].IneqGender = TempOdds / (1.0 + TempOdds);
		}
	}
}

void Pop::AssignAlcohol()
{
	int ic, ig, ia, ir, ih, ind, s;
	double Temp, a, b, p, q, x, bound;

	for (ic = 0; ic<InitPop; ic++){
		r[ic] = rg.Random();
		rprisk[ic] = rg.Random();
	}

	ind = 2;
	bound = 0.0;
	s = 0;
	a = 0.0;
	for (ic = 0; ic<InitPop; ic++){
		ig = Register[ic].SexInd;
		b = StdDevDrinkProb[ig] / (7.0 * AlcPropnReported);
		if (Register[ic].DrinkProbConstant == 0.0){
			p = r[ic];
			q = 1.0 - r[ic];
		}
		else{
			p = Register[ic].DrinkProbConstant;
			q = 1.0 - Register[ic].DrinkProbConstant;
		}
		cdfnor(&ind, &p, &q, &x, &a, &b, &s, &bound);
		Register[ic].DrinkProbConstant = x;
		b = StdDevDrinksPerDD[ig] / AlcPropnReported;
		// Check if we previously assigned confounding between alcohol and risk group
		if (Register[ic].DrinksPerDDconstant == 0.0){
			p = rprisk[ic];
			q = 1.0 - rprisk[ic];
		}
		else{
			p = Register[ic].DrinksPerDDconstant;
			q = 1.0 - Register[ic].DrinksPerDDconstant;
		}
		cdfnor(&ind, &p, &q, &x, &a, &b, &s, &bound);
		Register[ic].DrinksPerDDconstant = x;
		if (Register[ic].CurrAge>=15){ 
			Temp = BaseDrinkProb[ig];
			if (Register[ic].CurrUrban == 1){ Temp += UrbanEffectDrinkProb[ig]; }
			ir = Register[ic].Race;
			Temp += RaceEffectDrinkProb[ir][ig];
			ia = (Register[ic].AgeGroup - 3) / 2;
			if (ia >= 4){ ia = 4; }
			Temp += AgeEffectDrinkProb[ia][ig];
			ih = 0;
			if (Register[ic].HighestGrade >= 13){ ih = 3; }
			else if (Register[ic].HighestGrade == 12){ ih = 2; }
			else if (Register[ic].HighestGrade >= 8){ ih = 1; }
			Temp += EduEffectDrinkProb[ih][ig];
			if (Register[ic].MarriedInd == 1){ Temp += MarriedEffectDrinkProb[ig]; }
			Register[ic].DailyDrinkProb = exp(Temp) / (7.0 * AlcPropnReported) + 
				Register[ic].DrinkProbConstant;
			if (Register[ic].DailyDrinkProb < 0.0){ Register[ic].DailyDrinkProb = 0.0; }
			if (Register[ic].DailyDrinkProb > 1.0){ Register[ic].DailyDrinkProb = 1.0; }
			Temp = BaseDrinksPerDD[ig];
			if (Register[ic].CurrUrban == 1){ Temp += UrbanEffectDrinksPerDD[ig]; }
			Temp += RaceEffectDrinksPerDD[ir][ig];
			Temp += AgeEffectDrinksPerDD[ia][ig];
			Temp += EduEffectDrinksPerDD[ih][ig];
			if (Register[ic].Employed == 1){ Temp += EmployedEffectDrinksPerDD[ig]; }
			Temp = exp(Temp);
			Temp += ConscientiousEffectDrinksPerDD * Register[ic].Conscientiousness;
			if (ig == 0){ Temp += GenderIneqEffectDrinksPerDD * (Register[ic].IneqGender - 0.2); }
			Register[ic].DrinksPerDD = Temp / AlcPropnReported + Register[ic].DrinksPerDDconstant;
			if (Register[ic].DrinksPerDD < 1.0){ Register[ic].DrinksPerDD = 1.0; }
		}
		else{
			Register[ic].DailyDrinkProb = 0.0;
			Register[ic].DrinksPerDD = 1.0;
		}
	}
}

void Pop::AssignInitHousehold()
{
	AssignHeadship();
	AssignAdultsHH();
	AssignChildrenHH();
}

void Pop::AssignHeadship()
{
	// Determine whether adults are household heads, at the start of the simulation

	int ic, ig, ir, ih, ia, PID, HHID;
	double TempOdds, ProbHead;

	for (ic = 0; ic<InitPop; ic++){
		r[ic] = rg.Random();}

	for (ic = 0; ic < InitPop; ic++){
		if (Register[ic].CurrAge >= 15 && Register[ic].HouseholdID == 0){
			ig = Register[ic].SexInd;
			ia = Register[ic].AgeGroup - 2;
			if (ia > 10){ ia = 10; }
			ir = Register[ic].Race;
			ih = 0;
			if (Register[ic].HighestGrade > 0 && Register[ic].HighestGrade <= 7){ ih = 1; }
			if (Register[ic].HighestGrade > 7 && Register[ic].HighestGrade < 12){ ih = 2; }
			if (Register[ic].HighestGrade == 12){ ih = 3; }
			if (Register[ic].HighestGrade == 13){ ih = 4; }
			TempOdds = BaseOddsHead[ig];
			if (Register[ic].MarriedInd == 1){ TempOdds *= MarriedEffectHead[ig]; }
			TempOdds *= AgeEffectHead[ia][ig];
			TempOdds *= RaceEffectHead[ir][ig];
			if (Register[ic].Employed == 1){ TempOdds *= EmployedEffectHead[ig]; }
			if (Register[ic].InSchool == 1){ TempOdds *= InSchoolEffectHead[ig]; }
			TempOdds *= EduEffectHead[ih][ig];
			if (ig == 1 && Register[ic].DOLB > 0.0){ TempOdds *= BirthEffectHead; }
			ProbHead = TempOdds / (1.0 + TempOdds);
			if (r[ic] < ProbHead){
				NewHousehold(ic + 1);
				if (Register[ic].MarriedInd == 1){
					// The spouse can't be a HH head, otherwise the condition HouseholdID=0
					// would not have held when we tested for it earlier.
					PID = Register[ic].IDprimary;
					if (Register[PID - 1].CurrUrban == Register[ic].CurrUrban){
						HHID = Register[ic].HouseholdID;
						HHregister[HHID - 1].AddMember(PID);
						Register[PID - 1].HouseholdID = HHID;
					}
				}
			}
		}
	}
}

void Pop::AssignAdultsHH()
{
	// Determine whether adults who are NOT household heads are living with parents, 
	// at the start of the simulation, and assign them to households based on this.

	int ic, ii, ig, ir, ih, ia, PID, HHID, Siblings, SID, RandomStart, Children, CID;
	double TempOdds, ProbALWP;

	for (ic = 0; ic<InitPop; ic++){
		r[ic] = rg.Random();}
	for (ic = 0; ic<InitPop; ic++){
		r2[ic] = rg.Random();}

	// Check if person is living with parent
	for (ic = 0; ic < InitPop; ic++){
		if (Register[ic].CurrAge >= 15 && Register[ic].HouseholdID == 0){
			ig = Register[ic].SexInd;
			ia = Register[ic].AgeGroup - 2;
			if (ia > 7){ ia = 7; }
			ir = Register[ic].Race;
			ih = 0;
			if (Register[ic].HighestGrade > 0 && Register[ic].HighestGrade <= 7){ ih = 1; }
			if (Register[ic].HighestGrade > 7 && Register[ic].HighestGrade < 12){ ih = 2; }
			if (Register[ic].HighestGrade == 12){ ih = 3; }
			if (Register[ic].HighestGrade == 13){ ih = 4; }
			TempOdds = BaseOddsALWP[ig];
			if (Register[ic].MarriedInd == 1){ TempOdds *= MarriedEffectALWP[ig]; }
			TempOdds *= AgeEffectALWP[ia][ig];
			TempOdds *= RaceEffectALWP[ir][ig];
			if (Register[ic].Employed == 1){ TempOdds *= EmployedEffectALWP[ig]; }
			if (Register[ic].InSchool == 1){ TempOdds *= InSchoolEffectALWP[ig]; }
			TempOdds *= EduEffectALWP[ih][ig];
			ProbALWP = TempOdds / (1.0 + TempOdds);
			if (r[ic] < ProbALWP){
				if (Register[ic].ParentID[1] > 0){ // Mother is alive
					PID = Register[ic].ParentID[1];
					if (Register[PID - 1].HouseholdID > 0){
						HHID = Register[PID - 1].HouseholdID;
						if (HHregister[HHID - 1].IDhead == PID){ // Mother is a HH head
							if (HHregister[HHID - 1].Size < MaxHHsize){
								HHregister[HHID - 1].AddMember(ic + 1); 
								Register[ic].HouseholdID = HHID;
							}
						}
					}
				}
				if (Register[ic].HouseholdID == 0 && Register[ic].ParentID[0] > 0){ // Father is alive
					PID = Register[ic].ParentID[0];
					if (Register[PID - 1].HouseholdID > 0){
						HHID = Register[PID - 1].HouseholdID;
						if (HHregister[HHID - 1].IDhead == PID){ // Father is a HH head
							if (HHregister[HHID - 1].Size < MaxHHsize){
								HHregister[HHID - 1].AddMember(ic + 1);
								Register[ic].HouseholdID = HHID;
							}
						}
					}
				}
				if (Register[ic].MarriedInd == 1 && Register[ic].HouseholdID > 0){
					// The spouse can't be a HH head, otherwise the condition HouseholdID=0
					// would not have held when we tested for it earlier.
					PID = Register[ic].IDprimary;
					if (Register[PID - 1].CurrUrban == Register[ic].CurrUrban){
						HHID = Register[ic].HouseholdID;
						if (HHregister[HHID - 1].Size < MaxHHsize){
							HHregister[HHID - 1].AddMember(PID);
							Register[PID - 1].HouseholdID = HHID;
						}
					}
				}
			}
		}
	}
	
	// For those adults not yet assigned to a household, check if they have a maternal 
	// sibling they can live with
	for (ic = 0; ic < InitPop; ic++){
		if (Register[ic].CurrAge >= 15 && Register[ic].HouseholdID == 0){
			PID = Register[ic].ParentID[1];
			if (PID > 0){
				Siblings = 0;
				while (Register[PID - 1].ChildIDs[Siblings] > 0){ Siblings += 1; }
				if (Siblings > 0){
					RandomStart = r2[ic] * Siblings;
					for (ii = 0; ii < Siblings; ii++){
						SID = RandomStart + ii;
						if (SID >= Siblings){ SID = SID - Siblings; }
						SID = Register[PID - 1].ChildIDs[SID];
						if (Register[SID - 1].HouseholdID > 0){
							HHID = Register[SID - 1].HouseholdID;
							if (HHregister[HHID - 1].IDhead == SID && HHregister[HHID - 1].Size < MaxHHsize){
								HHregister[HHID - 1].AddMember(ic + 1);
								Register[ic].HouseholdID = HHID;
								break;
							}
						}
					}
				}
				if (Register[ic].MarriedInd == 1 && Register[ic].HouseholdID > 0){
					PID = Register[ic].IDprimary;
					if (Register[PID - 1].CurrUrban == Register[ic].CurrUrban){
						HHID = Register[ic].HouseholdID;
						if (HHregister[HHID - 1].Size < MaxHHsize){
							HHregister[HHID - 1].AddMember(PID);
							Register[PID - 1].HouseholdID = HHID;
						}
					}
				}
			}
		}
	}
	
	// Failing that, check if they have a paternal sibling they can live with
	for (ic = 0; ic < InitPop; ic++){
		if (Register[ic].CurrAge >= 15 && Register[ic].HouseholdID == 0){
			PID = Register[ic].ParentID[0];
			if (PID > 0){
				Siblings = 0;
				while (Register[PID - 1].ChildIDs[Siblings] > 0){ Siblings += 1; }
				if (Siblings > 0){
					RandomStart = r2[ic] * Siblings;
					for (ii = 0; ii < Siblings; ii++){
						SID = RandomStart + ii;
						if (SID >= Siblings){ SID = SID - Siblings; }
						SID = Register[PID - 1].ChildIDs[SID];
						if (Register[SID - 1].HouseholdID > 0){
							HHID = Register[SID - 1].HouseholdID;
							if (HHregister[HHID - 1].IDhead == SID && HHregister[HHID - 1].Size < MaxHHsize){
								HHregister[HHID - 1].AddMember(ic + 1);
								Register[ic].HouseholdID = HHID;
								break;
							}
						}
					}
				}
				if (Register[ic].MarriedInd == 1 && Register[ic].HouseholdID > 0){
					PID = Register[ic].IDprimary;
					if (Register[PID - 1].CurrUrban == Register[ic].CurrUrban){
						HHID = Register[ic].HouseholdID;
						if (HHregister[HHID - 1].Size < MaxHHsize){
							HHregister[HHID - 1].AddMember(PID);
							Register[PID - 1].HouseholdID = HHID;
						}
					}
				}
			}
		}
	}
	
	// If still not assigned to a household, check if they have a child to live with
	for (ic = 0; ic < InitPop; ic++){
		if (Register[ic].CurrAge >= 40 && Register[ic].HouseholdID == 0){
			// Assuming child would have to be aged at least 25 for them to be classified head
			Children = 0;
			while (Register[ic].ChildIDs[Children] > 0){ Children += 1; }
			if (Children > 0){
				RandomStart = r2[ic] * Children;
				for (ii = 0; ii < Children; ii++){
					CID = RandomStart + ii;
					if (CID >= Children){ CID = CID - Children; }
					CID = Register[ic].ChildIDs[CID];
					if (Register[CID - 1].HouseholdID > 0 && Register[CID - 1].CurrAge > 25){
						HHID = Register[CID - 1].HouseholdID;
						if (HHregister[HHID - 1].IDhead == CID && HHregister[HHID - 1].Size < MaxHHsize){
							HHregister[HHID - 1].AddMember(ic + 1);
							Register[ic].HouseholdID = HHID;
							break;
						}
					}
				}
			}
			if (Register[ic].MarriedInd == 1 && Register[ic].HouseholdID > 0){
				PID = Register[ic].IDprimary;
				if (Register[PID - 1].CurrUrban == Register[ic].CurrUrban){
					HHID = Register[ic].HouseholdID;
					if (HHregister[HHID - 1].Size < MaxHHsize){
						HHregister[HHID - 1].AddMember(PID);
						Register[PID - 1].HouseholdID = HHID;
					}
				}
			}
		}
	}
	
	// Failing that, assign to the household of a surviving parent, even if that parent
	// is not a household head
	for (ic = 0; ic < InitPop; ic++){
		if (Register[ic].CurrAge >= 15 && Register[ic].HouseholdID == 0){
			PID = Register[ic].ParentID[1];
			if (PID > 0){
				if (Register[PID - 1].HouseholdID > 0){
					HHID = Register[PID - 1].HouseholdID;
					if (HHregister[HHID - 1].Size < MaxHHsize){
						HHregister[HHID - 1].AddMember(ic + 1);
						Register[ic].HouseholdID = HHID;
					}
				}
			}
			PID = Register[ic].ParentID[0];
			if (PID > 0 && Register[ic].HouseholdID == 0){
				if (Register[PID - 1].HouseholdID > 0){
					HHID = Register[PID - 1].HouseholdID;
					if (HHregister[HHID - 1].Size < MaxHHsize){
						HHregister[HHID - 1].AddMember(ic + 1);
						Register[ic].HouseholdID = HHID;
					}
				}
			}
			if (Register[ic].MarriedInd == 1 && Register[ic].HouseholdID > 0){
				PID = Register[ic].IDprimary;
				if (Register[PID - 1].CurrUrban == Register[ic].CurrUrban){
					HHID = Register[ic].HouseholdID;
					if (HHregister[HHID - 1].Size < MaxHHsize){
						HHregister[HHID - 1].AddMember(PID);
						Register[PID - 1].HouseholdID = HHID;
					}
				}
			}
		}
	}
	
	// If all else fails, reassign the individual to be a household head.
	// Not ideal - there are a lot of unmarried young adults who are double orphans
	// at the start of the simulation, who will end up in this group?
	for (ic = 0; ic < InitPop; ic++){
		if (Register[ic].CurrAge >= 15 && Register[ic].HouseholdID == 0){
			NewHousehold(ic + 1);
			if (Register[ic].MarriedInd == 1 && Register[ic].HouseholdID > 0){
				PID = Register[ic].IDprimary;
				if (Register[PID - 1].CurrUrban == Register[ic].CurrUrban){
					HHID = Register[ic].HouseholdID;
					HHregister[HHID - 1].AddMember(PID);
					Register[PID - 1].HouseholdID = HHID;
				}
			}
		}
	}
}

void Pop::AssignChildrenHH()
{
	int ic, ii, ir, ia, PID, HHID, PID2;
	double TempOdds, ProbCLWP;

	for (ic = 0; ic<InitPop; ic++){
		r[ic] = rg.Random();}
	for (ic = 0; ic<InitPop; ic++){
		r2[ic] = rg.Random();}

	// Check if person is living with parent who is household head
	for (ic = 0; ic < InitPop; ic++){
		if (Register[ic].CurrAge < 15 && Register[ic].HouseholdID == 0){
			ia = Register[ic].AgeGroup;
			ir = Register[ic].Race;
			TempOdds = BaseOddsCLWP;
			TempOdds *= AgeEffectCLWP[ia];
			TempOdds *= RaceEffectCLWP[ir];
			if (Register[ic].InSchool == 1){ TempOdds *= InSchoolEffectCLWP; }
			ProbCLWP = TempOdds / (1.0 + TempOdds);
			if (r[ic] < ProbCLWP){
				if (Register[ic].ParentID[1] > 0){ // Mother is alive
					PID = Register[ic].ParentID[1];
					if (Register[PID - 1].HouseholdID > 0){
						HHID = Register[PID - 1].HouseholdID;
						if (HHregister[HHID - 1].IDhead == PID){ // Mother is a HH head
							HHregister[HHID - 1].AddMember(ic + 1);
							Register[ic].HouseholdID = HHID;
						}
					}
				}
				if (Register[ic].HouseholdID == 0 && Register[ic].ParentID[0] > 0){ // Father is alive
					PID = Register[ic].ParentID[0];
					if (Register[PID - 1].HouseholdID > 0){
						HHID = Register[PID - 1].HouseholdID;
						if (HHregister[HHID - 1].IDhead == PID){ // Father is a HH head
							HHregister[HHID - 1].AddMember(ic + 1);
							Register[ic].HouseholdID = HHID;
						}
					}
				}
			}
			if (Register[ic].HouseholdID == 0){
				// Check if they have a maternal grandparent they can live with
				if (Register[ic].ParentID[1] > 0){ // Mother is alive
					PID = Register[ic].ParentID[1];
					if (Register[PID - 1].ParentID[1] > 0){
						PID2 = Register[PID - 1].ParentID[1];
						HHID = Register[PID2 - 1].HouseholdID;
						if (HHregister[HHID - 1].IDhead == PID2){ // Grandma is a HH head
							HHregister[HHID - 1].AddMember(ic + 1);
							Register[ic].HouseholdID = HHID;
						}
					}
					if (Register[ic].HouseholdID == 0 && Register[PID - 1].ParentID[0] > 0){
						PID2 = Register[PID - 1].ParentID[0];
						HHID = Register[PID2 - 1].HouseholdID;
						if (HHregister[HHID - 1].IDhead == PID2){ // Grandpa is a HH head
							HHregister[HHID - 1].AddMember(ic + 1);
							Register[ic].HouseholdID = HHID;
						}
					}
				}
				// Failing that, check if they have a paternal grandparent they can live with
				if (Register[ic].HouseholdID == 0 && Register[ic].ParentID[0] > 0){ // Father is alive
					PID = Register[ic].ParentID[0];
					if (Register[PID - 1].ParentID[1] > 0){
						PID2 = Register[PID - 1].ParentID[1];
						HHID = Register[PID2 - 1].HouseholdID;
						if (HHregister[HHID - 1].IDhead == PID2){ // Grandma is a HH head
							HHregister[HHID - 1].AddMember(ic + 1);
							Register[ic].HouseholdID = HHID;
						}
					}
					if (Register[ic].HouseholdID == 0 && Register[PID - 1].ParentID[0] > 0){
						PID2 = Register[PID - 1].ParentID[0];
						HHID = Register[PID2 - 1].HouseholdID;
						if (HHregister[HHID - 1].IDhead == PID2){ // Grandpa is a HH head
							HHregister[HHID - 1].AddMember(ic + 1);
							Register[ic].HouseholdID = HHID;
						}
					}
				}
			}
			if (Register[ic].HouseholdID == 0){
				// As a last resort, assign them to the household of the parent, even if the parent is 
				// not a household head. There should be at least 1 surviving parent at the start of the
				// simulation, and all adults should be assigned to a household by this point.
				PID = Register[ic].ParentID[1];
				if (PID > 0){
					HHID = Register[PID - 1].HouseholdID;
					HHregister[HHID - 1].AddMember(ic + 1);
					Register[ic].HouseholdID = HHID;
				}
				PID = Register[ic].ParentID[0];
				if (PID > 0 && Register[ic].HouseholdID == 0){
					HHID = Register[PID - 1].HouseholdID;
					HHregister[HHID - 1].AddMember(ic + 1);
					Register[ic].HouseholdID = HHID;
				}
			}
		}
	}
}

void Pop::AssignContraception()
{
	int ic, ir, ia, il, ip, ih;
	double BaseContr, ORcontr, ProbContr;

	for (ic = 0; ic<InitPop; ic++){
		r[ic] = rg.Random();}
	for (ic = 0; ic<InitPop; ic++){
		r2[ic] = rg.Random();}

	for (ic = 0; ic < InitPop; ic++){
		if (Register[ic].AgeGroup>2 && Register[ic].AgeGroup < 10 &&
			Register[ic].SexInd==1 && Register[ic].DOLB < 1985.5){
			ir = Register[ic].Race;
			ih = Register[ic].HighestGrade;
			ia = 0;
			if (Register[ic].AgeGroup>4){ ia = 1; }
			if (Register[ic].AgeGroup>6){ ia = 2; }
			ip = 0;
			if (Register[ic].DOLB == 0){ ip = 1; }
			// First determine whether using contraception
			if (Register[ic].VirginInd == 1){ ProbContr = ContrVirgin; }
			else{ 
				BaseContr = InitContr[ia][ir]; 
				il = 1 - Register[ic].MarriedInd;
				ORcontr = ORcontrUnmarriedNeverPreg[il][ip][ir];
				ORcontr *= pow(ORcontrEdu, ih - 10);
				if (Register[ic].CurrPartners == 0 && Register[ic].FSWind == 0){
					ORcontr *= ORcontrAbstinent;}
				// Still need to add code to take account of absent partner, urban-rural dif
				ProbContr = 1.0 / (1.0 + ORcontr * (1.0 - BaseContr) / BaseContr);
			}
			if (r[ic] > ProbContr){ Register[ic].CurrContr = 0; }
			else{
				// Then determine what type of contraception
				BaseContr = InitInjectable[ia][ir];
				ORcontr = pow(ORinjectableEdu, ih - 10);
				if (ip == 1){ ORcontr *= ORinjectableNeverPreg; }
				ProbContr = 1.0 / (1.0 + ORcontr * (1.0 - BaseContr) / BaseContr);
				if (r2[ic] < ProbContr){ 
					Register[ic].CurrContr = 1; // Using injectable
					Register[ic].EverInjectable = 1;
				} 
				else{
					// Check if using pill or sterilized
					if (ip == 1){ Register[ic].CurrContr = 2; }
					else{
						r2[ic] = (r2[ic] - ProbContr) / (1.0 - ProbContr);
						ia = Register[ic].AgeGroup - 3;
						if (r2[ic] < SterilizationNotPill[ia]){ Register[ic].CurrContr = 3; }
						else{ 
							Register[ic].CurrContr = 2; 
							Register[ic].EverPill = 1;
						}
					}
				}
			}
		}
	}

	// Store outputs
	if (FixedUncertainty == 1){
		for (ic = 0; ic < 36; ic++){
			ContrPrev.out[CurrSim - 1][ic] = 0;}
		for (ic = 0; ic < InitPop; ic++){
			if (Register[ic].AgeGroup>2 && Register[ic].AgeGroup < 10 &&
				Register[ic].SexInd == 1 && Register[ic].DOLB < 1985.5 &&
				(Register[ic].DOLB>0 || Register[ic].MarriedInd==1)){
				ir = Register[ic].Race;
				ia = 0;
				if (Register[ic].AgeGroup>4){ ia = 1; }
				if (Register[ic].AgeGroup>6){ ia = 2; }
				ih = Register[ic].CurrContr;
				ContrPrev.out[CurrSim - 1][ir*12+ia*4+ih] += 1;
			}
		}
	}
}

void Pop::AssignPrison()
{
	int ic, ir, ia, ia2, ind;
	double NeverInPrison[76]; // % never in prison at baseline, ages 15-90
	double InPrison[76][2]; // % in prison at baseline, ages 15-90 (1st index), 1st time vs repeat
	double Released[76]; // % previously in prison at baseline, ages 15-90
	double ExitRate, SentenceDur, rtemp, a, b, p, q, x;

	ind = 2;
	ExitRate = 1.0 / ((1.0 - UnsentencedRelease) * (SentenceLength[0] / SentenceLength[1]) * 
		((1.0 - ProbParole) + ProbParole * 0.5) + MeanUnsentencedDur);
	NeverInPrison[0] = 1.0;
	InPrison[0][0] = 0.0;
	InPrison[0][1] = 0.0;
	Released[0] = 0.0;
	for (ia = 1; ia < 76; ia++){
		ia2 = (ia - 1) / 5;
		NeverInPrison[ia] = NeverInPrison[ia - 1] * exp(-PrisonEntryRate[ia2]);
		InPrison[ia][0] = (NeverInPrison[ia - 1] - NeverInPrison[ia]) * exp(-0.5 * ExitRate) +
			InPrison[ia - 1][0] * exp(-ExitRate);
		InPrison[ia][1] = Released[ia - 1] * (1.0 - exp(-PrisonEntryRate[ia2] * PrisonReentryAdj)) * 
			exp(-0.5 * ExitRate) + InPrison[ia - 1][1] * exp(-ExitRate);
		Released[ia] = 1.0 - NeverInPrison[ia] - InPrison[ia][0] - InPrison[ia][1];
	}

	for (ic = 0; ic<InitPop; ic++){
		r[ic] = rg.Random();}
	for (ic = 0; ic<InitPop; ic++){
		r2[ic] = rg.Random();}

	for (ic = 0; ic < InitPop; ic++){
		if (Register[ic].AgeGroup>2 && Register[ic].SexInd == 0 && Register[ic].CasualInd == 0 
			&& Register[ic].HetCasualInd == 0){
			ir = Register[ic].Race;
			ia = Register[ic].CurrAge - 15;
			// Determine whether in prison
			if (r[ic] < InPrison[ia][0] * PrisonEntryRace[ir]){
				Register[ic].Imprisoned = 2;
				Register[ic].PrevImprisoned = 0;
				SentenceDur = -log(r2[ic])/ExitRate;
				Register[ic].ReleaseDate = 1985.5 + SentenceDur;
			}
			else if (r[ic] < (InPrison[ia][0] + InPrison[ia][1]) * PrisonEntryRace[ir]){
				Register[ic].Imprisoned = 2;
				Register[ic].PrevImprisoned = 1;
				SentenceDur = -log(r2[ic]) / ExitRate;
				Register[ic].ReleaseDate = 1985.5 + SentenceDur;
			}
			else if (r[ic] < (InPrison[ia][0] + InPrison[ia][1] + Released[ia]) * PrisonEntryRace[ir]){
				Register[ic].Imprisoned = 0;
				Register[ic].PrevImprisoned = 1;
				// Arbitrarily set previous release date to start of 1985
				Register[ic].ReleaseDate = 1985.0; 
			}
			else{
				Register[ic].Imprisoned = 0;
				Register[ic].PrevImprisoned = 0;
			}
		}
	}
}

void Pop::SetEmployment()
{
	int ic, ia, ig, ih, ir, iy;
	double odds, testprob, newrand, EmployRate;

	for (ic = 0; ic<Register.size(); ic++){
		r2[ic] = rg.Random();
	}
	iy = CurrYear - StartYear;

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AgeGroup < 3 || Register[ic].InSchool == 1 || Register[ic].AliveInd == 0){
			Register[ic].Employed = 0;
		}
		else{
			odds = BaseEmployed;
			ig = Register[ic].SexInd;
			odds *= SexEffectEmployed[ig];
			if (Register[ic].CurrUrban == 1){ odds *= UrbanEffectEmployed[1]; }
			ir = Register[ic].Race;
			odds *= RaceEffectEmployed[ir];
			ia = Register[ic].AgeGroup - 3;
			if (ia > 10){ ia = 10; }
			odds *= AgeEffectEmployed[ia];
			ih = Register[ic].HighestGrade;
			odds *= EduEffectEmployed[ih];
			if (r2[ic] < (odds / (1.0 + odds))){ Register[ic].Employed = 1; }
			else{ Register[ic].Employed = 0; }
		}
	}
}

void Pop::ChooseLTpartner(int ID, double rnd1, double rnd2)
{
	// This function is used ONLY at baseline, to match married individuals.

	int ic, ia;
	int iage, page; // individual's age group - 2 (so that 0 ==> age 10-14), partner's age group - 2
	int igender, pgender; // individual's gender and partner's gender
	int irace; // individual's race
	int irisk, prisk; // individual's risk group -1 (so that 0 ==> high risk), partner's risk
	int iurban; // individual's urban/rural location
	int TotalRisk[2]; // total potential LT partners by risk group
	int EligibleByAge[16];
	int MaxMSM = 0;
	int igaym;
	double WeightsByAge[16];
	double Normalizer, TotalSelectionProb2;

	iage = Register[ID-1].AgeGroup - 2;
	igender = Register[ID-1].SexInd;
	irace = Register[ID - 1].Race;
	if (Register[ID - 1].SexInd == 1 || Register[ID - 1].MalePref <1.0){ pgender = 1 - igender; }
	else{ pgender = 0; }
	irisk = Register[ID-1].RiskGroup - 1;
	igaym = 0;
	if (Register[ID - 1].MalePref == 1.0 && igender == 0){ igaym = 1; }
	iurban = Register[ID-1].CurrUrban; 

	int tpp = Register.size();

	// (a) First choose risk group of LT partner
	// =========================================

	TotalRisk[0] = 0;
	TotalRisk[1] = 0;
	Normalizer = 0.0;
	for(ic=0; ic<tpp; ic++){
		if(Register[ic].MarriedInd==1 && Register[ic].IDprimary==0 && 
			Register[ic].SexInd==pgender && Register[ic].Race==irace &&
			(Register[ic].CurrUrban == iurban || Register[ic].Visiting==1)){
			if ((igender == 0 && igaym == 0) || // Straight/bi men
				(igender == 1 && Register[ic].MalePref < 1.0) || // Women
				(igaym == 1 && Register[ic].MalePref > 0.0 && ic!=(ID-1))){ // Gay men
					if (Register[ic].RiskGroup == 1){ TotalRisk[0] += 1; }
					else{ TotalRisk[1] += 1; }
			}
		}
	}
	if (igaym == 1){ // Seeking same-sex partner
		if (rnd1 < ActualPropnLTH_MSM[irisk]){
			if (TotalRisk[0] > 0){ prisk = 1; }
			else if (TotalRisk[1] > 0){ prisk = 2; }
			else{ MaxMSM = 1; }
		}
		else{
			if (TotalRisk[1] > 0){ prisk = 2; }
			else if (TotalRisk[0] > 0){ prisk = 1; }
			else{ MaxMSM = 1; }
		}
	}
	else{ // Seeking heterosexual partner
		if (rnd1 < ActualPropnLTH[irisk][igender]){
			if (TotalRisk[0] > 0){ prisk = 1; }
			else if (TotalRisk[1] > 0){ prisk = 2; }
			else{ MaxNewPartnerInd = 1; }
		}
		else{
			if (TotalRisk[1] > 0){ prisk = 2; }
			else if (TotalRisk[0] > 0){ prisk = 1; }
			else{ MaxNewPartnerInd = 1; }
		}
	}
	
	// (b) Next choose individual partner
	// ==================================

	// (b1) Calculate age weights
	if(MaxNewPartnerInd==0 && MaxMSM==0){
		for(ia=0; ia<16; ia++){
			EligibleByAge[ia] = 0;}
		for (ic = 0; ic < tpp; ic++){
			if (Register[ic].MarriedInd == 1 && Register[ic].IDprimary == 0 &&
				Register[ic].SexInd == pgender && Register[ic].RiskGroup == prisk &&
				Register[ic].Race == irace && (Register[ic].CurrUrban==iurban ||
				Register[ic].Visiting==1)){
				if ((igender == 0 && igaym == 0) || // Straight/bi men
					(igender == 1 && Register[ic].MalePref < 1.0) || // Women
					(igaym == 1 && Register[ic].MalePref > 0.0 && ic != (ID - 1))){ // Gay men
						ia = Register[ic].AgeGroup - 2; // -2 because age groups start at 10
						EligibleByAge[ia] += 1;
				}
			}
		}
		for(ia=0; ia<16; ia++){
			if(igender==0 && igaym==0 && EligibleByAge[ia] > 0){
				WeightsByAge[ia] = AgePrefM[iage][ia]/EligibleByAge[ia];
				Normalizer += AgePrefM[iage][ia];
			}
			else if(igender==1 && EligibleByAge[ia] > 0){
				WeightsByAge[ia] = AgePrefF[iage][ia]/EligibleByAge[ia];
				Normalizer += AgePrefF[iage][ia];
			}
			else if (igaym == 1 && EligibleByAge[ia] > 0){
				WeightsByAge[ia] = AgePrefMSM[iage][ia] / EligibleByAge[ia];
				Normalizer += AgePrefMSM[iage][ia];
			}
			else{
				WeightsByAge[ia] = 0.0;}
		}
		for(ia=0; ia<16; ia++){
			WeightsByAge[ia] = WeightsByAge[ia]/Normalizer;}
	}

	// (b2) Calculate individual weights
	if(MaxNewPartnerInd==0 && MaxMSM==0 && Normalizer>0.0){
		TotalSelectionProb2 = 0.0;
		for(ic=0; ic<tpp; ic++){
			if(Register[ic].MarriedInd==1 && Register[ic].IDprimary==0 && 
				Register[ic].SexInd==pgender && Register[ic].RiskGroup==prisk &&
				Register[ic].Race == irace && (Register[ic].CurrUrban == iurban ||
				Register[ic].Visiting == 1)){
				if ((igender == 0 && igaym == 0) || // Straight/bi men
					(igender == 1 && Register[ic].MalePref < 1.0) || // Women
					(igaym == 1 && Register[ic].MalePref > 0.0 && ic != (ID - 1))){ // Gay men
						page = Register[ic].AgeGroup - 2;
						TotalSelectionProb2 += WeightsByAge[page];
						Register[ic].CumSelectionProb = TotalSelectionProb2;
				}
			}
		}
		if(rnd2>TotalSelectionProb2){ // Due to rounding there is a very small chance
			// that rnd2 may be > TotalSelectionProb2
			rnd2 = TotalSelectionProb2;}
	}

	// (b3) Choose partner
	if(MaxNewPartnerInd==0){
		if(Normalizer==0 || MaxMSM==1){ // The only available partners are in ineligble age groups
										// or we have run out of potential MSM partners
			Register[ID-1].MarriedInd = 0;}
		else{
			for(ic=0; ic<tpp; ic++){
				if(Register[ic].MarriedInd==1 && Register[ic].IDprimary==0 && 
					Register[ic].SexInd==pgender && Register[ic].RiskGroup==prisk &&
					Register[ic].Race == irace && (Register[ic].CurrUrban == iurban ||
					Register[ic].Visiting == 1) && rnd2 <= Register[ic].CumSelectionProb){
					if ((igender == 0 && igaym == 0) || // Straight/bi men
						(igender == 1 && Register[ic].MalePref < 1.0) || // Women
						(igaym == 1 && Register[ic].MalePref > 0.0 && ic != (ID - 1))){ // Gay men
							Register[ID - 1].IDprimary = ic + 1;
							Register[ic].IDprimary = ID;
							break;
					}
				}
			}
		}
	}
	else{
		Register[ID-1].MarriedInd = 0;}
}

void Pop::ChooseSTpartner(int ID, double rnd1, double rnd2)
{
	// This function is used ONLY at baseline, to match individuals in short-term relationships.

	int ic, ia;
	int iage, page; // individual's age group - 2 (so that 0 ==> age 10-14), partner's age group - 2
	int igender, pgender; // individual's gender and partner's gender
	int irisk, prisk; // individual's risk group -1 (so that 0 ==> high risk), partner's risk
	int irace; // individual's race
	int iurban; // individual's urban-rural location
	int igaym; // 1 for men who are exclusively homosexual, 0 otherwise
	int TotalRisk[2]; // total potential ST partners by risk group
	int EligibleByAge[16];
	int MaxMSM = 0;
	double WeightsByAge[16];
	double Normalizer, TotalSelectionProb2;

	iage = Register[ID - 1].AgeGroup - 2;
	igender = Register[ID - 1].SexInd;
	irace = Register[ID - 1].Race;
	iurban = Register[ID - 1].CurrUrban;
	igaym = 0;
	if (Register[ID - 1].MalePref == 1.0 && igender == 0){ igaym = 1; }
	if (igaym == 0){ pgender = 1 - igender; }
	else{ pgender = 0; }
	irisk = Register[ID - 1].RiskGroup - 1;

	// (a) First choose risk group of ST partner
	// =========================================

	TotalRisk[0] = 0;
	TotalRisk[1] = 0;
	Normalizer = 0.0;
	for (ic = 0; ic<InitPop; ic++){
		if (Register[ic].CurrPartners>0 && Register[ic].IDprimary == 0 && 
			Register[ic].SexInd == pgender && Register[ic].Race == irace && 
			(Register[ic].CurrUrban == iurban || Register[ic].Visiting == 1)){
			if ((igender == 0 && igaym == 0) || // Straight/bi men
				(igender == 1 && Register[ic].MalePref < 1.0) || // Women
				(igaym == 1 && Register[ic].MalePref > 0.0 && ic != (ID - 1))){ // Gay men
				if (Register[ic].RiskGroup == 1){
					TotalRisk[0] += 1;}
				else{
					TotalRisk[1] += 1;}
				if (Register[ic].CurrPartners == 2 && Register[ic].ID2ndary == 0){
					TotalRisk[0] += 1;}
			}
		}
	}
	if (igaym == 0){
		if (rnd1 < ActualPropnSTH[irisk][igender]){
			if (TotalRisk[0] > 0){ prisk = 1; }
			else if (TotalRisk[1] > 0){ prisk = 2; }
			else{ MaxNewPartnerInd = 1; }
		}
		else{
			if (TotalRisk[1] > 0){ prisk = 2; }
			else if (TotalRisk[0] > 0){ prisk = 1; }
			else{ MaxNewPartnerInd = 1; }
		}
	}
	else{
		if (rnd1 < ActualPropnSTH_MSM[irisk]){
			if (TotalRisk[0] > 0){ prisk = 1; }
			else if (TotalRisk[1] > 0){ prisk = 2; }
			else{ MaxMSM = 1; }
		}
		else{
			if (TotalRisk[1] > 0){ prisk = 2; }
			else if (TotalRisk[0] > 0){ prisk = 1; }
			else{ MaxMSM = 1; }
		}
	}
	
	// (b) Next choose individual partner
	// ==================================

	// (b1) Calculate age weights
	if(MaxNewPartnerInd==0 && MaxMSM==0){
		for(ia=0; ia<16; ia++){
			EligibleByAge[ia] = 0;}
		for(ic=0; ic<InitPop; ic++){
			if(Register[ic].CurrPartners>0 && Register[ic].IDprimary==0 && 
				Register[ic].SexInd==pgender && Register[ic].RiskGroup==prisk &&
				Register[ic].Race == irace && (Register[ic].CurrUrban == iurban ||
				Register[ic].Visiting == 1)){
				if ((igender == 0 && igaym == 0) || // Straight/bi men
					(igender == 1 && Register[ic].MalePref < 1.0) || // Women
					(igaym == 1 && Register[ic].MalePref > 0.0 && ic != (ID - 1))){ // Gay men
						ia = Register[ic].AgeGroup - 2; // -2 because age groups start at 10
						EligibleByAge[ia] += 1;
				}
			}
			if(Register[ic].CurrPartners==2 && Register[ic].ID2ndary==0 && 
				Register[ic].SexInd==pgender && Register[ic].RiskGroup==prisk &&
				Register[ic].Race == irace && (Register[ic].CurrUrban == iurban ||
				Register[ic].Visiting == 1)){
				if ((igender == 0 && igaym == 0) || // Straight/bi men
					(igender == 1 && Register[ic].MalePref < 1.0) || // Women
					(igaym == 1 && Register[ic].MalePref > 0.0 && ic != (ID - 1))){ // Gay men
						ia = Register[ic].AgeGroup - 2;
						EligibleByAge[ia] += 1;
				}
			}
		}
		for(ia=0; ia<16; ia++){
			if(igender==0 && igaym==0 && EligibleByAge[ia] > 0){
				WeightsByAge[ia] = AgePrefM[iage][ia]/EligibleByAge[ia];
				Normalizer += AgePrefM[iage][ia];
			}
			else if (igender == 0 && igaym == 1 && EligibleByAge[ia] > 0){
				WeightsByAge[ia] = AgePrefMSM[iage][ia] / EligibleByAge[ia];
				Normalizer += AgePrefMSM[iage][ia];
			}
			else if(igender==1 && EligibleByAge[ia] > 0){
				WeightsByAge[ia] = AgePrefF[iage][ia]/EligibleByAge[ia];
				Normalizer += AgePrefF[iage][ia];
			}
			else{
				WeightsByAge[ia] = 0.0;}
		}
		for(ia=0; ia<16; ia++){
			WeightsByAge[ia] = WeightsByAge[ia]/Normalizer;}
	}

	// (b2) Calculate individual weights
	if(MaxNewPartnerInd==0 && MaxMSM==0 && Normalizer>0.0){
		TotalSelectionProb2 = 0.0;
		for(ic=0; ic<InitPop; ic++){
			if(Register[ic].CurrPartners>0 && Register[ic].SexInd==pgender && 
				Register[ic].RiskGroup == prisk && Register[ic].Race == irace && 
				(Register[ic].CurrUrban == iurban || Register[ic].Visiting == 1)){
				if ((igender == 0 && igaym == 0) || // Straight/bi men
					(igender == 1 && Register[ic].MalePref < 1.0) || // Women
					(igaym == 1 && Register[ic].MalePref > 0.0 && ic != (ID - 1))){ // Gay men
						page = Register[ic].AgeGroup - 2;
						double temp = TotalSelectionProb2;
						if (Register[ic].IDprimary == 0){
							TotalSelectionProb2 += WeightsByAge[page];}
						if (Register[ic].CurrPartners == 2 && Register[ic].ID2ndary == 0){
							TotalSelectionProb2 += WeightsByAge[page];}
						Register[ic].CumSelectionProb = TotalSelectionProb2;
				}
			}
		}
		if(rnd2>TotalSelectionProb2){ // Due to rounding there is a very small chance
			// that rnd2 may be > TotalSelectionProb2
			rnd2 = TotalSelectionProb2;}
	}

	// (b3) Choose partner
	if(MaxNewPartnerInd==0){
		if(Normalizer==0.0 || MaxMSM==1){ // The only available partners are in ineligble age groups
										  // or we have run out of potential MSM partners
			if(Register[ID-1].CurrPartners>0 && Register[ID-1].IDprimary==0){
				Register[ID-1].CurrPartners = 0;}
			else if(Register[ID-1].CurrPartners==2 && Register[ID-1].ID2ndary==0){
				Register[ID-1].CurrPartners = 1;}
			else{
				cout<<"Error 1:indiv shouldn't be assigned a new partner"<<endl;}
		}
		else{
			for(ic=0; ic<InitPop; ic++){
				if(Register[ic].CurrPartners>0 && Register[ic].SexInd==pgender && 
					Register[ic].RiskGroup == prisk && Register[ic].Race == irace && 
					(Register[ic].CurrUrban == iurban || Register[ic].Visiting == 1) &&
					rnd2 <= Register[ic].CumSelectionProb){
					if ((igender == 0 && igaym == 0) || // Straight/bi men
						(igender == 1 && Register[ic].MalePref < 1.0) || // Women
						(igaym == 1 && Register[ic].MalePref > 0.0 && ic != (ID - 1))){ // Gay men
							if (Register[ID - 1].IDprimary == 0){Register[ID - 1].IDprimary = ic + 1;}
							else{Register[ID - 1].ID2ndary = ic + 1;}
							if (Register[ic].IDprimary == 0){Register[ic].IDprimary = ID;}
							else{Register[ic].ID2ndary = ID;}
							break;
					}
				}
			}
		}
	}
	else{ // No partners available in either risk group
		if(Register[ID-1].CurrPartners>0 && Register[ID-1].IDprimary==0){
			Register[ID-1].CurrPartners = 0;}
		else if(Register[ID-1].CurrPartners==2 && Register[ID-1].ID2ndary==0){
			Register[ID-1].CurrPartners = 1;}
		else{
			cout<<"Error 2:indiv shouldn't be assigned a new partner"<<endl;}
	}
}

void Pop::AssignHIV()
{
	int ic, InitHigh15to49, InitHighMarried, InitHighUnmarried;
	int CumHighRisk1, CumHighRisk2, index1, nextID;

	// Note that this function hasn't been updated to be consistent with AssignHIV1990,
	// since it seems unlikely that we will need to use this function in future.

	//int seed = 6810; // Note that I've arbitrarily fixed the seeds for now.
	//if(FixedPeriod==0){
	//	seed += CurrSim;}
	//CRandomMersenne rg(seed);
	for(ic=0; ic<InitPop; ic++){ 
		r[ic] = rg.Random();
	}

	for(ic=0; ic<InitPop; ic++){ 
		if(Register[ic].RiskGroup==1 && Register[ic].VirginInd==0 && Register[ic].AgeGroup>=3 &&
			Register[ic].AgeGroup<10 && r[ic]<InitHIVprevHigh && ConstantInitialHIV==0){
				Register[ic].HIVstage = 2; // To be consistent with TSHISA
				Register[ic].DateInfect = StartYear; // Assuming they acquired HIV 6 months prior
													 // to the start of the simulation 
			}
		else{
			Register[ic].HIVstage = 0;}
	}
	if(ConstantInitialHIV==1){
		InitHigh15to49 = 0;
		InitHighMarried = 0;
		for(ic=0; ic<InitPop; ic++){
			if(Register[ic].RiskGroup==1 && Register[ic].VirginInd==0 && Register[ic].AgeGroup>=3 &&
				Register[ic].AgeGroup<10){
					InitHigh15to49 += 1;
					if(Register[ic].MarriedInd==1){
						InitHighMarried += 1;}
				}
		}
		InitHighUnmarried = InitHigh15to49-InitHighMarried;
		vector<int> indices1(InitHighUnmarried);
		vector<int> indices2(InitHighMarried);
		CumHighRisk1 = 0;
		CumHighRisk2 = 0;
		for(ic=0; ic<InitPop; ic++){
			if(Register[ic].RiskGroup==1 && Register[ic].VirginInd==0 && Register[ic].AgeGroup>=3 &&
				Register[ic].AgeGroup<10){
					// Define indices vectors to contain IDs of people in the high risk group
					if(Register[ic].MarriedInd==0){
						indices1[CumHighRisk1] = ic + 1;
						CumHighRisk1 += 1;
					}
					else{
						indices2[CumHighRisk2] = ic + 1;
						CumHighRisk2 += 1;
					}
				}
		}
		ic = 0;
		while(ic < (InitHIVprevHigh * InitHighUnmarried)){
			index1 = r[ic] * indices1.size();
			nextID = indices1[index1];
			Register[nextID-1].HIVstage = 2; // To be consistent with TSHISA
			Register[nextID-1].DateInfect = StartYear;
			indices1[index1] = indices1[indices1.size()-1];
			indices1.pop_back();
			ic += 1;
		}
		while(ic < (InitHIVprevHigh * InitHigh15to49)){
			index1 = r[ic] * indices2.size();
			nextID = indices2[index1];
			Register[nextID-1].HIVstage = 2; // To be consistent with TSHISA
			Register[nextID-1].DateInfect = StartYear;
			indices2[index1] = indices2[indices2.size()-1];
			indices2.pop_back();
			ic += 1;
		}
	}

	// Assign SuscepHIVadj
	int ind;
	double x, y, a, b, p, q;

	//int seed2 = 3814; // Note that I've arbitrarily fixed the seeds for now.
	//if(FixedPeriod==0){
	//	seed2 += CurrSim;}
	//CRandomMersenne rg2(seed2);
	for(ic=0; ic<InitPop; ic++){ 
		r[ic] = rg.Random();}
	for(ic=0; ic<InitPop; ic++){
		if(AllowHIVsuscepAdj==1){
			ind = 2;
			// Note that the following formulas for a and b apply only when the gamma mean is 1.
			a = 1.0/pow(SDsuscepHIVadj, 2.0);
			b = a;
			p = r[ic];
			q = 1 - r[ic];
			cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
			Register[ic].SuscepHIVadj = x;
		}
		else{
			Register[ic].SuscepHIVadj = 1.0;}
	}
}

void Pop::AssignHIV1990(int race)
{
	int ia, ig, ic, CumNewHIV, NewID;
	double ExpectedHIV, Temp, HighPrevByAge[7][2], MSMadj;

	//int seed = 6810; // Note that I've arbitrarily fixed the seeds for now.
	//if (FixedPeriod <= 5){
	//	seed += CurrSim;}
	//CRandomMersenne rg(seed);
	for (ic = 0; ic<Register.size(); ic++){
		r2[ic] = rg.Random();
	}
	for (ia = 0; ia < 7; ia++){
		for (ig = 0; ig < 2; ig++){
			HighPrevByAge[ia][ig] = HlabisaRatio[ia][ig] * InitHIVprevHigh * InitHIVprevRace[race];}
	}

	if (ConstantInitialHIV == 0){ // Allowing stochastic variation in # infections
		for (ic = 0; ic < Register.size(); ic++){
			if (Register[ic].RiskGroup == 1 && Register[ic].VirginInd == 0 && 
				Register[ic].AgeGroup >= 3 && Register[ic].AgeGroup < 10 && Register[ic].AliveInd==1 &&
				Register[ic].Race==race){
				// Note that we're not making an urban-rural adjustment if ConstantInitialHIV = 0
					ia = Register[ic].AgeGroup - 3;
					ig = Register[ic].SexInd;
					if (Register[ic].SexInd == 1 || Register[ic].MalePref == 0.0){ MSMadj = 1.0; }
					else{ MSMadj = RatioInitPrevMSM; }
					if (r2[ic] < HighPrevByAge[ia][ig] * MSMadj){
						Register[ic].HIVstage = 2; // To be consistent with TSHISA
						Register[ic].GetInitCD4andVL();
						Register[ic].DateInfect = 1990.0; // Assuming they acquired HIV 
							// 6 months prior to mid-1990 
					}
			}
		}
	}
	if (ConstantInitialHIV == 1){ // Initial # HIV infections is fixed
		ExpectedHIV = 0.0;
		for (ic = 0; ic<Register.size(); ic++){
			if (Register[ic].RiskGroup == 1 && Register[ic].VirginInd==0 && Register[ic].AgeGroup >= 3 &&
				Register[ic].AgeGroup<10 && Register[ic].AliveInd==1 && Register[ic].Race==race &&
				Register[ic].CurrUrban==1){
					ia = Register[ic].AgeGroup - 3;
					ig = Register[ic].SexInd;
					if (Register[ic].SexInd == 1 || Register[ic].MalePref == 0.0){ MSMadj = 1.0; }
					else{ MSMadj = RatioInitPrevMSM; }
					ExpectedHIV += HighPrevByAge[ia][ig] * MSMadj;
					Register[ic].CumSelectionProb = ExpectedHIV;
			}
		}
		CumNewHIV = 0;
		while (CumNewHIV < ExpectedHIV){
			for (ic = 0; ic < Register.size(); ic++){
				Temp = Register[ic].CumSelectionProb / ExpectedHIV;
				if (r2[CumNewHIV]<Temp && Register[ic].RiskGroup == 1 &&
					Register[ic].VirginInd==0 && Register[ic].AgeGroup >= 3 &&
					Register[ic].AgeGroup < 10 && Register[ic].AliveInd == 1 &&
					Register[ic].Race==race && Register[ic].CurrUrban==1){
						NewID = ic;
						break;
				}
			}
			CumNewHIV += 1;
			Register[NewID].HIVstage = 2; // To be consistent with TSHISA
			Register[NewID].GetInitCD4andVL();
			Register[NewID].DateInfect = 1990.0; // Assuming they acquired HIV 6 months prior to mid-1990
		}
	}
}

void Pop::AssignSTIs()
{
	int ic, id, group, offset, IDprimary, ID2ndary, ir;
	int RiskPrimary, Risk2ndary; // risk groups of primary and 2ndary partners

	//int seed = 781; // Note that I've arbitrarily fixed the seeds for now.
	//if(FixedPeriod==0){
	//	seed += CurrSim;}
	//CRandomMersenne rg(seed);
	for(ic=0; ic<InitPop; ic++){ 
		for(id=0; id<8; id++){
			rSTI[ic][id] = rg.Random();}
	}
	for(ic=0; ic<InitPop; ic++){ 
		// First determine the behaviour group
		if(Register[ic].VirginInd==1){
			if(Register[ic].RiskGroup==1){group = 0;}
			else{group = 1;}
		}
		else{
			if(Register[ic].RiskGroup==1){group = 2;}
			else{group = 14;}
			if(Register[ic].CurrPartners==1){
				group += 1;
				IDprimary = Register[ic].IDprimary;
				RiskPrimary = Register[IDprimary-1].RiskGroup;
				if(Register[ic].MarriedInd==1){
					group += 2;}
				group += RiskPrimary - 1;
			}
			if(Register[ic].CurrPartners==2){
				group += 5;
				IDprimary = Register[ic].IDprimary;
				ID2ndary = Register[ic].ID2ndary;
				RiskPrimary = Register[IDprimary-1].RiskGroup;
				Risk2ndary = Register[ID2ndary-1].RiskGroup;
				if(Register[ic].MarriedInd==0){
					group += RiskPrimary + Risk2ndary - 2;}
				else{
					group += (RiskPrimary - 1) * 2 + Risk2ndary + 2;}
			}
			if(Register[ic].FSWind==1){group = 19;}
		}
		// Determine offset
		offset = group * 16 + Register[ic].AgeGroup - 2;
		// Assign status for individual STIs
		if(HSVind==1){
			if(Register[ic].AgeGroup<2){
				Register[ic].HSVstage = 0;}
			else{
				ir = Register[ic].Race;
				if (rSTI[ic][0]>HSVprevRace[ir]){ Register[ic].HSVstage = 0; }
				else{
					rSTI[ic][0] = rSTI[ic][0] / HSVprevRace[ir];
					if (Register[ic].SexInd == 0){
						Register[ic].HSVstage = HSVtransitionM.GetSTDstage(offset, rSTI[ic][0]);}
					else{
						Register[ic].HSVstage = HSVtransitionF.GetSTDstage(offset, rSTI[ic][0]);}
				}
			}
		}
		if(TPind==1){
			if(Register[ic].AgeGroup<2){
				Register[ic].TPstage = 0;}
			else if(Register[ic].SexInd==0){
				Register[ic].TPstage = TPtransitionM.GetSTDstage(offset, rSTI[ic][1]);}
			else{
				Register[ic].TPstage = TPtransitionF.GetSTDstage(offset, rSTI[ic][1]);}
		}
		if(HDind==1){
			if(Register[ic].AgeGroup<2){
				Register[ic].HDstage = 0;}
			else if(Register[ic].SexInd==0){
				Register[ic].HDstage = HDtransitionM.GetSTDstage(offset, rSTI[ic][2]);}
			else{
				Register[ic].HDstage = HDtransitionF.GetSTDstage(offset, rSTI[ic][2]);}
		}
		if(NGind==1){
			if(Register[ic].AgeGroup<2){
				Register[ic].NGstage = 0;}
			else if(Register[ic].SexInd==0){
				Register[ic].NGstage = NGtransitionM.GetSTDstage(offset, rSTI[ic][3]);}
			else{
				Register[ic].NGstage = NGtransitionF.GetSTDstage(offset, rSTI[ic][3]);}
		}
		if(CTind==1){
			if(Register[ic].AgeGroup<2){
				Register[ic].CTstage = 0;}
			else if(Register[ic].SexInd==0){
				Register[ic].CTstage = CTtransitionM.GetSTDstage(offset, rSTI[ic][4]);}
			else{
				Register[ic].CTstage = CTtransitionF.GetSTDstage(offset, rSTI[ic][4]);}
		}
		if(TVind==1){
			if(Register[ic].AgeGroup<2){
				Register[ic].TVstage = 0;}
			else if(Register[ic].SexInd==0){
				Register[ic].TVstage = TVtransitionM.GetSTDstage(offset, rSTI[ic][5]);}
			else{
				Register[ic].TVstage = TVtransitionF.GetSTDstage(offset, rSTI[ic][5]);}
		}
		if(BVind==1){
			if(Register[ic].AgeGroup<2 || Register[ic].SexInd==0){
				Register[ic].BVstage = 0;}
			else{
				Register[ic].BVstage = BVtransitionF.GetSTDstage(offset, rSTI[ic][6]);}
		}
		if(VCind==1){
			if(Register[ic].AgeGroup<2 || Register[ic].SexInd==0){
				Register[ic].VCstage = 0;}
			else{
				Register[ic].VCstage = VCtransitionF.GetSTDstage(offset, rSTI[ic][7]);}
		}
	}
}

void Pop::ResetToBaselineStart()
{
	int ic, HHID;

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1){
			Register[ic].BaselineAge = Register[ic].CurrAge;
			// BaselinePoor is now calculated in the UpdateHHincome function
			/*HHID = Register[ic].HouseholdID;
			if (HHID > 0){
				if (log(HHregister[HHID - 1].PerCapitaIncomeAdj) < MeanLogIncome){ 
					Register[ic].BaselinePoor = 1; }
				else{ Register[ic].BaselinePoor = 0; }
			}
			else{ Register[ic].BaselinePoor = 1; }*/
			Register[ic].BaselineSchool = Register[ic].InSchool;
			if (Register[ic].DrinksPerDD < 5.0 || Register[ic].DailyDrinkProb <= 1.0 / 30.0){
				Register[ic].BaselineBinge = 0;
			}
			else{ Register[ic].BaselineBinge = 1; }
			//Register[ic].BaselineUnemployed = 1 - Register[ic].Employed;
			Register[ic].BaselineVirgin = Register[ic].VirginInd;
			//Register[ic].CumUnprotected = 0;
			Register[ic].CumCasual = 0;
			Register[ic].CumSTIs = 0;
			Register[ic].CumHSV2 = 0;
			Register[ic].CumTeenPreg = 0;
		}
	}
}

void Pop::GetPopPyramid()
{
	int ia, ig, ic;

	for(ia=0; ia<18; ia++){
		PopPyramid[ia][0] = 0;
		PopPyramid[ia][1] = 0;
	}

	for(ic=0; ic<Register.size(); ic++){
		if(Register[ic].AliveInd==1){
			ia = Register[ic].AgeGroup;
			ig = Register[ic].SexInd;
			PopPyramid[ia][ig] += 1;
		}
	}
}

void Pop::GetBehavPyramid()
{
	int ia, ig, is, ic;

	for(ia=0; ia<18; ia++){
		for(ig=0; ig<2; ig++){
			for(is=0; is<3; is++){
				BehavPyramid[ia][ig*3+is] = 0;}
		}
	}

	for(ic=0; ic<Register.size(); ic++){
		if(Register[ic].AliveInd==1){
			ia = Register[ic].AgeGroup;
			ig = Register[ic].SexInd;
			if(Register[ic].VirginInd==1){
				BehavPyramid[ia][ig*3] += 1;}
			else if(Register[ic].MarriedInd==1){
				BehavPyramid[ia][ig*3+2] += 1;}
			else{
				BehavPyramid[ia][ig*3+1] += 1;}
		}
	}
}

void Pop::GetSexuallyExp()
{
	int ia, ig, ic, PID;
	double TotalYouth[11][2], SexuallyExp[11][2];
	double SexuallyActive[11][2], ActiveYounger[11][2], ActiveOlder[11][2];
	double AgeDif;

	for (ia = 0; ia<11; ia++){
		for (ig = 0; ig<2; ig++){
			TotalYouth[ia][ig] = 0.0;
			SexuallyExp[ia][ig] = 0.0;
			SexuallyActive[ia][ig] = 0.0;
			ActiveYounger[ia][ig] = 0.0; // Younger partner
			ActiveOlder[ia][ig] = 0.0; // Older partner
		}
	}

	if (CurrYear == 2005){
		for (ic = 0; ic < Register.size(); ic++){
			// Get number sexually experienced
			if (Register[ic].AliveInd == 1 && Register[ic].CurrAge >= 14 &&
				Register[ic].CurrAge < 25){
				ia = Register[ic].CurrAge - 14;
				ig = Register[ic].SexInd;
				TotalYouth[ia][ig] += Register[ic].PopWeight;
				if (Register[ic].VirginInd == 0){ SexuallyExp[ia][ig] += Register[ic].PopWeight; }
			}
			// Get number in age-disparate relationships
			if (Register[ic].AliveInd == 1 && Register[ic].CurrAge >= 15 &&
				Register[ic].CurrPartners>0){
				ia = Register[ic].AgeGroup - 3;
				if (ia > 10){ ia = 10; }
				ig = Register[ic].SexInd;
				PID = Register[ic].IDprimary;
				AgeDif = Register[ic].CurrAge - Register[PID - 1].CurrAge;
				SexuallyActive[ia][ig] += Register[ic].PopWeight;
				if (AgeDif >= 5.0){ ActiveYounger[ia][ig] += Register[ic].PopWeight; }
				if (AgeDif <= -5.0){ ActiveOlder[ia][ig] += Register[ic].PopWeight; }
			}
		}

		for (ia = 0; ia < 11; ia++){
			for (ig = 0; ig < 2; ig++){
				SexuallyExperienced.out[CurrSim - 1][ia + ig * 11] = 1.0 *
					SexuallyExp[ia][ig] / TotalYouth[ia][ig];
			}
			PartnerAgeDifM.out[CurrSim - 1][ia] = 1.0 * ActiveYounger[ia][0] / SexuallyActive[ia][0];
			PartnerAgeDifF.out[CurrSim - 1][ia] = 1.0 * ActiveYounger[ia][1] / SexuallyActive[ia][1];
			PartnerAgeDifM.out[CurrSim - 1][ia + 11] = 1.0 * ActiveOlder[ia][0] / SexuallyActive[ia][0];
			PartnerAgeDifF.out[CurrSim - 1][ia + 11] = 1.0 * ActiveOlder[ia][1] / SexuallyActive[ia][1];
		}
	}

	// Get proportion sexually experienced among adolescents
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].CurrAge >= 10 &&
			Register[ic].CurrAge < 20){
			TotalYouth[0][0] += Register[ic].PopWeight;
			if (Register[ic].VirginInd == 0){ SexuallyExp[0][0] += Register[ic].PopWeight; }
		}
	}
	SexuallyExpAdol.out[CurrSim - 1][CurrYear - StartYear] = SexuallyExp[0][0] / TotalYouth[0][0];
}

void Pop::GetAgeDisparate()
{
	int ic, ip, ia, PID;
	double SexuallyActive[2][2]; // by # partners (1 or 2) and by age disparity (0 = no disparity)
	double AgeDif;

	SexuallyActive[0][0] = 0.0;
	SexuallyActive[0][1] = 0.0;
	SexuallyActive[1][0] = 0.0;
	SexuallyActive[1][1] = 0.0;

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].CurrAge >= 15 &&
			Register[ic].CurrAge < 25 && Register[ic].SexInd==1){
			if (Register[ic].CurrPartners == 1){
				PID = Register[ic].IDprimary;
				AgeDif = Register[PID-1].CurrAge - Register[ic].CurrAge;
				if (AgeDif >= 5){ SexuallyActive[0][1] += Register[ic].PopWeight; }
				else{ SexuallyActive[0][0] += Register[ic].PopWeight; }
			}
			if (Register[ic].CurrPartners == 2){
				PID = Register[ic].ID2ndary;
				AgeDif = Register[PID - 1].CurrAge - Register[ic].CurrAge;
				if (AgeDif >= 5){ SexuallyActive[1][1] += Register[ic].PopWeight; }
				else{ SexuallyActive[1][0] += Register[ic].PopWeight; }
			}
		}
	}

	AgeDisparate2000.out[CurrSim - 1][0] = 1.0 * SexuallyActive[0][1] / (SexuallyActive[0][0] + SexuallyActive[0][1]);
	AgeDisparate2000.out[CurrSim - 1][1] = 1.0 * SexuallyActive[1][1] / (SexuallyActive[1][0] + SexuallyActive[1][1]);
}

void Pop::GetPrefMatrix(int type)
{
	int ic;
	int iage, page; // indiv age group (-2) and partner age group (-2)
	int partnerID;

	for(iage=0; iage<16; iage++){
		for(page=0; page<16; page++){
			PrefMatrix[iage][page] = 0;}
	}
	for(ic=0; ic<Register.size(); ic++){
		if(Register[ic].AliveInd==1 && Register[ic].SexInd==1){
			iage = Register[ic].AgeGroup - 2;
			if(type==1 && Register[ic].MarriedInd==1){
				partnerID = Register[ic].IDprimary;
				page = Register[partnerID-1].AgeGroup - 2;
				PrefMatrix[iage][page] += 1;
			}
			if(type==2 && Register[ic].MarriedInd==0 &&
				Register[ic].IDprimary>0){
					partnerID = Register[ic].IDprimary;
					page = Register[partnerID-1].AgeGroup - 2;
					PrefMatrix[iage][page] += 1;
			}
			if(type==2 && Register[ic].ID2ndary>0){
				partnerID = Register[ic].ID2ndary;
				page = Register[partnerID-1].AgeGroup - 2;
				PrefMatrix[iage][page] += 1;
			}
		}
	}
}

void Pop::GetNumberPartners()
{
	int ia, ig, is, ic;

	if (CurrSim == 1){
		for (ia = 0; ia < 16; ia++){
			for (ig = 0; ig < 2; ig++){
				for (is = 0; is < 5; is++){
					NumberPartners[ia][ig * 5 + is] = 0;
				}
			}
		}
	}

	for(ic=0; ic<Register.size(); ic++){
		if(Register[ic].AliveInd==1 && Register[ic].AgeGroup>=2 &&
			Register[ic].VirginInd == 0){
				ia = Register[ic].AgeGroup - 2;
				ig = Register[ic].SexInd;
				is = Register[ic].CurrPartners;
				if (IncludePopWeights == 0){
					if (Register[ic].MarriedInd == 1){
						NumberPartners[ia][ig * 5 + 2 + is] += 1;}
					else{
						NumberPartners[ia][ig * 5 + is] += 1;}
				}
				else{
					if (Register[ic].MarriedInd == 1){
						NumberPartners[ia][ig * 5 + 2 + is] += Register[ic].PopWeight;}
					else{
						NumberPartners[ia][ig * 5 + is] += Register[ic].PopWeight;}
				}
		}
	}
}

void Pop::GetTotSex()
{
	int ia, ig, ii, il, ic, ir;
	int CSWID, CSWsexActsTot[MaxCSWs][2][3];

	for(ia=0; ia<16; ia++){
		for(ig=0; ig<2; ig++){
			for(il=0; il<4; il++){
				TotSex[ia][ig*4+il] = 0;}
		}
	}

	for(ic=0; ic<Register.size(); ic++){
		if(Register[ic].AliveInd==1 && Register[ic].AgeGroup>=2 &&
			Register[ic].VirginInd==0){
				ia = Register[ic].AgeGroup - 2;
				ig = Register[ic].SexInd;
				if(Register[ic].MarriedInd==1){
					TotSex[ia][ig*4+2] += Register[ic].UVIprimary;
					TotSex[ia][ig*4+3] += Register[ic].PVIprimary;
				}
				else{
					TotSex[ia][ig*4] += Register[ic].UVIprimary;
					TotSex[ia][ig*4+1] += Register[ic].PVIprimary;
				}
				TotSex[ia][ig*4] += Register[ic].UVI2ndary;
				TotSex[ia][ig*4+1] += Register[ic].PVI2ndary;
				if(ig==0 && Register[ic].RiskGroup==1){
					if(Register[ic].UVICSW + Register[ic].PVICSW >0){
						CSWID = Register[ic].IDofCSW;
						ir = Register[ic].Race;
						for(ii=0; ii<TotCurrFSW[ir]; ii++){
							if(CSWregister[ii][ir]==CSWID){
								CSWsexActsTot[ii][0][ir] += Register[ic].UVICSW;
								CSWsexActsTot[ii][1][ir] += Register[ic].PVICSW;
								break;
							}
						}
					}
				}
		}
	}

	// Commercial sex acts
	/*for(ic=0; ic<TotCurrFSW; ic++){
		cout<<"Unprotected sex acts with CSW "<<ic+1<<": "<<CSWsexActsTot[ic][0]<<endl;
		cout<<"Protected sex acts with CSW "<<ic+1<<": "<<CSWsexActsTot[ic][1]<<endl;
	}*/
}

void Pop::GetCondomUse()
{
	int ii, il, ic, ir, ia, ig;
	double Protected[6][2], Unprotected[6][2], Unprot[5][3][2], SexuallyActive[5][3][2];

	// Condom use by educational attainment
	if (CurrYear == 1998 || CurrYear == 2003 || CurrYear == 2016){
		for (ii = 0; ii < 6; ii++){
			for (il = 0; il < 2; il++){
				Protected[ii][il] = 0.0;
				Unprotected[ii][il] = 0.0;
			}
		}
		for (ic = 0; ic < Register.size(); ic++){
			if (Register[ic].AliveInd == 1 && Register[ic].AgeGroup > 2 &&
				Register[ic].AgeGroup < 10 && Register[ic].VirginInd == 0 &&
				Register[ic].SexInd == 1){
				if (Register[ic].HighestGrade == 0){ ii = 0; }
				else if (Register[ic].HighestGrade < 6){ ii = 1; }
				else if (Register[ic].HighestGrade < 8){ ii = 2; }
				else if (Register[ic].HighestGrade < 12){ ii = 3; }
				else if (Register[ic].HighestGrade == 12){ ii = 4; }
				else{ ii = 5; }
				if (Register[ic].CurrPartners > 0){
					if (Register[ic].MarriedInd == 1){
						Unprotected[ii][0] += (1.0 - Register[ic].CondomPrimary) * Register[ic].PopWeight;
						Protected[ii][0] += Register[ic].CondomPrimary * Register[ic].PopWeight;
					}
					else{
						Unprotected[ii][1] += (1.0 - Register[ic].CondomPrimary) * Register[ic].PopWeight;
						Protected[ii][1] += Register[ic].CondomPrimary * Register[ic].PopWeight;
					}
				}
			}
		}
		if (CurrYear == 1998){
			for (ii = 0; ii < 6; ii++){
				CondomUseByEdu.out[CurrSim - 1][ii] = 1.0 * Protected[ii][0] /
					(Protected[ii][0] + Unprotected[ii][0]);
			}
			for (ii = 0; ii < 6; ii++){
				CondomUseByEdu.out[CurrSim - 1][ii + 6] = 1.0 * Protected[ii][1] /
					(Protected[ii][1] + Unprotected[ii][1]);
			}
		}
		if (CurrYear == 2003){
			for (ii = 0; ii < 6; ii++){
				CondomUseByEdu.out[CurrSim - 1][ii + 12] = 1.0 * Protected[ii][0] /
					(Protected[ii][0] + Unprotected[ii][0]);
			}
			for (ii = 0; ii < 6; ii++){
				CondomUseByEdu.out[CurrSim - 1][ii + 18] = 1.0 * Protected[ii][1] /
					(Protected[ii][1] + Unprotected[ii][1]);
			}
		}
		if (CurrYear == 2016){
			for (ii = 0; ii < 6; ii++){
				CondomUseByEdu.out[CurrSim - 1][ii + 24] = 1.0 * Protected[ii][1] /
					(Protected[ii][1] + Unprotected[ii][1]);
			}
		}
	}

	// Condom use by race
	if (CurrYear == 2003){
		for (ir = 0; ir < 3; ir++){
			Protected[ir][0] = 0.0;
			Unprotected[ir][0] = 0.0;
		}
		for (ic = 0; ic<Register.size(); ic++){
			if (Register[ic].AliveInd == 1 && Register[ic].AgeGroup > 2 &&
				Register[ic].AgeGroup < 10 && Register[ic].VirginInd == 0){
				ir = Register[ic].Race;
				if (Register[ic].CurrPartners > 0){
					Unprotected[ir][0] += (1.0 - Register[ic].CondomPrimary) * Register[ic].PopWeight;
					Protected[ir][0] += Register[ic].CondomPrimary * Register[ic].PopWeight;
				}
			}
		}
		for (ir = 0; ir < 3; ir++){
			CondomUseByRace.out[CurrSim - 1][ir] = 1.0 * Protected[ir][0] /
				(Protected[ir][0] + Unprotected[ir][0]);
		}
	}

	// Condom use by HIV diagnosis
	//if (CurrYear == 2005){GetUnprotected();}
}

void Pop::GetCondomUse2()
{
	int ic, ia;
	double protectedsex[2], unprotectedsex[2];

	for (ia = 0; ia < 2; ia++){
		protectedsex[ia] = 0.0; 
		unprotectedsex[ia] = 0.0;
	}

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 1 &&
			Register[ic].CurrPartners>0){
			if (Register[ic].AgeGroup >= 3 && Register[ic].AgeGroup <= 4){
				protectedsex[0] += Register[ic].CondomPrimary * Register[ic].PopWeight;
				unprotectedsex[0] += (1.0 - Register[ic].CondomPrimary) * Register[ic].PopWeight;
			}
			if (Register[ic].AgeGroup >= 5 && Register[ic].AgeGroup <= 9){
				protectedsex[1] += Register[ic].CondomPrimary * Register[ic].PopWeight;
				unprotectedsex[1] += (1.0 - Register[ic].CondomPrimary) * Register[ic].PopWeight;
			}
		}
	}

	CondomUse15to24.out[CurrSim - 1][CurrYear - StartYear] = protectedsex[0] /
		(protectedsex[0] + unprotectedsex[0]);
	CondomUse25to49.out[CurrSim - 1][CurrYear - StartYear] = protectedsex[1] /
		(protectedsex[1] + unprotectedsex[1]);
	CondomUse15to49.out[CurrSim - 1][CurrYear - StartYear] = (protectedsex[0] + protectedsex[1]) /
		(protectedsex[0] + unprotectedsex[0] + protectedsex[1] + unprotectedsex[1]);
}

void Pop::GetPrEP()
{
	int ic;
	double tot;

	tot = 0.0;
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].OnPrEP>0.0){ tot += Register[ic].PopWeight; }
	}

	TakingPrEP.out[CurrSim - 1][CurrYear - StartYear] = tot;
}

void Pop::GetNumbersByHIVstage()
{
	int ic, iy, ig, is;
	double tot;

	iy = CurrYear - StartYear;

	tot = 0.0;
	for (ig = 0; ig < 2; ig++){
		for (is = 0; is < 7; is++){
			AdultHIVstageTrend[iy][ig * 7 + is] = 0.0;}
	}
	for(ic=0; ic<Register.size(); ic++){
		if(Register[ic].AliveInd==1 && Register[ic].AgeGroup>=3){
			// Note that this gives total pop at ages 15+
			ig = Register[ic].SexInd;
			is = Register[ic].HIVstage;
			AdultHIVstageTrend[iy][ig * 7 + is] += Register[ic].PopWeight;
		}
		if (Register[ic].AliveInd == 1){ tot += Register[ic].PopWeight; }
	}
	TotPop.out[CurrSim - 1][iy] = tot;
}

void Pop::GetInitHIVprevH()
{
	int Positive[2], Negative[2];
	int ic, ig;

	Positive[0] = 0;
	Negative[0] = 0;
	Positive[1] = 0;
	Negative[1] = 0;
	for(ic=0; ic<InitPop; ic++){
		if(Register[ic].RiskGroup==1 && Register[ic].AgeGroup>2 && 
			Register[ic].AgeGroup<10 && Register[ic].VirginInd==0){
				ig = Register[ic].SexInd;
				if(Register[ic].HIVstage==0){
					Negative[ig] += 1;}
				else{
					Positive[ig] += 1;}
		}
	}
	cout<<"HIV-positive males in high risk group, 15-49: "<<Positive[0]<<endl;
	cout<<"HIV-negative males in high risk group, 15-49: "<<Negative[0]<<endl;
	cout<<"HIV-positive females in high risk group, 15-49: "<<Positive[1]<<endl;
	cout<<"HIV-negative females in high risk group, 15-49: "<<Negative[1]<<endl;
}

void Pop::GetHHprevByAge()
{
	int ic, ia, ig;
	double Total[9][2], Positive[9][2];

	for (ia = 0; ia < 9; ia++){
		for (ig = 0; ig < 2; ig++){
			Total[ia][ig] = 0.0;
			Positive[ia][ig] = 0.0;
		}
	}
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AgeGroup>2 && Register[ic].AgeGroup < 12 && Register[ic].AliveInd == 1){
			ia = Register[ic].AgeGroup - 3;
			ig = Register[ic].SexInd;
			Total[ia][ig] += Register[ic].PopWeight;
			if (Register[ic].HIVstage>0){
				Positive[ia][ig] += Register[ic].PopWeight;
			}
		}
	}
	for (ia = 0; ia < 9; ia++){
		for (ig = 0; ig < 2; ig++){
			if (CurrYear == 2005){
				PrevHH2005.out[CurrSim - 1][ia + 9 * ig] = 1.0 * Positive[ia][ig] / Total[ia][ig];}
			if (CurrYear == 2008){
				PrevHH2008.out[CurrSim - 1][ia + 9 * ig] = 1.0 * Positive[ia][ig] / Total[ia][ig];}
			if (CurrYear == 2012){
				PrevHH2012.out[CurrSim - 1][ia + 9 * ig] = 1.0 * Positive[ia][ig] / Total[ia][ig];}
		}
	}
}

void Pop::GetANCprevByAge()
{
	int ia, ic, is, ir, iu;
	double Total[7], Positive[7], SexExpByStage[5][7], temp1, temp2, UrbanIndByHIV[2][2];
	double RaceHIV[3][2];

	for (ia = 0; ia < 7; ia++){
		Total[ia] = 0.0;
		Positive[ia] = 0.0;
	}
	for (iu = 0; iu < 2; iu++){
		for (is = 0; is < 2; is++){ UrbanIndByHIV[iu][is] = 0.0; }
	}
	for (ir = 0; ir < 3; ir++){
		for (is = 0; is < 2; is++){ RaceHIV[ir][is] = 0.0; }
	}
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AgeGroup>2 && Register[ic].AgeGroup < 10 && 
			Register[ic].AliveInd == 1 && Register[ic].SexInd==1 &&
			Register[ic].DOLB>(1.0*CurrYear+0.5)){
			ia = Register[ic].AgeGroup - 3;
			ir = Register[ic].Race;
			is = 0;
			if (Register[ic].HIVstage > 0){ is = 1; }
			Total[ia] += Register[ic].PopWeight * PublicANCuse[ir]/(1.0 - StillbirthRate[is]);
			if (Register[ic].HIVstage>0){
				Positive[ia] += Register[ic].PopWeight * PublicANCuse[ir] / (1.0 - StillbirthRate[is]);
			}
			if (FixedUncertainty == 1){
				iu = Register[ic].CurrUrban;
				UrbanIndByHIV[iu][is] += Register[ic].PopWeight * PublicANCuse[ir] / (1.0 - StillbirthRate[is]);
				RaceHIV[ir][is] += Register[ic].PopWeight / (1.0 - StillbirthRate[is]);
			}
		}
	}
	if (FixedUncertainty == 1){
		// Calculate totals for purpose of PrevFertile calculations
		for (ia = 0; ia < 5; ia++){
			for (is = 0; is < 7; is++){
				SexExpByStage[ia][is] = 0.0;}
		}
		for (ic = 0; ic < Register.size(); ic++){
			if (Register[ic].AgeGroup>2 && Register[ic].AgeGroup < 8 &&
				Register[ic].AliveInd == 1 && Register[ic].SexInd == 1 &&
				Register[ic].VirginInd == 0){
				ia = Register[ic].AgeGroup - 3;
				is = Register[ic].HIVstage;
				SexExpByStage[ia][is] += Register[ic].PopWeight;
			}
		}
	}
	// Calculate HIV prevalence (with continuity correction)
	// The continuity isn't necessary since we apply a minimum of 0.0001 when we calculate ModelPrev.
	if (Positive[0] > 0.0){ PrevPreg15.out[CurrSim - 1][CurrYear - 1990] = Positive[0] / Total[0]; }
	else{ PrevPreg15.out[CurrSim - 1][CurrYear - 1990] = 0.5 / Total[0]; }
	if (Positive[1] > 0.0){ PrevPreg20.out[CurrSim - 1][CurrYear - 1990] = Positive[1] / Total[1]; }
	else{ PrevPreg20.out[CurrSim - 1][CurrYear - 1990] = 0.5 / Total[1]; }
	if (Positive[2] > 0.0){ PrevPreg25.out[CurrSim - 1][CurrYear - 1990] = Positive[2] / Total[2]; }
	else{ PrevPreg25.out[CurrSim - 1][CurrYear - 1990] = 0.5 / Total[2]; }
	if (Positive[3] > 0.0){ PrevPreg30.out[CurrSim - 1][CurrYear - 1990] = Positive[3] / Total[3]; }
	else{ PrevPreg30.out[CurrSim - 1][CurrYear - 1990] = 0.5 / Total[3]; }
	if (Positive[4] > 0.0){ PrevPreg35.out[CurrSim - 1][CurrYear - 1990] = Positive[4] / Total[4]; }
	else{ PrevPreg35.out[CurrSim - 1][CurrYear - 1990] = 0.5 / Total[4]; }
	if ((Positive[5] + Positive[6]) > 0.0){ 
		PrevPreg40.out[CurrSim - 1][CurrYear - 1990] = (Positive[5] + Positive[6]) / (Total[5] + Total[6]); }
	else{ PrevPreg40.out[CurrSim - 1][CurrYear - 1990] = 0.5 / (Total[5] + Total[6]); }
	/*PrevPreg15.out[CurrSim - 1][CurrYear - 1990] = (Positive[0] / Total[0]) * 0.977 + 1.0 - 0.977;
	PrevPreg20.out[CurrSim - 1][CurrYear - 1990] = (Positive[1] / Total[1]) * 0.977 + 1.0 - 0.977;
	PrevPreg25.out[CurrSim - 1][CurrYear - 1990] = (Positive[2] / Total[2]) * 0.977 + 1.0 - 0.977;
	PrevPreg30.out[CurrSim - 1][CurrYear - 1990] = (Positive[3] / Total[3]) * 0.977 + 1.0 - 0.977;
	PrevPreg35.out[CurrSim - 1][CurrYear - 1990] = (Positive[4] / Total[4]) * 0.977 + 1.0 - 0.977;
	PrevPreg40.out[CurrSim - 1][CurrYear - 1990] = ((Positive[5] + Positive[6]) / (Total[5] + Total[6])) * 
		0.977 + 1.0 - 0.977;*/
	temp1 = 0.0;
	temp2 = 0.0;
	for (ia = 0; ia < 7; ia++){
		temp1 += Total[ia];
		temp2 += Positive[ia];
	}
	//PrevPregTotal.out[CurrSim - 1][CurrYear - 1990] = (temp2 / temp1) * 0.977 + 1.0 - 0.977;
	PrevPregTotal.out[CurrSim - 1][CurrYear - 1990] = temp2 / temp1;
	if (FixedUncertainty == 1){
		for (ia = 0; ia < 5; ia++){
			temp1 = 0.0;
			for (is = 1; is < 7; is++){
				temp1 += SexExpByStage[ia][is] * RelHIVfertility[is - 1];
			}
			if (ia == 0){ PrevFertile15.out[CurrSim - 1][CurrYear - 1990] = temp1 / (temp1 + SexExpByStage[ia][0]); }
			if (ia == 1){ PrevFertile20.out[CurrSim - 1][CurrYear - 1990] = temp1 / (temp1 + SexExpByStage[ia][0]); }
			if (ia == 2){ PrevFertile25.out[CurrSim - 1][CurrYear - 1990] = temp1 / (temp1 + SexExpByStage[ia][0]); }
			if (ia == 3){ PrevFertile30.out[CurrSim - 1][CurrYear - 1990] = temp1 / (temp1 + SexExpByStage[ia][0]); }
			if (ia == 4){ PrevFertile35.out[CurrSim - 1][CurrYear - 1990] = temp1 / (temp1 + SexExpByStage[ia][0]); }
		}
		/*PrevPregRural.out[CurrSim - 1][CurrYear - 1990] = (UrbanIndByHIV[0][1] / (UrbanIndByHIV[0][0] + 
			UrbanIndByHIV[0][1])) * 0.977 + 1.0 - 0.977;
		PrevPregUrban.out[CurrSim - 1][CurrYear - 1990] = (UrbanIndByHIV[1][1] / (UrbanIndByHIV[1][0] +
			UrbanIndByHIV[1][1])) * 0.977 + 1.0 - 0.977;*/
		PrevPregRural.out[CurrSim - 1][CurrYear - 1990] = (UrbanIndByHIV[0][1] / (UrbanIndByHIV[0][0] +
			UrbanIndByHIV[0][1]));
		PrevPregUrban.out[CurrSim - 1][CurrYear - 1990] = (UrbanIndByHIV[1][1] / (UrbanIndByHIV[1][0] +
			UrbanIndByHIV[1][1]));
		PrevPregB.out[CurrSim - 1][CurrYear - 1990] = (RaceHIV[0][1] / (RaceHIV[0][0] + RaceHIV[0][1])) *
			0.977 + 1.0 - 0.977;
		PrevPregC.out[CurrSim - 1][CurrYear - 1990] = (RaceHIV[1][1] / (RaceHIV[1][0] + RaceHIV[1][1])) *
			0.977 + 1.0 - 0.977;
		PrevPregW.out[CurrSim - 1][CurrYear - 1990] = (RaceHIV[2][1] / (RaceHIV[2][0] + RaceHIV[2][1])) *
			0.977 + 1.0 - 0.977;
	}
}

void Pop::GetHIVprevByEdu()
{
	int ii, ic;
	double Total[5], Positive[5], temp1, temp2;

	// Antenatal prevalence by education
	for (ii = 0; ii < 4; ii++){
		Total[ii] = 0.0;
		Positive[ii] = 0.0;
	}
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AgeGroup>2 && Register[ic].AgeGroup < 10 &&
			Register[ic].AliveInd == 1 && Register[ic].SexInd == 1 &&
			Register[ic].DOLB>(1.0*CurrYear+0.5)){
			if (Register[ic].HighestGrade == 0){ii = 0;}
			else if (Register[ic].HighestGrade <8){ ii = 1; }
			else if (Register[ic].HighestGrade <13){ ii = 2; }
			else { ii = 3; }
			Total[ii] += Register[ic].PopWeight;
			if (Register[ic].HIVstage>0){ Positive[ii] += Register[ic].PopWeight; }
		}
	}
	PrevPregNoEdu.out[CurrSim - 1][CurrYear - 1990] = Positive[0] / Total[0];
	PrevPregPrimary.out[CurrSim - 1][CurrYear - 1990] = Positive[1] / Total[1];
	PrevPreg2ndary.out[CurrSim - 1][CurrYear - 1990] = Positive[2] / Total[2];
	PrevPregTertiary.out[CurrSim - 1][CurrYear - 1990] = Positive[3] / Total[3];
	temp1 = 0.0;
	temp2 = 0.0;
	for (ii = 0; ii <4; ii++){
		temp1 += Total[ii];
		temp2 += Positive[ii];
	}
	PrevPregAllEdu.out[CurrSim - 1][CurrYear - 1990] = temp2 / temp1;
	
	// Household prevalence by education
	if (CurrYear == 2002){
		for (ii = 0; ii < 5; ii++){
			Total[ii] = 0.0;
			Positive[ii] = 0.0;
		}
		for (ic = 0; ic < Register.size(); ic++){
			if (Register[ic].AgeGroup>2 && Register[ic].AliveInd == 1){
				if (Register[ic].HighestGrade == 0){ ii = 0; }
				else if (Register[ic].HighestGrade <8){ ii = 1; }
				else if (Register[ic].HighestGrade <12){ ii = 2; }
				else if (Register[ic].HighestGrade ==12){ ii = 3; }
				else { ii = 4; }
				Total[ii] += Register[ic].PopWeight;
				if (Register[ic].HIVstage>0){ Positive[ii] += Register[ic].PopWeight; }
			}
		}
		for (ii = 0; ii < 5; ii++){
			PrevByEdu2002.out[CurrSim - 1][ii] = Positive[ii] / Total[ii];}
	}
}

void Pop::GetCD4profile()
{
	int ia, ic, is;
	double CD4stage[4], CD4stageDiag[4], total, totalage[5], diagnosed, treated;
	double totdur, totdurage[5];

	// Get CD4 profile
	for (is = 0; is<4; is++){ 
		CD4stage[is] = 0.0; 
		CD4stageDiag[is] = 0.0;
	}
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AgeGroup>2 && Register[ic].AliveInd == 1 &&
			Register[ic].HIVstage>0 && Register[ic].HIVstage<5){
			if (Register[ic].HIVstage == 4){ CD4stage[3] += Register[ic].PopWeight; }
			else if (Register[ic].HIVstage == 3){ CD4stage[2] += Register[ic].PopWeight; }
			else if (Register[ic].CD4 <500){ CD4stage[1] += Register[ic].PopWeight; }
			else { CD4stage[0] += Register[ic].PopWeight; }
			if (Register[ic].VCThistory == 2){
				if (Register[ic].HIVstage == 4){ CD4stageDiag[3] += Register[ic].PopWeight; }
				else if (Register[ic].HIVstage == 3){ CD4stageDiag[2] += Register[ic].PopWeight; }
				else if (Register[ic].CD4 <500){ CD4stageDiag[1] += Register[ic].PopWeight; }
				else { CD4stageDiag[0] += Register[ic].PopWeight; }
			}
		}
	}
	CD4500plus.out[CurrSim - 1][CurrYear - StartYear] = CD4stage[0];
	CD4350to499.out[CurrSim - 1][CurrYear - StartYear] = CD4stage[1];
	CD4200to349.out[CurrSim - 1][CurrYear - StartYear] = CD4stage[2];
	CD4under200.out[CurrSim - 1][CurrYear - StartYear] = CD4stage[3];
	DiagNoART500plus.out[CurrSim - 1][CurrYear - StartYear] = CD4stageDiag[0];
	DiagNoART500.out[CurrSim - 1][CurrYear - StartYear] = CD4stageDiag[1];
	DiagNoART350.out[CurrSim - 1][CurrYear - StartYear] = CD4stageDiag[2];
	DiagNoART200.out[CurrSim - 1][CurrYear - StartYear] = CD4stageDiag[3];

	// Get profile by level of care
	total = 0.0;
	diagnosed = 0.0;
	treated = 0.0;
	for (is = 0; is<4; is++){ CD4stage[is] = 0.0; }
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AgeGroup>2 && Register[ic].AliveInd == 1 &&
			Register[ic].HIVstage>0){
			total += Register[ic].PopWeight;
			if (Register[ic].VCThistory == 2){ diagnosed += Register[ic].PopWeight; }
			if (Register[ic].HIVstage == 5){ 
				treated += Register[ic].PopWeight;
				if (Register[ic].BaselineCD4 < 200){ CD4stage[3] += Register[ic].PopWeight; }
				else if (Register[ic].BaselineCD4 < 350){ CD4stage[2] += Register[ic].PopWeight; }
				else if (Register[ic].BaselineCD4 < 500){ CD4stage[1] += Register[ic].PopWeight; }
				else { CD4stage[0] += Register[ic].PopWeight; }
			}
		}
	}
	UndiagnosedAdults.out[CurrSim - 1][CurrYear - StartYear] = total - diagnosed;
	DiagnosedUntreated.out[CurrSim - 1][CurrYear - StartYear] = diagnosed - treated;
	TreatedAdults.out[CurrSim - 1][CurrYear - StartYear] = treated;
	AdultARTcoverage.out[CurrSim - 1][CurrYear - StartYear] = treated / total;
	TreatedAdults500plus.out[CurrSim - 1][CurrYear - StartYear] = CD4stage[0];
	TreatedAdults500.out[CurrSim - 1][CurrYear - StartYear] = CD4stage[1];
	TreatedAdults350.out[CurrSim - 1][CurrYear - StartYear] = CD4stage[2];
	TreatedAdults200.out[CurrSim - 1][CurrYear - StartYear] = CD4stage[3];

	// New outputs for Modelling Consortium HIV testing project
	double DiagSex[2], PosSex[2], DiagAge[3], PosAge[3], DiagMSM, PosMSM, DiagFSW, PosFSW;
	for (is = 0; is < 2; is++){
		DiagSex[is] = 0.0;
		PosSex[is] = 0.0;
	}
	for (is = 0; is < 3; is++){
		DiagAge[is] = 0.0;
		PosAge[is] = 0.0;
	}
	DiagMSM = 0.0; 
	PosMSM = 0.0; 
	DiagFSW = 0.0; 
	PosFSW = 0.0;
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AgeGroup>2 && Register[ic].AliveInd == 1 &&
			Register[ic].HIVstage>0){
			if (Register[ic].SexInd == 0){ PosSex[0] += Register[ic].PopWeight; }
			else{ PosSex[1] += Register[ic].PopWeight; }
			if (Register[ic].AgeGroup < 5){ PosAge[0] += Register[ic].PopWeight; }
			else if (Register[ic].AgeGroup < 10){ PosAge[1] += Register[ic].PopWeight; }
			else{ PosAge[2] += Register[ic].PopWeight; }
			if (Register[ic].FSWind == 1){ PosFSW += Register[ic].PopWeight; }
			if (Register[ic].RecentMSM > 0){ PosMSM += Register[ic].PopWeight; }
			if (Register[ic].VCThistory == 2){ 
				if (Register[ic].SexInd == 0){ DiagSex[0] += Register[ic].PopWeight; }
				else{ DiagSex[1] += Register[ic].PopWeight; }
				if (Register[ic].AgeGroup < 5){ DiagAge[0] += Register[ic].PopWeight; }
				else if (Register[ic].AgeGroup < 10){ DiagAge[1] += Register[ic].PopWeight; }
				else{ DiagAge[2] += Register[ic].PopWeight; }
				if (Register[ic].FSWind == 1){ DiagFSW += Register[ic].PopWeight; }
				if (Register[ic].RecentMSM > 0){ DiagMSM += Register[ic].PopWeight; }
			}
		}
	}
	MaleDiagnosed.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * DiagSex[0] / PosSex[0];
	FemDiagnosed.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * DiagSex[1] / PosSex[1];
	Diagnosed15to24.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * DiagAge[0] / PosAge[0];
	Diagnosed25to49.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * DiagAge[1] / PosAge[1];
	Diagnosed50plus.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * DiagAge[2] / PosAge[2];
	MSMdiagnosed.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * DiagMSM / PosMSM;
	FSWdiagnosed.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * DiagFSW / PosFSW;

	// Average duration of HIV infection in pregnant women
	total = 0.0;
	totdur = 0.0;
	for (ia = 0; ia < 5; ia++){
		totalage[ia] = 0.0;
		totdurage[ia] = 0.0;
	}
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].HIVstage>0 &&
			Register[ic].AgeGroup>2 && Register[ic].AgeGroup < 10 &&
			Register[ic].SexInd == 1 && Register[ic].VirginInd==0){ //&& Register[ic].DOLB>(1.0*CurrYear + 0.5)
			total += Register[ic].PopWeight;
			totdur += (1.0 * CurrYear + 0.5 - Register[ic].DateInfect) * Register[ic].PopWeight;
			ia = Register[ic].AgeGroup - 3;
			if (ia >= 0 && ia < 5){
				totalage[ia] += Register[ic].PopWeight;
				//totdurage[ia] += 1.0 * CurrYear + 0.5 - Register[ic].DateInfect;
				if (Register[ic].DateInfect > 1.0 * CurrYear - 0.5){ totdurage[ia] += Register[ic].PopWeight; }
			}
		}
	}
	AveHIVdurPreg.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * totdur / total;
	AveHIVdurPreg15.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * totdurage[0] / totalage[0];
	AveHIVdurPreg20.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * totdurage[1] / totalage[1];
	AveHIVdurPreg25.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * totdurage[2] / totalage[2];
	AveHIVdurPreg30.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * totdurage[3] / totalage[3];
	AveHIVdurPreg35.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * totdurage[4] / totalage[4];
}

void Pop::GetInitSTIprev()
{
	// Returns prevalence in 15-49 year old males and females
	
	int TotPop[2], TotHSV[2], TotTP[2], TotHD[2], TotNG[2];
	int TotCT[2], TotTV[2], TotBV, TotVC;
	int ic, ig;

	for(ig=0; ig<2; ig++){
		TotPop[ig] = 0;
		TotHSV[ig] = 0;
		TotTP[ig] = 0;
		TotHD[ig] = 0;
		TotNG[ig] = 0;
		TotCT[ig] = 0;
		TotTV[ig] = 0;
		TotBV = 0;
		TotVC = 0;
	}
	for(ic=0; ic<InitPop; ic++){
		if(Register[ic].AgeGroup>2 && Register[ic].AgeGroup<10){
			ig = Register[ic].SexInd;
			TotPop[ig] += 1;
			if(Register[ic].HSVstage>0){TotHSV[ig] += 1;}
			if(Register[ic].TPstage>0 && Register[ic].TPstage<5){TotTP[ig] += 1;}
			if(Register[ic].HDstage>0 && Register[ic].HDstage<3){TotHD[ig] += 1;}
			if(Register[ic].NGstage>0 && Register[ic].NGstage<3){TotNG[ig] += 1;}
			if(Register[ic].CTstage>0 && Register[ic].CTstage<3){TotCT[ig] += 1;}
			if(Register[ic].TVstage>0 && Register[ic].TVstage<3){TotTV[ig] += 1;}
			if(Register[ic].BVstage>1 && ig==1){TotBV += 1;}
			if(Register[ic].VCstage>0 && ig==1){TotVC += 1;}
		}
	}
	cout<<"Male HSV-2 prevalence, 15-49: "<<TotHSV[0]<<endl;
	cout<<"Female HSV-2 prevalence, 15-49: "<<TotHSV[1]<<endl;
	cout<<"Male TP prevalence, 15-49: "<<TotTP[0]<<endl;
	cout<<"Female TP prevalence, 15-49: "<<TotTP[1]<<endl;
	cout<<"Male HD prevalence, 15-49: "<<TotHD[0]<<endl;
	cout<<"Female HD prevalence, 15-49: "<<TotHD[1]<<endl;
	cout<<"Male NG prevalence, 15-49: "<<TotNG[0]<<endl;
	cout<<"Female NG prevalence, 15-49: "<<TotNG[1]<<endl;
	cout<<"Male CT prevalence, 15-49: "<<TotCT[0]<<endl;
	cout<<"Female CT prevalence, 15-49: "<<TotCT[1]<<endl;
	cout<<"Male TV prevalence, 15-49: "<<TotTV[0]<<endl;
	cout<<"Female TV prevalence, 15-49: "<<TotTV[1]<<endl;
	cout<<"Female BV prevalence, 15-49: "<<TotBV<<endl;
	cout<<"Female VC prevalence, 15-49: "<<TotVC<<endl;
	cout<<"Male pop 15-49: "<<TotPop[0]<<endl;
	cout<<"Female pop 15-49: "<<TotPop[1]<<endl;
}

void Pop::GetCurrSTIprev()
{
	// Calculates prevalence in 15-49 year old males and females
	// Very similar to GetInitSTIprev function.
	
	if (FixedUncertainty == 1){
		double TotPop[3], TotHSV[3], TotTP[2], TotHD[2], TotNG[2];
		double TotCT[2], TotTV[2], TotBV, TotVC, TotHIV[3];
		int ic, ig, ir;

		for (ig = 0; ig < 2; ig++){
			TotPop[ig] = 0.0;
			TotHIV[ig] = 0.0;
			TotHSV[ig] = 0.0;
			TotTP[ig] = 0.0;
			TotHD[ig] = 0.0;
			TotNG[ig] = 0.0;
			TotCT[ig] = 0.0;
			TotTV[ig] = 0.0;
			TotBV = 0.0;
			TotVC = 0.0;
		}
		for (ic = 0; ic<Register.size(); ic++){
			if (Register[ic].AgeGroup>2 && Register[ic].AgeGroup<10 && Register[ic].AliveInd == 1){
				ig = Register[ic].SexInd;
				TotPop[ig] += Register[ic].PopWeight;
				//if (ig==0 && Register[ic].AgeGroup>5){ TotPop[ig] += 1; }
				//else{ if (ig==1 && Register[ic].AgeGroup>4){ TotPop[ig] += 1; } }
				if (Register[ic].HIVstage>0){ TotHIV[ig] += Register[ic].PopWeight; }
				if (Register[ic].HSVstage > 0){ TotHSV[ig] += Register[ic].PopWeight; }
				//if (Register[ic].HSVstage>0){
				//	if (ig==0 && Register[ic].AgeGroup>5){ TotHSV[ig] += 1; }
				//	else{ if (ig==1 && Register[ic].AgeGroup>4){ TotHSV[ig] += 1; } }
				//}
				if (Register[ic].TPstage > 0 && Register[ic].TPstage<5){ TotTP[ig] += Register[ic].PopWeight; }
				//if(Register[ic].TPstage>1){ TotTP[ig] += 1; } // TP seroprevalence
				if (Register[ic].HDstage>0 && Register[ic].HDstage<3){ TotHD[ig] += Register[ic].PopWeight; }
				if (Register[ic].NGstage>0 && Register[ic].NGstage<3){ TotNG[ig] += Register[ic].PopWeight; }
				if (Register[ic].CTstage>0 && Register[ic].CTstage<3){ TotCT[ig] += Register[ic].PopWeight; }
				if (Register[ic].TVstage>0 && Register[ic].TVstage<3){ TotTV[ig] += Register[ic].PopWeight; }
				if (Register[ic].BVstage>1 && ig == 1){ TotBV += Register[ic].PopWeight; }
				if (Register[ic].VCstage > 0 && ig == 1){ TotVC += Register[ic].PopWeight; }
			}
		}
		HIVprev15to49M.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHIV[0] / TotPop[0];
		HIVprev15to49F.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHIV[1] / TotPop[1];
		HIVprev15to49.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * (TotHIV[0] + TotHIV[1]) / 
			(TotPop[0] + TotPop[1]);
		HSVprev15to49M.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHSV[0] / TotPop[0];
		HSVprev15to49F.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHSV[1] / TotPop[1];
		TPprev15to49M.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotTP[0] / TotPop[0];
		TPprev15to49F.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotTP[1] / TotPop[1];
		HDprev15to49M.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHD[0] / TotPop[0];
		HDprev15to49F.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHD[1] / TotPop[1];
		NGprev15to49M.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotNG[0] / TotPop[0];
		NGprev15to49F.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotNG[1] / TotPop[1];
		CTprev15to49M.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotCT[0] / TotPop[0];
		CTprev15to49F.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotCT[1] / TotPop[1];
		TVprev15to49M.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotTV[0] / TotPop[0];
		TVprev15to49F.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotTV[1] / TotPop[1];
		BVprev15to49F.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotBV / TotPop[1];
		VCprev15to49F.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotVC / TotPop[1];

		// STI prevalence by race
		for (ir = 0; ir<3; ir++){
			TotPop[ir] = 0.0;
			TotHIV[ir] = 0.0;
			TotHSV[ir] = 0.0;
		}
		for (ic = 0; ic<Register.size(); ic++){
			if (Register[ic].AgeGroup>2 && Register[ic].AgeGroup<10 && Register[ic].AliveInd == 1){
				ir = Register[ic].Race;
				TotPop[ir] += Register[ic].PopWeight;
				if (Register[ic].HIVstage>0){ TotHIV[ir] += Register[ic].PopWeight; }
				if (Register[ic].HSVstage>0){ TotHSV[ir] += Register[ic].PopWeight; }
			}
		}
		HIVprev15to49B.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHIV[0] / TotPop[0];
		HIVprev15to49C.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHIV[1] / TotPop[1];
		HIVprev15to49W.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHIV[2] / TotPop[2];
		HSVprev15to49B.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHSV[0] / TotPop[0];
		HSVprev15to49C.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHSV[1] / TotPop[1];
		HSVprev15to49W.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHSV[2] / TotPop[2];

		// HIV prevalence by urban-rural location
		for (ig = 0; ig < 2; ig++){
			TotPop[ig] = 0.0;
			TotHIV[ig] = 0.0;
		}
		for (ic = 0; ic<Register.size(); ic++){
			if (Register[ic].AgeGroup>2 && Register[ic].AgeGroup < 10 && Register[ic].AliveInd == 1){
				if (Register[ic].CurrUrban == 0){
					TotPop[0] += Register[ic].PopWeight;
					if (Register[ic].HIVstage > 0){ TotHIV[0] += Register[ic].PopWeight; }
				}
				else{
					TotPop[1] += Register[ic].PopWeight;
					if (Register[ic].HIVstage > 0){ TotHIV[1] += Register[ic].PopWeight; }
				}
			}
		}
		RuralHIV.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHIV[0] / TotPop[0];
		UrbanHIV.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHIV[1] / TotPop[1];
		TotPop15to49_R.out[CurrSim - 1][CurrYear - StartYear] = TotPop[0];
		TotPop15to49_U.out[CurrSim - 1][CurrYear - StartYear] = TotPop[1];
	}

	// Calculate the STI prevalence required for likelihood calculation purposes
	// (I haven't added code for HD, BV or VC as yet.)
	if (HIVcalib == 1 || (HIVind == 1 && FixedUncertainty == 1)){
		GetCSWprev(&HIVtransitionF, 0);
	}
	if (HSVcalib == 1){
		GetANCprev(&HSVtransitionF, 1);
		GetHHprev(&HSVtransitionF, 1);
		GetHHprev(&HSVtransitionM, 1);
		GetFPCprev(&HSVtransitionF, 1);
		GetCSWprev(&HSVtransitionF, 1);
	}
	if (TPcalib == 1){
		GetANCprev(&TPtransitionF, 2);
		GetHHprev(&TPtransitionF, 2);
		GetHHprev(&TPtransitionM, 2);
		GetFPCprev(&TPtransitionF, 2);
		GetCSWprev(&TPtransitionF, 2);
	}
	if (NGcalib == 1){
		GetANCprev(&NGtransitionF, 4);
		GetHHprev(&NGtransitionF, 4);
		GetHHprev(&NGtransitionM, 4);
		GetFPCprev(&NGtransitionF, 4);
		GetCSWprev(&NGtransitionF, 4);
	}
	if (CTcalib == 1){
		GetANCprev(&CTtransitionF, 5);
		GetHHprev(&CTtransitionF, 5);
		GetHHprev(&CTtransitionM, 5);
		GetFPCprev(&CTtransitionF, 5);
		GetCSWprev(&CTtransitionF, 5);
	}
	if (TVcalib == 1){
		GetANCprev(&TVtransitionF, 6);
		GetHHprev(&TVtransitionF, 6);
		GetHHprev(&TVtransitionM, 6);
		GetFPCprev(&TVtransitionF, 6);
		GetCSWprev(&TVtransitionF, 6);
	}
}

double Pop::GetQstatistic(int Cmatrix[3][3], int MatDim)
{
	int ir, ic;
	double RowTotal, trace, Qstatistic;

	trace = 0.0;
	for(ir=0; ir<MatDim; ir++){
		RowTotal = 0.0;
		for(ic=0; ic<MatDim; ic++){
			RowTotal += Cmatrix[ir][ic];}
		trace += Cmatrix[ir][ir]/RowTotal;
	}

	Qstatistic = (trace - 1.0)/(MatDim - 1.0);

	return Qstatistic;
}

void Pop::GetSTIconcordance()
{
	int ic, pID, pID2, iSTI, pSTI;
	int HIVmatrix[3][3], HSVmatrix[3][3], TPmatrix[3][3], HDmatrix[3][3];
	int NGmatrix[3][3], CTmatrix[3][3], TVmatrix[3][3];

	for(ic=0; ic<Register.size(); ic++){
		if(Register[ic].AgeGroup>2 && Register[ic].AliveInd==1 &&
			Register[ic].CurrPartners>0 && Register[ic].SexInd==1){
				pID = Register[ic].IDprimary;
				pID2 = Register[ic].ID2ndary;
				if(HIVind==1){
					if(Register[ic].HIVstage>0){iSTI=1;}
					else{iSTI=0;}
					if(Register[pID-1].HIVstage>0){pSTI=1;}
					else{pSTI=0;}
					HIVmatrix[iSTI][pSTI] += 1;
					if(Register[ic].CurrPartners>1){
						if(Register[pID2-1].HIVstage>0){pSTI=1;}
						else{pSTI=0;}
						HIVmatrix[iSTI][pSTI] += 1;
					}
				}
				if(HSVind==1){
					if(Register[ic].HSVstage>0){iSTI=1;}
					else{iSTI=0;}
					if(Register[pID-1].HSVstage>0){pSTI=1;}
					else{pSTI=0;}
					HSVmatrix[iSTI][pSTI] += 1;
					if(Register[ic].CurrPartners>1){
						if(Register[pID2-1].HSVstage>0){pSTI=1;}
						else{pSTI=0;}
						HSVmatrix[iSTI][pSTI] += 1;
					}
				}
				if(TPind==1){
					if(Register[ic].TPstage>0){iSTI=1;}
					else{iSTI=0;}
					if(Register[ic].TPstage>4){iSTI=0;} // immune
					if(Register[pID-1].TPstage>0){pSTI=1;}
					else{pSTI=0;}
					if(Register[pID-1].TPstage>4){pSTI=0;}
					TPmatrix[iSTI][pSTI] += 1;
					if(Register[ic].CurrPartners>1){
						if(Register[pID2-1].TPstage>0){pSTI=1;}
						else{pSTI=0;}
						if(Register[pID2-1].TPstage>4){pSTI=0;}
						TPmatrix[iSTI][pSTI] += 1;
					}
				}
				if(HDind==1){
					if(Register[ic].HDstage>0){iSTI=1;}
					else{iSTI=0;}
					if(Register[ic].HDstage>2){iSTI=0;} // immune
					if(Register[pID-1].HDstage>0){pSTI=1;}
					else{pSTI=0;}
					if(Register[pID-1].HDstage>2){pSTI=0;}
					HDmatrix[iSTI][pSTI] += 1;
					if(Register[ic].CurrPartners>1){
						if(Register[pID2-1].HDstage>0){pSTI=1;}
						else{pSTI=0;}
						if(Register[pID2-1].HDstage>2){pSTI=0;}
						HDmatrix[iSTI][pSTI] += 1;
					}
				}
				if(NGind==1){
					if(Register[ic].NGstage>0){iSTI=1;}
					else{iSTI=0;}
					if(Register[ic].NGstage>2){iSTI=0;} // immune
					if(Register[pID-1].NGstage>0){pSTI=1;}
					else{pSTI=0;}
					if(Register[pID-1].NGstage>2){pSTI=0;}
					NGmatrix[iSTI][pSTI] += 1;
					if(Register[ic].CurrPartners>1){
						if(Register[pID2-1].NGstage>0){pSTI=1;}
						else{pSTI=0;}
						if(Register[pID2-1].NGstage>2){pSTI=0;}
						NGmatrix[iSTI][pSTI] += 1;
					}
				}
				if(CTind==1){
					if(Register[ic].CTstage>0){iSTI=1;}
					else{iSTI=0;}
					if(Register[ic].CTstage>2){iSTI=0;} // immune
					if(Register[pID-1].CTstage>0){pSTI=1;}
					else{pSTI=0;}
					if(Register[pID-1].CTstage>2){pSTI=0;}
					CTmatrix[iSTI][pSTI] += 1;
					if(Register[ic].CurrPartners>1){
						if(Register[pID2-1].CTstage>0){pSTI=1;}
						else{pSTI=0;}
						if(Register[pID2-1].CTstage>2){pSTI=0;}
						CTmatrix[iSTI][pSTI] += 1;
					}
				}
				if(TVind==1){
					if(Register[ic].TVstage>0){iSTI=1;}
					else{iSTI=0;}
					if(Register[ic].TVstage>2){iSTI=0;} // immune
					if(Register[pID-1].TVstage>0){pSTI=1;}
					else{pSTI=0;}
					if(Register[pID-1].TVstage>2){pSTI=0;}
					TVmatrix[iSTI][pSTI] += 1;
					if(Register[ic].CurrPartners>1){
						if(Register[pID2-1].TVstage>0){pSTI=1;}
						else{pSTI=0;}
						if(Register[pID2-1].TVstage>2){pSTI=0;}
						TVmatrix[iSTI][pSTI] += 1;
					}
				}
			}
	}

	/*if(HIVind==1){
		HIVconcordance.out[CurrSim-1][CurrYear-StartYear] = GetQstatistic(HIVmatrix, 2);}
	if(HSVind==1){
		HSVconcordance.out[CurrSim-1][CurrYear-StartYear] = GetQstatistic(HSVmatrix, 2);}
	if(TPind==1){
		TPconcordance.out[CurrSim-1][CurrYear-StartYear] = GetQstatistic(TPmatrix, 2);}
	if(HDind==1){
		HDconcordance.out[CurrSim-1][CurrYear-StartYear] = GetQstatistic(HDmatrix, 2);}
	if(NGind==1){
		NGconcordance.out[CurrSim-1][CurrYear-StartYear] = GetQstatistic(NGmatrix, 2);}
	if(CTind==1){
		CTconcordance.out[CurrSim-1][CurrYear-StartYear] = GetQstatistic(CTmatrix, 2);}
	if(TVind==1){
		TVconcordance.out[CurrSim-1][CurrYear-StartYear] = GetQstatistic(TVmatrix, 2);}*/
}

void Pop::SavePopPyramid(const char* filout)
{
	int ia;
	ofstream file(filout);

	for(ia=0; ia<18; ia++){
		file<<right<<PopPyramid[ia][0]<<"	"<<PopPyramid[ia][1]<<endl;}
	file.close();
}

void Pop::SaveBehavPyramid(const char* filout)
{
	int ia, is;
	ofstream file(filout);

	for(ia=0; ia<18; ia++){
		for(is=0; is<6; is++){
			file<<right<<BehavPyramid[ia][is]<<"	";}
		file<<endl;
	}
	file.close();
}

void Pop::SavePrefMatrix(const char* filout)
{
	int ia, is;
	ofstream file(filout);

	for(ia=0; ia<16; ia++){
		for(is=0; is<16; is++){
			file<<right<<PrefMatrix[ia][is]<<"	";}
		file<<endl;
	}
	file.close();
}

void Pop::SaveNumberPartners(const char* filout)
{
	int ia, is;
	ofstream file(filout);

	for(ia=0; ia<16; ia++){
		for(is=0; is<10; is++){
			file<<right<<NumberPartners[ia][is]<<"	";}
		file<<endl;
	}
	file.close();
}

void Pop::SaveFSWcontacts2016(const char* filout)
{
	int ia;
	ofstream file(filout);

	for (ia = 0; ia<16; ia++){
		file << right << FSWcontacts2016[ia] << endl;
	}
	file.close();
}

void Pop::SaveLifetimePartners(const char* filout)
{
	int ic;
	double AgeNow;
	ofstream file(filout);

	for(ic=0; ic<Register.size(); ic++){
		if(Register[ic].AliveInd==1 && Register[ic].SexInd==0 && Register[ic].CurrPartners>0 &&
			Register[ic].AgeGroup>=3 && Register[ic].AgeGroup<10){
				AgeNow = 0.5 + CurrYear - Register[ic].DOB;
				file<<right<<ic<<"	"<<AgeNow<<"	"<<Register[ic].LifetimePartners<<endl;
			}
	}
	file.close();
}

void Pop::SaveTotSex(const char* filout)
{
	int ia, is;
	ofstream file(filout);

	for(ia=0; ia<16; ia++){
		for(is=0; is<8; is++){
			file<<right<<TotSex[ia][is]<<"	";}
		file<<endl;
	}
	file.close();
}

void Pop::SaveMarriageRates(const char* filout)
{
	/*int ia, is;
	ofstream file(filout);

	for(ia=0; ia<81; ia++){
		for(is=0; is<4; is++){
			file<<right<<AgeEffectMarriage[ia][is]<<"	";}
		file<<endl;
	}
	file.close();*/
}

void Pop::SaveAdultHIVstage(const char* filout)
{
	int iy, is;
	ofstream file(filout);

	for(iy=0; iy<=40; iy++){
		for(is=0; is<12; is++){
			file<<right<<AdultHIVstageTrend[iy][is]<<"	";}
		file<<endl;
	}
	file.close();
}

void Pop::GetHIVsurvivalTimes()
{
	int ic, count;

	count = 0;
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].HIVstage>0 && 
			(Register[ic].DateInfect - Register[ic].DOB > 15.0)){
			HIVsurvivalTimes.out[count][0] = ic + 1; //ID
			HIVsurvivalTimes.out[count][1] = Register[ic].SexInd; // Gender
			HIVsurvivalTimes.out[count][2] = Register[ic].DateInfect -
				Register[ic].DOB; // Age at seroconversion
			HIVsurvivalTimes.out[count][3] = Register[ic].SPVL;
			if (Register[ic].AliveInd == 1){
				HIVsurvivalTimes.out[count][4] = 0; // Dead indicator
				HIVsurvivalTimes.out[count][5] = 0.5 + CurrYear -
					Register[ic].DateInfect; // Assuming outputs calculated at start of yr
			}
			else{
				HIVsurvivalTimes.out[count][4] = 1; 
				HIVsurvivalTimes.out[count][5] = Register[ic].DOD -
					Register[ic].DateInfect; 
			}
			count += 1;
		}
	}
	HIVsurvivalTimes.TotIndivs = count;

	HIVsurvivalTimes.RecordSample("HIVsurvivalTimes.txt");
}

void Pop::GetVLs()
{
	int ic, count;

	count = 0;
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].HIVstage>0 && Register[ic].AliveInd==1 &&
			(Register[ic].DateInfect - Register[ic].DOB > 15.0)){
			ViralLoads.out[count][0] = ic + 1; //ID
			ViralLoads.out[count][1] = Register[ic].CD4; // CD4
			ViralLoads.out[count][2] = Register[ic].logVL; // VL
			count += 1;
		}
	}
	ViralLoads.TotIndivs = count;

	ViralLoads.RecordSample("ViralLoads.txt");
}

void Pop::GetUnprotected()
{
	int ic, count;

	if (CurrSim == 1){ count = 0; }
	else{ count = UnprotectedRisk.TotIndivs; }
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].CurrPartners>0 &&
			Register[ic].CurrAge>=15){
			UnprotectedRisk.out[count][0] = ic + 1; //ID
			UnprotectedRisk.out[count][1] = Register[ic].SexInd; // Sex
			UnprotectedRisk.out[count][2] = Register[ic].CurrAge; // Age
			UnprotectedRisk.out[count][3] = Register[ic].MarriedInd; 
			UnprotectedRisk.out[count][4] = Register[ic].CurrPartners; 
			UnprotectedRisk.out[count][5] = 0; 
			if (Register[ic].HIVstage>0){ UnprotectedRisk.out[count][5] = 1; }
			if (Register[ic].VCThistory==2){ UnprotectedRisk.out[count][5] = 2; }
			if (Register[ic].CondomPrimary<1.0 || (Register[ic].CurrPartners==2 &&
				Register[ic].Condom2ndary<1.0)){ 
				UnprotectedRisk.out[count][6] = 1; 
			}
			else{ UnprotectedRisk.out[count][6] = 0; }
			count += 1;
		}
	}
	UnprotectedRisk.TotIndivs = count;

	if (CurrSim == 3){ UnprotectedRisk.RecordSample("UnprotectedRisk.txt"); }
}

void Pop::GetMigHIV()
{
	int ic, count;

	for (ic = 0; ic < Register.size(); ic++){ r2[ic] = rg.Random(); }

	// count = 0;
	count = MigrationHIVassn.TotIndivs;
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].CurrAge >= 15 &&
			((Register[ic].CurrUrban == 0 && r2[ic]<0.1) || Register[ic].DateMig >= 1.0*CurrYear - 1.5)){
			MigrationHIVassn.out[count][0] = ic + 1; //ID
			MigrationHIVassn.out[count][1] = Register[ic].SexInd; // Sex
			MigrationHIVassn.out[count][2] = Register[ic].CurrAge; // Age
			MigrationHIVassn.out[count][3] = Register[ic].MarriedInd;
			MigrationHIVassn.out[count][4] = Register[ic].CurrPartners;
			MigrationHIVassn.out[count][5] = 0;
			if (Register[ic].HIVstage>0){ MigrationHIVassn.out[count][5] = 1; }
			if (Register[ic].DateMig >= 1.0*CurrYear - 1.5 && Register[ic].DateMig < 1.0*CurrYear + 0.5){
				MigrationHIVassn.out[count][6] = 1; }
			else{ MigrationHIVassn.out[count][6] = 0; }
			MigrationHIVassn.out[count][7] = Register[ic].Race; // Race
			MigrationHIVassn.out[count][8] = CurrYear; // Year
			count += 1;
		}
	}
	MigrationHIVassn.TotIndivs = count;

	if (CurrSim == 5 && CurrYear==2008){ MigrationHIVassn.RecordSample("MigrationHIVassn.txt"); }
}

void Pop::GetPrisons()
{
	int ic, ir;
	double TotM, TotPrison, TotYouth, TotHIV, TotSchool, TotPrev, TotSentenced;

	TotM = 0.0;
	TotPrison = 0.0;
	TotYouth = 0.0;
	TotHIV = 0.0;
	TotSchool = 0.0;
	TotPrev = 0.0;
	TotSentenced = 0.0;
	if (CurrYear == 2005){
		for (ir = 0; ir < 3; ir++){
			PrisonRace.out[CurrSim - 1][ir] = 0.0; }
	}

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].CurrAge >= 15 &&
			Register[ic].SexInd==0){ 
			TotM += Register[ic].PopWeight;
			if (Register[ic].Imprisoned>0){
				TotPrison += Register[ic].PopWeight;
				if (Register[ic].CurrAge < 26){ TotYouth += Register[ic].PopWeight; }
				if (Register[ic].HighestGrade >= 12){ TotSchool += Register[ic].PopWeight; }
				if (Register[ic].HIVstage > 0){ TotHIV += Register[ic].PopWeight; }
				if (Register[ic].PrevImprisoned > 0){ TotPrev += Register[ic].PopWeight; }
				if (Register[ic].Imprisoned == 2){ TotSentenced += Register[ic].PopWeight; }
				if (CurrYear == 2007){
					ir = Register[ic].Race;
					PrisonRace.out[CurrSim - 1][ir] += Register[ic].PopWeight;
				}
			}
		}
	}
	
	FractionInPrison.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotPrison / TotM;
	PrisonYouth.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotYouth / TotPrison;
	PrisonCompletedSchool.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotSchool / TotPrison;
	PrisonHIV.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotHIV / TotPrison;
	PrisonPrevious.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotPrev / TotPrison;
	PrisonSentenced.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotSentenced / TotPrison;
}

void Pop::GetEduProfile()
{
	int ic, iy, ir, ig, ii, ia, MID, PID;
	double Alive[3][2], Completed[3][2], Alive20to24[3], Tertiary[3];
	double CompletedYrs[5], TotPop[5], AveParentEdu;
	int Enrolled[14][3], Alive2[14][3];
	int AttainmentAge[8][2], AttainmentGrade[12], TotAge[8], TotGrade[12];
	int BirthCohort, CompletedGrade[12];

	iy = CurrYear - StartYear;

	// Annually-updated outputs
	for (ir = 0; ir < 3; ir++){
		for (ig = 0; ig < 2; ig++){
			Alive[ir][ig] = 0.0;
			Completed[ir][ig] = 0.0;
		}
		Alive20to24[ir] = 0.0;
		Tertiary[ir] = 0.0;
	}
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1){
			if (Register[ic].AgeGroup >= 3 && Register[ic].AgeGroup <= 6){
				ir = Register[ic].Race;
				ig = Register[ic].SexInd;
				if (Register[ic].AgeGroup < 5){
					Alive[ir][ig] += Register[ic].PopWeight;
					if (Register[ic].HighestGrade >= 12){ Completed[ir][ig] += Register[ic].PopWeight; }
				}
				if (Register[ic].AgeGroup == 4){ Alive20to24[ir] += Register[ic].PopWeight; }
				if (Register[ic].HighestGrade == 12 && Register[ic].InSchool == 1){
					Tertiary[ir] += Register[ic].PopWeight;
				}
			}
		}
	}
	CompletedMatricAM.out[CurrSim - 1][iy] = 1.0 * Completed[0][0] / Alive[0][0];
	CompletedMatricAF.out[CurrSim - 1][iy] = 1.0 * Completed[0][1] / Alive[0][1];
	CompletedMatricCM.out[CurrSim - 1][iy] = 1.0 * Completed[1][0] / Alive[1][0];
	CompletedMatricCF.out[CurrSim - 1][iy] = 1.0 * Completed[1][1] / Alive[1][1];
	CompletedMatricWM.out[CurrSim - 1][iy] = 1.0 * Completed[2][0] / Alive[2][0];
	CompletedMatricWF.out[CurrSim - 1][iy] = 1.0 * Completed[2][1] / Alive[2][1];
	CompletedMatric.out[CurrSim - 1][iy] = 1.0 * (Completed[0][0] + Completed[0][1] + Completed[1][0] + 
		Completed[1][1] + Completed[2][0] + Completed[2][1]) / (Alive[0][0] + Alive[0][1] +
		Alive[1][0] + Alive[1][1] + Alive[2][0] + Alive[2][1]);
	TertiaryEnrolRatioA.out[CurrSim - 1][iy] = 1.0 * Tertiary[0] / Alive20to24[0];
	TertiaryEnrolRatioC.out[CurrSim - 1][iy] = 1.0 * Tertiary[1] / Alive20to24[1];
	TertiaryEnrolRatioW.out[CurrSim - 1][iy] = 1.0 * Tertiary[2] / Alive20to24[2];

	// Outputs in 2001
	/*if (CurrYear == 2001){
		for (ia = 0; ia < 14; ia++){
			for (ir = 0; ir < 3; ir++){ 
				Enrolled[ia][ir] = 0; 
				Alive2[ia][ir] = 0;
			}
		}
		for (ic = 0; ic < Register.size(); ic++){
			if (Register[ic].AliveInd == 1 && Register[ic].CurrAge <= 20 &&
				Register[ic].CurrAge >6){
				ia = Register[ic].CurrAge - 7;
				ir = Register[ic].Race;
				Alive2[ia][ir] += 1;
				if (Register[ic].InSchool == 1 && Register[ic].HighestGrade < 12){
					Enrolled[ia][ir] += 1;}
			}
		}
		for (ia = 0; ia < 14; ia++){
			Enrolled2001A.out[CurrSim - 1][ia] = 1.0 * Enrolled[ia][0] / Alive2[ia][0];
			Enrolled2001C.out[CurrSim - 1][ia] = 1.0 * Enrolled[ia][1] / Alive2[ia][1];
			Enrolled2001W.out[CurrSim - 1][ia] = 1.0 * Enrolled[ia][2] / Alive2[ia][2];
		}
	}*/
	/*if(CurrYear == 2001){
		for (ii = 0; ii < 5; ii++){
			CompletedYrs[ii] = 0.0;
			TotPop[ii] = 0.0;
		}
		for (ic = 0; ic < Register.size(); ic++){
			if (Register[ic].AliveInd == 1 && Register[ic].Race == 0 && 
				Register[ic].CurrAge >= 16 && Register[ic].CurrAge <= 20){
				AveParentEdu = 0.0;
				MID = Register[ic].ParentID[1];
				PID = Register[ic].ParentID[0];
				if (MID > 0){
					if (PID > 0){ AveParentEdu = 0.5 * (Register[MID - 1].HighestGrade + Register[PID - 1].HighestGrade); }
					else{ AveParentEdu = 1.0 * Register[MID - 1].HighestGrade; }
				}
				else if (PID > 0){ AveParentEdu = 1.0 * Register[PID - 1].HighestGrade; }
				if (AveParentEdu < 1.0){ ii = 0; }
				else if (AveParentEdu < 7.0){ ii = 1; }
				else if (AveParentEdu < 11.0){ ii = 2; }
				else if (AveParentEdu <= 12.0){ ii = 3; }
				else { ii = 4; }
				CompletedYrs[ii] += Register[ic].HighestGrade;
				TotPop[ii] += 1.0;
			}
		}
		for (ii = 0; ii < 5; ii++){
			if (TotPop[ii]>0.0){ EduByParentEdu.out[CurrSim - 1][ii] = 1.0 * CompletedYrs[ii] / TotPop[ii]; }
		}
	}*/

	// Outputs in 2007
	/*if (CurrYear == 2007){
		for (ia = 0; ia < 8; ia++){
			TotAge[ia] = 0;
			AttainmentAge[ia][0] = 0;
			AttainmentAge[ia][1] = 0;
		}
		BirthCohort = 0;
		for (ii = 0; ii < 12; ii++){CompletedGrade[ii] = 0;}
		for (ic = 0; ic < Register.size(); ic++){
			if (Register[ic].AliveInd == 1 && Register[ic].AgeGroup>4 &&
				Register[ic].AgeGroup <= 12){
				ia = Register[ic].AgeGroup - 5;
				TotAge[ia] += 1;
				if (Register[ic].HighestGrade >= 7){ AttainmentAge[ia][0] += 1; }
				if (Register[ic].HighestGrade >= 12){ AttainmentAge[ia][1] += 1; }
			}
			if (Register[ic].AliveInd == 1 && Register[ic].CurrAge>=28 &&
				Register[ic].CurrAge <= 32){
				BirthCohort += 1;
				if (Register[ic].HighestGrade >= 12){ CompletedGrade[11] += 1; }
				else if (Register[ic].HighestGrade > 0){ 
					ii = Register[ic].HighestGrade - 1;
					CompletedGrade[ii] += 1; 
				}
			}
		}
		for (ia = 0; ia < 8; ia++){
			for (ii = 0; ii < 2; ii++){
				AttainmentByAge2007.out[CurrSim - 1][ii * 8 + ia] = 1.0 * 
					AttainmentAge[ia][ii] / TotAge[ia];
			}
		}
		for (ii = 0; ii < 12; ii++){
			AttainmentByGrade2007.out[CurrSim - 1][ii] = 1.0 *
				CompletedGrade[ii] / BirthCohort;
		}
	}*/
}

void Pop::GetAssortativeness()
{
	int ic, ig, ii, irace, prace, iedu, pedu, PID;
	double TotPop[10], TotMarried[10];

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].MarriedInd == 1 &&
			Register[ic].SexInd == 0 && Register[ic].CurrAge <= 40){
			PID = Register[ic].IDprimary;
			irace = Register[ic].Race;
			prace = Register[PID-1].Race;
			RaceConcordance.out[CurrSim - 1][irace * 3 + prace] += Register[ic].PopWeight;
			if (Register[ic].HighestGrade < 8){ iedu = 0; }
			else if (Register[ic].HighestGrade < 10){ iedu = 1; }
			else if (Register[ic].HighestGrade < 12){ iedu = 2; }
			else if (Register[ic].HighestGrade < 13){ iedu = 3; }
			else { iedu = 4; }
			if (Register[PID - 1].HighestGrade < 8){ pedu = 0; }
			else if (Register[PID - 1].HighestGrade < 10){ pedu = 1; }
			else if (Register[PID - 1].HighestGrade < 12){ pedu = 2; }
			else if (Register[PID - 1].HighestGrade < 13){ pedu = 3; }
			else { pedu = 4; }
			EduConcordance.out[CurrSim - 1][iedu * 5 + pedu] += Register[ic].PopWeight;
		}
	}

	// Calculate levels of marriage by educational attainment
	for (ii = 0; ii < 10; ii++){
		TotPop[ii] = 0.0;
		TotMarried[ii] = 0.0;
	}
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].CurrAge >= 20 &&
			Register[ic].CurrAge < 60){
			ig = Register[ic].SexInd;
			if (Register[ic].HighestGrade < 8){ iedu = 0; }
			else if (Register[ic].HighestGrade < 10){ iedu = 1; }
			else if (Register[ic].HighestGrade < 12){ iedu = 2; }
			else if (Register[ic].HighestGrade < 13){ iedu = 3; }
			else { iedu = 4; }
			TotPop[ig * 5 + iedu] += Register[ic].PopWeight;
			if (Register[ic].MarriedInd == 1){ TotMarried[ig * 5 + iedu] += Register[ic].PopWeight; }
		}
	}
	for (ii = 0; ii < 10; ii++){
		MarriedEdu2011.out[CurrSim - 1][ii] = 1.0 * TotMarried[ii] / TotPop[ii];}
}

void Pop::GetCohabitTemp()
{
	int ic, ir, ig, PID;
	double TempCohabit[3][2], Denominator[3][2], sum1, sum2;

	for (ir = 0; ir < 3; ir++){
		for (ig = 0; ig < 2; ig++){
			TempCohabit[ir][ig] = 0.0;
			Denominator[ir][ig] = 0.0;
		}
	}

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].MarriedInd == 1){
			ir = Register[ic].Race;
			ig = Register[ic].SexInd;
			Denominator[ir][ig] += Register[ic].PopWeight;
			//if (Register[ic].CurrUrban == 1){ Denominator[ir][ig] += 1; }
			PID = Register[ic].IDprimary;
			//if (Register[ic].CurrUrban != Register[PID - 1].CurrUrban && Register[ic].Visiting == 1){
			if (Register[ic].CurrUrban == 1 && Register[PID - 1].CurrUrban == 0){
				TempCohabit[ir][ig] += Register[ic].PopWeight;
			}
		}
	}

	if (CurrYear == 1996){
		for (ir = 0; ir < 3; ir++){
			for (ig = 0; ig < 2; ig++){
				CohabitTemp.out[CurrSim - 1][ir + ig * 3] = 1.0 * TempCohabit[ir][ig] / Denominator[ir][ig];
			}
		}
	}
	sum1 = 0;
	sum2 = 0;
	for (ir = 0; ir < 3; ir++){
		for (ig = 0; ig < 2; ig++){
			sum1 += TempCohabit[ir][ig];
			sum2 += Denominator[ir][ig];
		}
	}
	CohabitTempTrend.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * sum1 / sum2;
}

void Pop::GetConcurrency()
{
	int ic, ir, iy;
	double TotPop[3], Concurrent[3];

	// Prevalence of concurrency in men, by race
	if (CurrYear >= 2000 && CurrYear <= 2010){
		for (ir = 0; ir < 3; ir++){
			TotPop[ir] = 0.0;
			Concurrent[ir] = 0.0;
		}
		for (ic = 0; ic < Register.size(); ic++){
			if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 0 &&
				Register[ic].AgeGroup>2 && Register[ic].AgeGroup < 10){
				ir = Register[ic].Race;
				TotPop[ir] += Register[ic].PopWeight;
				if (Register[ic].CurrPartners > 1){ Concurrent[ir] += Register[ic].PopWeight; }
			}
		}

		iy = CurrYear - 2000;
		ConcurrencyPrev.out[CurrSim - 1][iy] = 1.0 * Concurrent[0] / TotPop[0];
		ConcurrencyPrev.out[CurrSim - 1][iy + 11] = 1.0 * Concurrent[1] / TotPop[1];
		ConcurrencyPrev.out[CurrSim - 1][iy + 22] = 1.0 * Concurrent[2] / TotPop[2];
	}

	// Prevalence of concurrency in men and women aged 15-49
	TotPop[0] = 0.0;
	Concurrent[0] = 0.0;
	iy = CurrYear - StartYear;
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].AgeGroup>2 && Register[ic].AgeGroup < 10){
			TotPop[0] += Register[ic].PopWeight;
			if (Register[ic].CurrPartners > 1){ Concurrent[0] += Register[ic].PopWeight; }
		}
	}
	ConcurrencyTot.out[CurrSim - 1][iy] = 1.0 * Concurrent[0] / TotPop[0];
}

void Pop::GetMarriage()
{
	int ic, ia, ir, ig;
	int MarriedM[10][3], MarriedF[10][3], TotM[10][3], TotF[10][3];
	double MarriedTemp[2], TotalTemp[2], MarriedYth, TotYouth;

	for (ig = 0; ig < 2; ig++){
		MarriedTemp[ig] = 0.0;
		TotalTemp[ig] = 0.0;
	}
	MarriedYth = 0.0;
	TotYouth = 0.0;
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].AgeGroup > 2){
			ig = Register[ic].SexInd;
			TotalTemp[ig] += Register[ic].PopWeight;
			if (Register[ic].AgeGroup < 5){ TotYouth += Register[ic].PopWeight; }
			if (Register[ic].MarriedInd == 1){ 
				MarriedTemp[ig] += Register[ic].PopWeight; 
				if (Register[ic].AgeGroup < 5){ MarriedYth += Register[ic].PopWeight; }
			}
		}
	}
	MarriedTrendM.out[CurrSim - 1][CurrYear - StartYear] = MarriedTemp[0] / TotalTemp[0];
	MarriedTrendF.out[CurrSim - 1][CurrYear - StartYear] = MarriedTemp[1] / TotalTemp[1];
	MarriedYouth.out[CurrSim - 1][CurrYear - StartYear] = MarriedYth / TotYouth;

	if (CurrYear == 2001 || CurrYear == 2011){
		for (ia = 0; ia < 10; ia++){
			for (ir = 0; ir < 3; ir++){
				MarriedM[ia][ir] = 0;
				MarriedF[ia][ir] = 0;
				TotM[ia][ir] = 0;
				TotF[ia][ir] = 0;
			}
		}
		for (ic = 0; ic < Register.size(); ic++){
			if (Register[ic].AliveInd == 1 && Register[ic].AgeGroup>2 &&
				Register[ic].AgeGroup < 13){
				ia = Register[ic].AgeGroup - 3;
				ir = Register[ic].Race;
				if (Register[ic].SexInd == 0){
					TotM[ia][ir] += 1;
					if (Register[ic].MarriedInd == 1){ MarriedM[ia][ir] += 1; }
				}
				else{
					TotF[ia][ir] += 1;
					if (Register[ic].MarriedInd == 1){ MarriedF[ia][ir] += 1; }
				}
			}
		}

		for (ir = 0; ir < 3; ir++){
			for (ia = 0; ia < 10; ia++){
				if (TotM[ia][ir] > 0){
					if (CurrYear == 2001){
						Married2001M.out[CurrSim - 1][ir * 10 + ia] = 1.0 * MarriedM[ia][ir] / TotM[ia][ir];
					}
					else{ Married2011M.out[CurrSim - 1][ir * 10 + ia] = 1.0 * MarriedM[ia][ir] / TotM[ia][ir]; }
				}
				if (TotF[ia][ir] > 0){
					if (CurrYear == 2001){
						Married2001F.out[CurrSim - 1][ir * 10 + ia] = 1.0 * MarriedF[ia][ir] / TotF[ia][ir];
					}
					else{ Married2011F.out[CurrSim - 1][ir * 10 + ia] = 1.0 * MarriedF[ia][ir] / TotF[ia][ir]; }
				}
			}
		}
	}
}

void Pop::GetMigrationOutput()
{
	int ic, ia, ir, ih, ig;
	double urban, rural, UrbAge[18][2], TotAge[18][2], UrbRace[3], TotRace[3], UrbEdu[6], TotEdu[6];

	urban = 0.0;
	rural = 0.0;

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1){
			if (Register[ic].CurrUrban == 0){ rural += Register[ic].PopWeight; }
			else { urban += Register[ic].PopWeight; }
		}
	}
	UrbanTrend.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * urban / (urban + rural);

	if (CurrYear == 2011 || CurrYear == 1996){
		for (ia = 0; ia < 18; ia++){
			for (ig = 0; ig < 2; ig++){
				UrbAge[ia][ig] = 0.0;
				TotAge[ia][ig] = 0.0;
			}
		}
		for (ir = 0; ir < 3; ir++){
			UrbRace[ir] = 0.0;
			TotRace[ir] = 0.0;
		}
		for (ih = 0; ih < 6; ih++){
			UrbEdu[ih] = 0.0;
			TotEdu[ih] = 0.0;
		}
		for (ic = 0; ic < Register.size(); ic++){
			if (Register[ic].AliveInd == 1){
				ia = Register[ic].AgeGroup;
				ig = Register[ic].SexInd;
				ir = Register[ic].Race;
				if (Register[ic].HighestGrade == 0){ ih = 0; }
				else if (Register[ic].HighestGrade <= 5){ ih = 1; }
				else if (Register[ic].HighestGrade <= 7){ ih = 2; }
				else if (Register[ic].HighestGrade <= 11){ ih = 3; }
				else if (Register[ic].HighestGrade == 12){ ih = 4; }
				else { ih = 5; }
				TotAge[ia][ig] += Register[ic].PopWeight;
				TotRace[ir] += Register[ic].PopWeight;
				TotEdu[ih] += Register[ic].PopWeight;
				if (Register[ic].CurrUrban == 1){
					UrbAge[ia][ig] += Register[ic].PopWeight;
					UrbRace[ir] += Register[ic].PopWeight;
					UrbEdu[ih] += Register[ic].PopWeight;
				}
			}
		}
		for (ia = 0; ia < 18; ia++){
			for (ig = 0; ig < 2; ig++){
				if (CurrYear == 2011){ 
					UrbanAge.out[CurrSim - 1][ia + ig * 18] = 1.0 * UrbAge[ia][ig] / TotAge[ia][ig]; }
				else{ UrbanAge1996.out[CurrSim - 1][ia + ig * 18] = 1.0 * UrbAge[ia][ig] / TotAge[ia][ig]; }
			}
		}
		for (ir = 0; ir < 3; ir++){
			if (CurrYear == 2011){
				UrbanRace.out[CurrSim - 1][ir] = 1.0 * UrbRace[ir] / TotRace[ir]; }
			else{ UrbanRace1996.out[CurrSim - 1][ir] = 1.0 * UrbRace[ir] / TotRace[ir]; }
		}
		for (ih = 0; ih < 6; ih++){
			UrbanEdu.out[CurrSim - 1][ih] = 1.0 * UrbEdu[ih] / TotEdu[ih];
		}
	}
}

void Pop::GetAlcoholOutput()
{
	int ic, ia, ig, iy;
	double Binge[2][2], Total[2][2]; // 1st index is age, 2nd is sex
	double AnyAlc[2], Temp[2], TotDrinks, Temp2;

	for (ig = 0; ig < 2; ig++){
		for (ia = 0; ia < 2; ia++){
			Binge[ia][ig] = 0.0;
			Total[ia][ig] = 0.0;
		}
	}

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].CurrAge >= 15 && Register[ic].CurrAge < 50){
			if (Register[ic].CurrAge < 25){ ia = 0; }
			else{ ia = 1; }
			ig = Register[ic].SexInd;
			Total[ia][ig] += Register[ic].PopWeight;
			if ((Register[ic].DrinksPerDD >= 5.0 && Register[ic].DailyDrinkProb > 1.0 / 30.0) ||
				(Register[ic].DrinksPerDD * RatioMaxDrinksToAve >= 5.0 &&
				Register[ic].DailyDrinkProb * ProbMaxDrinksPerDD > 1.0 / 30.0)){
					Binge[ia][ig] += Register[ic].PopWeight;
			}
		}
	}

	iy = CurrYear - StartYear;
	BingeDrinking15to24M.out[CurrSim - 1][iy] = Binge[0][0] / Total[0][0];
	BingeDrinking15to24F.out[CurrSim - 1][iy] = Binge[0][1] / Total[0][1];
	BingeDrinking15to49M.out[CurrSim - 1][iy] = (Binge[0][0] + Binge[1][0]) / (Total[0][0] + Total[1][0]);
	BingeDrinking15to49F.out[CurrSim - 1][iy] = (Binge[0][1] + Binge[1][1]) / (Total[0][1] + Total[1][1]);

	for (ig = 0; ig < 2; ig++){
		Total[0][ig] = 0.0;
		Binge[0][ig] = 0.0;
	}
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].HighestGrade >= 7 && 
			Register[ic].HighestGrade < 11 && Register[ic].InSchool == 1){
			ig = Register[ic].SexInd;
			Total[0][ig] += Register[ic].PopWeight;
			if ((Register[ic].DrinksPerDD >= 5.0 && Register[ic].DailyDrinkProb > 1.0 / 30.0) ||
				(Register[ic].DrinksPerDD * RatioMaxDrinksToAve >= 5.0 &&
				Register[ic].DailyDrinkProb * ProbMaxDrinksPerDD > 1.0 / 30.0)){
				Binge[0][ig] += Register[ic].PopWeight;
			}
		}
	}
	BingeDrinkingHSchoolM.out[CurrSim - 1][iy] = Binge[0][0] / Total[0][0];
	BingeDrinkingHSchoolF.out[CurrSim - 1][iy] = Binge[0][1] / Total[0][1];

	for (ig = 0; ig < 2; ig++){
		AnyAlc[ig] = 0.0;
		Temp[ig] = 0.0;
	}
	TotDrinks = 0.0;
	Temp2 = 0.0;

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].CurrAge >= 15){
			ig = Register[ic].SexInd;
			Temp[ig] += Register[ic].PopWeight;
			if (Register[ic].DailyDrinkProb > 1.0 / 30.0){
				AnyAlc[ig] += Register[ic].PopWeight;
			}
			TotDrinks += Register[ic].DailyDrinkProb * Register[ic].DrinksPerDD *
				Register[ic].PopWeight;
			Temp2 += 1.0;
		}
	}
	AlcoholPast12moM.out[CurrSim - 1][iy] = AnyAlc[0] / Temp[0];
	AlcoholPast12moF.out[CurrSim - 1][iy] = AnyAlc[1] / Temp[1];
	AveDrinksPerDay.out[CurrSim - 1][iy] = TotDrinks / (Temp[0] + Temp[1]);
	//AveDrinksPerDay.out[CurrSim - 1][iy] = TotDrinks / Temp2;
}

void Pop::GetGenderNormOutput()
{
	int ic, ii, iy;
	double Denominator[4], Numerator[4];

	for (ii = 0; ii < 4; ii++){
		Numerator[ii] = 0.0;
		Denominator[ii] = 0.0;
	}
	iy = CurrYear - StartYear;
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].SexInd == 0 && Register[ic].AliveInd == 1 && Register[ic].CurrAge >= 15){
			Denominator[0] += Register[ic].PopWeight;
			Numerator[0] += Register[ic].PopWeight * Register[ic].IneqGender;
			if (Register[ic].CurrAge < 25){ ii = 1; }
			else if (Register[ic].CurrAge < 35){ ii = 2; }
			else{ ii = 3; }
			Denominator[ii] += Register[ic].PopWeight;
			Numerator[ii] += Register[ic].PopWeight * Register[ic].IneqGender;
		}
	}
	IneqGenderAllM.out[CurrSim - 1][iy] = Numerator[0] / Denominator[0];
	IneqGender15to24M.out[CurrSim - 1][iy] = Numerator[1] / Denominator[1];
	IneqGender25to34M.out[CurrSim - 1][iy] = Numerator[2] / Denominator[2];
	IneqGender35plusM.out[CurrSim - 1][iy] = Numerator[3] / Denominator[3];
}

void Pop::GetHouseholdOutput()
{
	int ic, ih, ii, iy, ID, HID, MID, PID, MMID, MPID, PMID, PPID, SHID, SMID, SPID, FoundRel, ia;
	double TotPop, HHsizes[5], HHmembers[5], temp, RelToHead[16][5];
	double HomelessAge[2], TotPopAge[2], FemHomeless, DurHomeless, AlcHomeless, EmployedHomeless;

	TotPop = 0.0;
	for (ii = 0; ii < 5; ii++){
		HHsizes[ii] = 0.0;
		HHmembers[ii] = 0.0;
		for (ia = 0; ia < 16; ia++){ RelToHead[ia][ii] = 0.0; }
	}

	for (ih = 0; ih < HHregister.size(); ih++){
		if (HHregister[ih].Active == 1){
			temp = 0.0;
			HID = HHregister[ih].IDhead;
			SHID = Register[HID - 1].IDprimary; // Spouse of head
			for (ic = 0; ic < HHregister[ih].Size; ic++){
				ID = HHregister[ih].Members[ic];
				if (Register[ID - 1].AliveInd == 1){
					ia = Register[ID - 1].AgeGroup;
					if (ia > 15){ ia = 15; }
					if (CurrYear == 2018){ RelToHead[ia][4] += Register[ID - 1].PopWeight; }
					temp += Register[ID - 1].PopWeight;
					if (ID == HID){ 
						HHmembers[0] += Register[ID - 1].PopWeight; 
						if (CurrYear == 2018){ RelToHead[ia][0] += Register[ID - 1].PopWeight; }
					}
					else if (ID == SHID){ 
						HHmembers[1] += Register[ID - 1].PopWeight; 
						if (CurrYear == 2018){ RelToHead[ia][1] += Register[ID - 1].PopWeight; }
					}
					else if (Register[ID - 1].ParentID[0] == HID || Register[ID - 1].ParentID[1] == HID ||
						(SHID > 0 && (Register[ID - 1].ParentID[0] == SHID ||
						Register[ID - 1].ParentID[1] == SHID))){
						HHmembers[2] += Register[ID - 1].PopWeight;
						if (CurrYear == 2018){ RelToHead[ia][2] += Register[ID - 1].PopWeight; }
					}
					else{
						// Check if individual is grandchild of head
						FoundRel = 0;
						MID = Register[ID - 1].ParentID[1];
						if (MID > 0){
							MMID = Register[MID - 1].ParentID[1];
							MPID = Register[MID - 1].ParentID[0];
							SMID = Register[MID - 1].IDprimary;
							if (HID == MMID || HID == MPID || (SHID > 0 && (MMID == SHID || MPID == SHID))){
								HHmembers[3] += Register[ID - 1].PopWeight;
								FoundRel = 1;
							}
							else if(SMID > 0){
								if (Register[SMID - 1].ParentID[0] == HID || Register[SMID - 1].ParentID[1] == HID){
									HHmembers[3] += Register[ID - 1].PopWeight;
									FoundRel = 1;
								}
							}
						}
						PID = Register[ID - 1].ParentID[0];
						if (PID > 0 && FoundRel == 0){
							PMID = Register[PID - 1].ParentID[1];
							PPID = Register[PID - 1].ParentID[0];
							SPID = Register[PID - 1].IDprimary;
							if (HID == PMID || HID == PPID || (SHID > 0 && (PMID == SHID || PPID == SHID))){
								HHmembers[3] += Register[ID - 1].PopWeight;
								FoundRel = 1;
							}
							else if (SPID > 0){
								if (Register[SPID - 1].ParentID[0] == HID || Register[SPID - 1].ParentID[1] == HID){
									HHmembers[3] += Register[ID - 1].PopWeight;
									FoundRel = 1;
								}
							}
						}
						if (FoundRel == 0){ HHmembers[4] += Register[ID - 1].PopWeight; }
						if (CurrYear == 2018 && FoundRel==1){ RelToHead[ia][3] += Register[ID - 1].PopWeight; }
					}
				}
			}
			TotPop += temp;
			if (HHregister[ih].Size <= 2){ HHsizes[0] += temp; }
			else if (HHregister[ih].Size <= 4){ HHsizes[1] += temp; }
			else if (HHregister[ih].Size <= 6){ HHsizes[2] += temp; }
			else if (HHregister[ih].Size <= 9){ HHsizes[3] += temp; }
			else { HHsizes[4] += temp; }
		}
	}

	HomelessAge[0] = 0.0;
	HomelessAge[1] = 0.0;
	TotPopAge[0] = 0.0; 
	TotPopAge[1] = 0.0;
	FemHomeless = 0.0; 
	DurHomeless = 0.0; 
	AlcHomeless = 0.0; 
	EmployedHomeless = 0.0;
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1){
			if (Register[ic].CurrAge < 18){ TotPopAge[0] += Register[ic].PopWeight; }
			else{ TotPopAge[1] += Register[ic].PopWeight; }
		}
		if (Register[ic].AliveInd == 1 && Register[ic].HouseholdID == 0){
			if (Register[ic].CurrAge < 18){ HomelessAge[0] += Register[ic].PopWeight; }
			else{
				HomelessAge[1] += Register[ic].PopWeight;
				if (Register[ic].SexInd == 1){ FemHomeless += Register[ic].PopWeight; }
				if (Register[ic].DailyDrinkProb > 2.0 / 7.0){ AlcHomeless += Register[ic].PopWeight; }
				if (Register[ic].Employed == 1){ EmployedHomeless += Register[ic].PopWeight; }
				DurHomeless += (1.0 * CurrYear + 0.5 - Register[ic].DateHomeless) * 
					Register[ic].PopWeight;
			}
		}
	}

	iy = CurrYear - StartYear;
	if (TotPop > 0){
		IndivsInHHsized1or2.out[CurrSim - 1][iy] = HHsizes[0] / TotPop;
		IndivsInHHsized3or4.out[CurrSim - 1][iy] = HHsizes[1] / TotPop;
		IndivsInHHsized5or6.out[CurrSim - 1][iy] = HHsizes[2] / TotPop;
		IndivsInHHsized7to9.out[CurrSim - 1][iy] = HHsizes[3] / TotPop;
		IndivsInHHsized10plus.out[CurrSim - 1][iy] = HHsizes[4] / TotPop;
		HeadOfHH.out[CurrSim - 1][iy] = HHmembers[0] / TotPop;
		PartnerOfHead.out[CurrSim - 1][iy] = HHmembers[1] / TotPop;
		ChildOfHead.out[CurrSim - 1][iy] = HHmembers[2] / TotPop;
		GrandchildOfHead.out[CurrSim - 1][iy] = HHmembers[3] / TotPop;
		Homeless.out[CurrSim - 1][iy] = (HomelessAge[0] + HomelessAge[1]) / TotPop;
	}
	if (HomelessAge[1] > 0.0){
		HomelessAdult.out[CurrSim - 1][iy] = HomelessAge[1] / TotPopAge[1];
		HomelessChild.out[CurrSim - 1][iy] = HomelessAge[0] / TotPopAge[0];
		HomelessChildPropn.out[CurrSim - 1][iy] = HomelessAge[0] / (HomelessAge[0] + HomelessAge[1]);
		HomelessFem.out[CurrSim - 1][iy] = FemHomeless / HomelessAge[1];
		HomelessEmployed.out[CurrSim - 1][iy] = EmployedHomeless / HomelessAge[1];
		HomelessAveDur.out[CurrSim - 1][iy] = DurHomeless / HomelessAge[1];
		HomelessDrinkGT2perWeek.out[CurrSim - 1][iy] = AlcHomeless / HomelessAge[1];
	}
	if (CurrYear == 2018){
		for (ia = 0; ia < 16; ia++){
			for (ii = 0; ii < 4; ii++){
				HHmembersByAge.out[CurrSim - 1][ii * 16 + ia] = RelToHead[ia][ii] / RelToHead[ia][4];
			}
		}
	}
}

void Pop::GetDiscordance()
{
	int ic, ig, PID, num[2], denom[2], is1, is2;
	int MarriedHIV[2][2]; // # marriages, by HIV status of male (1st index) & female (2nd index)

	for (ig = 0; ig<2; ig++){
		num[ig] = 0;
		denom[ig] = 0;
		MarriedHIV[ig][0] = 0;
		MarriedHIV[ig][1] = 0;
	}

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].HIVstage>0){
			ig = Register[ic].SexInd;
			if (Register[ic].CurrPartners > 0){
				PID = Register[ic].IDprimary;
				if (Register[PID - 1].SexInd != ig){
					// We are limiting calculation of discordance to heterosexual couples
					denom[ig] += 1;
					if (Register[PID - 1].HIVstage > 0){ num[ig] += 1; }
				}
			}
			if (Register[ic].CurrPartners == 2){
				PID = Register[ic].ID2ndary;
				if (Register[PID - 1].SexInd != ig){
					// We are limiting calculation of discordance to heterosexual couples
					denom[ig] += 1;
					if (Register[PID - 1].HIVstage > 0){ num[ig] += 1; }
				}
			}
		}
		if (Register[ic].AliveInd == 1 && Register[ic].MarriedInd == 1 && Register[ic].SexInd == 1){
			is1 = 0;
			is2 = 0;
			if (Register[ic].HIVstage > 0){ is2 = 1; }
			PID = Register[ic].IDprimary;
			if (Register[PID-1].HIVstage > 0){ is1 = 1; }
			MarriedHIV[is1][is2] += 1;
		}
	}

	DiscordantPropnM.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * num[0] / denom[0];
	DiscordantPropnF.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * num[1] / denom[1];
	DiscordantPropn.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * (MarriedHIV[0][1] +
		MarriedHIV[1][0]) / (MarriedHIV[0][1] + MarriedHIV[1][0] + MarriedHIV[1][1]);
	DiscordantFtoM.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * MarriedHIV[0][1] /
		(MarriedHIV[0][1] + MarriedHIV[1][0]);
}

void Pop::GetDiscordance2()
{
	int ic, ig, PID, is1, is2;
	double MarriedHIV[2][2]; // # heterosexual marriages, by HIV status of male (1st index) & female (2nd index)
	double UnmarriedHIV[2][2]; // # hetero ST relationships, by HIV status of male (1st index) & female (2nd index)

	for (ig = 0; ig<2; ig++){
		MarriedHIV[ig][0] = 0.0;
		MarriedHIV[ig][1] = 0.0;
		UnmarriedHIV[ig][0] = 0.0;
		UnmarriedHIV[ig][1] = 0.0;
	}

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 1){
			if (Register[ic].HIVstage > 0){ is2 = 1; }
			else{ is2 = 0; }
			if (Register[ic].CurrPartners > 0){
				PID = Register[ic].IDprimary;
				if (Register[PID - 1].HIVstage > 0){ is1 = 1; }
				else{ is1 = 0; }
				if (Register[ic].MarriedInd == 1){ MarriedHIV[is1][is2] += Register[ic].PopWeight; }
				else{ UnmarriedHIV[is1][is2] += Register[ic].PopWeight; }
			}
			if (Register[ic].CurrPartners == 2){
				PID = Register[ic].ID2ndary;
				if (Register[PID - 1].HIVstage > 0){ is1 = 1; }
				else{ is1 = 0; }
				UnmarriedHIV[is1][is2] += Register[ic].PopWeight;
			}
		}
	}

	DiscordantMposLT.out[CurrSim - 1][CurrYear - StartYear] = MarriedHIV[1][0];
	DiscordantMposST.out[CurrSim - 1][CurrYear - StartYear] = UnmarriedHIV[1][0];
	DiscordantFposLT.out[CurrSim - 1][CurrYear - StartYear] = MarriedHIV[0][1];
	DiscordantFposST.out[CurrSim - 1][CurrYear - StartYear] = UnmarriedHIV[0][1];
	ConcordantPosLT.out[CurrSim - 1][CurrYear - StartYear] = MarriedHIV[1][1];
	ConcordantPosST.out[CurrSim - 1][CurrYear - StartYear] = UnmarriedHIV[1][1];
	ConcordantNegLT.out[CurrSim - 1][CurrYear - StartYear] = MarriedHIV[0][0];
	ConcordantNegST.out[CurrSim - 1][CurrYear - StartYear] = UnmarriedHIV[0][0];
}

void Pop::GetCasualCalib()
{
	int ic, ig, iy;
	double temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, Casual[2], TempPop[2];

	iy = CurrYear - StartYear;
	temp1 = 0.0;
	temp2 = 0.0;
	temp3 = 0.0;
	temp4 = 0.0;
	temp5 = 0.0;
	temp6 = 0.0;
	temp7 = 0.0;
	temp8 = 0.0;
	temp9 = 0.0;
	for (ig = 0; ig < 2; ig++){
		Casual[ig] = 0.0;
		TempPop[ig] = 0.0;
	}

	if (CurrYear == 2005){
		for (ic = 0; ic < Register.size(); ic++){
			if (Register[ic].SexInd == 0 && Register[ic].AliveInd == 1 &&
				Register[ic].CurrAge >= 25 && Register[ic].CurrAge < 45 &&
				(Register[ic].CasualLast3mo + Register[ic].CurrPartners >= 2)){
				temp1 += Register[ic].PopWeight;
				if (Register[ic].CasualLast3mo >= 1){ temp2 += Register[ic].PopWeight; }
				temp3 += Register[ic].PopWeight * (Register[ic].CasualLast3mo + Register[ic].CurrPartners);
				if (Register[ic].RiskGroup == 1){ temp4 += Register[ic].PopWeight; }
				if (Register[ic].Employed == 1){ temp8 += Register[ic].PopWeight; }
			}
			if (Register[ic].SexInd == 1 && Register[ic].AliveInd == 1 &&
				Register[ic].DOLB >(1.0 * CurrYear + 0.5)){
				temp5 += Register[ic].PopWeight;
				if (Register[ic].EverCasual == 1){
					temp6 += Register[ic].PopWeight;
					if (Register[ic].RiskGroup == 1){ temp7 += Register[ic].PopWeight; }
				}
				if (Register[ic].MarriedInd == 1){ temp9 += Register[ic].PopWeight; }
			}
			if (Register[ic].AliveInd == 1 && Register[ic].CurrAge >= 15 &&
				Register[ic].CurrAge < 50){
				ig = Register[ic].SexInd;
				TempPop[ig] += Register[ic].PopWeight;
				if (Register[ic].HetCasualInd == 1){ Casual[ig] += Register[ic].PopWeight; }
			}
		}

		CasualCalib2005.out[CurrSim - 1][0] = temp2 / temp1;
		CasualCalib2005.out[CurrSim - 1][1] = temp3 / temp1;
		CasualCalib2005.out[CurrSim - 1][2] = temp4 / temp1;
		CasualCalib2005.out[CurrSim - 1][3] = temp6 / temp5;
		CasualCalib2005.out[CurrSim - 1][4] = temp7 / temp6;
		CasualCalib2005.out[CurrSim - 1][5] = Casual[0] / TempPop[0];
		CasualCalib2005.out[CurrSim - 1][6] = Casual[1] / TempPop[1];
		CasualCalib2005.out[CurrSim - 1][7] = temp8 / temp1;
		CasualCalib2005.out[CurrSim - 1][8] = temp9 / temp5;
	}
	else{
		for (ic = 0; ic < Register.size(); ic++){
			if (Register[ic].AliveInd == 1 && Register[ic].CurrAge >= 15 &&
				Register[ic].CurrAge < 50){
				ig = Register[ic].SexInd;
				TempPop[ig] += Register[ic].PopWeight;
				if (Register[ic].HetCasualInd == 1){ Casual[ig] += Register[ic].PopWeight; }
			}
		}
	}
	CasualPrevM.out[CurrSim - 1][iy] = Casual[0] / TempPop[0];
	CasualPrevF.out[CurrSim - 1][iy] = Casual[1] / TempPop[1];
}

void Pop::OneYear()
{
	int ii, Month, ic, iy;
	double HalfYr, MonthExcess;

	HalfYr = 0.5 * CycleS;
	BehavCycleCount = 0;
	iy = CurrYear - StartYear;
	if (CurrYear == BaselineStart && StructuralRCTcalib == 1 && StructIntScenario == 0){
		SaveBaselinePop(); 
	}
	if (CurrYear == 2000){
		//int seed = 100;
		//rg.RandomInit(seed);
		/*RRcasualLogDropIncomeF = 1.0;
		RRdebutLogDropIncomeF = 1.0;
		RR_STemployedM = 1.0;
		ORcondomBingePW = 1.0;
		RRcasualBinge[0] = 1.0;
		RRcasualBinge[1] = 1.0;
		for (ic = 0; ic < Register.size(); ic++){
			if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 0){
				Register[ic].IneqGender = 0.0;
			}
		}*/
	}

	// Update age group and reset flow variables
	if (FixedUncertainty == 1){ ResetFlow(); }
	if(CurrYear>StartYear){
		UpdateAgeGroup();
		UpdateUrban();
		if (NGind == 1 || CTind == 1 || TVind == 1 || HDind == 1 || TPind == 1 || HSVind == 1){
			ExternalSTIintroductions();}
		UpdateIneqGender();
		UpdateAlcohol();
		UpdateEmployment();
		UpdateHouseholds();
	}
	UpdateHHincome();
	UpdateMaleCirc(); 
	UpdateMalePref();
	if (IncludePopWeights == 1){ UpdatePopWeights(); }
	if (CurrYear == BaselineStart && StructuralRCTcalib == 1){ 
		ResetToBaselineStart();
	}
	
	// Calculate all stock variable outputs at the start of the year
	GetNumbersByHIVstage();
	//GetSTIconcordance();
	if ((CurrYear == 2005 || CurrYear == 2008 || CurrYear == 2012) && HIVcalib == 1){
		GetHHprevByAge();}
	if (MSMcalib == 1 || MSMcalibHIV == 1 || FixedUncertainty == 1){ 
		GetMSMcalib(); 
		if (MSMcalib == 1 || MSMcalibHIV == 1){ GetMSMcalib2(); }
	}
	if (FixedUncertainty == 1){ GetSexuallyExp(); }
	if (CurrYear == 2005 && FixedUncertainty == 1){ GetNumberPartners(); }
	if (CurrYear == 2000 && FixedUncertainty == 1){ GetAgeDisparate(); }
	if (HIVind == 1 && FixedUncertainty == 1){ GetCD4profile(); }
	//if (CurrYear == 2010){ GetHIVsurvivalTimes(); }
	//if (CurrYear == 2005){ GetVLs(); }
	if (FixedUncertainty == 1){ GetEduProfile(); }
	//if (FixedUncertainty == 1 && (CurrYear >= 1993 && CurrYear <=2008) && CurrSim<=5){ GetMigHIV(); }
	if (FixedUncertainty == 1){ GetPrisons(); }
	if (FixedUncertainty == 1 && CurrYear==2011){ GetAssortativeness(); }
	if (FixedUncertainty == 1){ GetCohabitTemp(); }
	if (FixedUncertainty == 1){ GetMarriage(); }
	if (FixedUncertainty == 1){ GetConcurrency(); }
	if (FixedUncertainty == 1){ GetCondomUse2(); }
	if (FixedUncertainty == 1 && AnnPrEPuptake[CurrYear-StartYear]>0.0){ GetPrEP(); }

	CalcNonHIVmort();
	/*if (CurrYear == 1990){
		for (int ia = 0; ia < 16; ia++){
			for (int ig = 0; ig < 2; ig++){ AgeEffectPartners[ia][ig] *= 0.5; }
			AgeEffectFSWcontact[ia] *= 0.5;
		}
	}
	if (CurrYear == 1990){
		PartnerEffectNew[0][0] = 0.0;
		PartnerEffectNew[1][0] = 0.0;
		PartnerEffectNew[0][1] = 0.0;
		PartnerEffectNew[1][1] = 0.0;
		for (ii = 1; ii < 5; ii++){ PartnerEffectFSWcontact[ii] = 0.0; }
	}*/
	/*if (CurrYear == 1990){
		for (int ia = 0; ia < 4; ia++){
			SexualDebut[ia][0] *= 0.5;
			SexualDebut[ia][1] *= 0.5;
		}
	}*/
	UpdateFecundability();
	//if(CurrYear==StartYear){UpdateCondomUse();}
	UpdateCondomUse(); 
	UpdateCondomEffects();
	CalcInitFSWageDbn(iy);
	CalcCurrMarriageRates();
	CalcNonHIVfert();
	CalcContrOutput();
	if (FixedUncertainty == 1){ GetMigrationOutput(); }
	//if (FixedUncertainty == 1){ GetAlcoholOutput(); }
	if (FixedUncertainty == 1){ GetHouseholdOutput(); }
	if (FixedUncertainty == 1 && CurrYear>=1990){ GetDiscordance2(); }
	if (FixedUncertainty == 1){ GetCasualCalib(); }
	GetCurrSTIprev(); 
	if (FixedUncertainty == 1 && CurrYear >= 2019){ CalcStockCosts(); }
	if (CurrYear <= 2017 && CurrYear >= 1990 && HIVcalib == 1){
		GetANCprevByAge();}
	if (CurrYear <= 2012 && CurrYear >= 1990 && FixedUncertainty==1){
		GetHIVprevByEdu();}
	if (StructuralDriverCalib == 1){
		CalcStructuralAssns();
		ResetAnnPartners();
	}
	//if(CurrYear==StartYear){UpdateProbCure();}
	UpdateProbCure();
	UpdateSTDtransitionProbs();
	
	// Project changes in behaviour, disease & demography over year
	for(ii=0; ii<CycleS; ii++){
		Month = ii * 12 / CycleS;
		MonthExcess = 12.0 * ii / CycleS - Month;
		if (MonthExcess == 0.0){ 
			UpdateConception(Month); 
			UpdateContraception(Month);
			UpdateDisclosure();
			UpdateIncarceration(Month);
			if (AnnPrEPuptake[iy]>0.0){ UpdatePrEP(Month); }
			if (AnnVaccineUptake[iy]>0.0){ UpdateVaccination(Month); }
			if (StructuralRCTcalib == 1 && CurrYear >= BaselineStart && StructIntScenario > 0){
				StructInterventionWaning(Month);
			}
			if (Month == 0 && FixedUncertainty == 1){ 
				GetAlcoholOutput();
				GetGenderNormOutput();
			}
		}
		if (ii == HalfYr){ 
			UpdateEduProfile();
			if(MSMcalib == 1 || MSMcalibHIV == 1 || FixedUncertainty == 1){ ResetRecentMSM(); } 
		}
		OneBehavCycle();
	}
	if (FixedUncertainty == 1){ CalcHIVandSTIincidence(); }
	if (FixedUncertainty == 1 && CurrYear >= 2019){ CalcFlowCosts(); }
	CurrYear += 1;
}

void Pop::ResetFlow()
{
	int iy, ig;

	iy = CurrYear - StartYear;

	NewHIV.out[CurrSim - 1][iy] = 0;
	NewHIV_R.out[CurrSim - 1][iy] = 0;
	NewHIV_U.out[CurrSim - 1][iy] = 0;
	NewHIVexp.out[CurrSim - 1][iy] = 0.0;
	NewHSV.out[CurrSim - 1][iy] = 0;
	NewTP.out[CurrSim - 1][iy] = 0;
	NewNG.out[CurrSim - 1][iy] = 0;
	NewCT.out[CurrSim - 1][iy] = 0;
	NewTV.out[CurrSim - 1][iy] = 0;
	NewHIV_MSM.out[CurrSim - 1][iy] = 0;
	for (ig = 0; ig < 2; ig++){
		NewHIVbySex[ig] = 0.0;
		NewCTbySex[ig] = 0.0;
		NewNGbySex[ig] = 0.0;
		NewTVbySex[ig] = 0.0;
	}
	TotBirths.out[CurrSim - 1][iy] = 0;
	TeenBirths.out[CurrSim - 1][iy] = 0;
	SWsexActs.out[CurrSim - 1][iy] = 0;
	SWsexActsProt.out[CurrSim - 1][iy] = 0;

	// HCT outputs for Modelling Consortium project
	TotalTestsOI.out[CurrSim - 1][iy] = 0;
	PosTestsOI.out[CurrSim - 1][iy] = 0;
	NewDiagOI.out[CurrSim - 1][iy] = 0;
	TotalTestsPrEP.out[CurrSim - 1][iy] = 0;
	PosTestsPrEP.out[CurrSim - 1][iy] = 0;
	TotalTestsANC.out[CurrSim - 1][iy] = 0;
	PosTestsANC.out[CurrSim - 1][iy] = 0;
	NewDiagANC.out[CurrSim - 1][iy] = 0;
	TotalTestsGen.out[CurrSim - 1][iy] = 0;
	PosTestsGen.out[CurrSim - 1][iy] = 0;
	NewDiagGen.out[CurrSim - 1][iy] = 0;
	TotalTestsMMC.out[CurrSim - 1][iy] = 0;
	PosTestsMMC.out[CurrSim - 1][iy] = 0;
	TotalTestsPartner.out[CurrSim - 1][iy] = 0;
	PosTestsPartner.out[CurrSim - 1][iy] = 0;
	NewDiagPartner.out[CurrSim - 1][iy] = 0;
	TotalTestsHH_U.out[CurrSim - 1][iy] = 0;
	PosTestsHH_U.out[CurrSim - 1][iy] = 0;
	NewDiagHH_U.out[CurrSim - 1][iy] = 0;
	TotalTestsHH_R.out[CurrSim - 1][iy] = 0;
	PosTestsHH_R.out[CurrSim - 1][iy] = 0;
	NewDiagHH_R.out[CurrSim - 1][iy] = 0;
	TotalTestsMobile_U.out[CurrSim - 1][iy] = 0;
	PosTestsMobile_U.out[CurrSim - 1][iy] = 0;
	NewDiagMobile_U.out[CurrSim - 1][iy] = 0;
	TotalTestsMobile_R.out[CurrSim - 1][iy] = 0;
	PosTestsMobile_R.out[CurrSim - 1][iy] = 0;
	NewDiagMobile_R.out[CurrSim - 1][iy] = 0;
	TotalTestsFSW.out[CurrSim - 1][iy] = 0;
	PosTestsFSW.out[CurrSim - 1][iy] = 0;
	NewDiagFSW.out[CurrSim - 1][iy] = 0;
	TotalTestsMSM.out[CurrSim - 1][iy] = 0;
	PosTestsMSM.out[CurrSim - 1][iy] = 0;
	NewDiagMSM.out[CurrSim - 1][iy] = 0;
	TotalTestsSchool.out[CurrSim - 1][iy] = 0;
	PosTestsSchool.out[CurrSim - 1][iy] = 0;
	NewDiagSchool.out[CurrSim - 1][iy] = 0;
	TotalTestsANCpartner0.out[CurrSim - 1][iy] = 0;
	PosTestsANCpartner0.out[CurrSim - 1][iy] = 0;
	NewDiagANCpartner0.out[CurrSim - 1][iy] = 0;
	TotalTestsANCpartner1.out[CurrSim - 1][iy] = 0;
	PosTestsANCpartner1.out[CurrSim - 1][iy] = 0;
	NewDiagANCpartner1.out[CurrSim - 1][iy] = 0;
	TotalTestsPrison.out[CurrSim - 1][iy] = 0;
	PosTestsPrison.out[CurrSim - 1][iy] = 0;
	NewDiagPrison.out[CurrSim - 1][iy] = 0;
	TotalTestsSTI.out[CurrSim - 1][iy] = 0;
	PosTestsSTI.out[CurrSim - 1][iy] = 0;
	NewDiagSTI.out[CurrSim - 1][iy] = 0;
	TotalTestsWork.out[CurrSim - 1][iy] = 0;
	PosTestsWork.out[CurrSim - 1][iy] = 0;
	NewDiagWork.out[CurrSim - 1][iy] = 0;
	TotalTestsFPC.out[CurrSim - 1][iy] = 0;
	PosTestsFPC.out[CurrSim - 1][iy] = 0;
	NewDiagFPC.out[CurrSim - 1][iy] = 0;
	AcuteTests.out[CurrSim - 1][iy] = 0;
	LYsLost.out[CurrSim - 1][iy] = 0.0;
	LYsLostExp.out[CurrSim - 1][iy] = 0.0;
	ARTdeathsExp.out[CurrSim - 1][iy] = 0.0;
	NewART200.out[CurrSim - 1][iy] = 0;
	NewART350.out[CurrSim - 1][iy] = 0;
	NewART500.out[CurrSim - 1][iy] = 0;
	NewART500plus.out[CurrSim - 1][iy] = 0;
	NewARTexp.out[CurrSim - 1][iy] = 0.0;
	MMCoperations.out[CurrSim - 1][iy] = 0;
	ProtSexActs.out[CurrSim - 1][iy] = 0;
	TotBirthsHIV.out[CurrSim - 1][iy] = 0;
	TotBirthsART.out[CurrSim - 1][iy] = 0;
	AIDSdeaths.out[CurrSim - 1][iy] = 0;
	NonAIDSdeaths.out[CurrSim - 1][iy] = 0;
}

void Pop::ResetAnnPartners()
{
	int ic, ig;
	double MultPartners[2], Total[2];

	if (FixedUncertainty == 1){
		for (ig = 0; ig < 2; ig++){
			MultPartners[ig] = 0.0;
			Total[ig] = 0.0;
		}
		for (ic = 0; ic < Register.size(); ic++){
			if (Register[ic].AliveInd == 1 && Register[ic].CurrAge >= 15.0 && Register[ic].CurrAge < 50.0){
				ig = Register[ic].SexInd;
				Total[ig] += Register[ic].PopWeight;
				if (Register[ic].AnnPartners > 1){
					MultPartners[ig] += Register[ic].PopWeight;
				}
			}
		}
		MultPartnersM.out[CurrSim - 1][CurrYear - StartYear] = MultPartners[0] / Total[0];
		MultPartnersF.out[CurrSim - 1][CurrYear - StartYear] = MultPartners[1] / Total[1];
	}

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].CurrAge >= 10.0){
			Register[ic].AnnPartners = Register[ic].CurrPartners;
		}
	}
}

void Pop::UpdateAgeGroup()
{
	int ic;
	double AgeExact;

	for(ic=0; ic<Register.size(); ic++){
		if(Register[ic].AliveInd==1){
			AgeExact = 1.0 * CurrYear + 0.5 - Register[ic].DOB;
			if(AgeExact < 90.0){
				Register[ic].CurrAge = AgeExact;
				Register[ic].AgeGroup = AgeExact/5;}
			else{
				Register[ic].CurrAge = 90;
				Register[ic].AgeGroup = 17;}
		}
	}
}

void Pop::UpdateMaleCirc()
{
	int ic, ia, ir, M15to49, HHID;
	double circprob, MMCcircrate, UncircMen[16], denominator;
	// circprob is the base prob of MMC uptake in HIV-negative boys aged 10-14

	for (ic = 0; ic<Register.size(); ic++){
		r2[ic] = rg.Random(); }
	for (ia = 0; ia < 16; ia++){
		UncircMen[ia] = 0; }

	for (ic = 0; ic<Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 0 &&
			Register[ic].CircInd == 0 && Register[ic].CurrAge>0){
			ia = Register[ic].CurrAge;
			ir = Register[ic].Race;
			circprob = (MCprevBaseline[ia][ir] - MCprevBaseline[ia - 1][ir]) /
				(1.0 - MCprevBaseline[ia - 1][ir]);
			if (r2[ic] < circprob){ Register[ic].CircInd = 1; }
		}
	}

	if (MMCtoM_ratio[CurrYear - StartYear] > 0.0){
		for (ic = 0; ic<Register.size(); ic++){
			r2[ic] = rg.Random();}
		M15to49 = 0;
		for (ic = 0; ic < Register.size(); ic++){
			if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 0){
				if (Register[ic].AgeGroup>2 && Register[ic].AgeGroup < 10){ M15to49 += 1; }
				if (Register[ic].CircInd == 0 && Register[ic].HIVstage == 0 && Register[ic].AgeGroup>1){
					HHID = Register[ic].HouseholdID;
					if (HHID > 0){
						if (HHregister[HHID - 1].PerCapitaIncomeAdj > 0){
							Register[ic].CumSelectionProb = pow(RR_VMMClogIncome, 
								log(HHregister[HHID-1].PerCapitaIncomeAdj) - MeanLogIncome);
						}
						else{ Register[ic].CumSelectionProb = pow(RR_VMMClogIncome, -MeanLogIncome); }
					}
					else{ Register[ic].CumSelectionProb = pow(RR_VMMClogIncome, -MeanLogIncome); }
					ia = Register[ic].AgeGroup - 2;
					UncircMen[ia] += Register[ic].CumSelectionProb;
				}
			}
		}
		denominator = 0.0;
		for (ia = 0; ia < 16; ia++){
			denominator += MMCuptake[ia] * UncircMen[ia];}
		circprob = MMCtoM_ratio[CurrYear - StartYear] * M15to49 / denominator;
		if (CurrYear >= 2025){ circprob = 0.16; }
		for (ic = 0; ic<Register.size(); ic++){
			if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 0 &&
				Register[ic].CircInd == 0 && Register[ic].VCThistory<2){
					ia = Register[ic].AgeGroup - 2;
					//if (CurrYear >= 2000 && Register[ic].CumSelectionProb > 0.0 && Register[ic].CumSelectionProb < 1.0){
					//	Register[ic].CumSelectionProb = 1.0; }
					if (r2[ic] < circprob * MMCuptake[ia] * Register[ic].CumSelectionProb){
						TotalTestsMMC.out[CurrSim - 1][CurrYear - StartYear] += Register[ic].PopWeight;
						if (Register[ic].HIVstage < 2){ 
							Register[ic].CircInd = 1; 
							Register[ic].VCThistory = 1;
							if (FixedUncertainty == 1){ MMCoperations.out[CurrSim - 1][CurrYear - StartYear] += Register[ic].PopWeight; }
						}
						else{
							Register[ic].GetDiagnosed(ic + 1, 0);
							PosTestsMMC.out[CurrSim - 1][CurrYear - StartYear] += Register[ic].PopWeight;
						}
						if (Register[ic].HIVstage == 1){ AcuteTests.out[CurrSim - 1][CurrYear - StartYear] += Register[ic].PopWeight; }
					}
			}
		}
	}
}

void Pop::UpdateUrban()
{
	int ia, ic, ig, ir, ih, PID, HHID;
	double ProbMove;

	for (ic = 0; ic<Register.size(); ic++){
		r2[ic] = rg.Random();}

	for (ic = 0; ic < Register.size(); ic++){
		//if (CurrSim == 2 && CurrYear == 2007 && StructIntScenario == 7){ cout << "Updating location for indiv " << ic << endl; }
		if (Register[ic].AliveInd == 1  && Register[ic].CurrAge >= 15){
			ia = Register[ic].AgeGroup;
			ig = Register[ic].SexInd;
			ir = Register[ic].Race;
			if (Register[ic].HighestGrade == 0){ ih = 0; }
			else if (Register[ic].HighestGrade <= 5){ ih = 1; }
			else if (Register[ic].HighestGrade <= 7){ ih = 2; }
			else if (Register[ic].HighestGrade <= 11){ ih = 3; }
			else if (Register[ic].HighestGrade == 12){ ih = 4; }
			else { ih = 5; }
			if (Register[ic].CurrUrban == 0){
				ProbMove = RuralToUrban[ia][ig] * RuralAdjRace[ir] * RuralAdjEdu[ih]; }
			else{ ProbMove = UrbanToRural[ia][ig] * UrbanAdjRace[ir] * UrbanAdjEdu[ih]; }
			if (Register[ic].MarriedInd == 1){
				PID = Register[ic].IDprimary;
				if (Register[ic].CurrUrban == Register[PID - 1].CurrUrban){
					ProbMove *= MarriedMigAdj[ir]; }
				if (Register[ic].CurrUrban != Register[PID - 1].CurrUrban){
					ProbMove *= SeparatedMigAdj[ir];}
			}
			if (Register[ic].ParentID[0] > 0){
				PID = Register[ic].ParentID[0];
				if (Register[ic].CurrUrban != Register[PID - 1].CurrUrban &&
					Register[PID - 1].AliveInd == 1){
					ProbMove *= 2.0;
				}
			}
			if (Register[ic].ParentID[1] > 0){
				PID = Register[ic].ParentID[1];
				if (Register[ic].CurrUrban != Register[PID - 1].CurrUrban &&
					Register[PID - 1].AliveInd == 1){
					ProbMove *= 2.0;
				}
			}
			ProbMove = 1.0 - exp(-ProbMove); // Converting rate to probability
			if (r2[ic] < ProbMove){
				Register[ic].CurrUrban = 1 - Register[ic].CurrUrban;
				Register[ic].DateMig = 0.5 + CurrYear;
				HHID = Register[ic].HouseholdID;
				Register[ic].ChangeHHafterMigration(ic + 1);
				if (HHID > 0){ HHregister[HHID - 1].RemoveMember(ic + 1); }
			}
		}
	}
}

void Pop::ExternalSTIintroductions()
{
	// This function has been introduced to prevent stochastic extinction of STIs. Each year, one 
	// new case of each STI is introduced randomly into the FSW population, assumed to arise from sexual
	// contact with clients outside the modelled population.

	int ii, FSWposition, FSWsID;
	double rs[6];

	for (ii = 0; ii<6; ii++){
		rs[ii] = rg.Random();
	}
	for (ii = 0; ii < 6; ii++){
		FSWposition = rs[ii] * TotCurrFSW[0];
		FSWsID = CSWregister[FSWposition][0];
		if (ii == 0 && NGind == 1 && Register[FSWsID - 1].NGstage == 0){ Register[FSWsID - 1].NGstage = 2; }
		if (ii == 1 && CTind == 1 && Register[FSWsID - 1].CTstage == 0){ Register[FSWsID - 1].CTstage = 2; }
		if (ii == 2 && TVind == 1 && Register[FSWsID - 1].TVstage == 0){ Register[FSWsID - 1].TVstage = 2; }
		if (ii == 3 && HDind == 1 && Register[FSWsID - 1].HDstage == 0){ Register[FSWsID - 1].HDstage = 2; }
		if (ii == 4 && TPind == 1 && Register[FSWsID - 1].TPstage == 0){ Register[FSWsID - 1].TPstage = 3; }
		if (ii == 5 && HSVind == 1 && Register[FSWsID - 1].HSVstage == 0){ Register[FSWsID - 1].HSVstage = 2; }
	}
}

void Pop::UpdateEduProfile()
{
	int ic, ir, ia, ig, ii, ih, ip;
	double ProbRepetition, ProbDropout, residual, temp1;
	double AveParentEdu[3][2], TotParents[3][2]; // By race of child (1st index) and sex of parent (2nd index)
	int DropoutsF[2]; // Number of girls dropping out of school due to pregnancy/other reasons

	for (ic = 0; ic<Register.size(); ic++){
		r2[ic] = rg.Random();}
	DropoutsF[0] = 0;
	DropoutsF[1] = 0;

	// Calculate AveParentEdu
	for (ir = 0; ir < 3; ir++){
		for (ig = 0; ig < 2; ig++){
			AveParentEdu[ir][ig] = 0.0;
			TotParents[ir][ig] = 0.0;
		}
	}
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].CurrAge >= 5 &&
			Register[ic].HighestGrade < 12 && Register[ic].InSchool == 1){
			ir = Register[ic].Race;
			for (ig = 0; ig < 2; ig++){
				if (Register[ic].ParentID[ig] > 0){
					ii = Register[ic].ParentID[ig];
					ih = Register[ii-1].HighestGrade;
					if (ih == 13){ ih = 15; }
					AveParentEdu[ir][ig] += 1.0 * ih;
				}
				else{ AveParentEdu[ir][ig] += 8.0; }
				TotParents[ir][ig] += 1.0;
			}
		}
	}
	for (ir = 0; ir < 3; ir++){
		for (ig = 0; ig < 2; ig++){
			AveParentEdu[ir][ig] = AveParentEdu[ir][ig] / TotParents[ir][ig];
		}
	}

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].CurrAge >= 5 &&
			Register[ic].CurrAge < 35){
			ir = Register[ic].Race;
			ig = Register[ic].SexInd;
			ii = Register[ic].HighestGrade;
			ia = Register[ic].CurrAge;
			// First consider entry/re-entry into schooling
			if (ii == 0 && Register[ic].InSchool == 0){
				if (ia < 10){
					if (r2[ic] < Grade1entry[ia - 5][ir]){
						Register[ic].InSchool = 1;}
				}
			}
			if (Register[ic].InSchool == 0 && ii<12 && ia>=13 && ia<23){ 
				if (StructIntScenario == 4 && Register[ic].BaselinePoor == 1 && 
					CurrYear >= BaselineStart && r2[ic] < ProbReturnSupport){
					Register[ic].InSchool = 1;
				}
			}
			// Second consider progression/dropout
			else if (Register[ic].InSchool == 1){
				temp1 = 0.0;
				if (Register[ic].ParentID[0] > 0){
					ip = Register[ic].ParentID[0];
					temp1 += Register[ip - 1].HighestGrade;
				}
				else{ temp1 += 8.0; }
				if (Register[ic].ParentID[1] > 0){
					ip = Register[ic].ParentID[1];
					temp1 += Register[ip - 1].HighestGrade;
				}
				else{ temp1 += 8.0; }
				temp1 = temp1 - AveParentEdu[ir][0] - AveParentEdu[ir][1];
				if (ii >= 12 || Dropout[ii][ir][ig] == 0){ ProbDropout = Dropout[ii][ir][ig]; }
				else{ ProbDropout = 1.0 / (1.0 + ((1.0 - Dropout[ii][ir][ig]) / Dropout[ii][ir][ig]) / pow(ParentEduEffectDropout, temp1)); }
				if (ii<9){ ProbDropout *= pow(DropoutRedn[ir], CurrYear - 1984 - ii); }
				if (StructIntScenario == 4 && Register[ic].BaselinePoor == 1 && ii < 12 && ii >= 7 &&
					CurrYear >= BaselineStart){
					ProbDropout *= RRdropoutSupport * RRdropoutIncSupport;
				}
				if (ii < 12){
					//ProbRepetition = (1.0 - Dropout[ii][ir][ig]) / (1.0 + ((1.0 - GradeRepetition[ii][ir][ig]) /
					//	GradeRepetition[ii][ir][ig]) / pow(ConscientiousEffectGradeRep, Register[ic].Conscientiousness));
					ProbRepetition = (1.0 - ProbDropout) / (1.0 + ((1.0 - GradeRepetition[ii][ir][ig]) /
						GradeRepetition[ii][ir][ig]) / (pow(ConscientiousEffectGradeRep, Register[ic].Conscientiousness) *
						pow(ParentEduEffectGradeRep, temp1)));
				}
				else{ ProbRepetition = (1.0 - Dropout[ii][ir][ig]) * GradeRepetition[ii][ir][ig]; }
				if (ig == 1 && Register[ic].DOLB>CurrYear && ii<12){
					//ProbDropout = 1.0 - (1.0 - Dropout[ii][ir][ig])*(1.0 - PregDropout[ir]);
					ProbDropout = 1.0 - (1.0 - ProbDropout)*(1.0 - PregDropout[ir]);
					ProbRepetition = 1.0 - ProbDropout;
				}
				if (r2[ic] < ProbDropout){
					Register[ic].InSchool = 0;
					if (ig == 1){
						//if (r2[ic] < Dropout[ii][ir][ig]){ DropoutsF[1] += Register[ic].PopWeight; }
						//else{ DropoutsF[0] += Register[ic].PopWeight; }
						if (Register[ic].DOLB>CurrYear && ii<12 && r2[ic] < PregDropout[ir]){ DropoutsF[0] += Register[ic].PopWeight; }
						else{ DropoutsF[1] += Register[ic].PopWeight; }
					}
				}
				else if (r2[ic] >= ProbDropout + ProbRepetition){
					Register[ic].HighestGrade += 1;
					//Register[ic].CondomPref *= CondomEdu[ii + 1] / CondomEdu[ii];
					if (Register[ic].HighestGrade == 8 && Register[ic].SexInd == 0){
						Register[ic].IneqGender *= EduEffectIneqGender[1]; }
					if (Register[ic].HighestGrade == 12 && Register[ic].SexInd == 0){
						Register[ic].IneqGender *= EduEffectIneqGender[2]; }
					if (ii == 11){ 
						residual = (r2[ic] - ProbDropout - ProbRepetition) / (1.0 -
							ProbDropout - ProbRepetition);
						temp1 = 1.0 / (1.0 + ((1.0 - TertiaryEnrol[ir]) / TertiaryEnrol[ir]) /
							pow(ConscientiousEffectTertiaryEd, Register[ic].Conscientiousness));
						if (residual > temp1){
							Register[ic].InSchool = 0;}
					}
					if (ii == 12){ Register[ic].InSchool = 0; }
				}
			}
		}
		if (Register[ic].AliveInd == 1 && Register[ic].CurrAge == 35){
			Register[ic].InSchool = 0;}
	}
	
	if (FixedUncertainty == 1){
		DropoutDuePregnancy.out[CurrSim - 1][CurrYear - StartYear] =
			1.0 * DropoutsF[0] / (DropoutsF[0] + DropoutsF[1]);
	}
}

void Pop::UpdateIneqGender()
{
	int ic;
	double TempOdds;

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].SexInd == 0 && Register[ic].AliveInd == 1){
			TempOdds = Register[ic].IneqGender / (1.0 - Register[ic].IneqGender);
			if (Register[ic].CurrAge == 10){
				//TempOdds *= pow(AgeEffectIneqGender, -10);
				if (Register[ic].CurrUrban == 0){ TempOdds *= RuralEffectIneqGender; }
			}
			else if (Register[ic].CurrAge > 20){
				TempOdds *= AgeEffectIneqGender;
			}
			if (Register[ic].CurrAge >= 10 && TempOdds > 0.0){ 
				Register[ic].IneqGender = TempOdds / (1.0 + TempOdds); }
		}
	}
}

void Pop::UpdateAlcohol()
{
	int ic, ia, ih, ir, ig;
	double Temp;

	for (ic = 0; ic<Register.size(); ic++){
		ig = Register[ic].SexInd;
		Register[ic].TempCondomAdj = 1.0;
		if (Register[ic].CurrAge >= 15 && Register[ic].AliveInd == 1){
			Temp = BaseDrinkProb[ig];
			if (Register[ic].CurrUrban == 1){ Temp += UrbanEffectDrinkProb[ig]; }
			ir = Register[ic].Race;
			Temp += RaceEffectDrinkProb[ir][ig];
			ia = (Register[ic].AgeGroup - 3) / 2;
			if (ia >= 4){ ia = 4; }
			Temp += AgeEffectDrinkProb[ia][ig];
			ih = 0;
			if (Register[ic].HighestGrade >= 13){ ih = 3; }
			else if (Register[ic].HighestGrade == 12){ ih = 2; }
			else if (Register[ic].HighestGrade >= 8){ ih = 1; }
			Temp += EduEffectDrinkProb[ih][ig];
			if (Register[ic].MarriedInd == 1){ Temp += MarriedEffectDrinkProb[ig]; }
			Register[ic].DailyDrinkProb = exp(Temp) / (7.0 * AlcPropnReported) +
				Register[ic].DrinkProbConstant;
			if (Register[ic].DailyDrinkProb < 0.0){ Register[ic].DailyDrinkProb = 0.0; }
			if (Register[ic].DailyDrinkProb > 1.0){ Register[ic].DailyDrinkProb = 1.0; }
			Temp = BaseDrinksPerDD[ig];
			if (Register[ic].CurrUrban == 1){ Temp += UrbanEffectDrinksPerDD[ig]; }
			Temp += RaceEffectDrinksPerDD[ir][ig];
			Temp += AgeEffectDrinksPerDD[ia][ig];
			Temp += EduEffectDrinksPerDD[ih][ig];
			if (Register[ic].Employed == 1){ Temp += EmployedEffectDrinksPerDD[ig]; }
			Temp = exp(Temp);
			Temp += ConscientiousEffectDrinksPerDD * Register[ic].Conscientiousness;
			if (ig == 0){ Temp += GenderIneqEffectDrinksPerDD * (Register[ic].IneqGender - 0.2); }
			Register[ic].DrinksPerDD = Temp / AlcPropnReported + Register[ic].DrinksPerDDconstant;
			if (Register[ic].DrinksPerDD < 1.0){ Register[ic].DrinksPerDD = 1.0; }
			if (Register[ic].DrinksPerDD >= 5.0){
				Register[ic].TempCondomAdj = pow(ORcondomBingePW, Register[ic].DailyDrinkProb * 7.0);
			}
		}
	}
}

void Pop::UpdateEmployment()
{
	int ic, ia, ig, ih, ir, is, iy, HHID;
	double odds, testprob, newrand, EmployRate, TotEmployed[3][2], WorkingAge[3][2], temp1, temp2;
	int PrevEmployed, EmployedOut[11][2]; // Purely for checking calibration to LFS data

	for (ic = 0; ic<Register.size(); ic++){
		r2[ic] = rg.Random();}
	iy = CurrYear - StartYear;

	/*for (ia = 0; ia < 11; ia++){
		EmployedOut[ia][0] = 0;
		EmployedOut[ia][1] = 0;
	}*/
	
	for (ic = 0; ic < Register.size(); ic++){
		if (CurrYear == BaselineStart && StructuralRCTcalib == 1 && Register[ic].AliveInd == 1){
			Register[ic].BaselineUnemployed = 1 - Register[ic].Employed;
		}
		if (Register[ic].AgeGroup < 3 || Register[ic].AgeGroup > 12 || Register[ic].InSchool == 1 || 
			Register[ic].AliveInd == 0){
			Register[ic].Employed = 0;}
		else{
			odds = BaseEmployed2[iy];
			if (Register[ic].CurrUrban == 1){ odds *= UrbanEffectEmployed2[1]; }
			ir = Register[ic].Race;
			if (CurrYear <= 2008){ odds *= RaceEffectEmployed2[ir]; }
			else{ odds *= RaceEffectEmployed3[ir]; }
			ig = Register[ic].SexInd;
			if (ig == 1){ odds *= SexEffectEmployed2[ir]; }
			ia = Register[ic].AgeGroup - 3;
			//if (CurrYear == 2015){ EmployedOut[ia][0] += 1; }
			odds *= AgeEffectEmployed2[ia];
			ih = 0;
			if (Register[ic].HighestGrade > 0 && Register[ic].HighestGrade <= 7){ ih = 1; }
			if (Register[ic].HighestGrade > 7 && Register[ic].HighestGrade <= 11){ ih = 2; }
			if (Register[ic].HighestGrade == 12){ ih = 3; }
			if (Register[ic].HighestGrade == 13){ ih = 4; }
			odds *= EduEffectEmployed2[ih];
			if (Register[ic].HIVstage > 0){
				is = Register[ic].HIVstage - 1;
				odds *= HIVeffectEmployed[is];
			}
			if (Register[ic].SexInd == 1 && (0.5 + CurrYear - Register[ic].DOLB) < 15.0){
				odds *= ChildEffectEmployed;
			}
			PrevEmployed = Register[ic].Employed;
			if (PrevEmployed == 1){ odds *= PrevEmployedEffect; }
			else if (StructIntScenario == 5 && CurrYear >= BaselineStart && Register[ic].CurrAge < 50){
				if (Register[ic].BaselinePoor == 1){
					odds *= 1.0 / ORunemployedTraining;
				}
			}
			if (r2[ic] < (odds / (1.0 + odds))){
				Register[ic].Employed = 1;
				//if (CurrYear == 2015){ EmployedOut[ia][1] += 1; }
				// Test if individual gets tested
				/*if (TestingWorkplace[iy] > 0.0){
					newrand = r2[ic] / (odds / (1.0 + odds));
					testprob = TestingWorkplace[iy] * WorkTestUptake[ig] * EmployedReachable;
					if (Register[ic].VCThistory == 2){ testprob *= RetestAdjDiagnosed[10]; }
					if (Register[ic].HIVstage == 5){ testprob *= RetestAdjART[10]; }
					if (newrand < testprob){
						TotalTestsWork.out[CurrSim-1][iy] += 1;
						if (Register[ic].VCThistory == 2){ 
							Register[ic].GetRediagnosed(ic + 1, 3);
							PosTestsWork.out[CurrSim - 1][iy] += 1; 
						}
						else{ Register[ic].VCThistory = 1; }
						if (Register[ic].VCThistory == 1 && Register[ic].HIVstage > 1){
							Register[ic].GetDiagnosed(ic + 1, 3);
							PosTestsWork.out[CurrSim - 1][iy] += 1;
							NewDiagWork.out[CurrSim - 1][iy] += 1; 
						}
					}
				}*/
			}
			else{ 
				Register[ic].Employed = 0; 
				if (PrevEmployed == 1){ Register[ic].CheckReturnToNest(ic + 1); }
			}
		}
	}

	/*if (CurrYear == 2015 && CurrSim==1){
		for (ia = 0; ia < 11; ia++){
			EmployRate = 1.0 * EmployedOut[ia][1] / EmployedOut[ia][0];
			cout << "Employment rate at age " << ia << ": " << EmployRate << endl;
		}
	}*/

	if (FixedUncertainty == 1){
		for (ir = 0; ir < 3; ir++){
			for (ig = 0; ig < 2; ig++){
				TotEmployed[ir][ig] = 0.0;
				WorkingAge[ir][ig] = 0.0;
			}
		}
		for (ic = 0; ic < Register.size(); ic++){
			if (Register[ic].AgeGroup >= 3 && Register[ic].AgeGroup < 13 && Register[ic].AliveInd == 1){
				ig = Register[ic].SexInd;
				ir = Register[ic].Race;
				WorkingAge[ir][ig] += Register[ic].PopWeight;
				if (Register[ic].Employed == 1){ TotEmployed[ir][ig] += Register[ic].PopWeight; }
			}
		}
		temp1 = 0.0;
		temp2 = 0.0;
		for (ir = 0; ir < 3; ir++){
			for (ig = 0; ig < 2; ig++){
				temp1 += TotEmployed[ir][ig];
				temp2 += WorkingAge[ir][ig];
			}
		}
		EmployedPropn.out[CurrSim - 1][iy] = temp1 / temp2;
		EmployedPropnAM.out[CurrSim - 1][iy] = TotEmployed[0][0] / WorkingAge[0][0];
		EmployedPropnAF.out[CurrSim - 1][iy] = TotEmployed[0][1] / WorkingAge[0][1];
		EmployedPropnCM.out[CurrSim - 1][iy] = TotEmployed[1][0] / WorkingAge[1][0];
		EmployedPropnCF.out[CurrSim - 1][iy] = TotEmployed[1][1] / WorkingAge[1][1];
		EmployedPropnWM.out[CurrSim - 1][iy] = TotEmployed[2][0] / WorkingAge[2][0];
		EmployedPropnWF.out[CurrSim - 1][iy] = TotEmployed[2][1] / WorkingAge[2][1];
	}
}

void Pop::UpdateHHincome()
{
	int ia, ig, ir, ih, ii, ij, k, iy, ic, temp;
	double CumulativeIncome[20], Gini, zeroinc, zeroinc_denom, CashTransfer, GrantValue;

	iy = CurrYear - StartYear;

	// Calculate average monthly income in current year (on a log scale)
	if (CurrYear <= 2006){ ii = 1; }
	else if (CurrYear <= 2010){ ii = 2; }
	else if (CurrYear <= 2014){ ii = 3; }
	else{ ii = 4; }
	for (ia = 0; ia < 10; ia++){
		for (ig = 0; ig < 2; ig++){
			for (ir = 0; ir < 3; ir++){
				for (ih = 0; ih < 5; ih++){
					if (CurrYear <= 2002){
						CurrAveIncome[ia][ig][ir][ih] = BaseLogIncome[0] +
							EduEffectLogIncome[ih][0] + RaceEffectLogIncome[ir][0]; 
						if (ia > 0){ CurrAveIncome[ia][ig][ir][ih] += AgeEffectLogIncome[ia-1][0]; }
						if (ig > 0){ CurrAveIncome[ia][ig][ir][ih] += FemaleEffectLogIncome[0]; }
						CurrAveIncome[ia][ig][ir][ih] += AnnRealEarningsGrowth * (CurrYear - 2002) +
							log(ConsPriceIndex[iy] / ConsPriceIndex[17]);
					}
					else if (CurrYear <= 2018){
						CurrAveIncome[ia][ig][ir][ih] += 0.25 * (BaseLogIncome[ii] - BaseLogIncome[ii-1] +
							EduEffectLogIncome[ih][ii] - EduEffectLogIncome[ih][ii-1] +
							RaceEffectLogIncome[ir][ii] - RaceEffectLogIncome[ir][ii-1]);
						if (ia > 0){ 
							CurrAveIncome[ia][ig][ir][ih] += 0.25 * (
								AgeEffectLogIncome[ia - 1][ii] - AgeEffectLogIncome[ia - 1][ii-1]);
						}
						if (ig > 0){
							CurrAveIncome[ia][ig][ir][ih] += 0.25 * (
								FemaleEffectLogIncome[ii] - FemaleEffectLogIncome[ii-1]);
						}
					}
					else{
						CurrAveIncome[ia][ig][ir][ih] += AnnRealEarningsGrowth +
							log(ConsPriceIndex[iy] / ConsPriceIndex[iy-1]);
					}
				}
			}
		}
	}

	// Calculate monthly private pension
	CurrPrivatePension = PrivatePension2019 * exp(AnnRealPensionGrowth * (iy - 34)) * 
		ConsPriceIndex[iy] / ConsPriceIndex[34];

	// Call GetPerCapitaIncome
	for (ii = 0; ii < HHregister.size(); ii++){
		if (HHregister[ii].Active == 1){
			HHregister[ii].GetPerCapitaIncome(ii);
		}
	}

	// Calculate average per capita income on log scale 
	temp = 0.0;
	MeanLogIncome = 0.0;
	for (ii = 0; ii < Register.size(); ii++){
		if (Register[ii].AliveInd == 1 && Register[ii].HouseholdID > 0){
			ih = Register[ii].HouseholdID;
			temp += Register[ii].PopWeight;
			if (HHregister[ih-1].PerCapitaIncomeAdj > 0.0){
				MeanLogIncome += Register[ii].PopWeight *
					log(HHregister[ih-1].PerCapitaIncomeAdj);
			}
		}
		if (Register[ii].AliveInd == 1 && Register[ii].HouseholdID == 0){
			temp += Register[ii].PopWeight;
			// Assuming log(income) = 0 for homeless people
		}
	}
	zeroinc_denom = temp;
	MeanLogIncome = MeanLogIncome / temp;
	if (FixedUncertainty == 1){ MeanHHincomePC.out[CurrSim - 1][iy] = MeanLogIncome; }

	// Update BaselinePoor
	//if (StructuralRCTcalib == 1 && (StructIntScenario == 0 && FixedUncertainty == 1 && CurrYear == BaselineStart)){
	if (StructuralRCTcalib == 1 && ((StructIntScenario == 0 && FixedUncertainty == 0 && CurrYear == BaselineStart) ||
		(StructIntScenario >= 3 && StructIntScenario <= 5 && FixedUncertainty == 1 && CurrYear >= BaselineStart - 1))){
		for (ic = 0; ic < Register.size(); ic++){
			if (Register[ic].AliveInd == 1){
				ih = Register[ic].HouseholdID;
				if (ih > 0){
					if (log(HHregister[ih - 1].PerCapitaIncomeAdj) < MeanLogIncome){
						Register[ic].BaselinePoor = 1;
					}
					else{ Register[ic].BaselinePoor = 0; }
				}
				else{ Register[ic].BaselinePoor = 1; }
			}
		}
	}
	
	// Add cash transfers (only relevant when simulating intervention scenarios)
	if ((StructIntScenario == 3 || StructIntScenario == 4) && CurrYear >= BaselineStart){
		if (StructIntScenario == 3 && FixedUncertainty == 1){ GrantValue = 1600.0; }
		else{ GrantValue = 800.0; }
		GrantValue *= (ConsPriceIndex[iy] / ConsPriceIndex[20]) / 12.0;
		for (ii = 0; ii < HHregister.size(); ii++){
			if (StructIntScenario == 4){
				// Count the number of youth who qualify for the school support
				temp = 0;
				for (ic = 0; ic < HHregister[ii].Size; ic++){
					ij = HHregister[ii].Members[ic];
					if (Register[ij - 1].InSchool == 1 && Register[ij - 1].HighestGrade >= 7 &&
						Register[ij - 1].HighestGrade < 12){
						temp += 1;
					}
					if (Register[ij - 1].InSchool == 0 && Register[ij - 1].CurrAge >= 13 &&
						Register[ij - 1].CurrAge < 23 && Register[ij - 1].HighestGrade < 12){
						temp += 1;
					}
				}
				CashTransfer = GrantValue * temp;
			}
			else{ CashTransfer = GrantValue; }
			if (HHregister[ii].Active == 1 && log(HHregister[ii].PerCapitaIncomeAdj) < MeanLogIncome){
				HHregister[ii].PerCapitaIncomeAdj += CashTransfer / pow(HHregister[ii].Size -
					HHregister[ii].Kids * (1.0 - ChildWeightEquivScale), HHsizeAdjEquivScale);
				HHregister[ii].PerCapitaIncome += CashTransfer / HHregister[ii].Size;
			}
		}
	}

	// Calculate zero income proportion
	if (FixedUncertainty == 1){
		zeroinc = 0.0;
		for (ii = 0; ii < Register.size(); ii++){
			if (Register[ii].AliveInd == 1 && Register[ii].HouseholdID > 0){
				ih = Register[ii].HouseholdID;
				if (HHregister[ih - 1].PerCapitaIncome == 0.0){
					zeroinc += Register[ii].PopWeight;
				}
			}
		}
		ZeroIncomeHH.out[CurrSim - 1][iy] = zeroinc / zeroinc_denom;
	}

	// Calculate unsorted incomes (only relevant for Gini and median income output calcs)
	if (FixedUncertainty == 1){
		for (ir = 0; ir < 4; ir++){ TotalIncomes[ir] = 0; }
		for (ii = 0; ii < HHregister.size(); ii++){
			if (HHregister[ii].Active == 1){
				HHregister[ii].GetUnsortedIncomes();
			}
		}
		for (ii = 0; ii < Register.size(); ii++){
			if (Register[ii].AliveInd == 1 && Register[ii].HouseholdID == 0){
				ir = Register[ii].Race;
				temp = TotalIncomes[0];
				UnsortedIncomes[temp][0] = 1.0;
				TotalIncomes[0] += 1;
				temp = TotalIncomes[ir + 1];
				UnsortedIncomes[temp][ir + 1] = 1.0;
				TotalIncomes[ir + 1] += 1;
			}
		}
	}

	// Previously we adjusted condom use to take account of income effects at this point - now redundant.
	
	// Calculate sorted income measures
	if (FixedUncertainty == 1){
		for (ir = 0; ir < 4; ir++){
			// (a) Sort the incomes
			SortedIncomes[0] = UnsortedIncomes[0][ir];
			for (ii = 1; ii < TotalIncomes[ir]; ii++){
				if (UnsortedIncomes[ii][ir] < SortedIncomes[0]){
					for (ij = ii; ij > 0; ij--){ SortedIncomes[ij] = SortedIncomes[ij - 1]; }
					SortedIncomes[0] = UnsortedIncomes[ii][ir];
				}
				else if (UnsortedIncomes[ii][ir] >= SortedIncomes[ii - 1]){
					SortedIncomes[ii] = UnsortedIncomes[ii][ir];
				}
				else for (ij = 0; ij < ii - 1; ij++){
					if (UnsortedIncomes[ii][ir] >= SortedIncomes[ij] &&
						UnsortedIncomes[ii][ir] < SortedIncomes[ij + 1]){
						for (k = ii; k > ij + 1; k--){
							SortedIncomes[k] = SortedIncomes[k - 1];
						}
						SortedIncomes[ij + 1] = UnsortedIncomes[ii][ir];
						break;
					}
				}
			}
			// (b) Calculate the median incomes
			temp = 0.5 * TotalIncomes[ir];
			if (ir > 0){ MedianIncomeRace[ir - 1] = SortedIncomes[temp]; }
			if (ir == 0 && FixedUncertainty == 1){ MedianIncome.out[CurrSim - 1][iy] = SortedIncomes[temp]; }
			// (c) Calculate the income inequality outputs
			if (FixedUncertainty == 1){
				CumulativeIncome[0] = 0.0;
				Gini = 0.0;
				for (ii = 0; ii < 20; ii++){
					if (ii > 0){ CumulativeIncome[ii] = CumulativeIncome[ii - 1]; }
					Gini += CumulativeIncome[ii];
					for (ij = TotalIncomes[ir] * 0.05 * ii; ij < TotalIncomes[ir] * 0.05 * (ii + 1); ij++){
						CumulativeIncome[ii] += SortedIncomes[ij];
					}
					Gini += CumulativeIncome[ii];
				}
				Gini *= 0.05 / CumulativeIncome[19];
				if (ir == 0){ GiniIncome.out[CurrSim - 1][iy] = 1.0 - Gini; }
				if (ir == 1){ GiniIncomeB.out[CurrSim - 1][iy] = 1.0 - Gini; }
				if (ir == 2){ GiniIncomeC.out[CurrSim - 1][iy] = 1.0 - Gini; }
				if (ir == 3){ GiniIncomeW.out[CurrSim - 1][iy] = 1.0 - Gini; }
				if (ir == 0){
					PalmaRatio.out[CurrSim - 1][iy] =
						(CumulativeIncome[19] - CumulativeIncome[17]) / CumulativeIncome[7];
				}
			}
		}
	}
}

void Pop::UpdateHouseholds()
{
	int ic, NewHH;
	double odds, testprob, DrinksPerDay;

	// Leaving the nest

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].CurrAge >= 15 && Register[ic].AliveInd == 1){
			Register[ic].ChangeHHleavingNest(ic + 1);
		}
	}

	// Check if homeless individuals move back into households

	for (ic = 0; ic<Register.size(); ic++){
		r2[ic] = rg.Random();
	}

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].CurrAge >= 15 && Register[ic].HouseholdID == 0 && 
			Register[ic].AliveInd == 1){
			odds = BaseOddsNotHomeless;
			if (Register[ic].SexInd == 1){ odds *= ORnotHomelessFem; }
			if (Register[ic].Employed == 1){ odds *= ORnotHomelessEmployed; }
			testprob = odds / (1.0 + odds);
			if (r2[ic] < testprob){
				NewHH = Register[ic].FindRelativeToLiveWith(ic+1, 0);
				if (NewHH == 0){ NewHousehold(ic + 1); }
				Register[ic].DateHomeless = 0.0;
			}
		}
	}
}

void Pop::UpdatePopWeights()
{
	int ic, ia, ig, ir, iy, iy2, TempPop[18][2][3];
	double CurrPopWeights[18][2][3];

	iy = CurrYear - StartYear;
	if (CurrYear <= 2025){ iy2 = iy; }
	else{ iy2 = 40; }

	for (ia = 0; ia < 18; ia++){
		for (ig = 0; ig < 2; ig++){
			for (ir = 0; ir < 3; ir++){
				TempPop[ia][ig][ir] = 0;
			}
		}
	}
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1){
			ia = Register[ic].AgeGroup;
			ig = Register[ic].SexInd;
			ir = Register[ic].Race;
			TempPop[ia][ig][ir] += 1;
		}
	}

	for (ia = 0; ia < 18; ia++){
		for (ig = 0; ig < 2; ig++){
			for (ir = 0; ir < 3; ir++){
				if (TempPop[ia][ig][ir] > 0){
					CurrPopWeights[ia][ig][ir] = 1.0 * ThembisaTot[ia][ig][iy] *
						RaceWeights[ia][ig][iy2][ir] / TempPop[ia][ig][ir];
				}
				else{ CurrPopWeights[ia][ig][ir] = 1.0; }
			}
		}
	}

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1){
			ia = Register[ic].AgeGroup;
			ig = Register[ic].SexInd;
			ir = Register[ic].Race;
			Register[ic].PopWeight = CurrPopWeights[ia][ig][ir];
		}
	}
}

void Pop::CalcNonHIVmort()
{
	int ic, ia, ir;
	double AnnMortProb;

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1){
			// Get annual rate
			ia = Register[ic].CurrAge;
			ir = Register[ic].Race;
			if (Register[ic].SexInd == 0){
				Register[ic].NonHIVmortProb = NonAIDSmortM[ia][CurrYear - StartYear][ir];
			}
			else{
				Register[ic].NonHIVmortProb = NonAIDSmortF[ia][CurrYear - StartYear][ir];
			}
			// Convert annual rate into rate per behav cycle
			Register[ic].NonHIVmortProb = 1.0 - pow(1.0 - Register[ic].NonHIVmortProb,
				1.0 / CycleS);
		}
	}
}

void Pop::UpdateFecundability()
{
	int ic, ia, ind;
	double a, b, p, q, x;
	double TempProb;

	for (ic = 0; ic < Register.size(); ic++){
		r2[ic] = rg.Random();}

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 1){
			ia = Register[ic].CurrAge;
			if (ia >= 12 && ia < 20 && Register[ic].Fecundability == 0.0){
				TempProb = 1.0 - exp(-RateFecund[ia - 12]);
				if (r2[ic] < TempProb){ 
					ind = 2;
					// Note that the following formulas for a and b apply only when the gamma mean is 1.
					a = 1.0 / pow(SDfecundability, 2.0);
					b = a;
					p = r2[ic] / TempProb;
					q = 1.0 - r2[ic] / TempProb;
					cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
					Register[ic].Fecundability = x; 
				}
			}
			if (ia >= 35 && ia < 50 && Register[ic].Fecundability > 0.0){
				TempProb = 1.0 - exp(-RateInfecund[ia - 35]);
				if (r2[ic] < TempProb){ Register[ic].Fecundability = 0.0; }
			}
		}
	}
}

void Pop::CalcNonHIVfert()
{
	int ia, ic, ir, ii;
	double SexuallyActiveSum[35][3], TotalFemSum[35][3], TempFert[35][3];

	// First calculate SexuallyExpFert
	for(ia=0; ia<35; ia++){
		for (ir = 0; ir < 3; ir++){
			SexuallyActiveSum[ia][ir] = 0.0;
			TotalFemSum[ia][ir] = 0.0;
			// Rescaling to age at conception, not age at birth (start at age 14, not age 15).
			if (ia == 0){
				HIVnegFert[ia][ir] = 0.75 * FertilityTable[ia][CurrYear + 1 - StartYear][ir];}
			else{
				HIVnegFert[ia][ir] = 0.25 * FertilityTable[ia - 1][CurrYear - StartYear][ir] +
					0.75 * FertilityTable[ia][CurrYear + 1 - StartYear][ir];}
		}
	}
	for(ic=0; ic<Register.size(); ic++){
		if(Register[ic].AliveInd==1 && Register[ic].SexInd==1 &&
			Register[ic].CurrAge>=14 && Register[ic].CurrAge<49){
				ir = Register[ic].Race;
				ia = Register[ic].CurrAge - 14;
				TotalFemSum[ia][ir] += 1.0;
				Register[ic].NonHIVfertRate = 0.0;
				if (Register[ic].DOLB<(CurrYear + 0.5) && Register[ic].DOLW<(CurrYear + 0.5)){
					if (Register[ic].CurrPartners>0){
						Register[ic].NonHIVfertRate += 1.0 - Register[ic].CondomPrimary * CondomEffPreg;}
					if (Register[ic].CurrPartners==2){
						Register[ic].NonHIVfertRate += 1.0 - Register[ic].Condom2ndary * CondomEffPreg;}
					if (Register[ic].FSWind==1){
						Register[ic].NonHIVfertRate += 1.0 - CondomUseFSW * CondomEffPreg;}
					if (Register[ic].CurrContr > 0){
						ii = Register[ic].CurrContr;
						Register[ic].NonHIVfertRate *= (1.0 - ContrEffPreg[ii - 1]);
					}
					Register[ic].NonHIVfertRate *= Register[ic].Fecundability;
					SexuallyActiveSum[ia][ir] += Register[ic].NonHIVfertRate;
				}
			}
	}
	for(ia=0; ia<35; ia++){
		for (ir = 0; ir < 3; ir++){
			if (SexuallyActiveSum[ia][ir] > 0.0){
				SexuallyExpFert[ia][ir] = HIVnegFert[ia][ir] * TotalFemSum[ia][ir] / 
					SexuallyActiveSum[ia][ir];}
			else{ SexuallyExpFert[ia][ir] = HIVnegFert[ia][ir]; }
		}
	}

	// Then smooth SexuallyExpFert
	for (ia = 0; ia<35; ia++){
		for (ir = 0; ir < 3; ir++){
			TempFert[ia][ir] = 0.26 * SexuallyExpFert[ia][ir];
			if (ia>0){ TempFert[ia][ir] += 0.20 * SexuallyExpFert[ia - 1][ir]; }
			if (ia>1){ TempFert[ia][ir] += 0.12 * SexuallyExpFert[ia - 2][ir]; }
			if (ia>2){ TempFert[ia][ir] += 0.05 * SexuallyExpFert[ia - 3][ir]; }
			if (ia<34){ TempFert[ia][ir] += 0.20 * SexuallyExpFert[ia + 1][ir]; }
			if (ia<33){ TempFert[ia][ir] += 0.12 * SexuallyExpFert[ia + 2][ir]; }
			if (ia<32){ TempFert[ia][ir] += 0.05 * SexuallyExpFert[ia + 3][ir]; }
		}
	}
	for (ia = 0; ia<35; ia++){
		for (ir = 0; ir < 3; ir++){
			SexuallyExpFert[ia][ir] = TempFert[ia][ir]; }
	}

	// For simulating the impact of structural interventions
	if (StructuralRCTcalib == 1 && CurrYear >= BaselineStart){
		for (ia = 0; ia<6; ia++){ // Previously for ia < 35
			for (ir = 0; ir < 3; ir++){
				if ((StructIntScenario == 0 || FixedUncertainty == 1) && CurrYear == BaselineStart){
					//SexuallyExpFertBase[ia][ir][CurrYear - BaselineStart] = SexuallyExpFert[ia][ir];
					SexuallyExpFertBase[ia][ir] = SexuallyExpFert[ia][ir];
				}
				else{
					//SexuallyExpFert[ia][ir] = SexuallyExpFertBase[ia][ir][CurrYear - BaselineStart];
					SexuallyExpFert[ia][ir] = SexuallyExpFertBase[ia][ir];
				}
			}
		}
	}

	// Finally recalculate NonHIVfertRate (only relevant for purpose of calibration to ANC data)
	/*for (ic = 0; ic<Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 1 &&
			Register[ic].CurrAge >= 14 && Register[ic].CurrAge<49){
			ir = Register[ic].Race;
			ia = Register[ic].CurrAge - 14;
			Register[ic].NonHIVfertRate *= SexuallyExpFert[ia][ir]/12.0;
		}
	}*/
}

void Pop::CalcContrOutput()
{
	int ic, ii, ir, ia, ih, Tot, TotPill, TotInj;
	int TotRace[3][4], TotEdu[6][4], TotAge[7][4], AgeEver[4][4];

	// Set BirthIntervals and FertByMarriage to zero
	if (CurrYear == StartYear){
		for (ii = 0; ii < 15; ii++){ BirthIntervals.out[CurrSim - 1][ii] = 0; }
	}
	if (CurrYear == 1995){
		for (ii = 0; ii < 28; ii++){ FertByMarriage.out[CurrSim - 1][ii] = 0.0; }
	}
	
	// Calculate CurrInjectable, CurrPill
	Tot = 0;
	TotPill = 0;
	TotInj = 0;
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 1 &&
			Register[ic].AgeGroup>2 && Register[ic].AgeGroup<5 &&
			(Register[ic].CurrPartners>0 || Register[ic].FSWind == 1)){
			Tot += Register[ic].PopWeight;
			if (Register[ic].CurrContr == 1){ TotInj += Register[ic].PopWeight; }
			if (Register[ic].CurrContr == 2){ TotPill += Register[ic].PopWeight; }
		}
	}
	CurrInjectable.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotInj / Tot;
	CurrPill.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * TotPill / Tot;
	YouthHormonal.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * (TotInj + TotPill) / Tot;

	// Calculate DHS outputs
	if (CurrYear == 1998 || CurrYear == 2003){
		for (ii = 0; ii < 4; ii++){
			for (ir = 0; ir < 3; ir++){ TotRace[ir][ii] = 0; }
			for (ih = 0; ih < 6; ih++){ TotEdu[ih][ii] = 0; }
			for (ia = 0; ia < 7; ia++){ TotAge[ia][ii] = 0; }
			for (ia = 0; ia < 4; ia++){ AgeEver[ia][ii] = 0; }
		}
		for (ic = 0; ic < Register.size(); ic++){
			if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 1 &&
				Register[ic].AgeGroup>2 && Register[ic].AgeGroup < 10){
				ir = Register[ic].Race;
				ia = Register[ic].AgeGroup - 3;
				if (Register[ic].HighestGrade == 0){ ih = 0; }
				else if (Register[ic].HighestGrade <= 5){ ih = 1; }
				else if (Register[ic].HighestGrade <= 7){ ih = 2; }
				else if (Register[ic].HighestGrade <= 11){ ih = 3; }
				else if (Register[ic].HighestGrade == 12){ ih = 4; }
				else { ih = 5; }
				ii = Register[ic].CurrContr;
				if (Register[ic].CurrPartners>0 || Register[ic].FSWind == 1){
					TotRace[ir][0] += 1;
					if (ii > 0){ TotRace[ir][ii] += 1; }
					TotEdu[ih][0] += 1;
					if (ii > 0){ TotEdu[ih][ii] += 1; }
				}
				TotAge[ia][0] += 1;
				if (ii > 0){ TotAge[ia][ii] += 1; }
				if (ia < 4){
					AgeEver[ia][0] += 1;
					if (Register[ic].EverInjectable == 1){ AgeEver[ia][1] += 1; }
					if (Register[ic].EverPill == 1){ AgeEver[ia][2] += 1; }
				}
			}
		}
		for (ii = 0; ii < 3; ii++){
			for (ir = 0; ir < 3; ir++){
				if (CurrYear == 1998){
					CurrContrRace98.out[CurrSim - 1][ii * 3 + ir] = 1.0 * 
						TotRace[ir][ii + 1] / TotRace[ir][0];}
				else{
					CurrContrRace03.out[CurrSim - 1][ii * 3 + ir] = 1.0 *
						TotRace[ir][ii + 1] / TotRace[ir][0];}
			}
			for (ih = 0; ih < 6; ih++){
				if (CurrYear == 1998){
					CurrContrEdu98.out[CurrSim - 1][ii * 6 + ih] = 1.0 *
						TotEdu[ih][ii + 1] / TotEdu[ih][0];}
				else{
					CurrContrEdu03.out[CurrSim - 1][ii * 6 + ih] = 1.0 *
						TotEdu[ih][ii + 1] / TotEdu[ih][0];}
			}
			for (ia = 0; ia < 7; ia++){
				if (CurrYear == 1998){
					CurrUseContr98.out[CurrSim - 1][ii * 7 + ia] = 1.0 *
						TotAge[ia][ii + 1] / TotAge[ia][0];}
				else{
					CurrUseContr03.out[CurrSim - 1][ii * 7 + ia] = 1.0 *
						TotAge[ia][ii + 1] / TotAge[ia][0];}
			}
			for (ia = 0; ia < 4; ia++){
				if (CurrYear == 1998 && ii<2){
					EverUseContr98.out[CurrSim - 1][ii * 4 + ia] = 1.0 *
						AgeEver[ia][ii + 1] / AgeEver[ia][0];}
				else if (CurrYear == 2003 && ii<2){
					EverUseContr03.out[CurrSim - 1][ii * 4 + ia] = 1.0 *
						AgeEver[ia][ii + 1] / AgeEver[ia][0];}
			}
		}
	}
	if (CurrYear == 2016){
		for (ii = 0; ii < 4; ii++){
			for (ia = 0; ia < 7; ia++){ TotAge[ia][ii] = 0; }
		}
		for (ic = 0; ic < Register.size(); ic++){
			if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 1 &&
				Register[ic].AgeGroup>2 && Register[ic].AgeGroup < 10 &&
				(Register[ic].CurrPartners>0 || Register[ic].FSWind == 1)){
				// Note that here we're limiting analysis to sexually active women
				// (unlike in the analysis of the 1998 and 2003 DHS).
				ia = Register[ic].AgeGroup - 3;
				ii = Register[ic].CurrContr;
				TotAge[ia][0] += 1;
				if (ii > 0){ TotAge[ia][ii] += 1; }
			}
		}
		for (ii = 0; ii < 3; ii++){
			for (ia = 0; ia < 7; ia++){
				CurrUseContr16.out[CurrSim - 1][ii * 7 + ia] = 1.0 *
					TotAge[ia][ii + 1] / TotAge[ia][0];
			}
		}
	}

	// Calculate SexuallyExpFert
	if (FixedUncertainty == 1 && CurrYear==1998){
		for (ia = 0; ia < 35; ia++){
			SexuallyExpFertB.out[CurrSim - 1][ia] = SexuallyExpFert[ia][0];
			SexuallyExpFertC.out[CurrSim - 1][ia] = SexuallyExpFert[ia][1];
			SexuallyExpFertW.out[CurrSim - 1][ia] = SexuallyExpFert[ia][2];
		}
	}

	// Calculate denominators of FertByMarriage
	if (FixedUncertainty == 1 && CurrYear < 1998 && CurrYear >= 1995){
		for (ic = 0; ic < Register.size(); ic++){
			if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 1 &&
				Register[ic].AgeGroup>2 && Register[ic].AgeGroup < 10){
				ia = Register[ic].AgeGroup - 3;
				ii = Register[ic].MarriedInd;
				FertByMarriage.out[CurrSim - 1][14 + ia + ii * 7] += 1;
			}
		}
	}
}

void Pop::UpdateProbCure()
{
	/*MaleRxRate = InitMaleRxRate * (1.0 + RxPhaseIn[CurrYear-StartYear]);
	MaleTeenRxRate = InitMaleTeenRxRate * (1.0 + RxPhaseIn[CurrYear-StartYear]);
	FemRxRate = InitFemRxRate * (1.0 + RxPhaseIn[CurrYear-StartYear]);
	FemTeenRxRate = InitFemTeenRxRate * (1.0 + RxPhaseIn[CurrYear-StartYear]);*/
	/*FSWasympRxRate = InitFSWasympRxRate + RxPhaseIn[CurrYear-StartYear] * 0.5 *
		(0.25 - InitFSWasympRxRate);
	FSWasympCure = (InitFSWasympCure * (1.0 - RxPhaseIn[CurrYear-StartYear] * 0.5) *
		InitFSWasympRxRate + 1.0 * RxPhaseIn[CurrYear-StartYear] * 0.5 * 0.25)/
		((1.0 - RxPhaseIn[CurrYear-StartYear] * 0.5) * InitFSWasympRxRate +
		RxPhaseIn[CurrYear-StartYear] * 0.5 * 0.25);*/
	if(HSVind==1){
		/*HSVtransitionM.CorrectRxWithSM = InitCorrectRxHSV + (0.9 - InitCorrectRxHSV) * 
			RxPhaseIn[CurrYear-StartYear];
		HSVtransitionF.CorrectRxWithSM = HSVtransitionM.CorrectRxWithSM;*/
		HSVtransitionM.CalcProbCure();
		HSVtransitionF.CalcProbCure();}
	if(TPind==1){
		/*TPtransitionF.ANCpropnScreened = InitANCpropnScreened + RxPhaseIn[CurrYear-StartYear] *
			(1.0 - InitANCpropnScreened);
		TPtransitionF.ANCpropnTreated = InitANCpropnTreated + RxPhaseIn[CurrYear-StartYear] *
			(1.0 - InitANCpropnTreated);*/
		TPtransitionM.CalcProbCure();
		TPtransitionF.CalcProbCure();}
	if(HDind==1){
		HDtransitionM.CalcProbCure();
		HDtransitionF.CalcProbCure();}
	if(NGind==1){
		/*NGtransitionM.DrugEff = (1.0 - PropnCiproResistant[CurrYear-StartYear] *
			PropnCiproTreated[CurrYear-StartYear]) * InitDrugEffNG;
		NGtransitionF.DrugEff = NGtransitionM.DrugEff;*/
		NGtransitionM.CalcProbCure();
		NGtransitionF.CalcProbCure();}
	if(CTind==1){
		CTtransitionM.CalcProbCure();
		CTtransitionF.CalcProbCure();}
	if(TVind==1){
		//TVtransitionM.CorrectRxWithSM = InitCorrectRxTVM + (0.9 - InitCorrectRxTVM) * 
		//	RxPhaseIn[CurrYear-StartYear];
		TVtransitionM.CalcProbCure();
		TVtransitionF.CalcProbCure();}
	if(BVind==1){
		BVtransitionF.CalcProbCure();
		BVtransitionF.CalcProbPartialCure();}
	if(VCind==1){
		VCtransitionF.CalcProbCure();
		VCtransitionF.CalcProbPartialCure();}
}

void Pop::UpdateSTDtransitionProbs()
{
	if(HIVind==1){
		HIVtransitionM.CalcTransitionProbs();
		HIVtransitionF.CalcTransitionProbs();}
	if(HSVind==1){
		/*HSVtransitionM.RecurrenceRate = InitRecurrenceRateM * (1.0 - 0.8 * 
			RxPhaseIn[CurrYear-StartYear]);
		HSVtransitionF.RecurrenceRate = InitRecurrenceRateF * (1.0 - 0.8 * 
			RxPhaseIn[CurrYear-StartYear]);*/
		HSVtransitionM.CalcTransitionProbs();
		HSVtransitionF.CalcTransitionProbs();}
	if(TPind==1){
		TPtransitionM.CalcTransitionProbs();
		TPtransitionF.CalcTransitionProbs();}
	if(HDind==1){
		HDtransitionM.CalcTransitionProbs();
		HDtransitionF.CalcTransitionProbs();}
	if(NGind==1){
		NGtransitionM.CalcTransitionProbs();
		NGtransitionF.CalcTransitionProbs();}
	if(CTind==1){
		CTtransitionM.CalcTransitionProbs();
		CTtransitionF.CalcTransitionProbs();}
	if(TVind==1){
		TVtransitionM.CalcTransitionProbs();
		TVtransitionF.CalcTransitionProbs();}
	if(BVind==1){
		BVtransitionF.CalcTransitionProbs();}
	if(VCind==1){
		VCtransitionF.CalcTransitionProbs();}
}

void Pop::UpdateCondomUse()
{
	int ia, ib;
	double Rate15to19[3], x;

	// Note that in the revised model, the code for condom use in LT relationships is
	// redundant, i.e. changing RatioUltTo1998[1], ShapeBehavChange[1] and so on, will
	// have no effect on the results. The code and associated variables still need to 
	// be removed - future cleaning exercise.

	BaselineCondomUse = 0.08 + CondomScaling * (0.35 - 0.08);
	RatioUltTo1998[0] = 3.0 + CondomScaling * (10.0 - 3.0); // 20
	RatioUltTo1998[1] = 1.2 + CondomScaling * (7.5 - 1.2); // 15
	RatioUltTo1998[2] = 3.0 + CondomScaling * (12.0 - 3.0); // 24
	ShapeBehavChange[0] = 2.5 + CondomScaling * (4.0 - 2.5);
	ShapeBehavChange[1] = 3.6 + CondomScaling * (1.8 - 3.6);
	ShapeBehavChange[2] = 4.0 + CondomScaling * (4.0 - 4.0);
	RelEffectCondom[2] = 15.0 + CondomScaling * (10.0 - 15.0);

	for (ib = 0; ib<3; ib++){
		MedianToBehavChange[ib] = 13.0 * exp(-log(log(1.0 - log(RatioInitialTo1998[ib]) /
			log(RatioUltTo1998[ib])) / log(2.0)) / ShapeBehavChange[ib]);
		/*x = (BaselineCondomUse / (1.0 - BaselineCondomUse)) * RelEffectCondom[ib] *
			RatioInitialTo1998[ib] * pow(RatioUltTo1998[ib] / RatioInitialTo1998[ib], 1.0 -
			pow(0.5, pow((CurrYear - 1985) / MedianToBehavChange[ib], ShapeBehavChange[ib]))) *
			pow(RatioUltTo1998[ib], -0.25 * (1.0 - pow(0.5, pow((CurrYear - 1985)/
			26.0, 2.0 * ShapeBehavChange[ib]))));*/
		x = (BaselineCondomUse / (1.0 - BaselineCondomUse)) * RelEffectCondom[ib] *
			RatioInitialTo1998[ib] * pow(RatioUltTo1998[ib] / RatioInitialTo1998[ib], 1.0 -
			pow(0.5, pow((CurrYear - 1985) / MedianToBehavChange[ib], ShapeBehavChange[ib])));
		Rate15to19[ib] = x / (1.0 + x);
	}
	// Condom usage for females in ST & casual relationships
	for (ia = 0; ia < 16; ia++){
		x = (Rate15to19[0] / (1.0 - Rate15to19[0])) * exp(5.0 * (ia - 1) * AgeEffectCondom[0]);
		CondomUseST[ia][1] = x / (1.0 + x);
		CondomUseCasualHet[ia][1] = ORcondomCasual * CondomUseST[ia][1] /
			(1.0 + CondomUseST[ia][1] * (ORcondomCasual - 1.0));
	}
	// Condom usage for females in LT relationships
	for (ia = 0; ia < 16; ia++){
		x = (Rate15to19[1] / (1.0 - Rate15to19[1])) * exp(5.0 * (ia - 1) * AgeEffectCondom[1]);
		CondomUseLT[ia][1] = x / (1.0 + x);
	}
	// Condom use in males
	for (ia = 0; ia<16; ia++){
		CondomUseST[ia][0] = 0;
		CondomUseLT[ia][0] = 0;
		for (ib = 0; ib<16; ib++){
			CondomUseST[ia][0] += AgePrefM[ia][ib] * CondomUseST[ib][1];
			CondomUseLT[ia][0] += AgePrefM[ia][ib] * CondomUseLT[ib][1];
		}
		CondomUseCasual[ia] = ORcondomCasual * CondomUseST[ia][0] /
			(1.0 + CondomUseST[ia][0] * (ORcondomCasual - 1.0));
		CondomUseCasualHet[ia][0] = CondomUseCasual[ia];
	}
	// Condom use in FSW-client relationships
		CondomUseFSW = Rate15to19[2];
}

void Pop::UpdateCondomEffects()
{
	// This function adjusts levels of condom use to take account educational attainment and (in men)
	// inequitable gender norms. Effects of alcohol and income are accounted for in the 
	// UpdateAlcohol and UpdateHHincome functions respectively.

	int ic, ih, iy;
	double TempAdj, Median, BaseEduTrend, AltEduTrend, MedianAdjByEdu[16];

	iy = CurrYear - StartYear;
	BaseEduTrend = pow(1.0 * iy / MedianToBehavChange[0], ShapeBehavChange[0]);
	BaseEduTrend = pow(0.5, BaseEduTrend);
	for (ih = 0; ih < 16; ih++){
		Median = MedianToBehavChange[0] * pow(CondomEduMedian, (1.0 * ih - 10.0));
		AltEduTrend = pow(1.0 * iy / Median, ShapeBehavChange[0]);
		AltEduTrend = pow(0.5, AltEduTrend);
		MedianAdjByEdu[ih] = pow(RatioUltTo1998[0] / RatioInitialTo1998[0],
			BaseEduTrend - AltEduTrend);
	}

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].CurrAge >= 10){
			ih = Register[ic].HighestGrade;
			if (ih == 13){ ih = 15; }
			// First adjustment to take account of education effect that is constant over time
			TempAdj = pow(ORcondomPerYrSchool, (1.0 * ih - 10.0));
			//if (CurrYear < 2000){ TempAdj = pow(ORcondomPerYrSchool, (1.0 * ih - 10.0)); }
			//else{ TempAdj = pow(ORcondomPerYrSchool, 5.0); }
			// Second adjustment for effect of edu on median time to behaviour change
			// (Comment out the adjustment below if running a 'no change in condom use' scenario.)
			TempAdj *= MedianAdjByEdu[ih];
			// Third adjustment for inequitable gender norms
			if (Register[ic].SexInd == 0){
				TempAdj *= pow(ORcondomGenderIneq, Register[ic].IneqGender - 0.2);
			}
			Register[ic].TempCondomAdj *= TempAdj;
		}
	}
}

void Pop::CalcInitFSWageDbn(int iy)
{
	// Similar code to that in the SetInitSexActivity function in Thembisa
	int ia, ind;
	double Temp, Temp2, X, A, B, P;
	double FSWageLambda, FSWageAlpha;

	FSWageLambda = (FSWageMean[iy] - 10.0) / (FSWageSD[iy] * FSWageSD[iy]);
	FSWageAlpha = (FSWageMean[iy] - 10.0) * FSWageLambda;
	ind = 1;
	X = 0.0;
	A = FSWageAlpha;
	B = FSWageLambda;
	P = 5.0;
	cdfgam(&ind, &P, 0, &X, &A, &B, 0, 0);
	Temp = P;
	InitFSWageDbn[0] = Temp;
	for (ia = 1; ia<15; ia++){
		Temp2 = Temp;
		X = 5.0 * (ia + 1);
		P = 0.0;
		cdfgam(&ind, &P, 0, &X, &A, &B, 0, 0);
		Temp = P;
		InitFSWageDbn[ia] = Temp - Temp2;
	}
	InitFSWageDbn[15] = 1.0 - Temp;
}

void Pop::CalcCurrMarriageRates()
{
	int ia, ig, ir, iy;
	double PropnNeverMarried[77];

	for (ig = 0; ig < 2; ig++){
		for (ir = 0; ir < 3; ir++){
			for (ia = 0; ia < 5; ia++){ MarriageIncidence[ia][ig][ir] = 0.0; }
			for (ia = 0; ia <= 75; ia++){ // 0 corresponds to age 15 last birthday
				iy = CurrYear - (ia + 15) - StartYear; // Year of birth, in years after 1985
				if (ia < MarriageMin[ig] - 15){ MarriageIncidence[ia + 5][ig][ir] = 0.0; }
				else{
					PropnNeverMarried[ia] = 1.0 / (1.0 + pow((15.0 + ia - MarriageMin[ig]) *
						exp(-MarriageConstant[ig] - MarriageTrend[ig] * iy - RaceMarriage[ir]), 
						1.0 / MarriageShape[ir][ig]));
					PropnNeverMarried[ia + 1] = 1.0 / (1.0 + pow((15.0 + (ia + 1.0) - MarriageMin[ig]) *
						exp(-MarriageConstant[ig] - MarriageTrend[ig] * iy - RaceMarriage[ir]), 
						1.0 / MarriageShape[ir][ig]));
					MarriageIncidence[ia + 5][ig][ir] = 1.0 - PropnNeverMarried[ia + 1] / PropnNeverMarried[ia];
					MarriageIncidence[ia + 5][ig][ir] *= 1.0 / InSTrelationship[ia + 5][ig][ir];
				}
			}
		}
	}
}

void Pop::GetANCprev(STDtransition* a, int STDind)
{
	int ic, ii;
	double numerator, denominator;

	if (a->ANClogL.Observations>0){
		numerator = 0.0;
		denominator = 0.0;
		for (ii = 0; ii < Register.size(); ii++){
			if (Register[ii].AgeGroup>2 && Register[ii].AgeGroup < 10 &&
				Register[ii].AliveInd == 1 && Register[ii].SexInd == 1 &&
				Register[ii].DOLB>(1.0*CurrYear+0.5)){
				denominator += Register[ii].PopWeight;
				if (STDind == 1 && Register[ii].HSVstage>0){ numerator += Register[ii].PopWeight; }
				// Note that for calibration purposes we include the women who are TP-immune but exclude incubation
				if (STDind == 2 && Register[ii].TPstage>1){ numerator += Register[ii].PopWeight; }
				if (STDind == 3 && Register[ii].HDstage>0 && Register[ii].HDstage<3){ numerator += Register[ii].PopWeight; }
				if (STDind == 4 && Register[ii].NGstage>0 && Register[ii].NGstage<3){ numerator += Register[ii].PopWeight; }
				if (STDind == 5 && Register[ii].CTstage>0 && Register[ii].CTstage<3){ numerator += Register[ii].PopWeight; }
				if (STDind == 6 && Register[ii].TVstage>0 && Register[ii].TVstage<3){ numerator += Register[ii].PopWeight; }
				if (STDind == 7 && Register[ii].BVstage>1){ numerator += Register[ii].PopWeight; }
				if (STDind == 8 && Register[ii].VCstage>0){ numerator += Register[ii].PopWeight; }
			}
		}
		if (numerator > 0.0){ a->ANCprevUnsmoothed[CurrYear-StartYear] = numerator / denominator; }
		else{ a->ANCprevUnsmoothed[CurrYear - StartYear] = 0.5 / denominator; }
		// Add code for storing prevalence in current year?
		if (STDind == 1 && FixedUncertainty == 1 && HSVcalib == 1){
			HSVprevANC.out[CurrSim - 1][CurrYear - StartYear] = numerator / denominator;}
		if (STDind == 2 && FixedUncertainty == 1 && TPcalib == 1){
			TPprevANC.out[CurrSim - 1][CurrYear - StartYear] = numerator / denominator;}
		/*for (ic = 0; ic<a->ANClogL.Observations; ic++){
			if (a->ANClogL.StudyYear[ic] == CurrYear){
				a->ANClogL.ModelPrev[ic] = a->ANClogL.CurrentModelPrev;}
		}*/
	}
}

void Pop::GetHHprev(STDtransition* a, int STDind)
{
	int ic, ii, ia, id;
	double numerator, denominator;

	if (a->HouseholdLogL.Observations>0){
		for (ic = 0; ic<a->HouseholdLogL.Observations; ic++){
			if (a->HouseholdLogL.StudyYear[ic] == CurrYear){
				numerator = 0.0;
				denominator = 0.0;
				for (ii = 0; ii < Register.size(); ii++){
					ia = Register[ii].AgeGroup - 2;
					id = 0;
					if (a->HouseholdLogL.ExclVirgins[ic] == 1 && Register[ii].VirginInd == 1){ id = 1; }
					if (ia >(a->HouseholdLogL.AgeStart[ic] - 1) && ia < (a->HouseholdLogL.AgeEnd[ic] + 1) &&
						Register[ii].AliveInd == 1 && Register[ii].SexInd == a->SexInd && id == 0){
						denominator += Register[ii].PopWeight;
						if (STDind == 1 && Register[ii].HSVstage>0){ numerator += Register[ii].PopWeight; }
						// Note that for calibration purposes we include the people who are TP-immune but exclude incubation
						if (STDind == 2 && Register[ii].TPstage>1){ numerator += 1.0; }
						if (STDind == 3 && Register[ii].HDstage>0 && Register[ii].HDstage<3){ numerator += Register[ii].PopWeight; }
						if (STDind == 4 && Register[ii].NGstage>0 && Register[ii].NGstage<3){ numerator += Register[ii].PopWeight; }
						if (STDind == 5 && Register[ii].CTstage>0 && Register[ii].CTstage<3){ numerator += Register[ii].PopWeight; }
						if (STDind == 6 && Register[ii].TVstage>0 && Register[ii].TVstage<3){ numerator += Register[ii].PopWeight; }
						if (STDind == 7 && Register[ii].BVstage>1){ numerator += Register[ii].PopWeight; }
						if (STDind == 8 && Register[ii].VCstage>0){ numerator += Register[ii].PopWeight; }
					}
				}
				a->HouseholdLogL.ModelPrev[ic] = numerator / denominator;
			}
		}
	}
}

void Pop::GetFPCprev(STDtransition* a, int STDind)
{
	int ic, ii;
	double numerator, denominator, TempFP;

	if (a->FPClogL.Observations>0){
		numerator = 0.0;
		denominator = 0.0;
		for (ii = 0; ii < Register.size(); ii++){
			if (Register[ii].AgeGroup>1 && Register[ii].AliveInd == 1 && Register[ii].SexInd == 1 &&
				(Register[ii].CurrContr == 1 || Register[ii].CurrContr ==2)){
				TempFP = Register[ii].PopWeight;
				denominator += TempFP;
				if (STDind == 1 && Register[ii].HSVstage>0){ numerator += TempFP; }
				// Note that for calibration purposes we include the women who are TP-immune but exclude incubation
				if (STDind == 2 && Register[ii].TPstage>1){ numerator += TempFP; }
				if (STDind == 3 && Register[ii].HDstage>0 && Register[ii].HDstage<3){ numerator += TempFP; }
				if (STDind == 4 && Register[ii].NGstage>0 && Register[ii].NGstage<3){ numerator += TempFP; }
				if (STDind == 5 && Register[ii].CTstage>0 && Register[ii].CTstage<3){ numerator += TempFP; }
				if (STDind == 6 && Register[ii].TVstage>0 && Register[ii].TVstage<3){ numerator += TempFP; }
				if (STDind == 7 && Register[ii].BVstage>1){ numerator += TempFP; }
				if (STDind == 8 && Register[ii].VCstage>0){ numerator += TempFP; }
			}
		}
		a->FPClogL.CurrentModelPrev = numerator / denominator;
		// Add code for storing prevalence in current year?
		for (ic = 0; ic<a->FPClogL.Observations; ic++){
			if (a->FPClogL.StudyYear[ic] == CurrYear){
				a->FPClogL.ModelPrev[ic] = a->FPClogL.CurrentModelPrev;
			}
		}
	}
}

void Pop::GetCSWprev(STDtransition* a, int STDind)
{
	int ii, ia;
	double numerator, denominator;

	if (a->CSWlogL.Observations>0){
		numerator = 0.0;
		denominator = 0.0;
		for (ii = 0; ii < Register.size(); ii++){
			if (Register[ii].FSWind == 1 && Register[ii].AliveInd == 1 &&
				Register[ii].SexInd == 1){
				denominator += Register[ii].PopWeight;
				if (STDind == 0 && Register[ii].HIVstage>0){ numerator += Register[ii].PopWeight; }
				if (STDind == 1 && Register[ii].HSVstage>0){ numerator += Register[ii].PopWeight; }
				// Note that for calibration purposes we include the women who are TP-immune but exclude incubation
				if (STDind == 2 && Register[ii].TPstage>1){ numerator += Register[ii].PopWeight; }
				if (STDind == 3 && Register[ii].HDstage>0 && Register[ii].HDstage<3){ numerator += Register[ii].PopWeight; }
				if (STDind == 4 && Register[ii].NGstage>0 && Register[ii].NGstage<3){ numerator += Register[ii].PopWeight; }
				if (STDind == 5 && Register[ii].CTstage>0 && Register[ii].CTstage<3){ numerator += Register[ii].PopWeight; }
				if (STDind == 6 && Register[ii].TVstage>0 && Register[ii].TVstage<3){ numerator += Register[ii].PopWeight; }
				if (STDind == 7 && Register[ii].BVstage>1){ numerator += Register[ii].PopWeight; }
				if (STDind == 8 && Register[ii].VCstage>0){ numerator += Register[ii].PopWeight; }
			}
		}
		a->CSWprevUnsmoothed[CurrYear-StartYear] = numerator / denominator;
	}
}

void Pop::CalcHIVandSTIincidence()
{
	int iy;
	iy = CurrYear - StartYear;

	HIVinc15to49_U.out[CurrSim - 1][iy] = NewHIV_U.out[CurrSim - 1][iy] / (TotPop15to49_U.out[CurrSim - 1][iy] *
		(1.0 - UrbanHIV.out[CurrSim - 1][iy]));
	HIVinc15to49_R.out[CurrSim - 1][iy] = NewHIV_R.out[CurrSim - 1][iy] / (TotPop15to49_R.out[CurrSim - 1][iy] *
		(1.0 - RuralHIV.out[CurrSim - 1][iy]));
	HIVinc15to49.out[CurrSim - 1][iy] = (NewHIV_U.out[CurrSim - 1][iy] + NewHIV_R.out[CurrSim - 1][iy]) / 
		(TotPop15to49_U.out[CurrSim - 1][iy] * (1.0 - UrbanHIV.out[CurrSim - 1][iy]) +
		TotPop15to49_R.out[CurrSim - 1][iy] * (1.0 - RuralHIV.out[CurrSim - 1][iy]));
	NewHIVexpM.out[CurrSim - 1][iy] = NewHIVbySex[0];
	NewHIVexpF.out[CurrSim - 1][iy] = NewHIVbySex[1];

	NewCTNGTV.out[CurrSim - 1][iy] = NewCT.out[CurrSim - 1][iy] + NewNG.out[CurrSim - 1][iy] +
		NewTV.out[CurrSim - 1][iy];
	NewCTNGTV_M.out[CurrSim - 1][iy] = NewCTbySex[0] + NewNGbySex[0] + NewTVbySex[0];
	NewCTNGTV_F.out[CurrSim - 1][iy] = NewCTbySex[1] + NewNGbySex[1] + NewTVbySex[1];
}

void Pop::CalcStockCosts()
{
	int iy, it;
	double CumART, CurrART, CumDeaths, Treated;

	iy = CurrYear - StartYear;
	
	CumART = 0.0;
	CurrART = 0.0;
	CumDeaths = 0.0;
	for (it = 0; it < iy; it++){
		CumART += NewART200.out[CurrSim - 1][it] + NewART350.out[CurrSim - 1][it] +
			NewART500.out[CurrSim - 1][it] + NewART500plus.out[CurrSim - 1][it];
		CurrART += NewARTexp.out[CurrSim - 1][it] * (ARTresumption + ARTinterruption * exp(-
			(ARTresumption + ARTinterruption) * (iy - it + 0.5))) / (ARTresumption + ARTinterruption);
		CumDeaths += ARTdeathsExp.out[CurrSim - 1][it];
	}
	Treated = CurrART * (1.0 - CumDeaths / CumART);

	TotalCosts.out[CurrSim - 1][iy] = AnnARTcost * Treated;
	TotalCosts.out[CurrSim - 1][iy] += AnnCTXcost * DiagnosedUntreated.out[CurrSim - 1][iy] * 0.9;
	TotalCosts.out[CurrSim - 1][iy] += AnnPrEPcost * TakingPrEP.out[CurrSim - 1][iy];
}

void Pop::CalcFlowCosts()
{
	int iy;

	iy = CurrYear - StartYear;
	TotalCosts.out[CurrSim - 1][iy] += CondomCost * ProtSexActs.out[CurrSim - 1][iy] * 3.45; // Wastage factor
	TotalCosts.out[CurrSim - 1][iy] += MMCcost * MMCoperations.out[CurrSim - 1][iy];
	TotalCosts.out[CurrSim - 1][iy] += (PMTCTcost[0] + PMTCTcost[1]) * (TotBirthsHIV.out[CurrSim - 1][iy] -
		TotBirthsART.out[CurrSim - 1][iy]);
	TotalCosts.out[CurrSim - 1][iy] += (InfantTestCost[0] * 0.687 + InfantTestCost[1] * 0.92) * 
		TotBirthsHIV.out[CurrSim - 1][iy];
	TotalCosts.out[CurrSim - 1][iy] += PalliativeCareCost * AIDSdeaths.out[CurrSim - 1][iy];

	// Costs of baseline testing modalities
	TotalTestCosts.out[CurrSim - 1][iy] = BaseModalityCost[0][0] * PosTestsPartner.out[CurrSim - 1][iy] +
		BaseModalityCost[0][1] * (TotalTestsPartner.out[CurrSim - 1][iy] - PosTestsPartner.out[CurrSim - 1][iy]);
	TotalTestCosts.out[CurrSim - 1][iy] += BaseModalityCost[1][0] * PosTestsANC.out[CurrSim - 1][iy] +
		BaseModalityCost[1][1] * (TotalTestsANC.out[CurrSim - 1][iy] - PosTestsANC.out[CurrSim - 1][iy]) * (1.0 + ANCretestFreq[iy]);
	TotalTestCosts.out[CurrSim - 1][iy] += BaseModalityCost[2][0] * PosTestsOI.out[CurrSim - 1][iy] +
		BaseModalityCost[2][1] * (TotalTestsOI.out[CurrSim - 1][iy] - PosTestsOI.out[CurrSim - 1][iy]);
	TotalTestCosts.out[CurrSim - 1][iy] += BaseModalityCost[3][0] * PosTestsMMC.out[CurrSim - 1][iy] +
		BaseModalityCost[3][1] * (TotalTestsMMC.out[CurrSim - 1][iy] - PosTestsMMC.out[CurrSim - 1][iy]);
	TotalTestCosts.out[CurrSim - 1][iy] += BaseModalityCost[4][0] * PosTestsPrEP.out[CurrSim - 1][iy] +
		BaseModalityCost[4][1] * (TotalTestsPrEP.out[CurrSim - 1][iy] - PosTestsPrEP.out[CurrSim - 1][iy]);
	TotalTestCosts.out[CurrSim - 1][iy] += BaseModalityCost[5][0] * PosTestsPrison.out[CurrSim - 1][iy] +
		BaseModalityCost[5][1] * (TotalTestsPrison.out[CurrSim - 1][iy] - PosTestsPrison.out[CurrSim - 1][iy]);
	TotalTestCosts.out[CurrSim - 1][iy] += BaseModalityCost[6][0] * PosTestsGen.out[CurrSim - 1][iy] +
		BaseModalityCost[6][1] * (TotalTestsGen.out[CurrSim - 1][iy] - PosTestsGen.out[CurrSim - 1][iy]);
	TotalTestCosts.out[CurrSim - 1][iy] += BaseModalityCost[7][0] * PosTestsSTI.out[CurrSim - 1][iy] +
		BaseModalityCost[7][1] * (TotalTestsSTI.out[CurrSim - 1][iy] - PosTestsSTI.out[CurrSim - 1][iy]);

	// Costs of new testing modalities
	// 1. If there are community mobilization campaigns accompanying mobile testing, change the NewModalityCost arguments
	// from 2 and 3 to 7 in the 3rd and 4th commands respectively.
	// 2. If there is assisted partner notification, change the 1st command in the previous section to refer to
	// NewModalityCost[9] instead of BaseModalityCost[9], but ONLY if CurrYear >= intervention start year (2019).
	// 3. If there is an offer of self-testing accompanying home-based testing, change the cost assumptions in the first
	// 2 rows of the NewTestingModality array to the values in cells T16:T17 of the 'ST-HBT' sheet in Craig's workbook.
	// 4. If there is an offer of self-testing accompanying antenatal partner testing, change the cost assumptions in the 2nd-last
	// row of the BaseTestingModality array to the values in cells T16:T17 of the 'ST-APN-ANC' sheet in Craig's workbook.
	TotalTestCosts.out[CurrSim - 1][iy] += NewModalityCost[0][0] * PosTestsHH_U.out[CurrSim - 1][iy] +
		NewModalityCost[0][1] * (TotalTestsHH_U.out[CurrSim - 1][iy] - PosTestsHH_U.out[CurrSim - 1][iy]);
	TotalTestCosts.out[CurrSim - 1][iy] += NewModalityCost[1][0] * PosTestsHH_R.out[CurrSim - 1][iy] +
		NewModalityCost[1][1] * (TotalTestsHH_R.out[CurrSim - 1][iy] - PosTestsHH_R.out[CurrSim - 1][iy]);
	TotalTestCosts.out[CurrSim - 1][iy] += NewModalityCost[2][0] * PosTestsMobile_U.out[CurrSim - 1][iy] +
		NewModalityCost[2][1] * (TotalTestsMobile_U.out[CurrSim - 1][iy] - PosTestsMobile_U.out[CurrSim - 1][iy]);
	TotalTestCosts.out[CurrSim - 1][iy] += NewModalityCost[3][0] * PosTestsMobile_R.out[CurrSim - 1][iy] +
		NewModalityCost[3][1] * (TotalTestsMobile_R.out[CurrSim - 1][iy] - PosTestsMobile_R.out[CurrSim - 1][iy]);
	TotalTestCosts.out[CurrSim - 1][iy] += NewModalityCost[4][0] * PosTestsMSM.out[CurrSim - 1][iy] +
		NewModalityCost[4][1] * (TotalTestsMSM.out[CurrSim - 1][iy] - PosTestsMSM.out[CurrSim - 1][iy]);
	TotalTestCosts.out[CurrSim - 1][iy] += NewModalityCost[5][0] * PosTestsFSW.out[CurrSim - 1][iy] +
		NewModalityCost[5][1] * (TotalTestsFSW.out[CurrSim - 1][iy] - PosTestsFSW.out[CurrSim - 1][iy]);
	TotalTestCosts.out[CurrSim - 1][iy] += NewModalityCost[6][0] * PosTestsWork.out[CurrSim - 1][iy] +
		NewModalityCost[6][1] * (TotalTestsWork.out[CurrSim - 1][iy] - PosTestsWork.out[CurrSim - 1][iy]);
	TotalTestCosts.out[CurrSim - 1][iy] += NewModalityCost[8][0] * PosTestsSchool.out[CurrSim - 1][iy] +
		NewModalityCost[8][1] * (TotalTestsSchool.out[CurrSim - 1][iy] - PosTestsSchool.out[CurrSim - 1][iy]);
	TotalTestCosts.out[CurrSim - 1][iy] += NewModalityCost[10][0] * PosTestsANCpartner0.out[CurrSim - 1][iy] +
		NewModalityCost[10][1] * (TotalTestsANCpartner0.out[CurrSim - 1][iy] - PosTestsANCpartner0.out[CurrSim - 1][iy]);
	TotalTestCosts.out[CurrSim - 1][iy] += NewModalityCost[10][0] * PosTestsANCpartner1.out[CurrSim - 1][iy] +
		NewModalityCost[10][1] * (TotalTestsANCpartner1.out[CurrSim - 1][iy] - PosTestsANCpartner1.out[CurrSim - 1][iy]);
	TotalTestCosts.out[CurrSim - 1][iy] += NewModalityCost[11][0] * PosTestsFPC.out[CurrSim - 1][iy] +
		NewModalityCost[11][1] * (TotalTestsFPC.out[CurrSim - 1][iy] - PosTestsFPC.out[CurrSim - 1][iy]);

	TotalCosts.out[CurrSim - 1][iy] += TotalTestCosts.out[CurrSim - 1][iy];
}

void Pop::ResetStructOutcomes()
{
	int ic;
	double duration;

	duration = 1.0 * BehavCycleCount / CycleS;

	if (duration == 0.0 || duration == 0.5){
		for (ic = 0; ic < Register.size(); ic++){
			Register[ic].CumSTIs = 0;
		}
	}
	if (duration == 0.0){
		for (ic = 0; ic < Register.size(); ic++){
			Register[ic].CumCasual = 0;
		}
	}
}

void Pop::UpdateConception(int month)
{
	int ic, ir, ia, ind, is, ii, s, ID1;
	double gestation, x, y, a, b, p, q, bound;

	//int seed = 8821 + CurrYear * 33 + BehavCycleCount * 71;
	//if(CurrYear>=StartYear+FixedPeriod){
	//	seed += CurrSim;}
	//CRandomMersenne rg(seed);
	for (ic = 0; ic<MaxPop; ic++){
		rbirth[ic] = rg.Random();
		rconcep[ic] = rg.Random();
	}
	for (ic = 0; ic<Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 1 &&
			Register[ic].CurrAge >= 14 && Register[ic].CurrAge<49){
			ir = Register[ic].Race;
			ia = Register[ic].CurrAge - 14;
			Register[ic].NonHIVfertRate = 0.0;
			if (Register[ic].DOLB<(1.0*month / 12.0 + CurrYear + 0.5) && 
				Register[ic].DOLW<(1.0*month / 12.0 + CurrYear + 0.5)){
				if (Register[ic].CurrPartners>0){
					Register[ic].NonHIVfertRate += 1.0 - Register[ic].CondomPrimary * CondomEffPreg;
					ID1 = Register[ic].IDprimary;
				}
				if (Register[ic].CurrPartners == 2){
					Register[ic].NonHIVfertRate += 1.0 - Register[ic].Condom2ndary * CondomEffPreg;
					if (Register[ID1 - 1].ChildIDs[19] > 0){ ID1 = Register[ic].ID2ndary; }
				}
				if (Register[ic].FSWind == 1){
					Register[ic].NonHIVfertRate += 1.0 - CondomUseFSW * CondomEffPreg;}
				if (Register[ic].CurrContr>0){
					ii = Register[ic].CurrContr;
					Register[ic].NonHIVfertRate *= (1.0 - ContrEffPreg[ii - 1]);
				}
				Register[ic].NonHIVfertRate *= Register[ic].Fecundability * 
					SexuallyExpFert[ia][ir] / 12.0;
				is = Register[ic].HIVstage;
				if (is > 0){ Register[ic].NonHIVfertRate *= RelHIVfertility[is - 1]; }
				if (Register[ic].ChildIDs[19] > 0){
					Register[ic].NonHIVfertRate = 0.0; // Not allowing >20 children
				}
				if (Register[ic].CurrPartners > 0){
					if (Register[ID1 - 1].ChildIDs[19] > 0){
						Register[ic].NonHIVfertRate = 0.0; // Not allowing >20 children
					} 
				}
				if (rconcep[ic] < Register[ic].NonHIVfertRate){
					Register[ic].DOLC = 1.0 * CurrYear + 0.5 + 1.0 * month / 12.0 +
						(rconcep[ic] / Register[ic].NonHIVfertRate) / 12.0;
					if (StructuralRCTcalib == 1 && CurrYear >= BaselineStart && Register[ic].CurrAge < 20){
						Register[ic].CumTeenPreg = 1;
					}
					if (Register[ic].CurrPartners > 0){ Register[ic].FatherBaby = ID1; }
					else{ Register[ic].FatherBaby = 0; }
					Register[ic].PrevContr = Register[ic].CurrContr;
					if (Register[ic].CurrContr < 3){ Register[ic].CurrContr = 0; }
					bound = 0.0;
					ind = 2;
					a = MeanGestation;
					b = SDgestation;
					s = 0;
					p = rbirth[ic];
					q = 1.0 - rbirth[ic];
					cdfnor(&ind, &p, &q, &x, &a, &b, &s, &bound);
					gestation = x;
					double Interval = -Register[ic].DOLB;
					Register[ic].DOLB = Register[ic].DOLC + gestation;
					if (CurrYear >= 1993 && CurrYear < 1998 && FixedUncertainty == 1){
						Interval += Register[ic].DOLB;
						if (Interval < 1.5){ BirthIntervals.out[CurrSim - 1][ir * 5] += 1; }
						else if (Interval < 2.0){ BirthIntervals.out[CurrSim - 1][ir * 5 + 1] += 1; }
						else if (Interval < 3.0){ BirthIntervals.out[CurrSim - 1][ir * 5 + 2] += 1; }
						else if (Interval < 4.0){ BirthIntervals.out[CurrSim - 1][ir * 5 + 3] += 1; }
						else { BirthIntervals.out[CurrSim - 1][ir * 5 + 4] += 1; }
					}
				}
			}
		}
	}
}

void Pop::UpdateContraception(int month)
{
	int ic, ia, ir, iy, PID;
	double SterilProb, StartRate, BaseProb, AdjProb, CumOR;

	for (ic = 0; ic < Register.size(); ic++){ 
		r2[ic] = rg.Random(); 
		revent[ic] = rg.Random();
	}
	if (CurrYear <= 1996){ NoPrevUseContr = NoPrevUseContr1; }
	else if (CurrYear >= 2000){ NoPrevUseContr = NoPrevUseContr2; }
	else{ NoPrevUseContr = NoPrevUseContr1 + (NoPrevUseContr2 - 
		NoPrevUseContr1) * (CurrYear - 1996)/4.0; }
	iy = CurrYear - StartYear;

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 1 &&
			Register[ic].AgeGroup >= 2){
			if (Register[ic].CurrContr == 0){
				// Check if woman starts hormonal contraception postnatally
				// We check whether weaning occurred in the LAST month, since the date of
				// weaning will not yet have been calculated for births in the current month.
				if ((Register[ic].DOLW - Register[ic].DOLB) >= 0.5 &&
					(Register[ic].DOLB + 0.5) > 1.0 * (CurrYear + 0.5 + (month-1)/12.0) &&
					(Register[ic].DOLB + 0.5) < 1.0 * (CurrYear + 0.5 + month/12.0)){
					Register[ic].TestStartContr(3, 0);
				}
				else if ((Register[ic].DOLW - Register[ic].DOLB) < 0.5 &&
					Register[ic].DOLW > 1.0 * (CurrYear + 0.5 + (month-1)/12.0) &&
					Register[ic].DOLW < 1.0 * (CurrYear + 0.5 + month/12.0)){
					Register[ic].TestStartContr(3, 0);
				}
				// Test if woman starts hormonal contraception for other reasons
				if (Register[ic].DOLB < (CurrYear + 0.5 + month / 12.0) &&
					Register[ic].CurrContr == 0 && Register[ic].AgeGroup < 10){
					StartRate = StartContrOther * pow(Register[ic].Fecundability, 0.5) / 12.0;
						ia = Register[ic].AgeGroup;
						if (ia>3){ StartRate *= AgeEffectContr[ia - 3][0]; }
						AdjProb = 1.0 - exp(-StartRate);
						if (r2[ic] < AdjProb){
							// Start hormonal contraception - select method
							ir = Register[ic].Race;
							BaseProb = InjectablePref[ir];
							CumOR = pow(ORinjectableEdu, Register[ic].HighestGrade - 10);
							//if (CurrYear >= 2000){ CumOR = pow(ORinjectableEdu, 3.0); }
							CumOR *= pow(ORinjectableAge, 1.0*(Register[ic].CurrAge - 25.0)/5.0);
							if (CurrYear<1995){
								CumOR *= ORinjectableInit + (1.0 - ORinjectableInit) *
									(CurrYear - StartYear) / 10.0;
							}
							if (Register[ic].DOLB>0){ CumOR *= 1.0 / ORinjectableNeverPreg; }
							AdjProb = 1.0 / (1.0 - (1.0 - 1.0 / BaseProb) / CumOR);
							if (Register[ic].PrevContr == 1){
								AdjProb = 1.0 - (1.0 - AdjProb) * NewContrMethodWeight;}
							if (Register[ic].PrevContr == 2){
								AdjProb = AdjProb * NewContrMethodWeight;}
							if (revent[ic] < AdjProb){ Register[ic].CurrContr = 1; }
							else{ Register[ic].CurrContr = 2; }
						}
				}
			}
			else if (Register[ic].CurrContr<3){
				// Test if woman stops using hormonal contraception
				if (Register[ic].CurrPartners>0 || Register[ic].FSWind == 1){
					BaseProb = 1.0 - exp(-StopContr[1]);
					PID = Register[ic].IDprimary; 
					if (Register[ic].CurrPartners>0 && (Register[PID - 1].Imprisoned>0 ||
						Register[ic].CurrUrban != Register[PID - 1].CurrUrban)){
						if (Register[ic].CurrPartners == 2){ PID = Register[ic].ID2ndary; }
						if (Register[ic].CurrPartners == 1 || (Register[ic].CurrPartners == 2 &&
							(Register[ic].CurrUrban != Register[PID - 1].CurrUrban ||
							Register[PID - 1].Imprisoned>0))){
							BaseProb = 1.0 - exp(-StopContrAbsentPartner);
						}
					}
				}
				else{ BaseProb = 1.0 - exp(-StopContr[0]); }
				if (r2[ic] < BaseProb){
					Register[ic].PrevContr = Register[ic].CurrContr;
					Register[ic].CurrContr = 0;
				}
			}
			// Test if woman gets sterilized
			if (Register[ic].DOLB<(CurrYear + 0.5 + month / 12.0) &&
				Register[ic].DOLB > 0.0 && Register[ic].AgeGroup < 10 &&
				Register[ic].CurrContr<3 && Register[ic].Fecundability>0.0 &&
				(Register[ic].CurrPartners>0 || Register[ic].FSWind==1)){
				ia = Register[ic].AgeGroup;
				ir = Register[ic].Race;
				SterilProb = 1.0 - exp(-SterilizationRates[ia - 3] *
					RaceEffectSteril[ir] * RRsterilizationByYear[iy] / 12.0);
				if (r2[ic] < SterilProb){ Register[ic].CurrContr = 3; }
			}
			if (Register[ic].CurrContr == 1){ Register[ic].EverInjectable = 1; }
			if (Register[ic].CurrContr == 2){ Register[ic].EverPill = 1; }
		}
	}
}

void Pop::UpdateDisclosure()
{
	int ic;

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].VCThistory == 2){
			if (Register[ic].CurrPartners>0 && Register[ic].DisclosedPrimary == 0){
				Register[ic].TestDisclose(ic + 1, 1, 1);}
			if (Register[ic].CurrPartners==2 && Register[ic].Disclosed2ndary == 0){
				Register[ic].TestDisclose(ic + 1, 2, 1);}
		}
	}
}

void Pop::UpdateIncarceration(int month)
{
	int ic, ia, ir, ic2, ij, ind, iy, Homeless, HHID;
	double BaseProb, a, b, p, q, x, p2;

	for (ic = 0; ic < Register.size(); ic++){
		r2[ic] = rg.Random();
		revent[ic] = rg.Random();
	}
	iy = CurrYear - StartYear;

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 0 && Register[ic].AgeGroup>2){
			if (Register[ic].Imprisoned==0){
				ia = Register[ic].AgeGroup - 3;
				BaseProb = PrisonEntryRate[ia];
				if (Register[ic].HighestGrade >= 12){ BaseProb *= PrisonEntryEdu; }
				ir = Register[ic].Race;
				BaseProb *= PrisonEntryRace[ir];
				if (Register[ic].ReleaseDate > 0.0){ BaseProb *= PrisonReentryAdj; }
				BaseProb = 1.0 - exp(-BaseProb / 12.0);
				if (Register[ic].RiskGroup == 1 && Register[ic].VirginInd == 0){
					BaseProb *= PrisonEntryRisk[0];}
				else{ BaseProb *= PrisonEntryRisk[1]; }
				if (r2[ic] < BaseProb){
					Register[ic].Imprisoned = 1; 
					if (Register[ic].CasualInd == 1){
						Register[ic].CasualInd = 0;
						TotCurrCasual = TotCurrCasual - 1;
						// Update CasualRegister
						for (ic2 = 0; ic2<MaxCasual; ic2++){
							if (CasualRegister[ic2] == ic+1){
								for (ij = ic2; ij<MaxCasual - 1; ij++){
									CasualRegister[ij] = CasualRegister[ij + 1]; }
								CasualRegister[MaxCasual - 1] = 0;
							}
						}
					}
					if (Register[ic].HetCasualInd == 1){
						Register[ic].HetCasualInd = 0;
						TotCurrCasualHet[ir][0] = TotCurrCasualHet[ir][0] - 1;
						// Update CasualRegisterM
						for (ic2 = 0; ic2<MaxCasual; ic2++){
							if (CasualRegisterM[ic2][ir] == ic + 1){
								for (ij = ic2; ij<MaxCasual - 1; ij++){
									CasualRegisterM[ij][ir] = CasualRegisterM[ij + 1][ir];
								}
								CasualRegisterM[MaxCasual - 1][ir] = 0;
							}
						}
					}
				}
			}
			else if (Register[ic].Imprisoned == 1){
				BaseProb = 1.0 - exp(-1.0 / (MeanUnsentencedDur * 12.0));
				if (r2[ic] < BaseProb){
					if (revent[ic] < UnsentencedRelease){
						Register[ic].Imprisoned = 0;
						Register[ic].PrevImprisoned = 1;
						Register[ic].ReleaseDate = 1.0 * CurrYear + 0.5 + 1.0 * (month/12.0);
						HHID = Register[ic].HouseholdID;
						Homeless = Register[ic].TestHomeless(ic + 1);
						if (Homeless == 1 && HHID > 0){ HHregister[HHID - 1].RemoveMember(ic + 1); }
					}
					else{
						Register[ic].Imprisoned = 2;
						ind = 2;
						a = SentenceLength[0];
						b = SentenceLength[1];
						p = (revent[ic] - UnsentencedRelease) / (1.0 - UnsentencedRelease);
						q = 1 - p;
						cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
						p2 = r2[ic] / BaseProb;
						if (p2 < ProbParole){ x *= 0.5; }
						if (x < 1.0 / 12.0){ x = 1.0 / 12.0; }
						Register[ic].ReleaseDate = 1.0 * CurrYear + 0.5 + 1.0 * month / 12.0 + x;
					}
				}
			}
			else {
				if (Register[ic].ReleaseDate >= (1.0 * CurrYear + 0.5 + 1.0 * month / 12.0) &&
					Register[ic].ReleaseDate < (1.0 * CurrYear + 0.5 + 1.0 * (month + 1) / 12.0)){
					Register[ic].Imprisoned = 0;
					Register[ic].PrevImprisoned = 1;
					HHID = Register[ic].HouseholdID;
					Homeless = Register[ic].TestHomeless(ic + 1);
					if (Homeless == 1 && HHID > 0){ HHregister[HHID - 1].RemoveMember(ic + 1); }
				}
			}
		}
	}
}

void Pop::UpdatePrEP(int month)
{
	int ic, PID, PID2;
	double PrEPstart, PrEPstop, PrEPstartMSM;

	for (ic = 0; ic < Register.size(); ic++){
		r2[ic] = rg.Random();}
	PrEPstart = 1.0 - exp(-AnnPrEPuptake[CurrYear-StartYear]/12.0);
	PrEPstop = 1.0 - exp(-PrEPdiscontinue / 12.0);
	PrEPstartMSM = 1.0 - exp(-PrEPuptakeMSM[CurrYear - StartYear] / 12.0);

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].OnPrEP == 0.0){
			// Test if individual starts PrEP
			//if (Register[ic].CurrPartners>0 && Register[ic].DOLC<(CurrYear + 0.5 + 1.0 * month/12.0) && 
			//	(Register[ic].DOLB>(CurrYear + 0.5 + 1.0 * month / 12.0) || 
			//	Register[ic].DOLW>(CurrYear + 0.5 + 1.0 * month / 12.0))){
			//if (Register[ic].CurrContr == 1 && (Register[ic].CurrPartners>0 || Register[ic].FSWind == 1)){
			//if ((Register[ic].CurrPartners>0 || Register[ic].FSWind == 1 || Register[ic].CasualInd == 1) &&
			//	Register[ic].AgeGroup<5 && Register[ic].SexInd == 0){
			//if (Register[ic].InSchool == 0 && Register[ic].HighestGrade<12 && (Register[ic].CurrPartners>0 ||
			//	Register[ic].FSWind == 1 || Register[ic].CasualInd == 1) && Register[ic].AgeGroup<4){
			//if (Register[ic].InSchool == 1 && Register[ic].HighestGrade<12 && (Register[ic].CurrPartners>0 ||
			//	Register[ic].FSWind == 1 || Register[ic].CasualInd == 1)){
			//if (Register[ic].DateMig>(CurrYear - 1.5 + 1.0 * month/12.0)){
			//if (Register[ic].CurrPartners>0 && Register[PID-1].HIVstage>0 && Register[PID-1].DisclosedPrimary==1){
			if (Register[ic].FSWind == 1){
				if (r2[ic]<PrEPstart){ 
					TotalTestsPrEP.out[CurrSim - 1][CurrYear - StartYear] += Register[ic].PopWeight;
					if (Register[ic].HIVstage < 2){ // Assuming acute infection wouldn't be detectable
						Register[ic].OnPrEP = 1.0; 
						Register[ic].VCThistory = 1;
					}
					else if(Register[ic].VCThistory<2){ // Newly diagnosed
						Register[ic].GetDiagnosed(ic + 1, 0);
						PosTestsPrEP.out[CurrSim - 1][CurrYear - StartYear] += Register[ic].PopWeight;
					}
					if (Register[ic].HIVstage == 1){ AcuteTests.out[CurrSim - 1][CurrYear - StartYear] += Register[ic].PopWeight; }
				}
			}
			PID = Register[ic].IDprimary;
			PID2 = Register[ic].ID2ndary;
			if (Register[ic].SexInd == 0 && Register[ic].MalePref>0.0 && (Register[ic].CasualInd == 1 ||
				(Register[ic].CurrPartners > 0 && Register[PID - 1].SexInd == 0) ||
				(Register[ic].CurrPartners == 2 && Register[PID2 - 1].SexInd == 0))){
				if (r2[ic]<PrEPstartMSM){
					TotalTestsPrEP.out[CurrSim - 1][CurrYear - StartYear] += Register[ic].PopWeight;
					if (Register[ic].HIVstage < 2){ // Assuming acute infection wouldn't be detectable
						Register[ic].OnPrEP = 1.0;
						Register[ic].VCThistory = 1;
					}
					else if (Register[ic].VCThistory<2){ // Newly diagnosed
						Register[ic].GetDiagnosed(ic + 1, 0);
						PosTestsPrEP.out[CurrSim - 1][CurrYear - StartYear] += Register[ic].PopWeight;
					}
					if (Register[ic].HIVstage == 1){ AcuteTests.out[CurrSim - 1][CurrYear - StartYear] += Register[ic].PopWeight; }
				}
			}
		}
		else if (Register[ic].AliveInd == 1 && Register[ic].OnPrEP > 0.0){
			// Test if individual stops PrEP
			//if (Register[ic].CurrPartners>0 && Register[ic].DOLC<(CurrYear + 0.5 + 1.0 * month / 12.0) &&
			//	(Register[ic].DOLB>(CurrYear + 0.5 + 1.0 * month / 12.0) ||
			//	Register[ic].DOLW>(CurrYear + 0.5 + 1.0 * month / 12.0))){
			//if (Register[ic].CurrContr == 1 && (Register[ic].CurrPartners>0 || Register[ic].FSWind == 1)){
			//if ((Register[ic].CurrPartners>0 || Register[ic].FSWind == 1 || Register[ic].CasualInd == 1) && 
			//	Register[ic].AgeGroup<5 && Register[ic].SexInd==0){
			//if (Register[ic].InSchool == 0 && Register[ic].HighestGrade<12 && (Register[ic].CurrPartners>0 ||
			//	Register[ic].FSWind == 1 || Register[ic].CasualInd == 1) && Register[ic].AgeGroup<4){
			//if (Register[ic].InSchool == 1 && Register[ic].HighestGrade<12 && (Register[ic].CurrPartners>0 ||
			//	Register[ic].FSWind == 1 || Register[ic].CasualInd == 1)){
			//if (Register[ic].DateMig>(CurrYear - 1.5 + 1.0 * month / 12.0)){
			PID = Register[ic].IDprimary;
			PID2 = Register[ic].ID2ndary;
			//if (Register[ic].SexInd == 0 && Register[ic].MalePref>0.0 && (Register[ic].CasualInd == 1 ||
			//	(Register[ic].CurrPartners>0 && Register[PID - 1].SexInd == 0) ||
			//	(Register[ic].CurrPartners == 2 && Register[PID2 - 1].SexInd == 0))){
			//if (Register[ic].CurrPartners>0 && Register[PID - 1].HIVstage>0 && Register[PID - 1].DisclosedPrimary == 1){
			if (Register[ic].FSWind == 1 || (Register[ic].SexInd == 0 && Register[ic].MalePref>0.0 && 
				(Register[ic].CasualInd == 1 || (Register[ic].CurrPartners>0 && Register[PID - 1].SexInd == 0) ||
				(Register[ic].CurrPartners == 2 && Register[PID2 - 1].SexInd == 0)))){
				if (r2[ic]<PrEPstop){ Register[ic].OnPrEP = 0.0; }
				else if (Register[ic].HIVstage>1){ // Newly diagnosed
					Register[ic].OnPrEP = 0.0;
					Register[ic].GetDiagnosed(ic + 1, 0);
					PosTestsPrEP.out[CurrSim - 1][CurrYear - StartYear] += Register[ic].PopWeight;
				}
				else{ 
					TotalTestsPrEP.out[CurrSim - 1][CurrYear - StartYear] += 0.25 * Register[ic].PopWeight;
					if (Register[ic].HIVstage == 1){ AcuteTests.out[CurrSim - 1][CurrYear - StartYear] += Register[ic].PopWeight; }
				}
			}
			else{ Register[ic].OnPrEP = 0.0; }
		}
	}
}

void Pop::UpdateVaccination(int month)
{
	int ic, TotDoses, PID1, PID2;
	double VaccineStart, MonthsSince1stVacc;

	for (ic = 0; ic < Register.size(); ic++){
		r2[ic] = rg.Random();}
	VaccineStart = 1.0 - exp(-AnnVaccineUptake[CurrYear - StartYear] / 12.0);

	TotDoses = 0;
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].Date1stVacc == 0.0 && Register[ic].VCThistory<2){
			// Test if individual receives 1st vaccination
			if (Register[ic].AgeGroup>1 && Register[ic].AgeGroup<4){
			//if (Register[ic].InSchool == 0 && Register[ic].HighestGrade<12 && Register[ic].AgeGroup<4){
			//if (Register[ic].DateMig>(CurrYear - 1.5 + 1.0 * month / 12.0)){
			//if (Register[ic].DOLC<(CurrYear + 0.5 + 1.0 * month / 12.0) &&
			//	(Register[ic].DOLB>(CurrYear + 0.5 + 1.0 * month / 12.0) ||
			//	Register[ic].DOLW>(CurrYear + 0.5 + 1.0 * month / 12.0))){
			//if (Register[ic].SexInd == 1 && Register[ic].CurrContr==1){
			//if (Register[ic].InSchool == 1 && Register[ic].AgeGroup>2 && Register[ic].HighestGrade<12){
			//if (Register[ic].CurrPartners>0){ PID1 = Register[ic].IDprimary; }
			//if (Register[ic].CurrPartners==2){ PID2 = Register[ic].ID2ndary; }
			//if (Register[ic].SexInd == 0 && Register[ic].MalePref>0.0 && (Register[ic].CasualInd == 1 ||
			//	(Register[ic].CurrPartners>0 && Register[PID1 - 1].SexInd == 0) ||
			//	(Register[ic].CurrPartners == 2 && Register[PID2 - 1].SexInd == 0))){
			//if (Register[ic].CurrPartners>0 && Register[PID1-1].HIVstage>0 && Register[PID1-1].DisclosedPrimary == 1){
			//if (Register[ic].FSWind == 1){
			//if (Register[ic].AgeGroup >= 3 && Register[ic].AgeGroup<5 && Register[ic].SexInd==0){
				if (r2[ic]<VaccineStart){
					TotDoses += 1;
					Register[ic].VaccineEff = HIVvaccProt[0];
					Register[ic].Date1stVacc = 0.5 + CurrYear + 1.0 * month/12.0;
				}
			}
		}
		else if (Register[ic].AliveInd == 1 && Register[ic].VCThistory<2){
			Register[ic].VaccineEff = Register[ic].VaccineEff * exp(-VaccineWane);
			// Test if individual gets revaccinated or receives follow-up dose
			if (Register[ic].AgeGroup>1 && Register[ic].AgeGroup<4){
			//if (Register[ic].InSchool == 0 && Register[ic].HighestGrade<12 && Register[ic].AgeGroup<4){
			//if (Register[ic].DateMig>(CurrYear - 1.5 + 1.0 * month / 12.0)){
			//if (Register[ic].DOLC<(CurrYear + 0.5 + 1.0 * month / 12.0) &&
			//	(Register[ic].DOLB>(CurrYear + 0.5 + 1.0 * month / 12.0) ||
			//	Register[ic].DOLW>(CurrYear + 0.5 + 1.0 * month / 12.0))){
			//if (Register[ic].SexInd == 1 && Register[ic].CurrContr == 1){
			//if (Register[ic].InSchool == 1 && Register[ic].AgeGroup>2 && Register[ic].HighestGrade<12){
			//if (Register[ic].CurrPartners>0){ PID1 = Register[ic].IDprimary; }
			//if (Register[ic].CurrPartners==2){ PID2 = Register[ic].ID2ndary; }
			//if (Register[ic].SexInd == 0 && Register[ic].MalePref>0.0 && (Register[ic].CasualInd == 1 ||
			//	(Register[ic].CurrPartners>0 && Register[PID1 - 1].SexInd == 0) ||
			//	(Register[ic].CurrPartners == 2 && Register[PID2 - 1].SexInd == 0))){
			//if (Register[ic].CurrPartners>0 && Register[PID1 - 1].HIVstage>0 && Register[PID1 - 1].DisclosedPrimary == 1){
			//if (Register[ic].FSWind == 1){
			//if (Register[ic].AgeGroup >= 3 && Register[ic].AgeGroup < 5 && Register[ic].SexInd == 0){
				MonthsSince1stVacc = month + 12.0 * (0.5 + CurrYear -
					Register[ic].Date1stVacc);
				if (MonthsSince1stVacc == 1){
					if (r2[ic] < VaccAdherence){
						TotDoses += 1;
						Register[ic].VaccineEff = 1.0 - (1.0 - Register[ic].VaccineEff) *
							(1.0 - HIVvaccProt[0]);
					}
				}
				else if (MonthsSince1stVacc == 3 || MonthsSince1stVacc == 6 ||
					MonthsSince1stVacc == 12 || MonthsSince1stVacc == 36){
					// Add code if simulating vaccine impact for >5 years
					if (r2[ic] < VaccAdherence){
						TotDoses += 1;
						Register[ic].VaccineEff = 1.0 - (1.0 - Register[ic].VaccineEff) *
							(1.0 - HIVvaccProt[1]);
					}
				}
			}
		}
	}

	if (month == 0){ VaccineDoses.out[CurrSim - 1][CurrYear - StartYear] = TotDoses; }
	else{ VaccineDoses.out[CurrSim - 1][CurrYear - StartYear] += TotDoses; }
}

void Pop::StructInterventionWaning(int month)
{
	int ic, ig;
	double AlcAdj, NormAdj, TempAdj, TempAdj2, CurrTime, BaseDrinksPerDD, BaseDailyDrinkProb, BaseIneqGender;
	double RRcasualPostInt;

	CurrTime = 1.0 * month / 12.0 + CurrYear - BaselineStart;
	if (CurrTime == 0){
		// Check if people drop out of casual sex due to the intervention
		for (ic = 0; ic < Register.size(); ic++){
			r2[ic] = rg.Random();
		}
		if (StructIntScenario == 1){ AlcAdj = pow(RRalcoholSingle, 0.5); }
		else if (StructIntScenario == 2){ AlcAdj = pow(RRalcoholMultiple, 0.5); }
		else{ AlcAdj = 1.0; }
		if (StructIntScenario == 6){ NormAdj = RRgenderIneqComm; }
		else if (StructIntScenario == 7){ NormAdj = RRgenderIneqIndiv; }
		else{ NormAdj = 1.0; }
		for (ic = 0; ic < Register.size(); ic++){
			if ((Register[ic].HetCasualInd == 1 || Register[ic].CasualInd == 1) && Register[ic].AliveInd == 1){
				RRcasualPostInt = 1.0;
				if (Register[ic].SexInd == 0 && Register[ic].HetCasualInd == 1 && ((Register[ic].CurrAge >= 18 &&
					Register[ic].CurrAge < 50 && StructIntScenario == 6) || (Register[ic].CurrAge >= 15 &&
					Register[ic].CurrAge < 30 && StructIntScenario == 7)) && Register[ic].AliveInd == 1){
					// Inequitable norm adjustment
					RRcasualPostInt = pow(EffectIneqGenderCasual, Register[ic].IneqGender * (1.0 - NormAdj) / 0.1);
				}
				if (Register[ic].CurrAge >= 15 && Register[ic].DailyDrinkProb >= 1.0 / 30.0 &&
					Register[ic].DrinksPerDD >= 5.0){
					// Binge drinking adjustment
					ig = Register[ic].SexInd;
					if (Register[ic].DailyDrinkProb * pow(AlcAdj, 2.0) < 1.0 / 30.0){
						RRcasualPostInt *= 1.0 / RRcasualBinge[ig];
					}
					else if (Register[ic].DrinksPerDD * AlcAdj < 5.0){
						RRcasualPostInt *= 1.0 / RRcasualBinge[ig];
					}
				}
				if (RRcasualPostInt < 1.0){
					if (r2[ic] > RRcasualPostInt){
						Register[ic].HetCasualInd = 0;
						Register[ic].CasualInd = 0;
					}
				}
			}
		}
		// Then check if people cease unprotected sex due to the intervention
		for (ic = 0; ic < Register.size(); ic++){
			r2[ic] = rg.Random();
		}
		for (ic = 0; ic < Register.size(); ic++){
			if ((Register[ic].CondomPrimary == 0 || Register[ic].Condom2ndary == 0) &&
				Register[ic].AliveInd == 1 && Register[ic].CurrAge >= 15){
				Register[ic].TestCondomChange(ic + 1);
			}
		}
	}

	if (StructIntScenario == 1 || StructIntScenario == 2){
		if (StructIntScenario == 1){ AlcAdj = pow(RRalcoholSingle, 0.5); }
		if (StructIntScenario == 2){ AlcAdj = pow(RRalcoholMultiple, 0.5); }
		if (month > 0.0){
			TempAdj = (1.0 - (1.0 - AlcAdj) * pow(1.0 - ProbAlcoholReversion, CurrTime)) /
				(1.0 - (1.0 - AlcAdj) * pow(1.0 - ProbAlcoholReversion, CurrTime - 1.0/12.0));
			TempAdj2 = (1.0 - (1.0 - pow(AlcAdj, 2.0)) * pow(1.0 - ProbAlcoholReversion, CurrTime)) /
				(1.0 - (1.0 - pow(AlcAdj, 2.0)) * pow(1.0 - ProbAlcoholReversion, CurrTime - 1.0 / 12.0));
		}
		else if (CurrYear > BaselineStart){
			// In month 0 the function is being called after alcohol levels have just
			// been reset to baseline levels (UpdateAlcohol).
			TempAdj = (1.0 - (1.0 - AlcAdj) * pow(1.0 - ProbAlcoholReversion, CurrTime));
			TempAdj2 = (1.0 - (1.0 - pow(AlcAdj, 2.0)) * pow(1.0 - ProbAlcoholReversion, CurrTime));
		}
		else{ 
			TempAdj = AlcAdj; 
			TempAdj2 = pow(AlcAdj, 2.0);
		}
		for (ic = 0; ic < Register.size(); ic++){
			if (Register[ic].BaselineBinge == 1 && Register[ic].AliveInd == 1){
				BaseDailyDrinkProb = Register[ic].DailyDrinkProb;
				Register[ic].DailyDrinkProb *= TempAdj2;
				BaseDrinksPerDD = Register[ic].DrinksPerDD;
				Register[ic].DrinksPerDD *= TempAdj;
				if (Register[ic].DrinksPerDD >= 5.0 && BaseDrinksPerDD < 5.0){
					Register[ic].TempCondomAdj *= pow(ORcondomBingePW, Register[ic].DailyDrinkProb * 7.0);
				}
				else if (BaseDrinksPerDD >= 5.0 && Register[ic].DrinksPerDD < 5.0){
					Register[ic].TempCondomAdj *= 1.0 / pow(ORcondomBingePW, BaseDailyDrinkProb * 7.0);
				}
				else{
					Register[ic].TempCondomAdj *= pow(ORcondomBingePW, (Register[ic].DailyDrinkProb -
						BaseDailyDrinkProb) * 7.0);
				}
			}
		}
	}

	if (StructIntScenario == 6 || StructIntScenario == 7){
		if (StructIntScenario == 6){ NormAdj = RRgenderIneqComm; }
		if (StructIntScenario == 7){ NormAdj = RRgenderIneqIndiv; }
		if (month > 0.0 || CurrYear > BaselineStart){
			TempAdj = (1.0 - (1.0 - NormAdj) * pow(1.0 - ProbGenderIneqReversion, CurrTime)) /
				(1.0 - (1.0 - NormAdj) * pow(1.0 - ProbGenderIneqReversion, CurrTime - 1.0 / 12.0));
		}
		else{ TempAdj = NormAdj; }
		for (ic = 0; ic < Register.size(); ic++){
			if (Register[ic].SexInd == 0 && Register[ic].AliveInd == 1 && ((Register[ic].CurrAge >= 18 &&
				Register[ic].CurrAge < 50 && StructIntScenario == 6) || (Register[ic].CurrAge >= 15 &&
				Register[ic].CurrAge < 30 && StructIntScenario == 7))){
				// Adjust drinking level first, since it depends on previous gender norms
				BaseDrinksPerDD = Register[ic].DrinksPerDD;
				Register[ic].DrinksPerDD += Register[ic].IneqGender * (TempAdj - 1.0) * GenderIneqEffectDrinksPerDD / 
					AlcPropnReported;
				if (Register[ic].DrinksPerDD >= 5.0 && BaseDrinksPerDD < 5.0){
					Register[ic].TempCondomAdj *= pow(ORcondomBingePW, Register[ic].DailyDrinkProb * 7.0);
				}
				if (BaseDrinksPerDD >= 5.0 && Register[ic].DrinksPerDD < 5.0){
					Register[ic].TempCondomAdj *= 1.0 / pow(ORcondomBingePW, Register[ic].DailyDrinkProb * 7.0);
				}
				// Then update inequitable gender norm score
				BaseIneqGender = Register[ic].IneqGender;
				Register[ic].IneqGender *= TempAdj;
				Register[ic].TempCondomAdj *= pow(ORcondomGenderIneq, Register[ic].IneqGender - BaseIneqGender);
			}
		}
	}
}

void Pop::OneBehavCycle()
{
	int ic, STDcyclesPerBehavCycle;

	if (StructuralRCTcalib == 1 && CurrYear >= BaselineStart){ ResetStructOutcomes(); }

	BehavCycleCount += 1;
	for(ic=0; ic<Register.size(); ic++){
		if(Register[ic].AliveInd==1 && Register[ic].AgeGroup>=2){
			Register[ic].NewStatus = 0;}
	}
	UpdateBirths();

	// Call the OneSTDcycle the appropriate number of times
	STDcyclesPerBehavCycle = CycleD/CycleS;
	STDcycleCount = 0;
	for(ic=0; ic<STDcyclesPerBehavCycle; ic++){
		OneSTDcycle();}

	// Calculate non-HIV mortality
	UpdateNonHIVmort();
	
	// Calculate movements between sexual behaviour classes
	UpdateFSWprofile();
	BalanceSexualPartners();
	CalcPartnerTransitions();
	UpdateMSMcasual();
	UpdateHetCasual();
	UpdateVisitStatus();
	
	if (StructuralRCTcalib == 1 && CurrYear>=BaselineStart){ GetStructRCTests(); }
}

void Pop::UpdateBirths()
{
	int ic, ia, ii;

	for(ic=0; ic<Register.size(); ic++){
		if(Register[ic].AliveInd==1 && Register[ic].SexInd==1){
			if (Register[ic].DOLB>(1.0 * CurrYear + 0.5 + 1.0 * (BehavCycleCount - 1.0)/CycleS) && 
				Register[ic].DOLB<(1.0 * CurrYear + 0.5 + 1.0 * BehavCycleCount / CycleS)){
				NewBirth(ic+1);
				if (FixedUncertainty == 1 && CurrYear >= 1995 && CurrYear < 1998){
					ia = Register[ic].AgeGroup - 3;
					ii = Register[ic].MarriedInd;
					if (ia>=0){ FertByMarriage.out[CurrSim - 1][ia + ii * 7] += 1; }
				}
			}
		}
	}
}

void Pop::UpdateNonHIVmort()
{
	int ic;
	double rmort[MaxPop], tt;

	//int seed = 5020 + CurrYear * 27 + BehavCycleCount * 85;
	//if(CurrYear>=StartYear+FixedPeriod){
	//	seed += CurrSim;}
	//CRandomMersenne rg(seed);
	for(ic=0; ic<MaxPop; ic++){
		rmort[ic] = rg.Random();}
	for(ic=0; ic<Register.size(); ic++){
		if(Register[ic].AliveInd==1){
			if(Register[ic].AgeGroup==0){ 
				if(Register[ic].DOB > (CurrYear + 0.5 + (BehavCycleCount - 1.0)/CycleS)){
					// Time from start of behav cycle to DOB (in years)
					tt = Register[ic].DOB - (CurrYear + 0.5 + (BehavCycleCount - 1.0)/CycleS);
					// Time from DOB to end of behav cycle (in behav cycles)
					tt = 1.0 - (tt * CycleS);
					Register[ic].NonHIVmortProb = 1.0 - pow(1.0 - 
						Register[ic].NonHIVmortProb, tt);
				}
			}
			if(rmort[ic]<Register[ic].NonHIVmortProb){
				SetToDead(ic+1, 0);}
			if (FixedUncertainty == 1 && Register[ic].AgeGroup > 2 && Register[ic].HIVstage > 4){
				ARTdeathsExp.out[CurrSim - 1][CurrYear - StartYear] += Register[ic].NonHIVmortProb * Register[ic].PopWeight;
			}
		}
	}
}

void Pop::UpdateFSWprofile()
{
	int ia, ic, ii, ij, is, ir;
	double rexit[MaxCSWs], rentry[MaxCSWs];
	double ExitProb, EntryProb, temp, normalizer[3], DoubleIntDif, IntPart;
	double SingleHRfem[16][3];
	int CSWID, offset, ActualNewCSW, tempID[MaxCSWs], TotFSWbyAge[16][3], EntryInd;

	FSWexit = 1.0 / DurFSW[CurrYear - StartYear];

	// Firstly randomly determine which sex workers retire
	for (ir = 0; ir < 3; ir++){
		for (ic = 0; ic < TotCurrFSW[ir]; ic++){
			rexit[ic] = rg.Random();}
		for (ic = 0; ic < MaxCSWs; ic++){
			tempID[ic] = CSWregister[ic][ir];}
		offset = 0;
		for (ic = 0; ic<MaxCSWs; ic++){
			CSWID = CSWregister[ic][ir];
			if (CSWID == 0){
				break;}
			ia = Register[CSWID - 1].AgeGroup - 2;
			is = Register[CSWID - 1].HIVstage;
			ExitProb = exp(-FSWexit / CycleS);
			if (is>0){
				ExitProb = pow(ExitProb, HIVeffectFSWexit[is - 1]);}
			ExitProb = 1.0 - ExitProb;
			if (Register[CSWID - 1].Employed == 1){ ExitProb = 1.0; }
			if (rexit[ic] < ExitProb){
				TotCurrFSW[ir] = TotCurrFSW[ir] - 1;
				Register[CSWID - 1].FSWind = 0;
				for (ij = 0; ij < MaxCSWs - ic - 1 - offset; ij++){
					tempID[ic + ij - offset] = tempID[ic + ij + 1 - offset];
				}
				offset += 1;
			}
			Register[CSWID - 1].NewStatus = 1;
		}
		for (ic = 0; ic < MaxCSWs; ic++){
			CSWregister[ic][ir] = tempID[ic];}
	}

	// Secondly get rates of entry into CSW and age-specific weights
	for (ir = 0; ir < 3; ir++){
		TotCurrFSW[ir] = 0;
		DesiredFSWcontacts[ir] = 0.0;
		for (ia = 0; ia < 16; ia++){
			FSWentry[ia][ir] = 0.0;
			TotFSWbyAge[ia][ir] = 0;
			SingleHRfem[ia][ir] = 0.0;
		}
	}
	if (FixedUncertainty == 1 && CurrYear == 2016 && BehavCycleCount == 1 && CurrSim == 1){
		for (ii = 0; ii < 16; ii++){ FSWcontacts2016[ii] = 0.0; }
	}
	for(ic=0; ic<Register.size(); ic++){
		if(Register[ic].AliveInd==1 && Register[ic].AgeGroup>=2 && Register[ic].RiskGroup==1){
			if(Register[ic].SexInd==1 && Register[ic].CurrPartners==0 &&
				Register[ic].VirginInd == 0 && Register[ic].Employed == 0){
					ia = Register[ic].AgeGroup - 2;
					ir = Register[ic].Race;
					if(Register[ic].FSWind==0){
						is = Register[ic].HIVstage;
						if(is==0){
							SingleHRfem[ia][ir] += 1.0;}
						else{
							if (Register[ic].VCThistory < 2){ SingleHRfem[ia][ir] += HIVeffectFSWentry[is - 1]; }
							else{ SingleHRfem[ia][ir] += HIVeffectFSWentry[is - 1] * RR_FSWentryDiagnosed; }
						}
					}
					else{
						TotFSWbyAge[ia][ir] += 1;
						TotCurrFSW[ir] += 1;
					}
				}
			if(Register[ic].SexInd==0 && Register[ic].VirginInd==0 && Register[ic].Imprisoned==0){
				ia = Register[ic].CurrAge - 10;
				ir = Register[ic].Race;
				is = Register[ic].HIVstage;
				temp = FSWcontactConstant * AgeEffectFSWcontact[ia];
				temp *= (1.0 - Register[ic].MalePref);
				if (is > 0){ temp *= HIVeffectPartners[is - 1]; }
				if (Register[ic].CurrUrban == 1){ temp *= SWurban[0]; }
				else{ temp *= SWurban[1]; }
				if (Register[ic].Employed == 1){ temp *= RR_FSWcontactEmployedM; }
				//if (Register[ic].Employed == 1 || CurrYear >= 2000){ temp *= RR_FSWcontactEmployedM; }
				if(Register[ic].CurrPartners==0){ temp *= PartnerEffectFSWcontact[0]; }
				else{ temp *= RaceEffectNew[ir]; }
				if(Register[ic].CurrPartners==1 && Register[ic].MarriedInd==0){ temp *= PartnerEffectFSWcontact[1]; }
				if(Register[ic].CurrPartners==1 && Register[ic].MarriedInd==1){ temp *= PartnerEffectFSWcontact[2]; }
				if(Register[ic].CurrPartners==2 && Register[ic].MarriedInd==0){ temp *= PartnerEffectFSWcontact[3]; }
				if(Register[ic].CurrPartners==2 && Register[ic].MarriedInd==1){ temp *= PartnerEffectFSWcontact[4]; }
				DesiredFSWcontacts[ir] += temp;
				// Outputs for age pattern of male contacts with FSWs (optional)
				if (FixedUncertainty == 1 && CurrYear == 2016 && BehavCycleCount == 1){
					ii = Register[ic].AgeGroup - 2;
					FSWcontacts2016[ii] += temp * Register[ic].PopWeight;
				}
			}
		}
	}
	for (ir = 0; ir < 3; ir++){
		if (TotCurrFSW[ir] < DesiredFSWcontacts[ir] / AnnNumberClients){
			normalizer[ir] = 0.0;
			for (ia = 0; ia<16; ia++){
				if (SingleHRfem[ia][ir]>0.0){
					FSWentry[ia][ir] = ((DesiredFSWcontacts[ir] / AnnNumberClients) *
						InitFSWageDbn[ia] - TotFSWbyAge[ia][ir]) / SingleHRfem[ia][ir];
				}
				else{
					FSWentry[ia][ir] = 0.0;}
				if (FSWentry[ia][ir] < 0.0){
					FSWentry[ia][ir] = 0.0;}
				else{
					normalizer[ir] += FSWentry[ia][ir] * SingleHRfem[ia][ir];
				}
			}
		}
		else{
			for (ia = 0; ia < 16; ia++){
				FSWentry[ia][ir] = 0.0;}
		}
	}

	// Lastly determine which women become sex workers
	for (ir = 0; ir < 3; ir++){
		if (TotCurrFSW[ir] < DesiredFSWcontacts[ir] / AnnNumberClients){
			for (ic = 0; ic < MaxCSWs; ic++){
				rentry[ic] = rg.Random();}

			// (a) Determine the number of new sex workers
			RequiredNewFSW = (DesiredFSWcontacts[ir] / AnnNumberClients) - TotCurrFSW[ir];
			DoubleIntDif = modf(RequiredNewFSW, &IntPart);
			ActualNewCSW = RequiredNewFSW - DoubleIntDif;
			if (rentry[0] < DoubleIntDif){
				ActualNewCSW += 1;
			}

			// (b) Determine which women become sex workers
			temp = 0.0;
			for (ic = 0; ic<Register.size(); ic++){
				if (Register[ic].AliveInd == 1 && Register[ic].AgeGroup >= 2 && Register[ic].RiskGroup == 1 &&
					Register[ic].SexInd == 1 && Register[ic].VirginInd == 0 && Register[ic].Race == ir &&
					Register[ic].CurrPartners == 0 && Register[ic].FSWind == 0 && Register[ic].Employed == 0){
					ia = Register[ic].AgeGroup - 2;
					is = Register[ic].HIVstage;
					EntryInd = 0;
					EntryProb = FSWentry[ia][ir] / normalizer[ir];
					if (is>0){
						if (Register[ic].VCThistory<2){ EntryProb *= HIVeffectFSWentry[is - 1]; }
						else{ EntryProb *= HIVeffectFSWentry[is - 1] * RR_FSWentryDiagnosed; }
					}
					for (ij = 1; ij <= ActualNewCSW; ij++){
						if (rentry[ij] > temp && rentry[ij] <= (temp + EntryProb)){
							EntryInd = 1;}
					}
					if (EntryInd == 1){
						Register[ic].NewStatus = 1;
						// Not nec to call SetNewStatusTo1 since EligibleByAge only gets 
						// calculated later (in the CalcPartnerTransition function).
						Register[ic].FSWind = 1;
						if (Register[ic].CurrContr == 0){ 
							Register[ic].TestStartContr(1, 0); }
						CSWregister[TotCurrFSW[ir]][ir] = ic + 1;
						TotCurrFSW[ir] += 1;
					}
					temp += EntryProb;
				}
			}
		}
	}
}

void Pop::UpdateVisitStatus()
{
	int ic, PID1, PID2, Eligible;
	double AveLength, ProbReturn, ProbVisit;

	for (ic = 0; ic < Register.size(); ic++){
		r2[ic] = rg.Random();}

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1){
			AveLength = MaxVisitLength / pow(Register[ic].VisitFreq, 0.5);
			if (Register[ic].Visiting == 1 && AveLength < 365.0 / CycleS){
				Register[ic].Visiting = 0; }
			else if (Register[ic].Visiting == 1 && AveLength >= 365.0 / CycleS){
				ProbReturn = 365.0/(CycleS * AveLength); 
				if (r2[ic] < ProbReturn){ Register[ic].Visiting = 0; }
			}
			else{
				Eligible = 0;
				if (Register[ic].CurrPartners>0){
					PID1 = Register[ic].IDprimary;
					if (Register[ic].CurrUrban != Register[PID1 - 1].CurrUrban &&
						Register[PID1 - 1].Visiting == 0){
						Eligible = 1;
					}
				}
				if (Register[ic].CurrPartners==2 && Eligible==0){
					PID2 = Register[ic].ID2ndary;
					if (Register[ic].CurrUrban != Register[PID2 - 1].CurrUrban &&
						Register[PID2 - 1].Visiting == 0){
						Eligible = 1;
					}
				}
				if (Eligible == 1){
					ProbVisit = 0.5 * Register[ic].VisitFreq / CycleS;
					if (r2[ic] < ProbVisit){ Register[ic].Visiting = 1; }
				}
			}
		}
	}
}

void Pop::BalanceSexualPartners()
{
	int ia, ic, ig, ii, ij, ip, ir, is;
	int ID1, ID2, PR1, PR2, HHID;
	int TotFSWbyAge[16];
	double temp, ConcurPenalty, lambda;

	for(ir=0; ir<2; ir++){
		for(ig=0; ig<3; ig++){
			DesiredSTpartners[ir][ig] = 0.0;}
	}

	// First calculate the desired rates at which new partnerships are formed
	for(ic=0; ic<Register.size(); ic++){
		if(Register[ic].AliveInd==1 && Register[ic].AgeGroup>=2){
			if (Register[ic].Imprisoned==0 && (Register[ic].CurrPartners == 0 || 
				(Register[ic].CurrPartners == 1 && Register[ic].RiskGroup == 1))){
				ia = Register[ic].CurrAge - 10;
				ig = Register[ic].SexInd;
				ip = Register[ic].Race;
				ir = Register[ic].RiskGroup - 1;
				is = Register[ic].HIVstage;
				if(Register[ic].VirginInd==0){
					Register[ic].DesiredNewPartners = PartnershipFormation[ir][ig] * 
						AgeEffectPartners[ia][ig] * Register[ic].PartnerRateAdj;
					if (Register[ic].CurrPartners == 1){
						ConcurPenalty = RaceEffectNew[ip];
						if (ig == 0){ 
							ConcurPenalty *= (1.0 + EffectIneqGenderConcurrency * (Register[ic].IneqGender - 0.2)); }
						if (Register[ic].MarriedInd == 0){ ConcurPenalty *= PartnerEffectNew[0][ig]; }
						else{ ConcurPenalty *= PartnerEffectNew[1][ig]; }
						ID1 = Register[ic].IDprimary;
						if (Register[ic].CurrUrban != Register[ID1 - 1].CurrUrban ||
							Register[ID1 - 1].Imprisoned>0){
							lambda = pow(AbsentPartnerConcurAdj, (365.0 - MaxVisitLength * pow(
								Register[ic].VisitFreq, 0.5)) / 365.0);
							if (ConcurPenalty > 0.0){ ConcurPenalty = pow(ConcurPenalty, lambda); }
						}
						if (ConcurPenalty > 1.0){ ConcurPenalty = 1.0; }
						if (ConcurPenalty < 0.0){ ConcurPenalty = 0.0; }
						Register[ic].DesiredNewPartners *= ConcurPenalty;
					}
					if (ig == 0 && Register[ic].Employed==1){
						Register[ic].DesiredNewPartners *= RR_STemployedM;}
					if(is>0){
						Register[ic].DesiredNewPartners *= HIVeffectPartners[is-1];}
					if(Register[ic].FSWind==1){
						Register[ic].DesiredNewPartners = 0.0;}
					if (ig == 1 || Register[ic].MalePref == 0.0){
						DesiredSTpartners[ir][ig] += Register[ic].DesiredNewPartners;} 
					else{
						DesiredSTpartners[ir][0] += Register[ic].DesiredNewPartners * (1.0 - Register[ic].MalePref);
						DesiredSTpartners[ir][2] += Register[ic].DesiredNewPartners * Register[ic].MalePref;
					}
				}
				else{
					Register[ic].DesiredNewPartners = SexualDebut[ia][ig];
					if(ir==1){
						Register[ic].DesiredNewPartners *= DebutAdjLow[ig];}
					Register[ic].DesiredNewPartners *= DebutAdjRace[ip];
					//if (Register[ic].InSchool == 1 || (CurrYear >= 2000 && Register[ic].CurrAge <= 21)){
					//	Register[ic].DesiredNewPartners *= RRdebutInSchool[ig]; }
					if (Register[ic].InSchool == 1){ Register[ic].DesiredNewPartners *= RRdebutInSchool[ig];}
					if (ig == 1){
						HHID = Register[ic].HouseholdID;
						if (HHID > 0){
							if (HHregister[HHID - 1].PerCapitaIncomeAdj < MeanLogIncome){
								if (HHregister[HHID - 1].PerCapitaIncomeAdj > 0){
									Register[ic].DesiredNewPartners *= pow(RRdebutLogDropIncomeF,
										MeanLogIncome - log(HHregister[HHID - 1].PerCapitaIncomeAdj));
								}
								else{
									Register[ic].DesiredNewPartners *= pow(RRdebutLogDropIncomeF,
										MeanLogIncome);
								}
							}
						}
						if (HHID == 0 && Register[ic].CurrAge >= 15){
							Register[ic].DesiredNewPartners *= pow(RRdebutLogDropIncomeF,
								MeanLogIncome);
						}
					}
				}
			}
			else{
				Register[ic].DesiredNewPartners = 0.0;}
		}
	}

	// Then calculate the balancing factors for the rates of ST partner acquisition (heterosexual)
	DesiredPartnerRiskM[0][0] = (1.0 - AssortativeM) + AssortativeM * 
		DesiredSTpartners[0][1] / (DesiredSTpartners[0][1] + DesiredSTpartners[1][1]);
	DesiredPartnerRiskM[0][1] = 1.0 - DesiredPartnerRiskM[0][0];
	DesiredPartnerRiskM[1][1] = (1.0 - AssortativeM) + AssortativeM * 
		DesiredSTpartners[1][1] / (DesiredSTpartners[0][1] + DesiredSTpartners[1][1]);
	DesiredPartnerRiskM[1][0] = 1.0 - DesiredPartnerRiskM[1][1];

	DesiredPartnerRiskF[0][0] = (1.0 - AssortativeF) + AssortativeF * 
		DesiredSTpartners[0][0] / (DesiredSTpartners[0][0] + DesiredSTpartners[1][0]);
	DesiredPartnerRiskF[0][1] = 1.0 - DesiredPartnerRiskF[0][0];
	DesiredPartnerRiskF[1][1] = (1.0 - AssortativeF) + AssortativeF * 
		DesiredSTpartners[1][0] / (DesiredSTpartners[0][0] + DesiredSTpartners[1][0]);
	DesiredPartnerRiskF[1][0] = 1.0 - DesiredPartnerRiskF[1][1];

	DesiredPartnerRiskMSM[0][0] = (1.0 - AssortativeM) + AssortativeM *
		DesiredSTpartners[0][2] / (DesiredSTpartners[0][2] + DesiredSTpartners[1][2]);
	DesiredPartnerRiskMSM[0][1] = 1.0 - DesiredPartnerRiskMSM[0][0];
	DesiredPartnerRiskMSM[1][1] = (1.0 - AssortativeM) + AssortativeM *
		DesiredSTpartners[1][2] / (DesiredSTpartners[0][2] + DesiredSTpartners[1][2]);
	DesiredPartnerRiskMSM[1][0] = 1.0 - DesiredPartnerRiskMSM[1][1];

	if(AllowBalancing==1){
		for(ii=0; ii<2; ii++){
			for(ij=0; ij<2; ij++){
				AdjSTrateM[ii][ij] = (GenderEquality * DesiredSTpartners[ij][1] * 
					DesiredPartnerRiskF[ij][ii] + (1.0 - GenderEquality) * 
					DesiredSTpartners[ii][0] * DesiredPartnerRiskM[ii][ij])/
					(DesiredSTpartners[ii][0] * DesiredPartnerRiskM[ii][ij]);
				AdjSTrateF[ii][ij] = (GenderEquality * DesiredSTpartners[ii][1] * 
					DesiredPartnerRiskF[ii][ij] + (1.0 - GenderEquality) * 
					DesiredSTpartners[ij][0] * DesiredPartnerRiskM[ij][ii])/
					(DesiredSTpartners[ii][1] * DesiredPartnerRiskF[ii][ij]);
				AdjSTrateMSM[ii][ij] = (0.5 * DesiredSTpartners[ii][2] *
					DesiredPartnerRiskMSM[ii][ij] + 0.5 *
					DesiredSTpartners[ij][2] * DesiredPartnerRiskMSM[ij][ii]) /
					(DesiredSTpartners[ii][2] * DesiredPartnerRiskMSM[ii][ij]);
			}
		}
	}
	else{
		for(ii=0; ii<2; ii++){
			for(ij=0; ij<2; ij++){
				AdjSTrateM[ii][ij] = 1.0;
				AdjSTrateF[ii][ij] = 1.0;
				AdjSTrateMSM[ii][ij] = 1.0;
			}
		}
	}

	// Secondly calculate the desired rates at which marriages are formed
	/*for(ii=0; ii<2; ii++){
		for(ij=0; ij<2; ij++){
			DesiredMarriagesM[ii][ij] = 0.0;
			DesiredMarriagesF[ii][ij] = 0.0;
		}
	}
	for(ic=0; ic<Register.size(); ic++){
		if(Register[ic].AliveInd==1 && Register[ic].AgeGroup>=2){
			if(Register[ic].CurrPartners>0 && Register[ic].MarriedInd==0){
				ia = Register[ic].CurrAge - 10;
				ig = Register[ic].SexInd;
				ip = Register[ic].Race;
				ir = Register[ic].RiskGroup - 1;
				ID1 = Register[ic].IDprimary;
				if(Register[ID1-1].MarriedInd==0){ // This condition isn't included in THISA
					PR1 = Register[ID1-1].RiskGroup - 1;
					if(ig==0){
						DesiredMarriagesM[ir][PR1] += MarriageRate[PR1][ir] *
							AgeEffectMarriage[ia][ir] * RaceMarriage[ip];}
					else{
						DesiredMarriagesF[ir][PR1] += MarriageRate[PR1][2+ir] *
							AgeEffectMarriage[ia][2 + ir] * RaceMarriage[ip];}
				}
				if(Register[ic].CurrPartners==2){
					ID2 = Register[ic].ID2ndary;
					if(Register[ID2-1].MarriedInd==0){ 
						PR2 = Register[ID2-1].RiskGroup - 1;
						if(ig==0){
							DesiredMarriagesM[ir][PR2] += MarriageRate[PR2][ir] *
								AgeEffectMarriage[ia][ir] * RaceMarriage[ip];}
						else{
							DesiredMarriagesF[ir][PR2] += MarriageRate[PR2][2+ir] *
								AgeEffectMarriage[ia][2 + ir] * RaceMarriage[ip];}
					}
				}
			}
		}
	}*/

	// Then calculate the balancing factors for the rates of marriage
	if(AllowBalancing==1){
		/*for(ii=0; ii<2; ii++){
			for(ij=0; ij<2; ij++){
				AdjLTrateM[ii][ij] = (GenderEquality * DesiredMarriagesF[ij][ii] +
					(1.0 - GenderEquality) * DesiredMarriagesM[ii][ij])/DesiredMarriagesM[ii][ij];
				AdjLTrateF[ii][ij] = (GenderEquality * DesiredMarriagesF[ii][ij] + 
					(1.0 - GenderEquality) * DesiredMarriagesM[ij][ii])/DesiredMarriagesF[ii][ij];
			}
		}*/
	}
	else{
		for(ii=0; ii<2; ii++){
			for(ij=0; ij<2; ij++){
				AdjLTrateM[ii][ij] = 1.0;
				AdjLTrateF[ii][ij] = 1.0;
			}
		}
	}
}

void Pop::CalcPartnerTransitions()
{
	int ia, ic, ii, ir;
	int TotUnassigned, nextID, EventType, pID, pID2, prisk, index1, sexpref;
	double TotalSelectionProb, rr, TempCondom;
	vector<int> indices(Register.size());

	for(ic=0; ic<MaxPop; ic++){
		r2[ic] = rg.Random();}
	for(ic=0; ic<MaxPop; ic++){
		revent[ic] = rg.Random();}
	// New addition for sampling partner age;
	for(ic=0; ic<MaxPop; ic++){
		rpAge[ic] = rg.Random();} 
	for (ic = 0; ic<MaxPop; ic++){
		rpref[ic] = rg.Random();}
	
	TotUnassigned = 0; // Number who haven't yet been assigned a status
	for(ia=0; ia<16; ia++){
		MalePartners[ia][0].Reset();
		MalePartners[ia][1].Reset();
		FemPartners[ia][0].Reset();
		FemPartners[ia][1].Reset();
		MSMpartners[ia][0].Reset();
		MSMpartners[ia][1].Reset();
	}
	for(ic=0; ic<Register.size(); ic++){
		if(Register[ic].AliveInd==1 && Register[ic].AgeGroup>=2 && Register[ic].NewStatus==0){
			TotUnassigned += 1;
			if(Register[ic].DesiredNewPartners>0 && Register[ic].VirginInd==0){
				ia = Register[ic].AgeGroup - 2; // -2 because age groups start at 10
				ir = Register[ic].RiskGroup - 1;
				if(Register[ic].SexInd==0 && Register[ic].MalePref<1.0){
					MalePartners[ia][ir].AddMember(ic+1, 1);}
				if (Register[ic].SexInd == 0 && Register[ic].MalePref>0.0){
					MSMpartners[ia][ir].AddMember(ic + 1, 0);}
				if (Register[ic].SexInd == 1){
					FemPartners[ia][ir].AddMember(ic+1, 0);}
			}
		}
	}
	for(ii=0; ii<indices.size(); ii++){
		indices[ii] = ii;}
	ii = 0;
	while(TotUnassigned>0){
		nextID = 0;
		while (nextID == 0){
			index1 = r2[ii] * indices.size();
			nextID = indices[index1] + 1;
			if (Register[nextID - 1].AliveInd == 0 || Register[nextID - 1].AgeGroup < 2 ||
				Register[nextID - 1].NewStatus == 1){
				nextID = 0;
				ii += 1;
			}
			indices[index1] = indices[indices.size() - 1];
			indices.pop_back();
		}
		rr = revent[ii];
		if (Register[nextID - 1].SexInd == 1 || Register[nextID - 1].MalePref == 0){ sexpref = 0; }
		else{ // Gay/bisexual men
			if (rpref[ii] < Register[nextID - 1].MalePref){ sexpref = 1; }
			else{ sexpref = 0; } // Choose opposite sex partner
		}
		EventType = Register[nextID - 1].SelectEvent(nextID, rr, sexpref); 
		SetNewStatusTo1(nextID);
		TotUnassigned = TotUnassigned - 1;
		if(EventType==0){ // No change in behav status
			if(Register[nextID-1].CurrPartners>0){
				pID = Register[nextID-1].IDprimary;
				prisk = Register[pID-1].RiskGroup;
				if(prisk==2 && Register[pID-1].NewStatus==0){
					SetNewStatusTo1(pID);
					TotUnassigned = TotUnassigned - 1;
				}
			}
			if(Register[nextID-1].CurrPartners==2){
				pID = Register[nextID-1].ID2ndary;
				prisk = Register[pID-1].RiskGroup;
				if(prisk==2 && Register[pID-1].NewStatus==0){
					SetNewStatusTo1(pID);
					TotUnassigned = TotUnassigned - 1;
				}
			}
		}
		if(EventType==1){ // Acquire ST partner in high risk group
			if (Register[nextID - 1].SexInd == 1 ||
				(Register[nextID - 1].SexInd == 0 && Register[nextID - 1].MalePref == 0.0)){
					// Heterosexual default
					if (Register[nextID - 1].CurrPartners == 1){
						if (Register[nextID - 1].SexInd == 0){ rpAge[ii] = pow(rpAge[ii], AgePrefAdj2ndary); }
						else{ rpAge[ii] = pow(rpAge[ii], 1.0/AgePrefAdj2ndary); }
					}
					GetNewPartner(nextID, rpAge[ii], 1);
					if (MaxNewPartnerInd == 1){
						GetNewPartner(nextID, rpAge[ii], 2);
					}
			}
			else{ // MSM
				if (sexpref==1){
					GetNewMSMpartner(nextID, rpAge[ii], 1);
					if (MaxNewPartnerInd == 1){
						GetNewMSMpartner(nextID, rpAge[ii], 2);
					}
				}
				else{ // MSM selects a female partner
					GetNewPartner(nextID, rpAge[ii], 1);
					if (MaxNewPartnerInd == 1){
						GetNewPartner(nextID, rpAge[ii], 2);
					}
				}
			}
			if(MaxNewPartnerInd!=1){ // Indiv has been assigned a partner
				TotUnassigned = TotUnassigned - 1;}
			if(MaxNewPartnerInd==1){ // Will only happen if there are no 
				// potential partners in either high or low risk group. In this
				// case, assume there is no change in behavioural status.
				if(Register[nextID-1].CurrPartners>0){
					pID = Register[nextID-1].IDprimary;
					prisk = Register[pID-1].RiskGroup;
					if(prisk==2 && Register[pID-1].NewStatus==0){
						SetNewStatusTo1(pID);
						TotUnassigned = TotUnassigned - 1;
					}
				}
			}
		}
		if(EventType==2){ // Acquire ST partner in low risk group
			if (Register[nextID - 1].SexInd == 1 ||
				(Register[nextID - 1].SexInd == 0 && Register[nextID - 1].MalePref == 0.0)){
					// Heterosexual default
					if (Register[nextID - 1].CurrPartners == 1){
						if (Register[nextID - 1].SexInd == 0){ rpAge[ii] = pow(rpAge[ii], AgePrefAdj2ndary); }
						else{ rpAge[ii] = pow(rpAge[ii], 1.0 / AgePrefAdj2ndary); }
					}
					GetNewPartner(nextID, rpAge[ii], 2);
					if (MaxNewPartnerInd == 1){
						GetNewPartner(nextID, rpAge[ii], 1);}
			}
			else{ // MSM
				if (sexpref==1){
					GetNewMSMpartner(nextID, rpAge[ii], 2);
					if (MaxNewPartnerInd == 1){
						GetNewMSMpartner(nextID, rpAge[ii], 1);}
				}
				else{ // MSM selects a female partner
					GetNewPartner(nextID, rpAge[ii], 2);
					if (MaxNewPartnerInd == 1){
						GetNewPartner(nextID, rpAge[ii], 1);
					}
				}
			}
			if(MaxNewPartnerInd!=1){ // Indiv has been assigned a partner
				TotUnassigned = TotUnassigned - 1;}
			if(MaxNewPartnerInd==1){ // Will only happen if there are no 
				// potential partners in either high or low risk group. In this
				// case, assume there is no change in behavioural status.
				if(Register[nextID-1].CurrPartners>0){
					pID = Register[nextID-1].IDprimary;
					prisk = Register[pID-1].RiskGroup;
					if(prisk==2 && Register[pID-1].NewStatus==0){
						SetNewStatusTo1(pID);
						TotUnassigned = TotUnassigned - 1;
					}
				}
			}
		}
		if(EventType==3){ // Marry primary partner
			pID = Register[nextID-1].IDprimary;
			Register[nextID-1].MarriedInd = 1;
			Register[pID-1].MarriedInd = 1;
			ChangeHHafterMarriage(nextID, pID);
			SetNewStatusTo1(pID);
			TotUnassigned = TotUnassigned - 1;
			Register[nextID - 1].AssignCondom(1, nextID, pID);
			if(Register[pID-1].ID2ndary==nextID){
				TempCondom = Register[pID - 1].CondomPrimary;
				Register[pID-1].ID2ndary = Register[pID-1].IDprimary;
				Register[pID-1].IDprimary = nextID;
				Register[pID - 1].CondomPrimary = Register[pID - 1].Condom2ndary;
				Register[pID - 1].Condom2ndary = TempCondom;
			}
			if(Register[pID-1].CurrPartners==2){
				pID2 = Register[pID-1].ID2ndary;
				prisk = Register[pID2-1].RiskGroup;
				if(prisk==2 && Register[pID2-1].NewStatus==0){
					SetNewStatusTo1(pID2);
					TotUnassigned = TotUnassigned - 1;
				}
			}
			if(Register[nextID-1].CurrPartners==2){
				pID = Register[nextID-1].ID2ndary;
				prisk = Register[pID-1].RiskGroup;
				if(prisk==2 && Register[pID-1].NewStatus==0){
					SetNewStatusTo1(pID);
					TotUnassigned = TotUnassigned - 1;
				}
			}
		}
		if(EventType==4){ // Marry 2ndary partner
			pID = Register[nextID-1].ID2ndary;
			Register[nextID-1].MarriedInd = 1;
			Register[pID-1].MarriedInd = 1;
			ChangeHHafterMarriage(nextID, pID);
			SetNewStatusTo1(pID);
			TotUnassigned = TotUnassigned - 1;
			Register[nextID - 1].AssignCondom(1, nextID, pID);
			TempCondom = Register[nextID - 1].CondomPrimary;
			Register[nextID-1].ID2ndary = Register[nextID-1].IDprimary;
			Register[nextID-1].IDprimary = pID;
			Register[nextID - 1].CondomPrimary = Register[nextID - 1].Condom2ndary;
			Register[nextID - 1].Condom2ndary = TempCondom;
			if(Register[pID-1].ID2ndary==nextID){
				TempCondom = Register[pID - 1].CondomPrimary;
				Register[pID-1].ID2ndary = Register[pID-1].IDprimary;
				Register[pID-1].IDprimary = nextID;
				Register[pID - 1].CondomPrimary = Register[pID - 1].Condom2ndary;
				Register[pID - 1].Condom2ndary = TempCondom;
			}
			if(Register[pID-1].CurrPartners==2){
				pID2 = Register[pID-1].ID2ndary;
				prisk = Register[pID2-1].RiskGroup;
				if(prisk==2 && Register[pID2-1].NewStatus==0){
					SetNewStatusTo1(pID2);
					TotUnassigned = TotUnassigned - 1;
				}
			}
			pID = Register[nextID-1].ID2ndary; // Not the same as before
			prisk = Register[pID-1].RiskGroup;
			if(prisk==2 && Register[pID-1].NewStatus==0){
				SetNewStatusTo1(pID);
				TotUnassigned = TotUnassigned - 1;
			}
		}
		if(EventType==5){ // Divorce
			pID = Register[nextID-1].IDprimary;
			Register[nextID-1].MarriedInd = 0;
			Register[pID-1].MarriedInd = 0;
			Register[nextID-1].DOUD = 1.0 * (CurrYear + 0.5 + 1.0 * BehavCycleCount / CycleS);
			Register[pID-1].DOUD = 1.0 * (CurrYear + 0.5 + 1.0 * BehavCycleCount / CycleS);
			ChangeHHafterDivorce(nextID, pID);
			SetNewStatusTo1(pID);
			TotUnassigned = TotUnassigned - 1;
			if(Register[nextID-1].CurrPartners==1){
				Register[nextID-1].CurrPartners = 0;
				Register[nextID-1].IDprimary = 0;
				Register[nextID-1].DisclosedPrimary = 0;
			}
			else{
				Register[nextID-1].CurrPartners = 1;
				Register[nextID-1].IDprimary = Register[nextID-1].ID2ndary;
				Register[nextID-1].CondomPrimary = Register[nextID-1].Condom2ndary;
				Register[nextID-1].DisclosedPrimary = Register[nextID-1].Disclosed2ndary;
				Register[nextID-1].ID2ndary = 0;
				Register[nextID-1].Disclosed2ndary = 0;
			}
			if(Register[pID-1].CurrPartners==1){
				Register[pID-1].CurrPartners = 0;
				Register[pID-1].IDprimary = 0;
				Register[pID-1].DisclosedPrimary = 0;
			}
			else{
				Register[pID-1].CurrPartners = 1;
				Register[pID-1].IDprimary = Register[pID-1].ID2ndary;
				Register[pID-1].CondomPrimary = Register[pID-1].Condom2ndary;
				Register[pID-1].DisclosedPrimary = Register[pID-1].Disclosed2ndary;
				Register[pID-1].ID2ndary = 0;
				Register[pID-1].Disclosed2ndary = 0;
			}
		}
		if(EventType==6){ // End primary ST relationship
			pID = Register[nextID-1].IDprimary;
			if(Register[nextID-1].CurrPartners==1){
				Register[nextID-1].IDprimary = 0;
				Register[nextID-1].CurrPartners = 0;
				Register[nextID-1].DisclosedPrimary = 0;
			}
			else{
				Register[nextID-1].IDprimary = Register[nextID-1].ID2ndary;
				Register[nextID-1].CondomPrimary = Register[nextID-1].Condom2ndary;
				Register[nextID-1].DisclosedPrimary = Register[nextID-1].Disclosed2ndary;
				Register[nextID-1].ID2ndary = 0;
				Register[nextID-1].CurrPartners = 1;
				Register[nextID-1].Disclosed2ndary = 0;
			}
			if(Register[pID-1].CurrPartners==1){
				Register[pID-1].IDprimary = 0;
				Register[pID-1].CurrPartners = 0;
				Register[pID-1].DisclosedPrimary = 0;
			}
			else{
				if(Register[pID-1].IDprimary==nextID){
					Register[pID-1].IDprimary = Register[pID-1].ID2ndary;
					Register[pID-1].CondomPrimary = Register[pID-1].Condom2ndary;
					Register[pID-1].DisclosedPrimary = Register[pID-1].Disclosed2ndary;
				}
				Register[pID-1].ID2ndary = 0;
				Register[pID-1].CurrPartners = 1;
				Register[pID-1].Disclosed2ndary = 0;
			}
			SetNewStatusTo1(pID);
			TotUnassigned = TotUnassigned - 1;
		}
		if(EventType==7){ // End secondary ST relationship
			pID = Register[nextID-1].ID2ndary;
			Register[nextID-1].ID2ndary = 0;
			Register[nextID-1].CurrPartners = 1;
			Register[nextID-1].Disclosed2ndary = 0;
			if(Register[pID-1].CurrPartners==1){
				Register[pID-1].IDprimary = 0;
				Register[pID-1].CurrPartners = 0;
				Register[pID-1].DisclosedPrimary = 0;
			}
			else{
				if(Register[pID-1].IDprimary==nextID){
					Register[pID-1].IDprimary = Register[pID-1].ID2ndary;
					Register[pID-1].CondomPrimary = Register[pID-1].Condom2ndary;
					Register[pID-1].DisclosedPrimary = Register[pID-1].Disclosed2ndary;
				}
				Register[pID-1].ID2ndary = 0;
				Register[pID-1].CurrPartners = 1;
				Register[pID-1].Disclosed2ndary = 0;
			}
			SetNewStatusTo1(pID);
			TotUnassigned = TotUnassigned - 1;
		}
		/*TotUnassigned = 0;
		for(ic=0; ic<Register.size(); ic++){
			if(Register[ic].AliveInd==1 && Register[ic].AgeGroup>=2 && Register[ic].NewStatus==0){
				TotUnassigned += 1;}
		}*/
		ii += 1;
	}
}

void Pop::UpdateMSMcasual()
{
	int ic;
	double EntryRate, ExitRate, ExactAge, TempProb;

	for (ic = 0; ic<MaxPop; ic++){
		r2[ic] = rg.Random();}

	TotCurrCasual = 0;
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 0 && 
			Register[ic].MalePref>0.0 && Register[ic].VirginInd==0 &&
			Register[ic].Imprisoned==0){
			if (Register[ic].CasualInd == 0 && Register[ic].HetCasualInd == 0){
				// Model entry into casual sex
				ExactAge = 1.0 * (CurrYear + 0.5 + (BehavCycleCount - 1.0) /
					CycleS - Register[ic].DOB);
				EntryRate = CasualEntry * pow(CasualAgeAdj, ExactAge - 20.0) *
					Register[ic].MalePref;
				if (Register[ic].RiskGroup == 2){ 
					if (Register[ic].CurrPartners > 0){
						EntryRate = 0;}
					else{ EntryRate *= CasualLowAdj; }
				}
				else{
					if (Register[ic].CurrPartners > 0){
						EntryRate *= CasualHighAdj;}
				}
				if (Register[ic].CurrAge >= 15 && Register[ic].DailyDrinkProb >= 1.0 / 30.0 && 
					Register[ic].DrinksPerDD >= 5.0){
					EntryRate *= RRcasualBinge[0];
				}
				TempProb = 1.0 - exp(-EntryRate / CycleS);
				if (r2[ic] < TempProb){ 
					Register[ic].CasualInd = 1; 
					TotCurrCasual += 1;
					CasualRegister[TotCurrCasual - 1] = ic+1;
				}
			}
			else if (Register[ic].CasualInd == 1){
				if (MSMcalib == 1 || MSMcalibHIV == 1 || FixedUncertainty == 1){
					Register[ic].EverMSM = 1;
					Register[ic].RecentMSM += 1;
				}
				// Model exit from casual sex
				ExitRate = CasualExit / Register[ic].MalePref;
				if (Register[ic].RiskGroup == 2 && Register[ic].CurrPartners>0){
					TempProb = 1.0;}
				else{ TempProb = 1.0 - exp(-ExitRate / CycleS); }
				if (r2[ic] < TempProb){ Register[ic].CasualInd = 0; }
				else{ 
					TotCurrCasual += 1; 
					CasualRegister[TotCurrCasual - 1] = ic+1;
				}
			}
		}
	}
}

void Pop::UpdateHetCasual()
{
	int ic, ig, ir, TempTot, HHID, ExactAge;
	double EntryRate, ExitRate, TempProb, ConcurPenalty;

	for (ic = 0; ic<MaxPop; ic++){
		r2[ic] = rg.Random();
	}

	for (ir = 0; ir < 3; ir++){
		for (ig = 0; ig < 2; ig++){ TotCurrCasualHet[ir][ig] = 0; }
	}
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].VirginInd == 0 &&
			Register[ic].Imprisoned == 0){
			ig = Register[ic].SexInd;
			ir = Register[ic].Race;
			if (Register[ic].CasualInd == 0 && Register[ic].HetCasualInd == 0){
				// Model entry into casual sex
				//ExactAge = 1.0 * (CurrYear + 0.5 + (BehavCycleCount - 1.0) /
				//	CycleS - Register[ic].DOB);
				//if (ExactAge > 90){ ExactAge = 90; }
				//EntryRate = CasualEntryHet[ig] * pow(CasualAgeAdjHet[ig], ExactAge - 20.0);
				ExactAge = Register[ic].CurrAge;
				EntryRate = CasualEntryHet[ig] * AgeEffectPartners[ExactAge - 10][ig] / 
					BasePartnerAcqH[ig];
				if (ig == 0){ EntryRate *= (1.0 - Register[ic].MalePref); }
				if (Register[ic].RiskGroup == 2){
					if (Register[ic].CurrPartners > 0){
						EntryRate = 0.0;
					}
					else{ EntryRate *= CasualLowAdjHet[ig]; }
				}
				else{
					if (Register[ic].CurrPartners > 0){
						ConcurPenalty = CasualHighAdjHet;
						if (ig == 0){ 
							ConcurPenalty *= (1.0 + EffectIneqGenderConcurrency * (Register[ic].IneqGender - 0.2)); 
							if (ConcurPenalty > 1.0){ ConcurPenalty = 1.0; }
							if (ConcurPenalty < 0.0){ ConcurPenalty = 0.0; }
						}
						EntryRate *= ConcurPenalty;
					}
				}
				if (ig == 0){
					EntryRate *= pow(EffectIneqGenderCasual, (0.2 - Register[ic].IneqGender) / 0.1);
				}
				if (Register[ic].CurrAge >=15 && Register[ic].DailyDrinkProb >= 1.0 / 30.0 && 
					Register[ic].DrinksPerDD >= 5.0){
					EntryRate *= RRcasualBinge[ig];
				}
				//if (Register[ic].SexInd == 0 && (Register[ic].Employed == 1 || CurrYear>=2000)){
				if (Register[ic].SexInd == 0 && Register[ic].Employed == 1){
					EntryRate *= RRcasualEmployedM;
				}
				if (Register[ic].SexInd == 1){
					HHID = Register[ic].HouseholdID;
					if (HHID > 0){
						if (HHregister[HHID - 1].PerCapitaIncomeAdj < MeanLogIncome){
							if (HHregister[HHID - 1].PerCapitaIncomeAdj > 0){
								EntryRate *= pow(RRcasualLogDropIncomeF, MeanLogIncome -
									log(HHregister[HHID - 1].PerCapitaIncomeAdj));
							}
							else{
								EntryRate *= pow(RRcasualLogDropIncomeF, MeanLogIncome);
							}
						}
					}
					if (HHID == 0){
						EntryRate *= pow(RRcasualLogDropIncomeF, MeanLogIncome);
					}
				}
				TempProb = 1.0 - exp(-EntryRate / CycleS);
				if (r2[ic] < TempProb){
					Register[ic].HetCasualInd = 1;
					TotCurrCasualHet[ir][ig] += 1;
					TempTot = TotCurrCasualHet[ir][ig];
					if (ig == 0){ CasualRegisterM[TempTot - 1][ir] = ic + 1; }
					else{ CasualRegisterF[TempTot - 1][ir] = ic + 1; }
					Register[ic].EverCasual = 1;
				}
			}
			else if (Register[ic].HetCasualInd == 1){
				// Model exit from casual sex
				ExitRate = CasualExit;
				if (ig == 0 && Register[ic].MalePref > 0.0 && Register[ic].MalePref < 1.0){
					ExitRate = CasualExit / (1.0 - Register[ic].MalePref);
				}
				if ((Register[ic].RiskGroup == 2 && Register[ic].CurrPartners>0) ||
					(ig == 0 && Register[ic].MalePref == 1.0)){
						TempProb = 1.0;
				}
				else{ TempProb = 1.0 - exp(-ExitRate / CycleS); }
				if (r2[ic] < TempProb){ Register[ic].HetCasualInd = 0; }
				else{
					TotCurrCasualHet[ir][ig] += 1;
					TempTot = TotCurrCasualHet[ir][ig];
					if (ig == 0){ CasualRegisterM[TempTot - 1][ir] = ic + 1; }
					else{ CasualRegisterF[TempTot - 1][ir] = ic + 1; }
				}
			}
		}
	}
}

void Pop::GetNewPartner(int ID, double rnd, int rsk)
{
	// Similar to the ChooseSTpartner function.

	int ic, ia, ib;
	int iage, page; // individual's age group - 2 (so that 0 ==> age 10-14), partner's age group - 2
	int igender, pgender; // individual's gender and partner's gender
	int irisk; // individual's risk group -1 (so that 0 ==> high risk)
	//double EligibleByAge[16];
	double rnd3, WeightsByAge[16];
	double Normalizer, TotalSelectionProb2;

	iage = Register[ID-1].AgeGroup - 2;
	igender = Register[ID-1].SexInd;
	pgender = 1 - igender;
	irisk = Register[ID-1].RiskGroup - 1;

	MaxNewPartnerInd = 0;

	// (1) Select partner age
	if(igender==0){
		for(ia=0; ia<16; ia++){
			//if(rnd<CumAgePrefM[iage][ia] && FemPartners[ia][rsk-1].TotalDesire>0.0){
			if (rnd<CumAgePrefM[iage][ia] && FemPartners[ia][rsk - 1].Pool.size()>0){
				page = ia;
				break;
			}
			//if(rnd<CumAgePrefM[iage][ia] && FemPartners[ia][rsk-1].TotalDesire==0.0){
			if (rnd<CumAgePrefM[iage][ia] && FemPartners[ia][rsk - 1].Pool.size()==0){
				// Sample another age group
				if(ia==0){
					rnd3 = rnd/AgePrefM[iage][0];}
				else{
					rnd3 = (rnd - CumAgePrefM[iage][ia-1])/AgePrefM[iage][ia];}
				for(ib=0; ib<16; ib++){
					if(rnd3<CumAgePrefM[iage][ib]){
						//if(FemPartners[ib][rsk-1].TotalDesire>0.0){
						if (FemPartners[ib][rsk - 1].Pool.size()>0){
							page = ib;}
						else{
							// No luck after 2 tries, so give up.
							MaxNewPartnerInd = 1;}
						break;
					}
				}
				break;
			}
		}
	}
	else{
		for(ia=0; ia<16; ia++){
			//if(rnd<CumAgePrefF[iage][ia] && MalePartners[ia][rsk-1].TotalDesire>0.0){
			if (rnd<CumAgePrefF[iage][ia] && MalePartners[ia][rsk - 1].Pool.size()>0){
				page = ia;
				break;
			}
			//if(rnd<CumAgePrefF[iage][ia] && MalePartners[ia][rsk-1].TotalDesire==0.0){
			if (rnd<CumAgePrefF[iage][ia] && MalePartners[ia][rsk - 1].Pool.size()==0){
				// Sample another age group
				if(ia==0){
					rnd3 = rnd/AgePrefF[iage][0];}
				else{
					rnd3 = (rnd - CumAgePrefF[iage][ia-1])/AgePrefF[iage][ia];}
				for(ib=0; ib<16; ib++){
					if(rnd3<CumAgePrefF[iage][ib]){
						//if(MalePartners[ib][rsk-1].TotalDesire>0.0){
						if (MalePartners[ib][rsk - 1].Pool.size()>0){
							page = ib;}
						else{
							// No luck after 2 tries, so give up.
							MaxNewPartnerInd = 1;}
						break;
					}
				}
				break;
			}
		}
	}

	// (2) Select partner ID
	if(MaxNewPartnerInd==0){
		if(igender==0){
			ic = FemPartners[page][rsk-1].SamplePool2(ID);}
		else{
			ic = MalePartners[page][rsk-1].SamplePool2(ID);}
		// Because of rounding errors, there's a small chance that an individual
		// may be selected even if they're not eligible - so check eligibility.
		if(Register[ic-1].NewStatus==1){
			MaxNewPartnerInd = 1;}
		if (MaxNewPartnerInd == 0){
			if(Register[ID-1].IDprimary==0){
				Register[ID-1].IDprimary = ic;}
			else{
				Register[ID-1].ID2ndary = ic;}
			if(Register[ic-1].IDprimary==0){
				Register[ic-1].IDprimary = ID;}
			else{
				Register[ic-1].ID2ndary = ID;}
			Register[ID-1].AssignCondom(0, ID, ic);
			Register[ID-1].CurrPartners += 1;
			Register[ic-1].CurrPartners += 1;
			Register[ID-1].LifetimePartners += 1;
			Register[ic-1].LifetimePartners += 1;
			if (StructuralDriverCalib == 1 || StructuralRCTcalib == 1){ 
				Register[ID - 1].AnnPartners += 1;
				Register[ic - 1].AnnPartners += 1;
			}
			if (MSMcalib == 1 || MSMcalibHIV == 1){
				if (igender == 0 && Register[ID - 1].MalePref > 0.0){
					Register[ID - 1].EverBi = 1;
					Register[ID - 1].RecentBi += 1;
				}
				if (pgender == 0 && Register[ic - 1].MalePref > 0.0){
					Register[ic - 1].EverBi = 1;
					Register[ic - 1].RecentBi += 1;
				}
			}
			if(Register[ID-1].VirginInd==1){
				Register[ID-1].VirginInd = 0;}
			if(Register[ic-1].VirginInd==1){
				Register[ic-1].VirginInd = 0;}
			SetNewStatusTo1(ic);
		}
	}
}

void Pop::GetNewMSMpartner(int ID, double rnd, int rsk)
{
	// Similar to the ChooseSTpartner function.

	int ic, ia, ib;
	int iage, page; // individual's age group - 2 (so that 0 ==> age 10-14), partner's age group - 2
	int irisk; // individual's risk group -1 (so that 0 ==> high risk)
	//double EligibleByAge[16];
	double rnd3, WeightsByAge[16], temp;
	double Normalizer, TotalSelectionProb2;

	iage = Register[ID - 1].AgeGroup - 2;
	irisk = Register[ID - 1].RiskGroup - 1;

	MaxNewPartnerInd = 0;

	// (1) Select partner age
		for (ia = 0; ia<16; ia++){
			if (rnd<CumAgePrefMSM[iage][ia] && MSMpartners[ia][rsk - 1].Pool.size()>0 &&
				// This is to avoid MSM having sex with themselves
				(iage != ia || irisk!=(rsk-1) || MSMpartners[ia][rsk - 1].Pool.size()>1)){ 
				page = ia;
				break;
			}
			if (rnd<CumAgePrefMSM[iage][ia] && (MSMpartners[ia][rsk - 1].Pool.size() == 0 ||
				// This is to avoid MSM having sex with themselves
				(MSMpartners[ia][rsk - 1].Pool.size() == 1 && ia == iage && irisk == (rsk - 1)))){ 
				// Sample another age group
				if (ia == 0){
					rnd3 = rnd / AgePrefMSM[iage][0];
				}
				else{
					rnd3 = (rnd - CumAgePrefMSM[iage][ia - 1]) / AgePrefMSM[iage][ia];
				}
				for (ib = 0; ib<16; ib++){
					if (rnd3<CumAgePrefMSM[iage][ib]){
						if (MSMpartners[ib][rsk - 1].Pool.size()>0 && 
							(iage != ib || irisk!=(rsk-1) || MSMpartners[ib][rsk - 1].Pool.size()>1)){
							page = ib;
						}
						else{
							// No luck after 2 tries, so give up.
							MaxNewPartnerInd = 1;
						}
						break;
					}
				}
				break;
			}
		}
	

	// (2) Select partner ID
	if (MaxNewPartnerInd == 0){
		ic = MSMpartners[page][rsk - 1].SamplePool2(ID);
		// Because of rounding errors, there's a small chance that an individual
		// may be selected even if they're not eligible - so check eligibility.
		if (Register[ic - 1].NewStatus == 1 || ic==ID){
			MaxNewPartnerInd = 1;
		}
		if(MaxNewPartnerInd == 0){
			if (Register[ID - 1].IDprimary == 0){
				Register[ID - 1].IDprimary = ic;}
			else{
				Register[ID - 1].ID2ndary = ic;}
			if (Register[ic - 1].IDprimary == 0){
				Register[ic - 1].IDprimary = ID;}
			else{
				Register[ic - 1].ID2ndary = ID;}
			Register[ID - 1].AssignCondom(0, ID, ic);
			Register[ID - 1].CurrPartners += 1;
			Register[ic - 1].CurrPartners += 1;
			Register[ID - 1].LifetimePartners += 1;
			Register[ic - 1].LifetimePartners += 1;
			if (StructuralDriverCalib == 1 || StructuralRCTcalib == 1){
				Register[ID - 1].AnnPartners += 1;
				Register[ic - 1].AnnPartners += 1;
			}
			if (MSMcalib == 1 || MSMcalibHIV == 1 || FixedUncertainty == 1){
				Register[ID - 1].EverMSM = 1;
				Register[ic - 1].EverMSM = 1;
				Register[ID - 1].RecentMSM = +1;
				Register[ic - 1].RecentMSM = +1;
			}
			if (Register[ID - 1].VirginInd == 1){
				Register[ID - 1].VirginInd = 0;
			}
			if (Register[ic - 1].VirginInd == 1){
				Register[ic - 1].VirginInd = 0;
			}
			SetNewStatusTo1(ic);
		}
	}
}

void Pop::SetNewStatusTo1(int ID)
{
	int ia, ir, ii;
	double Decrement;
	
	Register[ID-1].NewStatus = 1;
	Decrement = Register[ID-1].DesiredNewPartners;
	if(Decrement>0.0 && Register[ID-1].VirginInd==0){
		ia = Register[ID-1].AgeGroup - 2;
		ir = Register[ID-1].RiskGroup - 1;
		if(Register[ID-1].SexInd==0){
			MalePartners[ia][ir].TotalDesire = 
				MalePartners[ia][ir].TotalDesire - Decrement *
				(1.0 - Register[ID-1].MalePref);
			MSMpartners[ia][ir].TotalDesire =
				MSMpartners[ia][ir].TotalDesire - Decrement *
				Register[ID - 1].MalePref;
		}
		else{
			FemPartners[ia][ir].TotalDesire = 
				FemPartners[ia][ir].TotalDesire - Decrement;
		}
	}
}

void Pop::SetToDead(int ID, int Cause)
{
	// Differs from SetNewStatusTo 1 in that (a) we have to adjust partner profile(s)
	// but (b) we don't have to correct EligibleByAge (since it hasn't yet been calculated).

	int ic, ij, pID1, pID2, HHID, ir, ia, ig;

	Register[ID-1].AliveInd = 0;
	Register[ID - 1].DOD = 1.0 * (CurrYear + 0.5 + 1.0 * BehavCycleCount / CycleS);
	Register[ID-1].NewStatus = 1;
	//HHID = Register[ID - 1].HouseholdID;
	//if (HHID > 0){ HHregister[HHID - 1].RemoveMember(ID); }
	if(Register[ID-1].IDprimary>0){
		pID1 = Register[ID-1].IDprimary;
		if(Register[pID1-1].CurrPartners==1){
			Register[pID1-1].IDprimary = 0;
			Register[pID1-1].CurrPartners = 0;
		}
		else{
			if(Register[pID1-1].IDprimary==ID){
				Register[pID1-1].IDprimary = Register[pID1-1].ID2ndary;}
			Register[pID1-1].ID2ndary = 0;
			Register[pID1-1].CurrPartners = 1;
		}
		Register[pID1-1].NewStatus = 1;
		if(Register[ID-1].MarriedInd==1){
			Register[pID1-1].MarriedInd = 0;
			Register[pID1-1].DOUD = Register[ID - 1].DOD;
		}
	}
	if(Register[ID-1].ID2ndary>0){
		pID2 = Register[ID-1].ID2ndary;
		if(Register[pID2-1].CurrPartners==1){
			Register[pID2-1].IDprimary = 0;
			Register[pID2-1].CurrPartners = 0;
		}
		else{
			if(Register[pID2-1].IDprimary==ID){
				Register[pID2-1].IDprimary = Register[pID2-1].ID2ndary;}
			Register[pID2-1].ID2ndary = 0;
			Register[pID2-1].CurrPartners = 1;
		}
		Register[pID2-1].NewStatus = 1;
	}
	if(Register[ID-1].FSWind==1){
		ir = Register[ID - 1].Race;
		TotCurrFSW[ir] = TotCurrFSW[ir] - 1;
		// Update CSWregister
		for(ic=0; ic<MaxCSWs; ic++){
			if (CSWregister[ic][ir] == ID){
				for(ij=ic; ij<MaxCSWs-1; ij++){
					CSWregister[ij][ir] = CSWregister[ij + 1][ir];}
				CSWregister[MaxCSWs - 1][ir] = 0;
			}
		}
	}
	if (Register[ID - 1].CasualInd == 1){
		TotCurrCasual = TotCurrCasual - 1;
		// Update CasualRegister
		for (ic = 0; ic<MaxCasual; ic++){
			if (CasualRegister[ic] == ID){
				for (ij = ic; ij<MaxCasual - 1; ij++){
					CasualRegister[ij] = CasualRegister[ij + 1];
				}
				CasualRegister[MaxCasual-1] = 0;
			}
		}
	}
	if (Register[ID - 1].HetCasualInd == 1){
		ig = Register[ID - 1].SexInd;
		ir = Register[ID - 1].Race;
		TotCurrCasualHet[ir][ig] = TotCurrCasualHet[ir][ig] - 1;
		// Update CasualRegister
		for (ic = 0; ic<MaxCasual; ic++){
			if (CasualRegisterM[ic][ir] == ID){
				for (ij = ic; ij<MaxCasual - 1; ij++){
					CasualRegisterM[ij][ir] = CasualRegisterM[ij + 1][ir];
				}
				CasualRegisterM[MaxCasual - 1][ir] = 0;
			}
			if (CasualRegisterF[ic][ir] == ID){
				for (ij = ic; ij<MaxCasual - 1; ij++){
					CasualRegisterF[ij][ir] = CasualRegisterF[ij + 1][ir];
				}
				CasualRegisterF[MaxCasual - 1][ir] = 0;
			}
		}
	}
	HHID = Register[ID - 1].HouseholdID;
	if (HHID > 0){ HHregister[HHID - 1].RemoveMember(ID); }
	if (FixedUncertainty == 1 && Cause==1){
		AIDSdeaths.out[CurrSim - 1][CurrYear - StartYear] += Register[ID-1].PopWeight;
		ia = Register[ID - 1].CurrAge;
		ig = Register[ID - 1].SexInd;
		if (ia >= 15){
			LYsLost.out[CurrSim - 1][CurrYear - StartYear] += LE_West26[ia][ig] * Register[ID - 1].PopWeight;
		}
	}
	if (FixedUncertainty == 1 && Cause == 0){
		NonAIDSdeaths.out[CurrSim - 1][CurrYear - StartYear] += Register[ID - 1].PopWeight;
	}
}

void Pop::NewBirth(int ID)
{
	int ic, popsize, ind, MatAgeGrp, s, iy, ir, PID1, PID2, ig, FatherKids;
	Indiv KidA;
	double r[25];
	double x, y, a, b, p, q, bound, DurBF, ScreeningProb, HRprob, TempOdds, AveRisk;

	memset(&KidA, 0, sizeof(KidA));
	//int seed = 4475 + CurrYear * 62 + BehavCycleCount * 98 + ID * 23;
	//if(CurrYear>=StartYear+FixedPeriod){
	//	seed += CurrSim;}
	//CRandomMersenne rg(seed);
	for(ic=0; ic<25; ic++){
		r[ic] = rg.Random();}
	/* The random numbers are used as follows:
	r[0] determines sex
	r[1] determines high/low sexual risk
	r[2] determines duration of breastfeeding
	r[3] determines if child is vertically infected
	r[4] determines relative rate of partner acquisition (N/A in most simulations)
	r[5] determines HIV susceptibility adjustment 
	r[6] determines same-sex preference
	r[7] determines extent of change in same-sex preference with age
	r[8] determines role preference in MSM relationships
	r[9] determines whether mother tests for HIV
	r[10] determines male circumcision (for infants circumcised)
	r[11] determines maternal testing for syphilis
	r[12] determines condom preference
	r[13] determines urban/rural location at birth
	r[14] determines desired visit frequency
	r[15] determines whether mother's partner tests for HIV
	r[16] determines conscientiousness
	r[17] determines inequitable gender norms
	r[18] determines daily drinking prob
	r[19] determines number of drinks per drinking day
	r[20] determines log income difference
	r[21] determines whether we use r[1] or r[19] to assign drinks per DD*/
	iy = CurrYear - StartYear;
	TotBirths.out[CurrSim - 1][iy] += Register[ID - 1].PopWeight;
	if (FixedUncertainty == 1 && Register[ID - 1].DOB > (1.0 * CurrYear + 0.5 + 1.0 * BehavCycleCount/CycleS - 20.0)){
		TeenBirths.out[CurrSim - 1][iy] += Register[ID - 1].PopWeight; }

	KidA.AliveInd = 1;
	KidA.PopWeight = 1.0;
	KidA.VirginInd = 1;
	KidA.MarriedInd = 0;
	KidA.CurrPartners = 0;
	KidA.FSWind = 0;
	KidA.IDprimary = 0;
	KidA.ID2ndary = 0;
	KidA.NonHIVfertRate = 0.0;
	KidA.Fecundability = 0.0;
	KidA.Visiting = 0;
	KidA.Imprisoned = 0;
	KidA.ReleaseDate = 0.0; 
	KidA.PrevImprisoned = 0;
	KidA.Employed = 0;
	KidA.DateMig = 0.0;
	KidA.OnPrEP = 0.0;
	KidA.VaccineEff = 0.0;
	KidA.Date1stVacc = 0.0;
	if (Register[ID - 1].ChildIDs[19] == 0){ KidA.ParentID[1] = ID; }
	else{ KidA.ParentID[1] = 0; }
	KidA.ParentID[0] = 0; // Default - then check if there's a father
	PID1 = Register[ID - 1].FatherBaby;
	if (PID1 > 0){
		FatherKids = 0;
		while (Register[PID1 - 1].ChildIDs[FatherKids] > 0){ FatherKids += 1; }
		if (FatherKids < 20){ KidA.ParentID[0] = PID1; }
	}
	KidA.DailyDrinkProb = 0.0;
	KidA.DrinksPerDD = 1.0;

	// Assign race
	if (PID1 == 0){ KidA.Race = Register[ID - 1].Race; }
	else if (Register[ID - 1].Race == Register[PID1 - 1].Race){
		KidA.Race = Register[ID - 1].Race; 
	}
	else{ KidA.Race = 1; } // Mixed race/coloured
	ir = KidA.Race;
	
	// Assign personality
	ind = 2;
	a = 0.0;
	b = 1.0;
	p = r[16];
	q = 1.0 - r[16];
	s = 0;
	bound = 0.0;
	cdfnor(&ind, &p, &q, &x, &a, &b, &s, &bound);
	KidA.Conscientiousness = x;
	
	// Assign sex, risk group, NonHIVmortProb and initial STD states
	if(r[0]<MaleBirthPropn[ir]){
		KidA.SexInd = 0;
		HRprob = 1.0 / (1.0 + ((1.0 - HighPropnM) / HighPropnM) / pow(ConscientiousEffectSex, 
			KidA.Conscientiousness));
		if(r[1]<HRprob){KidA.RiskGroup = 1;}
		else{KidA.RiskGroup = 2;}
		KidA.NonHIVmortProb = 1.0 - pow(1.0 - 
			InfantMort1st6mM[iy][ir], 2.0 / CycleS);
	}
	else{
		KidA.SexInd = 1;
		HRprob = 1.0 / (1.0 + ((1.0 - HighPropnF) / HighPropnF) / pow(ConscientiousEffectSex, 
			KidA.Conscientiousness));
		if(r[1]<HighPropnF){KidA.RiskGroup = 1;}
		else{KidA.RiskGroup = 2;}
		KidA.NonHIVmortProb = 1.0 - pow(1.0 - 
			InfantMort1st6mF[iy][ir], 2.0/CycleS);
	}
	ig = KidA.SexInd;
	KidA.CTstage = 0;
	KidA.HDstage = 0;
	KidA.HSVstage = 0;
	KidA.NGstage = 0;
	KidA.TPstage = 0;
	KidA.TVstage = 0;
	KidA.BVstage = 0;
	KidA.VCstage = 0;
	
	// Assign exact date of birth and age group
	KidA.DOB = Register[ID - 1].DOLB;
	KidA.AgeGroup = 0;
	KidA.CurrAge = 0;
	
	// Handle situation in which mother is HIV-positive
	MatAgeGrp = Register[ID - 1].AgeGroup - 3;
	ScreeningProb = PMTCTaccess[iy] * AcceptScreening;
	if (Register[ID - 1].VCThistory == 2){ ScreeningProb *= RetestAdjDiagnosed[2]; }
	if (Register[ID - 1].HIVstage == 5){ ScreeningProb *= RetestAdjART[2]; }
	if (r[9] < ScreeningProb){
		TotalTestsANC.out[CurrSim - 1][CurrYear - StartYear] += Register[ID - 1].PopWeight;
		if (Register[ID - 1].VCThistory < 2){ Register[ID - 1].VCThistory = 1; }
		if (Register[ID - 1].VCThistory == 2){ 
			Register[ID - 1].GetRediagnosed(ID, 2);
			PosTestsANC.out[CurrSim - 1][CurrYear - StartYear] += Register[ID - 1].PopWeight;
		}
		if (Register[ID - 1].HIVstage == 1){ AcuteTests.out[CurrSim - 1][iy] += Register[ID - 1].PopWeight; }
	}
	if (Register[ID - 1].HIVstage>0){
		if (FixedUncertainty == 1){ 
			TotBirthsHIV.out[CurrSim - 1][CurrYear - StartYear] += Register[ID - 1].PopWeight;
			if (Register[ID - 1].HIVstage == 5){ TotBirthsART.out[CurrSim - 1][CurrYear - StartYear] += Register[ID - 1].PopWeight; }
		}
		if (r[9] < PMTCTaccess[iy] * AcceptScreening){
			// As in Thembisa, the prob is independent of prior HIV diagnosis. But all of this code 
			// needs to be updated to take into account regimens other than sd NVP.
			CurrPerinatal = PropnInfectedAtBirth * (1.0 - AcceptNVP * RednNVP);
			CurrPostnatal = (1.0 - CurrPerinatal) * PropnInfectedAfterBirth *
				(1.0 - RednFF);
			if (Register[ID - 1].VCThistory == 1 && Register[ID - 1].HIVstage > 1){ 
				Register[ID - 1].GetDiagnosed(ID, 2);
				PosTestsANC.out[CurrSim - 1][CurrYear - StartYear] += Register[ID - 1].PopWeight;
				NewDiagANC.out[CurrSim - 1][CurrYear - StartYear] += Register[ID - 1].PopWeight;
			}
		}
		else{
			CurrPerinatal = PropnInfectedAtBirth;
			CurrPostnatal = (1.0 - CurrPerinatal) * PropnInfectedAfterBirth;
		}
		if(r[3]<CurrPerinatal){
			KidA.HIVstage = 2;
			KidA.CD4 = 300.0;
			KidA.logVL = 4.8;
			KidA.DateInfect = KidA.DOB;
		}
		else if(r[3]<CurrPerinatal+CurrPostnatal){
			KidA.HIVstage = 2;
			KidA.CD4 = 300.0;
			KidA.logVL = 4.8;
			// We arbitrarily assume age of postnatal infection is 6 months.
			// At this stage this is only used for determining if the transmission
			// was perinatal or postnatal, so the exact value is not important.
			KidA.DateInfect = KidA.DOB + 0.5;
		}
		else{
			KidA.HIVstage = 0;}
	}
	// Test whether male partner of pregnant woman is tested
	if (r[9] < ScreeningProb && CurrSim!=52 && Register.size() != 33042){
		if (Register[ID - 1].CurrPartners > 0){
			// Test whether primary partner is tested
			ir = Register[ID - 1].MarriedInd;
			PID1 = Register[ID - 1].IDprimary;
			if (Register[ID - 1].VCThistory == 1){ // Woman was diagnosed negative
				ScreeningProb = ANCpartnersInvited[iy][0] * ANCpartnerTested[0][ir];
				if (Register[PID1 - 1].VCThistory == 2){ ScreeningProb *= RetestAdjDiagnosed[11]; }
				if (Register[PID1 - 1].HIVstage == 5){ ScreeningProb *= RetestAdjART[11]; }
				if (r[15] < ScreeningProb){
					TotalTestsANCpartner0.out[CurrSim - 1][iy] += Register[PID1 - 1].PopWeight;
					if (Register[PID1 - 1].VCThistory == 2){ 
						Register[PID1 - 1].GetRediagnosed(PID1, 0);
						PosTestsANCpartner0.out[CurrSim - 1][iy] += Register[PID1 - 1].PopWeight;
					}
					else{ Register[PID1 - 1].VCThistory = 1; }
					if (Register[PID1 - 1].HIVstage > 1 && Register[PID1 - 1].VCThistory == 1){
						Register[PID1 - 1].GetDiagnosed(PID1, 0); 
						PosTestsANCpartner0.out[CurrSim - 1][iy] += Register[PID1 - 1].PopWeight;
						NewDiagANCpartner0.out[CurrSim - 1][iy] += Register[PID1 - 1].PopWeight;
					}
					if (Register[PID1 - 1].HIVstage == 1){ AcuteTests.out[CurrSim - 1][iy] += Register[PID1 - 1].PopWeight; }
				}
			}
			else{ // Woman was diagnosed positive
				ScreeningProb = ANCpartnersInvited[iy][1] * ANCpartnerTested[1][ir];
				if (Register[PID1 - 1].VCThistory == 2){ ScreeningProb *= RetestAdjDiagnosed[11]; }
				if (Register[PID1 - 1].HIVstage == 5){ ScreeningProb *= RetestAdjART[11]; }
				if (r[15] < ScreeningProb){
					TotalTestsANCpartner1.out[CurrSim - 1][iy] += Register[PID1 - 1].PopWeight;
					if (Register[PID1 - 1].VCThistory == 2){ 
						Register[PID1 - 1].GetRediagnosed(PID1, 0);
						PosTestsANCpartner1.out[CurrSim - 1][iy] += Register[PID1 - 1].PopWeight;
					}
					else{ Register[PID1 - 1].VCThistory = 1; }
					if (Register[PID1 - 1].HIVstage > 1 && Register[PID1 - 1].VCThistory == 1){
						Register[PID1 - 1].GetDiagnosed(PID1, 0); 
						PosTestsANCpartner1.out[CurrSim - 1][iy] += Register[PID1 - 1].PopWeight;
						NewDiagANCpartner1.out[CurrSim - 1][iy] += Register[PID1 - 1].PopWeight;
					}
					if (Register[PID1 - 1].HIVstage == 1){ AcuteTests.out[CurrSim - 1][iy] += Register[PID1 - 1].PopWeight; }
				}
			}
		}
		if (Register[ID - 1].CurrPartners == 2){
			// Test whether 2ndary partner is tested
			PID2 = Register[ID - 1].ID2ndary;
			if (Register[ID - 1].VCThistory == 1){ // Woman was diagnosed negative
				ScreeningProb = ANCpartnersInvited[iy][0] * ANCpartnerTested[0][0];
				if (Register[PID2 - 1].VCThistory == 2){ ScreeningProb *= RetestAdjDiagnosed[11]; }
				if (Register[PID2 - 1].HIVstage == 5){ ScreeningProb *= RetestAdjART[11]; }
				if (r[15] < ScreeningProb){
					TotalTestsANCpartner0.out[CurrSim - 1][iy] += Register[PID2 - 1].PopWeight;
					if (Register[PID2 - 1].VCThistory == 2){ 
						Register[PID2 - 1].GetRediagnosed(PID2, 0);
						PosTestsANCpartner0.out[CurrSim - 1][iy] += Register[PID2 - 1].PopWeight;
					}
					else{ Register[PID2 - 1].VCThistory = 1; }
					if (Register[PID2 - 1].HIVstage > 1 && Register[PID2 - 1].VCThistory==1){
						Register[PID2 - 1].GetDiagnosed(PID2, 0); 
						PosTestsANCpartner0.out[CurrSim - 1][iy] += Register[PID2 - 1].PopWeight;
						NewDiagANCpartner0.out[CurrSim - 1][iy] += Register[PID2 - 1].PopWeight;
					}
					if (Register[PID2 - 1].HIVstage == 1){ AcuteTests.out[CurrSim - 1][iy] += Register[PID2 - 1].PopWeight; }
				}
			}
			else{ // Woman was diagnosed positive
				ScreeningProb = ANCpartnersInvited[iy][1] * ANCpartnerTested[1][0];
				if (Register[PID2 - 1].VCThistory == 2){ ScreeningProb *= RetestAdjDiagnosed[11]; }
				if (Register[PID2 - 1].HIVstage == 5){ ScreeningProb *= RetestAdjART[11]; }
				if (r[15] < ScreeningProb){
					TotalTestsANCpartner1.out[CurrSim - 1][iy] += Register[PID2 - 1].PopWeight;
					if (Register[PID2 - 1].VCThistory == 2){ 
						Register[PID2 - 1].GetRediagnosed(PID2, 0);
						PosTestsANCpartner1.out[CurrSim - 1][iy] += Register[PID2 - 1].PopWeight;
					}
					else{ Register[PID2 - 1].VCThistory = 1; }
					if (Register[PID2 - 1].HIVstage > 1 && Register[PID2 - 1].VCThistory == 1){
						Register[PID2 - 1].GetDiagnosed(PID2, 0); 
						PosTestsANCpartner1.out[CurrSim - 1][iy] += Register[PID2 - 1].PopWeight;
						NewDiagANCpartner1.out[CurrSim - 1][iy] += Register[PID2 - 1].PopWeight;
					}
					if (Register[PID2 - 1].HIVstage == 1){ AcuteTests.out[CurrSim - 1][iy] += Register[PID2 - 1].PopWeight; }
				}
			}
		}
	}
	
	// Set date of weaning
	// Note that we previously calculated the date of weaning before maternal HIV diagnosis was modelled,
	// but since maternal diagnosis determines the rate of weaning, the order of calculation has changed.
	if (Register[ID - 1].VCThistory < 2){
		if (r[2] > EverFeed){ DurBF = 0.0; }
		else{
			DurBF = MedianFeed[0] * pow(log(r[2] / EverFeed) / log(0.5), 1.0 / ShapeFeed[0]) / 12.0;}
	}
	else{
		if (r[2] > (1.0 - HIVdiagNoBF[iy])){ DurBF = 0.0; }
		else{
			DurBF = MedianFeed[1] * pow(log(r[2] / (1.0 - HIVdiagNoBF[iy])) / log(0.5), 1.0 / ShapeFeed[1]) / 12.0;}
	}
	Register[ID - 1].DOLW = Register[ID - 1].DOLB + DurBF;
	
	// Assign PartnerRateAdj
	if(AllowPartnerRateAdj==1){
		ind = 2;
		// Note that the following formulas for a and b apply only when the gamma mean is 1.
		a = 1.0/pow(SDpartnerRateAdj, 2.0);
		b = a;
		p = r[4];
		q = 1 - r[4];
		cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
		KidA.PartnerRateAdj = x;
	}
	else{
		KidA.PartnerRateAdj = 1.0;}
	
	// Assign SuscepHIVadj
	if(AllowHIVsuscepAdj==1){
		ind = 2;
		// Note that the following formulas for a and b apply only when the gamma mean is 1.
		a = 1.0/pow(SDsuscepHIVadj, 2.0);
		b = a;
		p = r[5];
		q = 1 - r[5];
		cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
		KidA.SuscepHIVadj = x;
	}
	else{
		KidA.SuscepHIVadj = 1.0;}
	
	// Assign MSM characteristics
	KidA.CasualInd = 0;
	KidA.HetCasualInd = 0;
	KidA.InsertivePref = 0.0;
	if (KidA.SexInd == 0){
		if (r[6]>MSMfraction){ KidA.MalePref = 0.0; }
		else{
			if (r[6] > (MSMfraction*BiFraction)){ 
				KidA.MalePref = 1.0; 
				if (r[8] < RolePrefHomo[0]){ KidA.InsertivePref = 0.0; }
				else if (r[8] < RolePrefHomo[0] + RolePrefHomo[1]){
					KidA.InsertivePref = 0.5;
				}
				else{ KidA.InsertivePref = 1.0; }
			}
			else{
				// Assign initial male pref
				a = InitMalePrefBeta[0];
				b = InitMalePrefBeta[1];
				p = r[6] / (MSMfraction*BiFraction);
				q = 1 - p;
				bound = 0.0;
				ind = 2;
				cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
				KidA.MalePref = x;
				// Assign annual change in male pref
				a = AnnChangeMalePref[0];
				b = AnnChangeMalePref[1];
				p = r[7];
				q = 1.0 - r[7];
				s = 0;
				cdfnor(&ind, &p, &q, &x, &a, &b, &s, &bound);
				KidA.ChangeMalePref = x;
				// Assign role preference
				if (r[8] < RolePrefBi[0]){ KidA.InsertivePref = 0.0; }
				else if (r[8] < RolePrefBi[0] + RolePrefBi[1]){
					KidA.InsertivePref = 0.5;
				}
				else{ KidA.InsertivePref = 1.0; }
			}
		}
	}
	else{ KidA.MalePref = 1.0; }
	
	// Assign circumcision status
	if (KidA.SexInd == 0){
		if (r[10] < MCprevBaseline[0][KidA.Race]){KidA.CircInd = 1;}
		else{ KidA.CircInd = 0; }
	}
	
	// Handle situation in which mother has syphilis
	if (TPind == 1 && Register[ID - 1].TPstage>1 && Register[ID - 1].TPstage < 5){
		if (Register[ID - 1].TPstage == 4 && r[11] < TPtransitionF.ANCpropnCured){
			Register[ID - 1].TPstage = 6;}
		if (Register[ID - 1].TPstage == 3 && r[11] < TPtransitionF.ANCpropnCured){
			Register[ID - 1].TPstage = 5;}
		if (Register[ID - 1].TPstage == 2){
			if (r[11] < TPtransitionF.ANCpropnCured * TPtransitionF.PropnSuscepAfterRx){
				Register[ID - 1].TPstage = 0;}
			else if (r[11] < TPtransitionF.ANCpropnCured){
				Register[ID - 1].TPstage = 5;}
		}
	}

	// Assign CondomPref
	ind = 2;
	// Note that the following formulas for a and b apply only when the gamma mean is 1.
	a = 1.0 / pow(SDcondomPref, 2.0);
	b = a;
	p = r[12];
	q = 1 - r[12];
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
	KidA.CondomPref = x;
	ir = KidA.Race;
	//KidA.CondomPref *= CondomRace[ir] * CondomEdu[0];
	KidA.CondomPref *= CondomRace[ir];
	
	// Assign urban/rural indicator
	// This is now redundant because the AssignChildToHH function called at the end will
	// re-assign them to an urban/rural location depending on which household they join.
	KidA.CurrUrban = Register[ID - 1].CurrUrban;
	if (KidA.CurrUrban == 1 && KidA.Race == 0){
		if (r[13] < 0.1){ KidA.CurrUrban = 0; }
	}
	
	// Assign desired visit frequency
	a = VisitFreq[0];
	b = VisitFreq[1];
	p = r[14];
	q = 1 - r[14];
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
	KidA.VisitFreq = x;
	
	// r[15] is used for testing whether male partner(s) get tested.

	// Assign inequitable gender norm probability
	if (KidA.SexInd == 0){
		TempOdds = exp(-r[17]) * RaceEffectIneqGender[ir];
		AveRisk = HighPropnM + (1.0 - HighPropnM) * 2.0;
		TempOdds *= pow(HighRiskEffectIneqGender, AveRisk - KidA.RiskGroup);
		KidA.IneqGender = TempOdds / (1.0 + TempOdds);
	}
	else{ KidA.IneqGender = 0.0; } // Not defined for females
	
	// Assign constant component of daily drinking probability
	a = 0;
	b = StdDevDrinkProb[ig]/(7.0 * AlcPropnReported);
	if (r[21] < ConfoundingAlcSex){
		p = 1.0 - r[1]; // Note we reverse the order because low r implies high risk
		q = r[1];
	}
	else{
		p = r[18];
		q = 1.0 - r[18];
	}
	s = 0;
	cdfnor(&ind, &p, &q, &x, &a, &b, &s, &bound);
	KidA.DrinkProbConstant = x;
	
	// Assign constant component of drinks per drinking day
	a = 0;
	b = StdDevDrinksPerDD[ig] / AlcPropnReported;
	if (r[21] < ConfoundingAlcSex){
		p = 1.0 - r[1]; // Note we reverse the order because low r implies high risk
		q = r[1];
	}
	else{
		p = r[19];
		q = 1.0 - r[19];
	}
	s = 0;
	cdfnor(&ind, &p, &q, &x, &a, &b, &s, &bound);
	KidA.DrinksPerDDconstant = x;
	
	// Assign income potential
	ind = 2;
	a = 0.0;
	b = 1.0;
	p = r[20];
	q = 1.0 - r[20];
	s = 0;
	bound = 0.0;
	cdfnor(&ind, &p, &q, &x, &a, &b, &s, &bound);
	KidA.LogIncomeDif = x * StdDevLogIncome;
	KidA.LogIncomeDif += KidA.Conscientiousness * ConscientiousEffectLogIncome;
	
	Register.push_back(KidA);
	popsize = Register.size(); // the ID for the new infant
	if (Register[ID-1].ChildIDs[19] == 0){ Register[ID - 1].UpdateChildren(popsize); }
	PID1 = Register[ID - 1].FatherBaby;
	if (PID1>0 && FatherKids < 20){ Register[PID1 - 1].UpdateChildren(popsize); }
	Register[popsize - 1].AssignChildToHH(popsize);
}

void Pop::NewHousehold(int ID)
{
	HouseholdGroup NewHH(ID); 
	// The constructor automatically assigns initial values based on ID.
	int NumberHHs;

	HHregister.push_back(NewHH);
	NumberHHs = HHregister.size();
	Register[ID - 1].HouseholdID = NumberHHs;
}

void Pop::ChangeHHafterMarriage(int ID1, int ID2)
{
	int ii, ia, ig, ir, ih;
	int HHID, OldHH[2], PrevHH, NewHH, CurrHead[2], TempID;
	double OddsHeadship[2];

	// Determine household headship just prior to marriage
	// If an individual's household is already at its maximum size, and they are head,
	// we treat them as if they are no longer head, on the assumption they move out.
	CurrHead[0] = 0;
	CurrHead[1] = 0;
	HHID = Register[ID1 - 1].HouseholdID;
	OldHH[0] = HHID;
	if (HHID > 0){
		if (HHregister[HHID - 1].IDhead == ID1 && HHregister[HHID - 1].Size < MaxHHsize){ 
			CurrHead[0] = 1; 
		}
	}
	HHID = Register[ID2 - 1].HouseholdID;
	OldHH[1] = HHID;
	if (HHID > 0){
		if (HHregister[HHID - 1].IDhead == ID2 && HHregister[HHID - 1].Size < MaxHHsize){ 
			CurrHead[1] = 1; 
		}
	}

	// ID2 moves to household of ID1
	if (CurrHead[0] == 1 && CurrHead[1] == 0){
		NewHH = OldHH[0];
		HHregister[NewHH - 1].AddMember(ID2);
		Register[ID2 - 1].HouseholdID = NewHH;
		if (OldHH[1] > 0){
			PrevHH = OldHH[1];
			if (PrevHH > 0){ HHregister[PrevHH - 1].RemoveMember(ID2); }
		}
	}

	// ID1 moves to household of ID2
	else if (CurrHead[0] == 0 && CurrHead[1] == 1){
		NewHH = OldHH[1];
		HHregister[NewHH - 1].AddMember(ID1);
		Register[ID1 - 1].HouseholdID = NewHH;
		if (OldHH[0] > 0){
			PrevHH = OldHH[0];
			if (PrevHH > 0){ HHregister[PrevHH - 1].RemoveMember(ID1); }
		}
	}

	// Otherwise determine which member is head of the new household
	else{
		for (ii = 0; ii < 2; ii++){
			if (ii == 0){ TempID = ID1; }
			else{ TempID = ID2; }
			ig = Register[TempID - 1].SexInd;
			ia = Register[TempID - 1].AgeGroup - 2;
			if (ia > 10){ ia = 10; }
			ir = Register[TempID - 1].Race;
			ih = 0;
			if (Register[TempID - 1].HighestGrade > 0 && Register[TempID - 1].HighestGrade <= 7){ ih = 1; }
			if (Register[TempID - 1].HighestGrade > 7 && Register[TempID - 1].HighestGrade < 12){ ih = 2; }
			if (Register[TempID - 1].HighestGrade == 12){ ih = 3; }
			if (Register[TempID - 1].HighestGrade == 13){ ih = 4; }
			OddsHeadship[ii] = BaseOddsHead[ig] * MarriedEffectHead[ig];
			OddsHeadship[ii] *= AgeEffectHead[ia][ig];
			OddsHeadship[ii] *= RaceEffectHead[ir][ig];
			if (Register[TempID - 1].Employed == 1){ OddsHeadship[ii] *= EmployedEffectHead[ig]; }
			if (Register[TempID - 1].InSchool == 1){ OddsHeadship[ii] *= InSchoolEffectHead[ig]; }
			OddsHeadship[ii] *= EduEffectHead[ih][ig];
			if (ig == 1 && Register[TempID - 1].DOLB > 0.0){ OddsHeadship[ii] *= BirthEffectHead; }
		}
		if (CurrHead[0] == 0 && CurrHead[1] == 0){
			// A new household is formed
			if (OddsHeadship[0] >= OddsHeadship[1]){ 
				NewHousehold(ID1); 
				if (OldHH[0] > 0){
					PrevHH = OldHH[0];
					if (PrevHH > 0){ HHregister[PrevHH - 1].RemoveMember(ID1); }
				}
				NewHH = Register[ID1-1].HouseholdID;
				HHregister[NewHH - 1].AddMember(ID2);
				Register[ID2 - 1].HouseholdID = NewHH;
				if (OldHH[1] > 0){
					PrevHH = OldHH[1];
					if (PrevHH > 0){ HHregister[PrevHH - 1].RemoveMember(ID2); }
				}
			}
			else{
				NewHousehold(ID2);
				if (OldHH[1] > 0){
					PrevHH = OldHH[1];
					if (PrevHH > 0){ HHregister[PrevHH - 1].RemoveMember(ID2); }
				}
				NewHH = Register[ID2 - 1].HouseholdID;
				HHregister[NewHH - 1].AddMember(ID1);
				Register[ID1 - 1].HouseholdID = NewHH;
				if (OldHH[0] > 0){
					PrevHH = OldHH[0];
					if (PrevHH > 0){ HHregister[PrevHH - 1].RemoveMember(ID1); }
				}
			}
		}
		else{
			// CurrHead[0] = 1 and CurrHead[1] = 1
			if (OddsHeadship[0] >= OddsHeadship[1]){
				NewHH = OldHH[0];
				HHregister[NewHH - 1].AddMember(ID2);
				Register[ID2 - 1].HouseholdID = NewHH;			
				PrevHH = OldHH[1];
				if (PrevHH > 0){ HHregister[PrevHH - 1].RemoveMember(ID2); }
			}
			else{
				NewHH = OldHH[1];
				HHregister[NewHH - 1].AddMember(ID1);
				Register[ID1 - 1].HouseholdID = NewHH;			
				PrevHH = OldHH[0];
				if (PrevHH > 0){ HHregister[PrevHH - 1].RemoveMember(ID1); }
			}
		}
	}
}

void Pop::ChangeHHafterDivorce(int ID1, int ID2)
{
	int ir, ia, MAlive, FAlive;
	int LeavingID, CurrHH, CurrHead, NewHH, TempID;
	double OddsNewHH, ProbNewHH, RandHH;

	// Code only needs to be called if household members were already in the same home
	// (might not be the case for migrant couples or where one/both are homeless).
	
	CurrHH = Register[ID1 - 1].HouseholdID;
	RandHH = rg.Random();

	if (CurrHH == Register[ID2 - 1].HouseholdID & CurrHH > 0){
		// First identify which household member leaves the household.
		CurrHead = HHregister[CurrHH - 1].IDhead;
		if (CurrHead == ID1){ LeavingID = ID2; }
		else if (CurrHead == ID2){ LeavingID = ID1; }
		else if (Register[ID1 - 1].ParentID[0] == CurrHead || Register[ID1 - 1].ParentID[1] == CurrHead){
				LeavingID = ID2;
		}
		else{ LeavingID = ID1; }

		// Secondly determine whether the departing partner forms a new household.
		NewHH = 0;
		if (Register[LeavingID - 1].SexInd == 0){ NewHH = 1; }
		else{
			OddsNewHH = BaseOddsNHADS;
			MAlive = 0;
			FAlive = 0;
			if (Register[LeavingID - 1].ParentID[0] > 0){
				TempID = Register[LeavingID - 1].ParentID[0];
				if (Register[TempID - 1].AliveInd == 1){ FAlive = 1; }
			}
			if (Register[LeavingID - 1].ParentID[1] > 0){
				TempID = Register[LeavingID - 1].ParentID[1];
				if (Register[TempID - 1].AliveInd == 1){ MAlive = 1; }
			}
			if (MAlive == 1 || FAlive == 1){ OddsNewHH *= ParentAliveNHADS; }
			if (Register[LeavingID - 1].Employed == 1){ OddsNewHH *= EmployedEffectNHADS; }
			ir = Register[LeavingID - 1].Race;
			OddsNewHH *= RaceEffectNHADS[ir];
			ia = Register[LeavingID - 1].CurrAge;
			OddsNewHH *= pow(AgeEffectNHADS, ia) * pow(Age2EffectNHADS, ia * ia);
			ProbNewHH = OddsNewHH / (1.0 + OddsNewHH);
			if (RandHH < ProbNewHH){ NewHH = 1; }
		}

		// Thirdly move the departing partner to a different household.
		if (NewHH == 1){ NewHousehold(LeavingID); }
		else{ 
			NewHH = Register[LeavingID - 1].FindRelativeToLiveWith(LeavingID, 0); 
			// If there's no relative they can stay with, test if they become homeless or
			// form a new household.
			if (NewHH == 0){ NewHH = Register[LeavingID-1].TestHomeless(LeavingID); }
			if (NewHH == 0){ NewHousehold(LeavingID); }
		}

		// Fourthly remove the departing partner from the former household.
		HHregister[CurrHH - 1].RemoveMember(LeavingID);
	}
}

void Pop::GetStructRCTests()
{
	if (StructIntScenario == 0 || StructIntScenario == 1){ SingleSessionAlcoholCounselling.GetCurrModelEsts(); }
	if (StructIntScenario == 0 || StructIntScenario == 2){ MultiSessionAlcoholCounselling.GetCurrModelEsts(); }
	if (StructIntScenario == 0 || StructIntScenario == 3){ CashTransfers.GetCurrModelEsts(); }
	if (StructIntScenario == 0 || StructIntScenario == 4){ SchoolSupport.GetCurrModelEsts(); }
	if (StructIntScenario == 0 || StructIntScenario == 5){ VocationalTraining.GetCurrModelEsts(); }
	if (StructIntScenario == 0 || StructIntScenario == 6){ GenderTransformCommun.GetCurrModelEsts(); }
	if (StructIntScenario == 0 || StructIntScenario == 7){ GenderTransformIndiv.GetCurrModelEsts(); }
}

void Pop::OneSTDcycle()
{
	int ic;

	STDcycleCount += 1;
	GetSexActs();
	if (FixedUncertainty == 1 && STDcycleCount == 1 && BehavCycleCount == 1 &&
		(CurrYear == 1998 || CurrYear == 2003 || CurrYear == 2005 || CurrYear == 2016)){
		GetCondomUse();}
	if(HIVind==1){
		CalcHIVtesting();
		GetHIVtransitions();
	}
	if(HSVind==1 || TPind==1 || HDind==1 || NGind==1 || CTind==1 || 
		TVind==1 || BVind==1 || VCind==1){
			GetSTDtransitions();}

	// Update disease states
	for(ic=0; ic<Register.size(); ic++){
		if(Register[ic].AliveInd==1){
			if(HIVind==1){
				Register[ic].HIVstage = Register[ic].HIVstageE;}
			if(HSVind==1){
				Register[ic].HSVstage = Register[ic].HSVstageE;}
			if(TPind==1){
				Register[ic].TPstage = Register[ic].TPstageE;}
			if(HDind==1){
				Register[ic].HDstage = Register[ic].HDstageE;}
			if(NGind==1){
				Register[ic].NGstage = Register[ic].NGstageE;}
			if(CTind==1){
				Register[ic].CTstage = Register[ic].CTstageE;}
			if(TVind==1){
				Register[ic].TVstage = Register[ic].TVstageE;}
			if(Register[ic].SexInd==1){
				if(BVind==1 && Register[ic].AgeGroup>=2){
					Register[ic].BVstage = Register[ic].BVstageE;}
				if(VCind==1){
					Register[ic].VCstage = Register[ic].VCstageE;}
			}
		}
	}
}

void Pop::GetSexActs()
{
	int ic, ir, CSWID, startID;

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && (Register[ic].CasualInd == 1 || Register[ic].HetCasualInd == 1)){
			Register[ic].IDofCasual = 0;
			Register[ic].PAIcasual = 0;
			Register[ic].UAIcasual = 0;
		}
	}

	// Determine random starting point in Register
	r2[0] = rg.Random();
	startID = Register.size() * r2[0];

	for (ic = startID; ic<Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].VirginInd == 0 && Register[ic].SexInd == 0){
			Register[ic].SimulateSexActs(ic + 1);
		}
	}
	for (ic = 0; ic<startID; ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].VirginInd == 0 && Register[ic].SexInd == 0){
			Register[ic].SimulateSexActs(ic + 1);
		}
	}

	for(ic=0; ic<Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].VirginInd == 0 && Register[ic].SexInd == 1){
			if(Register[ic].CurrPartners==0){
				Register[ic].UVIprimary = 0;
				Register[ic].PVIprimary = 0;
			}
			if(Register[ic].CurrPartners<2){
				Register[ic].UVI2ndary = 0;
				Register[ic].PVI2ndary = 0;
			}
			if (Register[ic].HetCasualInd == 0){
				Register[ic].UAIcasual = 0;
				Register[ic].PAIcasual = 0;
			}
		}
	}
}

void Pop::CalcHIVtesting()
{
	int ic, ig, ia, is, iy, ih, iu, PID1, PID2;
	double VCTprob, OItestProb, HBCTprob, MobileProb, FSWprob, MSMprob, SchoolProb, STIprob, FPCprob, WorkProb;

	for (ic = 0; ic<Register.size(); ic++){
		r2[ic] = rg.Random();}
	for (ic = 0; ic<Register.size(); ic++){
		revent[ic] = rg.Random();}
	for (ic = 0; ic<Register.size(); ic++){
		rpAge[ic] = rg.Random();}
	iy = CurrYear - StartYear;

	for (ic = 0; ic<Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].AgeGroup>1 && Register[ic].Imprisoned == 0){
			ig = Register[ic].SexInd;
			ia = Register[ic].CurrAge-10;
			VCTprob = HCT1stTime[iy][ig] * VCTageEffect[ia][ig] / CycleD;
			if (Register[ic].VirginInd == 1){ VCTprob = 0.0; }
			ih = Register[ic].HighestGrade;
			VCTprob *= pow(VCTadjEdu, ih - 11);
			//if (CurrYear < 2000){ VCTprob *= pow(VCTadjEdu, ih - 11); }
			//else{ VCTprob *= pow(VCTadjEdu, 2.0); }
			iu = Register[ic].CurrUrban;
			HBCTprob = HBCTfreq[iy][iu] * HBCTuptake[ig] / CycleD;
			if (PropnST_HBCT[iy]>0.0){
				HBCTprob = (HBCTfreq[iy][iu] / CycleD) * (PropnST_HBCT[iy] / (1.0 + (1.0 - HBCTuptake[ig]) /
					(HBCTuptake[ig] * OR_STuptake)) + (1.0 - PropnST_HBCT[iy]) * HBCTuptake[ig]);
			}
			MobileProb = MobileTestCoverage[iy][iu] * (CommunityMobilization[iy] * MobileTestUptake[1] +
				(1.0 - CommunityMobilization[iy]) * MobileTestUptake[0]) / CycleD;
			is = Register[ic].HIVstage;
			OItestProb = OIincidence[is] * OIsDiagnosed[iy] / CycleD;
			if (Register[ic].Employed == 1){ WorkProb = TestingWorkplace[iy] * WorkTestUptake[ig] * EmployedReachable / CycleD; }
			else{ WorkProb = 0.0; }
			STIprob = 0.0;
			if (Register[ic].NGstage == 1 || Register[ic].CTstage == 1 || Register[ic].TVstage == 1 ||
				Register[ic].TPstage == 2 || Register[ic].HSVstage == 1 || Register[ic].HSVstage == 3 ||
				Register[ic].BVstage == 2 || Register[ic].VCstage == 2){
				if (Register[ic].FSWind == 1){ STIprob = 1.0 - exp(-FSWRxRate * PropnTreatedPublicF * 52.0 / CycleD); }
				else if (ig == 0){
					if (Register[ic].AgeGroup>=4){ STIprob = 1.0 - exp(-MaleRxRate * PropnTreatedPublicM * 52.0 / CycleD); }
					else{ STIprob = 1.0 - exp(-MaleTeenRxRate * PropnTreatedPublicM * 52.0 / CycleD); }
				}
				else {
					if (Register[ic].AgeGroup >= 4){ STIprob = 1.0 - exp(-FemRxRate * PropnTreatedPublicF * 52.0 / CycleD); }
					else{ STIprob = 1.0 - exp(-FemTeenRxRate * PropnTreatedPublicF * 52.0 / CycleD); }
				}
				STIprob *= TestingHIVinSTIs[iy] * 0.055;
				// 0.055 is to bring total STI cases in line with the numbers reported
			}
			FPCprob = 0.0;
			if (Register[ic].SexInd == 1 && (Register[ic].CurrContr == 1 || Register[ic].CurrContr == 2)){
				FPCprob = TestingFPC[iy] * FPCtestUptake / CycleD;}
			if (Register[ic].VCThistory == 1){ VCTprob *= RetestAdj; }
			if (Register[ic].VCThistory == 2){ 
				VCTprob *= RetestAdjDiagnosed[0];
				OItestProb *= RetestAdjDiagnosed[1];
				HBCTprob *= RetestAdjDiagnosed[3]; 
				MobileProb *= RetestAdjDiagnosed[4];
				STIprob *= RetestAdjDiagnosed[5];
				FPCprob *= RetestAdjDiagnosed[6];
				WorkProb *= RetestAdjDiagnosed[10];
			}
			if (Register[ic].HIVstage == 5){
				VCTprob *= RetestAdjART[0];
				OItestProb *= RetestAdjART[1];
				HBCTprob *= RetestAdjART[3];
				MobileProb *= RetestAdjART[4];
				STIprob *= RetestAdjART[5];
				FPCprob *= RetestAdjART[6];
				WorkProb *= RetestAdjART[10];
			}
			if (r2[ic] < (1.0 - exp(-VCTprob))){ // General HIV testing
				TotalTestsGen.out[CurrSim - 1][iy] += Register[ic].PopWeight;
				if (Register[ic].VCThistory == 2){ 
					Register[ic].GetRediagnosed(ic + 1, 0);
					PosTestsGen.out[CurrSim - 1][iy] += Register[ic].PopWeight;
				}
				else{ Register[ic].VCThistory = 1; }
				if (Register[ic].HIVstage>1 && Register[ic].VCThistory == 1){
					Register[ic].GetDiagnosed(ic + 1, 0);
					PosTestsGen.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					NewDiagGen.out[CurrSim - 1][iy] += Register[ic].PopWeight;
				}
				if (Register[ic].HIVstage == 1){ AcuteTests.out[CurrSim - 1][iy] += Register[ic].PopWeight; }
			}
			else if (r2[ic] < (1.0 - exp(-VCTprob - HBCTprob))){ // Tested through HBCT
				if (Register[ic].CurrUrban == 1){
					TotalTestsHH_U.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					if (Register[ic].VCThistory == 2){ 
						Register[ic].GetRediagnosed(ic + 1, 3);
						PosTestsHH_U.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					}
					else{ Register[ic].VCThistory = 1; }
					if (Register[ic].HIVstage>1 && Register[ic].VCThistory<2){
						Register[ic].GetDiagnosed(ic + 1, 3); 
						PosTestsHH_U.out[CurrSim - 1][iy] += Register[ic].PopWeight;
						NewDiagHH_U.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					}
					if (Register[ic].HIVstage == 1){ AcuteTests.out[CurrSim - 1][iy] += Register[ic].PopWeight; }
				}
				else{
					TotalTestsHH_R.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					if (Register[ic].VCThistory == 2){ 
						Register[ic].GetRediagnosed(ic + 1, 3);
						PosTestsHH_R.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					}
					else{ Register[ic].VCThistory = 1; }
					if (Register[ic].HIVstage>1 && Register[ic].VCThistory<2){
						Register[ic].GetDiagnosed(ic + 1, 3);
						PosTestsHH_R.out[CurrSim - 1][iy] += Register[ic].PopWeight;
						NewDiagHH_R.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					}
					if (Register[ic].HIVstage == 1){ AcuteTests.out[CurrSim - 1][iy] += Register[ic].PopWeight; }
				}
			}
			else if (r2[ic] < (1.0 - exp(-VCTprob - HBCTprob - MobileProb))){ // Mobile tested
				if (Register[ic].CurrUrban == 1){
					TotalTestsMobile_U.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					if (Register[ic].VCThistory == 2){ 
						Register[ic].GetRediagnosed(ic + 1, 3);
						PosTestsMobile_U.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					}
					else{ Register[ic].VCThistory = 1; }
					if (Register[ic].HIVstage>1 && Register[ic].VCThistory<2){
						Register[ic].GetDiagnosed(ic + 1, 3);
						PosTestsMobile_U.out[CurrSim - 1][iy] += Register[ic].PopWeight;
						NewDiagMobile_U.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					}
					if (Register[ic].HIVstage == 1){ AcuteTests.out[CurrSim - 1][iy] += Register[ic].PopWeight; }
				}
				else{
					TotalTestsMobile_R.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					if (Register[ic].VCThistory == 2){ 
						Register[ic].GetRediagnosed(ic + 1, 3);
						PosTestsMobile_R.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					}
					else{ Register[ic].VCThistory = 1; }
					if (Register[ic].HIVstage>1 && Register[ic].VCThistory<2){
						Register[ic].GetDiagnosed(ic + 1, 3);
						PosTestsMobile_R.out[CurrSim - 1][iy] += Register[ic].PopWeight;
						NewDiagMobile_R.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					}
					if (Register[ic].HIVstage == 1){ AcuteTests.out[CurrSim - 1][iy] += Register[ic].PopWeight; }
				}
			}
			else if (r2[ic] < (1.0 - exp(-VCTprob - HBCTprob - MobileProb - STIprob))){ // Tested in STI clinic
				TotalTestsSTI.out[CurrSim - 1][iy] += Register[ic].PopWeight;
				if (Register[ic].VCThistory == 2){ 
					Register[ic].GetRediagnosed(ic + 1, 0);
					PosTestsSTI.out[CurrSim - 1][iy] += Register[ic].PopWeight;
				}
				else{ Register[ic].VCThistory = 1; }
				if (Register[ic].HIVstage>1 && Register[ic].VCThistory == 1){
					Register[ic].GetDiagnosed(ic + 1, 0);
					PosTestsSTI.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					NewDiagSTI.out[CurrSim - 1][iy] += Register[ic].PopWeight;
				}
				if (Register[ic].HIVstage == 1){ AcuteTests.out[CurrSim - 1][iy] += Register[ic].PopWeight; }
			}
			else if (r2[ic] < (1.0 - exp(-VCTprob - HBCTprob - MobileProb - STIprob - FPCprob))){ // Tested in FP clinic
				TotalTestsFPC.out[CurrSim - 1][iy] += Register[ic].PopWeight;
				if (Register[ic].VCThistory == 2){ 
					Register[ic].GetRediagnosed(ic + 1, 0);
					PosTestsFPC.out[CurrSim - 1][iy] += Register[ic].PopWeight;
				}
				else{ Register[ic].VCThistory = 1; }
				if (Register[ic].HIVstage>1 && Register[ic].VCThistory == 1){
					Register[ic].GetDiagnosed(ic + 1, 0);
					PosTestsFPC.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					NewDiagFPC.out[CurrSim - 1][iy] += Register[ic].PopWeight;
				}
				if (Register[ic].HIVstage == 1){ AcuteTests.out[CurrSim - 1][iy] += Register[ic].PopWeight; }
			}
			else if (r2[ic] < (1.0 - exp(-VCTprob - HBCTprob - MobileProb - STIprob - FPCprob - WorkProb))){ // Tested at work
				TotalTestsWork.out[CurrSim - 1][iy] += Register[ic].PopWeight;
				if (Register[ic].VCThistory == 2){
					Register[ic].GetRediagnosed(ic + 1, 3);
					PosTestsWork.out[CurrSim - 1][iy] += Register[ic].PopWeight;
				}
				else{ Register[ic].VCThistory = 1; }
				if (Register[ic].HIVstage>1 && Register[ic].VCThistory == 1){
					Register[ic].GetDiagnosed(ic + 1, 3);
					PosTestsWork.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					NewDiagWork.out[CurrSim - 1][iy] += Register[ic].PopWeight;
				}
				if (Register[ic].HIVstage == 1){ AcuteTests.out[CurrSim - 1][iy] += Register[ic].PopWeight; }
			}
			else{ // Tested as an OI patient
				if (r2[ic] < (1.0 - exp(-VCTprob - HBCTprob - MobileProb - STIprob - FPCprob - WorkProb - OItestProb))){
					TotalTestsOI.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					if (Register[ic].VCThistory == 2){ 
						Register[ic].GetRediagnosed(ic + 1, 1);
						PosTestsOI.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					}
					else{ Register[ic].VCThistory = 1; }
					if (Register[ic].HIVstage>1 && Register[ic].VCThistory == 1){
						Register[ic].GetDiagnosed(ic + 1, 1);
						PosTestsOI.out[CurrSim - 1][iy] += Register[ic].PopWeight;
						NewDiagOI.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					}
					if (Register[ic].HIVstage == 1){ AcuteTests.out[CurrSim - 1][iy] += Register[ic].PopWeight; }
				}
			}
			if (Register[ic].FSWind == 1){ // Sex worker testing
				FSWprob = FSWtestUptake[iy] / CycleD;
				if (Register[ic].VCThistory == 2){ FSWprob *= RetestAdjDiagnosed[7]; }
				if (Register[ic].HIVstage == 5){ FSWprob *= RetestAdjART[7]; }
				FSWprob = 1.0 - exp(-FSWprob);
				if (revent[ic] < FSWprob){
					TotalTestsFSW.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					if (Register[ic].VCThistory == 2){ 
						Register[ic].GetRediagnosed(ic + 1, 3);
						PosTestsFSW.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					}
					else{ Register[ic].VCThistory = 1; }
					if (Register[ic].HIVstage>1 && Register[ic].VCThistory == 1){
						Register[ic].GetDiagnosed(ic + 1, 3);
						PosTestsFSW.out[CurrSim - 1][iy] += Register[ic].PopWeight;
						NewDiagFSW.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					}
					if (Register[ic].HIVstage == 1){ AcuteTests.out[CurrSim - 1][iy] += Register[ic].PopWeight; }
				}
			}
			// MSM testing
			if (Register[ic].SexInd == 0 && Register[ic].MalePref>0.0){
				MSMprob = MSMtestUptake[iy] / CycleD;
				if (Register[ic].VCThistory == 2){ MSMprob *= RetestAdjDiagnosed[8]; }
				if (Register[ic].HIVstage == 5){ MSMprob *= RetestAdjART[8]; }
				PID1 = Register[ic].IDprimary;
				PID2 = Register[ic].ID2ndary;
				// Note we can use the same random number as for FSW since the 2 groups are non-overlapping.
				if (revent[ic] < MSMprob && (Register[ic].CasualInd==1 || (Register[ic].CurrPartners>0 &&
					Register[PID1 - 1].SexInd == 0) || (Register[ic].CurrPartners==2 &&
					Register[PID2 - 1].SexInd == 0))){
					TotalTestsMSM.out[CurrSim - 1][iy] += Register[ic].PopWeight;
						if (Register[ic].VCThistory == 2){ 
							Register[ic].GetRediagnosed(ic + 1, 0);
							PosTestsMSM.out[CurrSim - 1][iy] += Register[ic].PopWeight;
						}
						else{ Register[ic].VCThistory = 1; }
						if (Register[ic].HIVstage>1 && Register[ic].VCThistory == 1){
							Register[ic].GetDiagnosed(ic + 1, 0);
							PosTestsMSM.out[CurrSim - 1][iy] += Register[ic].PopWeight;
							NewDiagMSM.out[CurrSim - 1][iy] += Register[ic].PopWeight;
						}
						if (Register[ic].HIVstage == 1){ AcuteTests.out[CurrSim - 1][iy] += Register[ic].PopWeight; }
				}
			}
			// Testing in schools
			if (Register[ic].InSchool == 1 && Register[ic].HighestGrade<12 && Register[ic].HighestGrade>6){ 
				SchoolProb = SchoolTestFreq[iy] * SchoolTestUptake[ig] / CycleD;
				if (Register[ic].VirginInd == 1){ SchoolProb *= SchoolTestVirginRR; }
				if (Register[ic].VCThistory == 2){ SchoolProb *= RetestAdjDiagnosed[9]; }
				if (Register[ic].HIVstage == 5){ SchoolProb *= RetestAdjART[9]; }
				if (rpAge[ic] < SchoolProb){
					TotalTestsSchool.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					if (Register[ic].VCThistory == 2){ 
						Register[ic].GetRediagnosed(ic + 1, 3);
						PosTestsSchool.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					}
					else{ Register[ic].VCThistory = 1; }
					if (Register[ic].HIVstage>1 && Register[ic].VCThistory == 1){
						Register[ic].GetDiagnosed(ic + 1, 3);
						PosTestsSchool.out[CurrSim - 1][iy] += Register[ic].PopWeight;
						NewDiagSchool.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					}
					if (Register[ic].HIVstage == 1){ AcuteTests.out[CurrSim - 1][iy] += Register[ic].PopWeight; }
				}
			}
		}
		// Testing in prisons
		if (Register[ic].AliveInd == 1 && Register[ic].Imprisoned > 0){
			VCTprob = TestingPrisons[iy];
			//VCTprob = TestingPrisons[iy] * 1.17/0.34;
			if (Register[ic].VCThistory == 2){VCTprob *= RetestAdjDiagnosed[12]; }
			if (Register[ic].HIVstage == 5){ VCTprob *= RetestAdjART[12]; }
			VCTprob = 1.0 - exp(-VCTprob / CycleD);
			if (revent[ic] < VCTprob){
				TotalTestsPrison.out[CurrSim - 1][iy] += Register[ic].PopWeight;
				if (Register[ic].VCThistory == 2){ 
					Register[ic].GetRediagnosed(ic + 1, 0);
					PosTestsPrison.out[CurrSim - 1][iy] += Register[ic].PopWeight;
				}
				else{ Register[ic].VCThistory = 1; }
				if (Register[ic].HIVstage > 1 && Register[ic].VCThistory == 1){
					Register[ic].GetDiagnosed(ic + 1, 0);
					PosTestsPrison.out[CurrSim - 1][iy] += Register[ic].PopWeight;
					NewDiagPrison.out[CurrSim - 1][iy] += Register[ic].PopWeight;
				}
				if (Register[ic].HIVstage == 1){ AcuteTests.out[CurrSim - 1][iy] += Register[ic].PopWeight; }
			}
		}
	}
}

void Pop::GetHIVtransitions()
{
	int ic;

	//int seed = 7819 + CurrYear * 35 + BehavCycleCount * 73 + STDcycleCount * 22;
	//if(CurrYear>=StartYear+FixedPeriod){
	//	seed += CurrSim;}
	//CRandomMersenne rg(seed);
	for(ic=0; ic<Register.size(); ic++){
		r2[ic] = rg.Random();}

	for(ic=0; ic<Register.size(); ic++){
		if(Register[ic].AliveInd==1 && Register[ic].HIVstage>0 &&
			Register[ic].AgeGroup<3 && (Register[ic].DateInfect -
			Register[ic].DOB) < 3.0){
				// Vertically-infected children 
				if(Register[ic].SexInd==0 && 
					Register[ic].DOB==Register[ic].DateInfect){
						MaleChild.Perinatal.GetNewStage(ic+1, r2[ic]);}
				else if(Register[ic].SexInd==0){
					MaleChild.Breastmilk.GetNewStage(ic+1, r2[ic]);}
				if(Register[ic].SexInd==1 && 
					Register[ic].DOB==Register[ic].DateInfect){
						FemChild.Perinatal.GetNewStage(ic+1, r2[ic]);}
				else if(Register[ic].SexInd==1){
					FemChild.Breastmilk.GetNewStage(ic+1, r2[ic]);}
		}
		else if(Register[ic].AliveInd==1 && (Register[ic].HIVstage>0 ||
			Register[ic].VirginInd==0)){
				Register[ic].GetNewHIVstate(ic+1, r2[ic]);}
	}
}

void Pop::GetSTDtransitions()
{
	int ic, id;

	//int seed = 3705 + CurrYear * 37 + BehavCycleCount * 81 + STDcycleCount * 53;
	//if(CurrYear>=StartYear+FixedPeriod){
	//	seed += CurrSim;}
	//CRandomMersenne rg(seed);
	for(ic=0; ic<Register.size(); ic++){
		for(id=0; id<8; id++){
			rSTI[ic][id] = rg.Random();}
	}

	for(ic=0; ic<Register.size(); ic++){
		if(Register[ic].AliveInd==1 && Register[ic].AgeGroup>=2){
			if(Register[ic].VirginInd==0){
				if(HSVind==1){Register[ic].GetNewHSVstate(ic+1, rSTI[ic][0]);}
				if(TPind==1){Register[ic].GetNewTPstate(ic+1, rSTI[ic][1]);}
				if(HDind==1){Register[ic].GetNewHDstate(ic+1, rSTI[ic][2]);}
				if(NGind==1){Register[ic].GetNewNGstate(ic+1, rSTI[ic][3]);}
				if(CTind==1){Register[ic].GetNewCTstate(ic+1, rSTI[ic][4]);}
				if(TVind==1){Register[ic].GetNewTVstate(ic+1, rSTI[ic][5]);}
			}
			if(Register[ic].SexInd==1){
				if(BVind==1){BVtransitionF.GetNewStage(ic+1, rSTI[ic][6]);}
				if(VCind==1){VCtransitionF.GetNewStage(ic+1, rSTI[ic][7]);}
			}
		}
	}
}

Partner::Partner(){}

PartnerCohort::PartnerCohort()
{	
	TotalDesire = 0.0;
}

void PartnerCohort::Reset()
{
	TotalDesire = 0.0;
	Pool.clear();
}

void PartnerCohort::AddMember(int PID, int Pref)
{
	Partner A;
	A.ID = PID;
	A.DesiredNewPartners = Register[PID-1].DesiredNewPartners;
	if (Pref == 0){ A.DesiredNewPartners *= Register[PID - 1].MalePref; }
	else{ A.DesiredNewPartners *= (1.0 - Register[PID - 1].MalePref); }

	Pool.push_back(A);
	TotalDesire += A.DesiredNewPartners;
}

int PartnerCohort::SamplePool(double rand1)
{
	int ReachedEnd; // Indicator of whether search algorithm has finished
	int ic, PartnerID;
	double CumProb;

	rand1 *= TotalDesire;
	ic = 0;
	CumProb = 0.0;
	ReachedEnd = 0;
	while(ReachedEnd==0){
		PartnerID = Pool[ic].ID;
		if(Register[PartnerID-1].NewStatus==0){
			CumProb += Pool[ic].DesiredNewPartners;
			if(rand1 < CumProb){
				ReachedEnd = 1;}
		}
		if(ic==Pool.size()-1 && ReachedEnd==0){
			// Prevent infinite loops that may arise due to rounding errors
			ReachedEnd = 1;}
		else{
			ic += 1;}
	}

	return PartnerID;
}

int PartnerCohort::SamplePool2(int ID1)
{
	int ii, PartnerID, FoundPartner, MaxSearches;
	int race1, edu1, edu2, urban1;
	double rpartner[100][2], MultAdj;

	race1 = Register[ID1 - 1].Race;
	if (race1 == 0){ MaxSearches = 10; }
	else{ MaxSearches = 100; }
	urban1 = Register[ID1 - 1].CurrUrban;

	for (ii = 0; ii < MaxSearches; ii++){
		rpartner[ii][0] = rg.Random();
		rpartner[ii][1] = rg.Random();
	}

	FoundPartner = 0;
	ii = 0;
	while (FoundPartner == 0){
		PartnerID = SamplePool(rpartner[ii][0]);
		if (PartnerID == 0){ cout << "Problem with SamplePool function" << endl; }
		if (Register[PartnerID - 1].NewStatus == 1){
			// This is an indicator that there aren't any available people in the risk
			// pool, so abort the search.
			FoundPartner = 1;
			MaxNewPartnerInd = 1;
		}
		else{
			MultAdj = 1.0;
			if (race1 != Register[PartnerID - 1].Race){ MultAdj *= RaceMixing[race1]; }
			edu1 = Register[ID1 - 1].HighestGrade;
			edu2 = Register[PartnerID - 1].HighestGrade;
			MultAdj *= EduMixing[edu1][edu2]; 
			if (Register[ID1 - 1].HIVstage == 5 && Register[PartnerID - 1].HIVstage < 5){
				MultAdj *= ARTdiscordance; }
			if (Register[PartnerID - 1].CurrUrban != urban1){ MultAdj *= LocationMixing; }
			// Prevent situation where MSM selects himself as a sexual partner
			if (PartnerID == ID1){ MultAdj = 0.0; }
			if (rpartner[ii][1] < MultAdj){ FoundPartner = 1; }
		}
		ii += 1;
		// If no luck after 10 trials, give up
		if (ii == MaxSearches & FoundPartner == 0){
			FoundPartner = 1;
			MaxNewPartnerInd = 1;
		}
	}

	return PartnerID;
}

void ReadSexAssumps(const char* input)
{
	int ia, ib, ig, is;
	ifstream file;

	file.open(input);
	if (file.fail()) {
		cerr << "Could not open input file.txt\n";
		exit(1);
	}
	file.ignore(255,'\n');
	file>>HighPropnM>>HighPropnF;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file >> ConscientiousEffectSex;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file>>AssortativeM>>AssortativeF;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>GenderEquality;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>AnnNumberClients;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>DebutMedian[0]>>DebutMedian[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> DebutShape[0] >> DebutShape[1];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>DebutAdjLow[0]>>DebutAdjLow[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> DebutAdjRace[0] >> DebutAdjRace[1] >> DebutAdjRace[2];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RRdebutLogDropIncomeF;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RRdebutInSchool[0] >> RRdebutInSchool[1];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>PartnershipFormation[0][0]>>PartnershipFormation[1][0]>>
		PartnershipFormation[0][1]>>PartnershipFormation[1][1];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>BasePartnerAcqH[0]>>BasePartnerAcqH[1];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>GammaMeanST[0]>>GammaMeanST[1];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>GammaStdDevST[0]>>GammaStdDevST[1];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>PartnerEffectNew[0][0]>>PartnerEffectNew[0][1]>>PartnerEffectNew[1][0]>>
		PartnerEffectNew[1][1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RaceEffectNew[0] >> RaceEffectNew[1] >> RaceEffectNew[2];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RaceMixing[0] >> RaceMixing[1] >> RaceMixing[2];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(is=0; is<6; is++){
		file>>HIVeffectPartners[is];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RR_STemployedM;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> MarriageConstant[0] >> MarriageConstant[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> MarriageTrend[0] >> MarriageTrend[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia < 3; ia++){
		file >> MarriageShape[ia][0] >> MarriageShape[ia][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> MarriageMin[0] >> MarriageMin[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> DivorceAdj;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> DivorceTrend;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ORremarriage[0] >> ORremarriage[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ORmarriage2ndaryEdu[0] >> ORmarriage2ndaryEdu[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ORmarriageTertiary[0] >> ORmarriageTertiary[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ORmarriageInSchool;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RaceMarriage[0] >> RaceMarriage[1] >> RaceMarriage[2];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> BaseOddsMarriage[0] >> BaseOddsMarriage[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<3; ia++){
		for (ib = 0; ib<2; ib++){
			file >> OR_MarriedAge[ia][ib];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<3; ia++){
		for (ib = 0; ib<2; ib++){
			file >> OR_MarriedRace[ia][ib];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> OR_MarriedEmployed[0] >> OR_MarriedEmployed[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<5; ia++){
		for (ib = 0; ib<2; ib++){
			file >> OR_MarriedEdu[ia][ib];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> OR_MarriedAgeEmployed[0] >> OR_MarriedAgeEmployed[1];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>MeanFSWcontacts>>GammaMeanFSW>>GammaStdDevFSW;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(ib=0; ib<5; ib++){
		file>>PartnerEffectFSWcontact[ib];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> SWurban[0] >> SWurban[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RR_FSWcontactEmployedM;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(is=0; is<6; is++){
		file>>HIVeffectFSWentry[is];}
	for(is=0; is<6; is++){
		file>>HIVeffectFSWexit[is];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file >> RR_FSWentryDiagnosed;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file>>MeanDurSTrel[0][0]>>MeanDurSTrel[0][1]>>MeanDurSTrel[1][0]>>MeanDurSTrel[1][1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RRbreakupBinge;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(ia=0; ia<16; ia++){
		file>>LTseparation[ia][0]>>LTseparation[ia][1];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(ia=0; ia<16; ia++){
		for(ib=0; ib<16; ib++){
			file>>AgePrefF[ia][ib];}
	}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(ia=0; ia<16; ia++){
		for(ib=0; ib<16; ib++){
			file>>AgePrefM[ia][ib];}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> AgePrefAdj2ndary;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(ia=0; ia<16; ia++){
		file>>FreqSexST[ia][1];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(ia=0; ia<16; ia++){
		file>>FreqSexLT[ia][1];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>BaselineCondomSvy;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(ib=0; ib<3; ib++){
		file>>RelEffectCondom[ib];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(ib=0; ib<3; ib++){
		file>>AgeEffectCondom[ib];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(ib=0; ib<3; ib++){
		file>>RatioInitialTo1998[ib];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(ib=0; ib<3; ib++){
		file>>RatioUltTo1998[ib];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(ib=0; ib<3; ib++){
		file>>MedianToBehavChange[ib];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ib = 0; ib<14; ib++){
		file >> CondomEdu[ib];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> CondomEduMedian;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ib = 0; ib<3; ib++){
		file >> CondomRace[ib];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>SDpartnerRateAdj;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> SDcondomPref;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> CondomScaling;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ORcondomBingePW;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ORcondomPerYrSchool;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> CasualEntryHet[0] >> CasualEntryHet[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> CasualAgeAdjHet[0] >> CasualAgeAdjHet[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> CasualLowAdjHet[0] >> CasualLowAdjHet[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> CasualHighAdjHet;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> CasualSexFreqHet;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RRcasualBinge[0] >> RRcasualBinge[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RRcasualEmployedM;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RRcasualLogDropIncomeF;
	file.close();

	// Calculate frequency of sex in men
	double sumxy;
	for(ia=0; ia<16; ia++){
		sumxy = 0;
		for(ib=0; ib<16; ib++){
			sumxy += AgePrefM[ia][ib] * FreqSexST[ib][1];}
		FreqSexST[ia][0] = sumxy;
	}
	for(ia=0; ia<16; ia++){
		sumxy = 0;
		for(ib=0; ib<16; ib++){
			sumxy += AgePrefM[ia][ib] * FreqSexLT[ib][1];}
		FreqSexLT[ia][0] = sumxy;
	}

	// Set BaselineCondomUse
	BaselineCondomUse = BaselineCondomSvy;

	// Calculate Weibull shape parameters for pace of behaviour change
	for(ib=0; ib<3; ib++){
		ShapeBehavChange[ib] = log(log(1.0 - log(RatioInitialTo1998[ib])/
			log(RatioUltTo1998[ib]))/log(2.0))/log(13.0/MedianToBehavChange[ib]);
	}

	// Calculate CondomEdu
	for (ib = 0; ib < 13; ib++){
		CondomEdu[ib] = pow(ORcondomPerYrSchool, ib - 10);}

	// Calculate CumAgePrefM and CumAgePrefF
	for(ia=0; ia<16; ia++){
		CumAgePrefM[ia][0] = AgePrefM[ia][0];
		CumAgePrefF[ia][0] = AgePrefF[ia][0];
		for(ib=1; ib<15; ib++){
			CumAgePrefM[ia][ib] = CumAgePrefM[ia][ib-1] + AgePrefM[ia][ib];
			CumAgePrefF[ia][ib] = CumAgePrefF[ia][ib-1] + AgePrefF[ia][ib];
		}
		CumAgePrefM[ia][15] = 1.0;
		CumAgePrefF[ia][15] = 1.0;
	}
}

void ReadSTDepi(const char* input)
{
	int ia, is, iz, ig;
	double alpha, beta, temp;
	ifstream file;

	file.open(input);
	if (file.fail()) {
		cerr << "Could not open input file.txt\n";
		exit(1);
	}
	file.ignore(255,'\n');
	file>>HSVtransitionM.AveDuration[0]>>HSVtransitionM.AveDuration[1]>>
		HSVtransitionM.AveDuration[2]>>HSVtransitionF.AveDuration[0]>>
		HSVtransitionF.AveDuration[1]>>HSVtransitionF.AveDuration[2];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(iz=0; iz<6; iz++){
		file>>TPtransitionM.AveDuration[iz];}
	for(iz=0; iz<6; iz++){
		file>>TPtransitionF.AveDuration[iz];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>HDtransitionM.AveDuration[0]>>HDtransitionM.AveDuration[1]>>HDtransitionM.AveDuration[2]>>
		HDtransitionF.AveDuration[0]>>HDtransitionF.AveDuration[1]>>HDtransitionF.AveDuration[2];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>NGtransitionM.AveDuration[0]>>NGtransitionM.AveDuration[1]>>NGtransitionM.AveDuration[2]>>
		NGtransitionF.AveDuration[0]>>NGtransitionF.AveDuration[1]>>NGtransitionF.AveDuration[2];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>CTtransitionM.AveDuration[0]>>CTtransitionM.AveDuration[1]>>CTtransitionM.AveDuration[2]>>
		CTtransitionF.AveDuration[0]>>CTtransitionF.AveDuration[1]>>CTtransitionF.AveDuration[2];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>TVtransitionM.AveDuration[0]>>TVtransitionM.AveDuration[1]>>TVtransitionM.AveDuration[2]>>
		TVtransitionF.AveDuration[0]>>TVtransitionF.AveDuration[1]>>TVtransitionF.AveDuration[2];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(is=0; is<5; is++){
		file>>HIVtransitionM.AveDuration[is];}
	for(is=0; is<5; is++){
		file>>HIVtransitionF.AveDuration[is];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>VCtransitionF.AveDuration[0]>>VCtransitionF.AveDuration[1];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(iz=1; iz<4; iz++){
		file>>BVtransitionF.CtsTransition[iz][0];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>BVtransitionF.CtsTransition[0][1]>>BVtransitionF.CtsTransition[2][1]>>
		BVtransitionF.CtsTransition[3][1];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>HSVtransitionM.RecurrenceRate>>HSVtransitionF.RecurrenceRate>>
		VCtransitionF.RecurrenceRate;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>HDtransitionM.PropnImmuneAfterRx>>NGtransitionM.PropnImmuneAfterRx>>
		CTtransitionM.PropnImmuneAfterRx>>TVtransitionM.PropnImmuneAfterRx>>
		HDtransitionF.PropnImmuneAfterRx>>NGtransitionF.PropnImmuneAfterRx>>
		CTtransitionF.PropnImmuneAfterRx>>TVtransitionF.PropnImmuneAfterRx;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>HDtransitionM.PropnImmuneAfterSR>>NGtransitionM.PropnImmuneAfterSR>>
		CTtransitionM.PropnImmuneAfterSR>>TVtransitionM.PropnImmuneAfterSR>>
		HDtransitionF.PropnImmuneAfterSR>>NGtransitionF.PropnImmuneAfterSR>>
		CTtransitionF.PropnImmuneAfterSR>>TVtransitionF.PropnImmuneAfterSR;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>HSVtransitionM.SymptomaticPropn>>HDtransitionM.SymptomaticPropn>>
		NGtransitionM.SymptomaticPropn>>CTtransitionM.SymptomaticPropn>>
		TVtransitionM.SymptomaticPropn>>HSVtransitionF.SymptomaticPropn>>
		HDtransitionF.SymptomaticPropn>>NGtransitionF.SymptomaticPropn>>
		CTtransitionF.SymptomaticPropn>>TVtransitionF.SymptomaticPropn>>
		BVtransitionF.SymptomaticPropn;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>HSVtransitionM.TransmProb>>TPtransitionM.TransmProb>>HDtransitionM.TransmProb>>
		NGtransitionM.TransmProb>>CTtransitionM.TransmProb>>TVtransitionM.TransmProb>>
		HSVtransitionF.TransmProb>>TPtransitionF.TransmProb>>HDtransitionF.TransmProb>>
		NGtransitionF.TransmProb>>CTtransitionF.TransmProb>>TVtransitionF.TransmProb>>
		HSVtransitionM.TransmProbSW>>TPtransitionM.TransmProbSW>>HDtransitionM.TransmProbSW>>
		NGtransitionM.TransmProbSW>>CTtransitionM.TransmProbSW>>TVtransitionM.TransmProbSW;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>HSVtransitionM.RelTransmLT>>TPtransitionM.RelTransmLT>>HDtransitionM.RelTransmLT>>
		NGtransitionM.RelTransmLT>>CTtransitionM.RelTransmLT>>TVtransitionM.RelTransmLT>>
		HIVtransitionM.RelTransmLT>>HSVtransitionF.RelTransmLT>>TPtransitionF.RelTransmLT>>
		HDtransitionF.RelTransmLT>>NGtransitionF.RelTransmLT>>CTtransitionF.RelTransmLT>>
		TVtransitionF.RelTransmLT>>HIVtransitionF.RelTransmLT;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(is=0; is<3; is++){
		file>>InitHIVtransm[is][0];}
	for(is=0; is<3; is++){
		file>>InitHIVtransm[is][1];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	// Note that in the next line we are reading the male parameters into the female
	// arrays and the female parameters into the male arrays. This is deliberate; it makes
	// things a lot simpler when calculating the InfectProb arrays from the TransProb arrays
	// (in the STDtransition class).
	file>>HSVtransitionF.CondomEff>>TPtransitionF.CondomEff>>HDtransitionF.CondomEff>>
		NGtransitionF.CondomEff>>CTtransitionF.CondomEff>>TVtransitionF.CondomEff>>
		HIVtransitionF.CondomEff>>HSVtransitionM.CondomEff>>TPtransitionM.CondomEff>>
		HDtransitionM.CondomEff>>NGtransitionM.CondomEff>>CTtransitionM.CondomEff>>
		TVtransitionM.CondomEff>>HIVtransitionM.CondomEff;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file >> HSVtransitionF.CircEff >> TPtransitionF.CircEff >> HDtransitionF.CircEff >>
		NGtransitionF.CircEff >> CTtransitionF.CircEff >> TVtransitionF.CircEff >>
		HIVtransitionF.CircEff >> HSVtransitionM.CircEff >> TPtransitionM.CircEff >>
		HDtransitionM.CircEff >> NGtransitionM.CircEff >> CTtransitionM.CircEff >>
		TVtransitionM.CircEff >> HIVtransitionM.CircEff;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for(iz=0; iz<4; iz++){
		file>>HSVtransitionM.HIVinfecIncrease[iz];}
	for(iz=0; iz<4; iz++){
		file>>HSVtransitionF.HIVinfecIncrease[iz];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(iz=0; iz<6; iz++){
		file>>TPtransitionM.HIVinfecIncrease[iz];}
	for(iz=0; iz<6; iz++){
		file>>TPtransitionF.HIVinfecIncrease[iz];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>HDtransitionM.HIVinfecIncrease[0]>>HDtransitionM.HIVinfecIncrease[1]>>
		HDtransitionF.HIVinfecIncrease[0]>>HDtransitionF.HIVinfecIncrease[1];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>NGtransitionM.HIVinfecIncrease[0]>>NGtransitionM.HIVinfecIncrease[1]>>
		NGtransitionF.HIVinfecIncrease[0]>>NGtransitionF.HIVinfecIncrease[1];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>CTtransitionM.HIVinfecIncrease[0]>>CTtransitionM.HIVinfecIncrease[1]>>
		CTtransitionF.HIVinfecIncrease[0]>>CTtransitionF.HIVinfecIncrease[1];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>TVtransitionM.HIVinfecIncrease[0]>>TVtransitionM.HIVinfecIncrease[1]>>
		TVtransitionF.HIVinfecIncrease[0]>>TVtransitionF.HIVinfecIncrease[1];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(is=0; is<6; is++){
		file>>HIVtransitionM.HIVinfecIncrease[is];}
	for(is=0; is<6; is++){
		file>>HIVtransitionF.HIVinfecIncrease[is];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>VCtransitionF.HIVinfecIncrease[0]>>VCtransitionF.HIVinfecIncrease[1];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(iz=0; iz<3; iz++){
		file>>BVtransitionF.HIVinfecIncrease[iz];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(iz=0; iz<4; iz++){
		file>>HSVtransitionM.HIVsuscepIncrease[iz];}
	for(iz=0; iz<4; iz++){
		file>>HSVtransitionF.HIVsuscepIncrease[iz];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(iz=0; iz<6; iz++){
		file>>TPtransitionM.HIVsuscepIncrease[iz];}
	for(iz=0; iz<6; iz++){
		file>>TPtransitionF.HIVsuscepIncrease[iz];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>HDtransitionM.HIVsuscepIncrease[0]>>HDtransitionM.HIVsuscepIncrease[1]>>
		HDtransitionF.HIVsuscepIncrease[0]>>HDtransitionF.HIVsuscepIncrease[1];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>NGtransitionM.HIVsuscepIncrease[0]>>NGtransitionM.HIVsuscepIncrease[1]>>
		NGtransitionF.HIVsuscepIncrease[0]>>NGtransitionF.HIVsuscepIncrease[1];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>CTtransitionM.HIVsuscepIncrease[0]>>CTtransitionM.HIVsuscepIncrease[1]>>
		CTtransitionF.HIVsuscepIncrease[0]>>CTtransitionF.HIVsuscepIncrease[1];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>TVtransitionM.HIVsuscepIncrease[0]>>TVtransitionM.HIVsuscepIncrease[1]>>
		TVtransitionF.HIVsuscepIncrease[0]>>TVtransitionF.HIVsuscepIncrease[1];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>VCtransitionF.HIVsuscepIncrease[0]>>VCtransitionF.HIVsuscepIncrease[1];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(iz=0; iz<3; iz++){
		file>>BVtransitionF.HIVsuscepIncrease[iz];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	// Note that in the next few lines we are reading the male parameters into the female
	// arrays and the female parameters into the male arrays. This is deliberate; it makes
	// things a lot simpler when calculating the InfectProb arrays from the TransProb arrays
	// (in the STDtransition class).
	for(ia=0; ia<16; ia++){
		file>>HSVtransitionF.SuscepIncrease[ia]>>TPtransitionF.SuscepIncrease[ia]>>
			HDtransitionF.SuscepIncrease[ia]>>NGtransitionF.SuscepIncrease[ia]>>
			CTtransitionF.SuscepIncrease[ia]>>TVtransitionF.SuscepIncrease[ia]>>
			HIVtransitionF.SuscepIncrease[ia];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(ia=0; ia<16; ia++){
		file>>HSVtransitionM.SuscepIncrease[ia]>>TPtransitionM.SuscepIncrease[ia]>>
			HDtransitionM.SuscepIncrease[ia]>>NGtransitionM.SuscepIncrease[ia]>>
			CTtransitionM.SuscepIncrease[ia]>>TVtransitionM.SuscepIncrease[ia]>>
			HIVtransitionM.SuscepIncrease[ia];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (is = 0; is<7; is++){
		file >> OIincidence[is];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(is=0; is<6; is++){
		file>>HSVsheddingIncrease[is];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(is=0; is<6; is++){
		file>>HSVrecurrenceIncrease[is];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(is=0; is<6; is++){
		file>>VCtransitionF.IncidenceIncrease[is];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>VCtransitionF.Incidence>>BVtransitionF.Incidence1;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>BVtransitionF.IncidenceMultTwoPartners>>BVtransitionF.IncidenceMultNoPartners;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>HSVsymptomInfecIncrease;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(iz=0; iz<3; iz++){
		file>>InfecIncreaseSyndrome[iz][0];}
	for(iz=0; iz<3; iz++){
		file>>InfecIncreaseSyndrome[iz][1];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(iz=0; iz<3; iz++){
		file>>SuscepIncreaseSyndrome[iz][0];}
	for(iz=0; iz<3; iz++){
		file>>SuscepIncreaseSyndrome[iz][1];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(is=0; is<6; is++){
		file>>RelHIVfertility[is];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>PropnInfectedAtBirth>>PropnInfectedAfterBirth;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>MaleChild.Perinatal.PreAIDSmedian>>FemChild.Perinatal.PreAIDSmedian>>
		MaleChild.Breastmilk.PreAIDSmedian>>FemChild.Breastmilk.PreAIDSmedian;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>MaleChild.Perinatal.PreAIDSshape>>FemChild.Perinatal.PreAIDSshape>>
		MaleChild.Breastmilk.PreAIDSshape>>FemChild.Breastmilk.PreAIDSshape;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>MaleChild.Perinatal.MeanAIDSsurvival>>FemChild.Perinatal.MeanAIDSsurvival>>
		MaleChild.Breastmilk.MeanAIDSsurvival>>FemChild.Breastmilk.MeanAIDSsurvival;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>MaleTeenRxRate>>MaleRxRate>>FemTeenRxRate>>FemRxRate>>FSWRxRate;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>PropnTreatedPublicM>>PropnTreatedPublicF>>PropnTreatedPrivateM>>
		PropnTreatedPrivateF;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ProbOIprecedesAIDSmort;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ARTinterruption >> ARTresumption;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file >> RRinterruptionM >> ARTinterruption20mult;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file>>HSVtransitionF.CorrectRxPreSM>>TPtransitionF.CorrectRxPreSM>>
		HDtransitionF.CorrectRxPreSM>>NGtransitionF.CorrectRxPreSM>>
		CTtransitionF.CorrectRxPreSM>>TVtransitionF.CorrectRxPreSM>>
		BVtransitionF.CorrectRxPreSM>>VCtransitionF.CorrectRxPreSM;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>HSVtransitionM.CorrectRxWithSM>>TPtransitionM.CorrectRxWithSM>>
		HDtransitionM.CorrectRxWithSM>>NGtransitionM.CorrectRxWithSM>>
		CTtransitionM.CorrectRxWithSM>>TVtransitionM.CorrectRxWithSM>>
		HSVtransitionF.CorrectRxWithSM>>TPtransitionF.CorrectRxWithSM>>
		HDtransitionF.CorrectRxWithSM>>NGtransitionF.CorrectRxWithSM>>
		CTtransitionF.CorrectRxWithSM>>TVtransitionF.CorrectRxWithSM>>
		BVtransitionF.CorrectRxWithSM;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>HSVtransitionF.DrugEff>>TPtransitionF.DrugEff>>HDtransitionF.DrugEff>>
		NGtransitionF.DrugEff>>CTtransitionF.DrugEff>>TVtransitionF.DrugEff>>
		BVtransitionF.DrugEff>>VCtransitionF.DrugEff;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>BVtransitionF.DrugPartialEff>>VCtransitionF.DrugPartialEff;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>HSVtransitionF.TradnalEff>>TPtransitionF.TradnalEff>>HDtransitionF.TradnalEff>>
		NGtransitionF.TradnalEff>>CTtransitionF.TradnalEff>>TVtransitionF.TradnalEff>>
		BVtransitionF.TradnalEff>>VCtransitionF.TradnalEff;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>TPtransitionF.ANCpropnScreened>>TPtransitionF.ANCpropnTreated;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file >> VCTageMean[0] >> VCTageMean[1] >> VCTageSD[0] >> VCTageSD[1] >> RetestAdj >> VCTadjEdu;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RR_ARTstartPer100CD4;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file>>AcceptScreening>>AcceptNVP>>RednNVP>>RednFF;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>SecondaryRxMult>>SecondaryCureMult>>TPtransitionM.PropnSuscepAfterRx>>
		TPtransitionF.PropnSuscepAfterRx;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>FSWasympRxRate>>FSWasympCure;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>InitHIVprevHigh;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> InitHIVprevRace[0] >> InitHIVprevRace[1] >> InitHIVprevRace[2];
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>RatioUltToInitHIVtransm;
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	file>>SDsuscepHIVadj;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HIVtransitionM.MeanInitCD4 >> HIVtransitionF.MeanInitCD4;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HIVtransitionM.SDinitCD4 >> HIVtransitionF.SDinitCD4;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HIVtransitionM.MeanInitVL >> HIVtransitionF.MeanInitVL;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HIVtransitionM.SDinitVL >> HIVtransitionF.SDinitVL;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HIVtransitionM.AgeEffectInitVL >> HIVtransitionF.AgeEffectInitVL;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HIVtransitionM.ChangeCD4 >> HIVtransitionF.ChangeCD4;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HIVtransitionM.VLeffectCD4 >> HIVtransitionF.VLeffectCD4;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HIVtransitionM.AgeEffectCD4 >> HIVtransitionF.AgeEffectCD4;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HIVtransitionM.ChangeVL >> HIVtransitionF.ChangeVL;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HIVtransitionM.VLeffectVL >> HIVtransitionF.VLeffectVL;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HIVtransitionM.CD4effectVL >> HIVtransitionF.CD4effectVL;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HIVtransitionM.AgeEffectVL >> HIVtransitionF.AgeEffectVL;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HIVtransitionM.MortZeroCD4 >> HIVtransitionF.MortZeroCD4;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HIVtransitionM.RednMortCD4 >> HIVtransitionF.RednMortCD4;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HIVtransitionM.ARTmortRedn >> HIVtransitionF.ARTmortRedn;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HIVtransitionM.ARTmortConstant >> HIVtransitionF.ARTmortConstant;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> VLeffectTransm[0] >> VLeffectTransm[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HIVsusceptInjectable;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HIVinfectHormonal;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> PregEffectHIVtransm[0] >> PregEffectHIVtransm[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HSVprevRace[0] >> HSVprevRace[1] >> HSVprevRace[2];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> PrEPefficacy;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> PrEPefficacyMSM;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> PrEPdiscontinue;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> VaccineWane;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HIVvaccProt[0] >> HIVvaccProt[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> VaccAdherence;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia < 16; ia++){
		file >> MMCuptake[ia];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RR_VMMClogIncome;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> PublicANCuse[0] >> PublicANCuse[1] >> PublicANCuse[2];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> StillbirthRate[0] >> StillbirthRate[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HBCTuptake[0] >> HBCTuptake[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> OR_STuptake;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> MobileTestUptake[0] >> MobileTestUptake[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> SchoolTestUptake[0] >> SchoolTestUptake[1] >> SchoolTestVirginRR;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ANCpartnerTested[0][0] >> ANCpartnerTested[1][0] >> ANCpartnerTested[0][1] >> ANCpartnerTested[1][1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> WorkTestUptake[0] >> WorkTestUptake[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> EmployedReachable;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> FPCtestUptake;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iz = 0; iz < 15; iz++){
		file >> RetestAdjDiagnosed[iz];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iz = 0; iz < 15; iz++){
		file >> RetestAdjART[iz];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RR_ARTstartCommunity;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RR_ARTstartST;
	file.close();

	if (VaryParameters == 1){
		SimulateParameters();}

	// Set the effectiveness of SM for VC
	VCtransitionF.CorrectRxWithSM = VCtransitionF.CorrectRxPreSM;

	// Specify parameter values that apply to both sexes.
	HSVtransitionM.CorrectRxPreSM = HSVtransitionF.CorrectRxPreSM;
	TPtransitionM.CorrectRxPreSM = TPtransitionF.CorrectRxPreSM;
	HDtransitionM.CorrectRxPreSM = HDtransitionF.CorrectRxPreSM;
	NGtransitionM.CorrectRxPreSM = NGtransitionF.CorrectRxPreSM;
	CTtransitionM.CorrectRxPreSM = CTtransitionF.CorrectRxPreSM;
	TVtransitionM.CorrectRxPreSM = TVtransitionF.CorrectRxPreSM;

	HSVtransitionM.DrugEff = HSVtransitionF.DrugEff;
	TPtransitionM.DrugEff = TPtransitionF.DrugEff;
	HDtransitionM.DrugEff = HDtransitionF.DrugEff;
	NGtransitionM.DrugEff = NGtransitionF.DrugEff;
	CTtransitionM.DrugEff = CTtransitionF.DrugEff;
	TVtransitionM.DrugEff = TVtransitionF.DrugEff;

	HSVtransitionM.TradnalEff = HSVtransitionF.TradnalEff;
	TPtransitionM.TradnalEff = TPtransitionF.TradnalEff;
	HDtransitionM.TradnalEff = HDtransitionF.TradnalEff;
	NGtransitionM.TradnalEff = NGtransitionF.TradnalEff;
	CTtransitionM.TradnalEff = CTtransitionF.TradnalEff;
	TVtransitionM.TradnalEff = TVtransitionF.TradnalEff;

	// Convert annualized recurrence and incidence rates into weekly rates
	//HSVtransitionM.RecurrenceRate = HSVtransitionM.RecurrenceRate/52.0;
	//HSVtransitionF.RecurrenceRate = HSVtransitionF.RecurrenceRate/52.0;
	VCtransitionF.RecurrenceRate = VCtransitionF.RecurrenceRate/52.0;
	VCtransitionF.Incidence = VCtransitionF.Incidence/52.0;

	// Calculate remaining elements of the CtsTransition matrix and AveDuration for BV
	BVtransitionF.CtsTransition[1][2] = BVtransitionF.Incidence1 * 
		BVtransitionF.SymptomaticPropn;
	BVtransitionF.CtsTransition[1][3] = BVtransitionF.Incidence1 * 
		(1.0 - BVtransitionF.SymptomaticPropn);
	BVtransitionF.AveDuration[0] = 1.0/BVtransitionF.CtsTransition[0][1];
	BVtransitionF.AveDuration[1] = 1.0/(BVtransitionF.CtsTransition[1][0] + 
		BVtransitionF.CtsTransition[1][2] + BVtransitionF.CtsTransition[1][3]);
	BVtransitionF.AveDuration[2] = 1.0/(BVtransitionF.CtsTransition[2][0] + 
		BVtransitionF.CtsTransition[2][1]);
	BVtransitionF.AveDuration[3] = 1.0/(BVtransitionF.CtsTransition[3][0] + 
		BVtransitionF.CtsTransition[3][1]);

	// Calculate average ART survival in kids
	MaleChild.Perinatal.AveYrsOnART = HIVtransitionM.AveDuration[4]/52.0;
	MaleChild.Breastmilk.AveYrsOnART = HIVtransitionM.AveDuration[4]/52.0;
	FemChild.Perinatal.AveYrsOnART = HIVtransitionF.AveDuration[4]/52.0;
	FemChild.Breastmilk.AveYrsOnART = HIVtransitionF.AveDuration[4]/52.0;

	// Set the HIV transmission probs in current year to their initial values
	HIVtransitionM.TransmProb[0] = InitHIVtransm[0][0];
	HIVtransitionF.TransmProb[0] = InitHIVtransm[0][1];
	for(is=0; is<3; is++){
		HIVtransitionM.TransmProb[is+1] = InitHIVtransm[1][0];
		HIVtransitionM.TransmProb[is+4] = InitHIVtransm[2][0];
		HIVtransitionF.TransmProb[is+1] = InitHIVtransm[1][1];
		HIVtransitionF.TransmProb[is+4] = InitHIVtransm[2][1];
	}
	if(CofactorType==0){
		HIVtransitionM.TransmProb[1] *= 2.0;
		HIVtransitionM.TransmProb[3] *= 0.5;
		HIVtransitionM.TransmProb[4] *= 2.0;
		HIVtransitionM.TransmProb[6] *= 0.5;
		HIVtransitionF.TransmProb[1] *= 2.0;
		HIVtransitionF.TransmProb[3] *= 0.5;
		HIVtransitionF.TransmProb[4] *= 2.0;
		HIVtransitionF.TransmProb[6] *= 0.5;
	}

	// Set the RelTransmCSW values
	HSVtransitionM.RelTransmCSW = HSVtransitionM.TransmProbSW/HSVtransitionM.TransmProb;
	TPtransitionM.RelTransmCSW = TPtransitionM.TransmProbSW/TPtransitionM.TransmProb;
	HDtransitionM.RelTransmCSW = HDtransitionM.TransmProbSW/HDtransitionM.TransmProb;
	NGtransitionM.RelTransmCSW = NGtransitionM.TransmProbSW/NGtransitionM.TransmProb;
	CTtransitionM.RelTransmCSW = CTtransitionM.TransmProbSW/CTtransitionM.TransmProb;
	TVtransitionM.RelTransmCSW = TVtransitionM.TransmProbSW/TVtransitionM.TransmProb;

	// Set the initial STD treatment parameters
	InitDrugEffNG = NGtransitionM.DrugEff;
	InitMaleRxRate = MaleRxRate;
	InitMaleTeenRxRate = MaleTeenRxRate;
	InitFemRxRate = FemRxRate;
	InitFemTeenRxRate = FemTeenRxRate;
	InitFSWasympRxRate = FSWasympRxRate;
	InitFSWasympCure = FSWasympCure;
	InitCorrectRxHSV = HSVtransitionM.CorrectRxWithSM;
	InitCorrectRxTVM = TVtransitionM.CorrectRxWithSM;
	InitANCpropnScreened = TPtransitionF.ANCpropnScreened;
	InitANCpropnTreated = TPtransitionF.ANCpropnTreated;
	//InitRecurrenceRateM = HSVtransitionM.RecurrenceRate;
	//InitRecurrenceRateF = HSVtransitionF.RecurrenceRate;

	// Calculate relative rates of HIV testing by age
	for (ig = 0; ig < 2; ig++){
		beta = pow(VCTageSD[ig], 2.0) / VCTageMean[ig];
		alpha = VCTageMean[ig] / beta;
		for (ia = 0; ia < 81; ia++){
			VCTageEffect[ia][ig] = pow((10.0 + ia) / 25.0, alpha - 1.0) * exp((15.0 - ia) / beta);
		}
	}

	// Calculate the relative rates of ART interruption by age and sex
	// The parameters here are explained in section 3.5 of the Thembisa 4.7 national report.
	for (ia = 0; ia < 25; ia++){
		AgeAdjInterrupt[ia][1] = 1.0 + (exp(0.4 * (ia - 10)) / (1.0 + exp(0.4 * (ia - 10)))) *
			(2.0 * ARTinterruption20mult * exp(-0.3 * (ia - 10)) / (1.0 + exp(-0.3 * (ia - 10))));
		if (ia >= 10){
			AgeAdjInterrupt[ia][0] = AgeAdjInterrupt[ia][1] * RRinterruptionM;
		}
		else{
			temp = AgeAdjInterrupt[ia][1];
			AgeAdjInterrupt[ia][1] = 1.5 + (ARTinterruption20mult - 1.5) *
				(temp - 1.0) / (ARTinterruption20mult - 1.0);
			AgeAdjInterrupt[ia][0] = 1.5 + (ARTinterruption20mult * RRinterruptionM -
				1.5) * (temp - 1.0) / (ARTinterruption20mult - 1.0);
		}
	}
	for (ia = 25; ia < 81; ia++){
		AgeAdjInterrupt[ia][0] = RRinterruptionM;
		AgeAdjInterrupt[ia][1] = 1.0;
	}
}

void ReadConception()
{
	int ir, ia, ic;
	ifstream file;

	file.open("Conception.txt");
	if (file.fail()) {
		cerr << "Could not open Conception.txt\n";
		exit(1);
	}
	file.ignore(255, '\n');
	file >> MeanGestation;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> SDgestation;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> SDfecundability;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia < 8; ia++){
		file >> RateFecund[ia];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia < 15; ia++){
		file >> RateInfecund[ia];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia < 3; ia++){
		for (ir = 0; ir < 3; ir++){
			file >> InitContr[ia][ir];}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia < 3; ia++){
		for (ir = 0; ir < 3; ir++){
			file >> InitInjectable[ia][ir];}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ir = 0; ir < 3; ir++){
		file >> ORcontrUnmarriedNeverPreg[0][0][ir] >> ORcontrUnmarriedNeverPreg[1][0][ir] >>
			ORcontrUnmarriedNeverPreg[0][1][ir] >> ORcontrUnmarriedNeverPreg[1][1][ir];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ORcontrEdu;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ORcontrAbstinent;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ContrVirgin;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ORinjectableEdu;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ORinjectableNeverPreg;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia < 7; ia++){ file >> SterilizationNotPill[ia]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia < 7; ia++){ file >> SterilizationRates[ia]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RaceEffectSteril[0] >> RaceEffectSteril[1] >> RaceEffectSteril[2];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> StartContrNewRel[0] >> StartContrNewRel[1] >> StartContrNewRel[2];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> StartContrPostnatal[0] >> StartContrPostnatal[1] >> StartContrPostnatal[2];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> StartContrOther;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> PrevBirthEffectContr[0] >> PrevBirthEffectContr[1] >> PrevBirthEffectContr[2];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ir = 0; ir < 2; ir++){
		for (ia = 0; ia < 7; ia++){ file >> AgeEffectContr[ia][ir]; }
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> EduEffectContr;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HIVeffectContr;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> CondomEffectContr;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> NoPrevUseContr1;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> NoPrevUseContr2;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> InjectablePref[0] >> InjectablePref[1] >> InjectablePref[2];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ORinjectableAge;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ORinjectableInit;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> NewContrMethodWeight;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> StopContr[0] >> StopContr[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> StopContrAbsentPartner;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> CondomEffPreg;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ContrEffPreg[0] >> ContrEffPreg[1] >> ContrEffPreg[2];
	file.close();
}

void ReadBreastfeeding()
{
	int ii;
	ifstream file;

	file.open("Breastfeeding.txt");
	if (file.fail()) {
		cerr << "Could not open Breastfeeding.txt\n";
		exit(1);
	}
	file.ignore(255, '\n');
	file >> EverFeed;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii < 2; ii++){
		file >> MedianFeed[ii];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii < 2; ii++){
		file >> ShapeFeed[ii];}
	file.close();
}

void ReadMigration()
{
	int ii, ia, ir;
	ifstream file;

	file.open("MigrationAssumps.txt");
	if (file.fail()) {
		cerr << "Could not open MigrationAssumps.txt\n";
		exit(1);
	}
	file.ignore(255, '\n');
	for (ia = 0; ia < 18; ia++){
		for (ii = 0; ii < 8; ii++){ file >> InitUrbanPropn[ia][ii][0]; }
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia < 18; ia++){
		for (ii = 0; ii < 8; ii++){ file >> InitUrbanPropn[ia][ii][1]; }
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia < 18; ia++){
		file >> UrbanToRural[ia][0] >> UrbanToRural[ia][1];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia < 18; ia++){
		file >> RuralToUrban[ia][0] >> RuralToUrban[ia][1];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> UrbanAdjRace[0] >> UrbanAdjRace[1] >> UrbanAdjRace[2];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RuralAdjRace[0] >> RuralAdjRace[1] >> RuralAdjRace[2];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii < 6; ii++){
		file >> UrbanAdjEdu[ii];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii < 6; ii++){
		file >> RuralAdjEdu[ii];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> VisitFreq[0] >> VisitFreq[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> LocationMixing;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> MaxVisitLength;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> MaxCoitalFreqWeight;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> AbsentPartnerConcurAdj;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ir = 0; ir < 3; ir++){
		file >> MarriedMigAdj[ir]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ir = 0; ir < 3; ir++){
		file >> SeparatedMigAdj[ir];}
	file.close();

	for (ia = 0; ia < 18; ia++){
		for (ii = 0; ii < 6; ii++){
			if (ia < 3){
				InitUrbanPropn[ia][ii][1] *= 0.82;
				InitUrbanPropn[ia][ii][1] *= 0.82;
			}
			else{
				InitUrbanPropn[ia][ii][1] *= 0.97;
				InitUrbanPropn[ia][ii][1] *= 0.82;
			}
		}
	}
}

void ReadDisclosure()
{
	ifstream file;

	file.open("Disclosure.txt");
	if (file.fail()) {
		cerr << "Could not open Disclosure.txt\n";
		exit(1);
	}
	file.ignore(255, '\n');
	file >> DiscloseProb;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> DiscloseRate;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> DiscloseMale;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> DiscloseART;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> DiscloseMarried;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> DiscloseEffectRelDur;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> DiscloseEffectCondom;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ARTdiscordance;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ProbReferral[0] >> ProbReferral[1];
	file.close();
}

void ReadGenderNormAssumps()
{
	ifstream file;

	file.open("GenderNormAssumps.txt");
	if (file.fail()) {
		cerr << "Could not open GenderNormAssumps.txt\n";
		exit(1);
	}
	file.ignore(255, '\n');
	file >> AgeEffectIneqGender;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> EduEffectIneqGender[0] >> EduEffectIneqGender[1] >> EduEffectIneqGender[2];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RuralEffectIneqGender;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RaceEffectIneqGender[0] >> RaceEffectIneqGender[1] >> RaceEffectIneqGender[2];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HighRiskEffectIneqGender;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> EffectIneqGenderConcurrency;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> EffectIneqGenderCasual;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ORcondomGenderIneq;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RRgenderIneqIndiv;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RRgenderIneqComm;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ProbGenderIneqReversion;
	file.close();
}

void ReadAlcoholAssumps()
{
	int ii;
	ifstream file;

	file.open("AlcoholAssumps.txt");
	if (file.fail()) {
		cerr << "Could not open AlcoholAssumps.txt\n";
		exit(1);
	}
	file.ignore(255, '\n');
	file >> BaseDrinkProb[0] >> BaseDrinkProb[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> UrbanEffectDrinkProb[0] >> UrbanEffectDrinkProb[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii < 5; ii++){
		file >> AgeEffectDrinkProb[ii][0] >> AgeEffectDrinkProb[ii][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> MarriedEffectDrinkProb[0] >> MarriedEffectDrinkProb[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii < 4; ii++){
		file >> EduEffectDrinkProb[ii][0] >> EduEffectDrinkProb[ii][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii < 3; ii++){
		file >> RaceEffectDrinkProb[ii][0] >> RaceEffectDrinkProb[ii][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> StdDevDrinkProb[0] >> StdDevDrinkProb[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> BaseDrinksPerDD[0] >> BaseDrinksPerDD[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> UrbanEffectDrinksPerDD[0] >> UrbanEffectDrinksPerDD[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii < 5; ii++){
		file >> AgeEffectDrinksPerDD[ii][0] >> AgeEffectDrinksPerDD[ii][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> EmployedEffectDrinksPerDD[0] >> EmployedEffectDrinksPerDD[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii < 4; ii++){
		file >> EduEffectDrinksPerDD[ii][0] >> EduEffectDrinksPerDD[ii][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii < 3; ii++){
		file >> RaceEffectDrinksPerDD[ii][0] >> RaceEffectDrinksPerDD[ii][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ConscientiousEffectDrinksPerDD;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> GenderIneqEffectDrinksPerDD;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> StdDevDrinksPerDD[0] >> StdDevDrinksPerDD[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ConfoundingAlcSex;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> AlcPropnReported;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RatioMaxDrinksToAve;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ProbMaxDrinksPerDD;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RRalcoholSingle;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RRalcoholMultiple;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ProbAlcoholReversion;
	file.close();
}

void ReadHouseholdAssumps()
{
	int ia, ir, ih;
	double OddsLeaveNest[36]; // Odds of having left nest (ages 15-50)
	ifstream file;

	file.open("HouseholdAssumps.txt");
	if (file.fail()) {
		cerr << "Could not open HouseholdAssumps.txt\n";
		exit(1);
	}
	file.ignore(255, '\n');
	file >> BaseOddsHead[0] >> BaseOddsHead[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> MarriedEffectHead[0] >> MarriedEffectHead[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia < 11; ia++){
		file >> AgeEffectHead[ia][0] >> AgeEffectHead[ia][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ir = 0; ir < 3; ir++){
		file >> RaceEffectHead[ir][0] >> RaceEffectHead[ir][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> EmployedEffectHead[0] >> EmployedEffectHead[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> InSchoolEffectHead[0] >> InSchoolEffectHead[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ih = 0; ih < 5; ih++){
		file >> EduEffectHead[ih][0] >> EduEffectHead[ih][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> BirthEffectHead;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> BaseOddsALWP[0] >> BaseOddsALWP[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> MarriedEffectALWP[0] >> MarriedEffectALWP[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia < 8; ia++){
		file >> AgeEffectALWP[ia][0] >> AgeEffectALWP[ia][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ir = 0; ir < 3; ir++){
		file >> RaceEffectALWP[ir][0] >> RaceEffectALWP[ir][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> EmployedEffectALWP[0] >> EmployedEffectALWP[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> InSchoolEffectALWP[0] >> InSchoolEffectALWP[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ih = 0; ih < 5; ih++){
		file >> EduEffectALWP[ih][0] >> EduEffectALWP[ih][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> BaseOddsCLWP;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> AgeEffectCLWP[0] >> AgeEffectCLWP[1] >> AgeEffectCLWP[2];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RaceEffectCLWP[0] >> RaceEffectCLWP[1] >> RaceEffectCLWP[2];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> InSchoolEffectCLWP;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> BaseOddsCLSP[0] >> BaseOddsCLSP[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> AgeEffectCLSP[0] >> AgeEffectCLSP[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> Age2EffectCLSP[0] >> Age2EffectCLSP[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ir = 0; ir < 3; ir++){
		file >> RaceEffectCLSP[ir][0] >> RaceEffectCLSP[ir][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> SchoolingEffectCLSP[0] >> SchoolingEffectCLSP[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> BaseOddsNewHH[0] >> BaseOddsNewHH[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> MarriedEffectNewHH[0] >> MarriedEffectNewHH[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia < 7; ia++){
		file >> AgeEffectNewHH[ia][0] >> AgeEffectNewHH[ia][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ir = 0; ir < 3; ir++){
		file >> RaceEffectNewHH[ir][0] >> RaceEffectNewHH[ir][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> EmployedEffectNewHH[0] >> EmployedEffectNewHH[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> InSchoolEffectNewHH[0] >> InSchoolEffectNewHH[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ih = 0; ih < 5; ih++){
		file >> EduEffectNewHH[ih][0] >> EduEffectNewHH[ih][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> UrbanEffectNewHH[0] >> UrbanEffectNewHH[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> BaseOddsNHADS;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> AgeEffectNHADS;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> Age2EffectNHADS;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RaceEffectNHADS[0] >> RaceEffectNHADS[1] >> RaceEffectNHADS[2];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> EmployedEffectNHADS;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ParentAliveNHADS;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> BaseOddsLeaveNest;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ORleaveNestAge;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ORleaveNestAge2;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RRleaveNestFemale;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RRleaveNestRace[0] >> RRleaveNestRace[1] >> RRleaveNestRace[2];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ORleaveNestEmployed;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RRleaveNestInSchool;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RRleaveNestHH10plus;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> BaseOddsHHALN;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> OR_HHALNfemale;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> OR_HHALNage;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> OR_HHALNage2;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> OR_HHALNrace[0] >> OR_HHALNrace[1] >> OR_HHALNrace[2];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> OR_HHALNemployed;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> OR_HHALNinSchool;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ProbHomeless;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> BaseOddsNotHomeless;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ORnotHomelessFem;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ORnotHomelessEmployed;
	file.close();

	for (ia = 0; ia < 36; ia++){
		OddsLeaveNest[ia] = BaseOddsLeaveNest * pow(ORleaveNestAge, ia + 15) *
			pow(ORleaveNestAge2, (ia + 15) * (ia + 15)) * ORleaveNestEmployed;
		// Convert odds to probability
		OddsLeaveNest[ia] = OddsLeaveNest[ia] / (1.0 + OddsLeaveNest[ia]);
	}
	for (ia = 0; ia < 35; ia++){
		AnnProbLeaveNest[ia] = (OddsLeaveNest[ia + 1] - OddsLeaveNest[ia]) / 
			(1.0 - OddsLeaveNest[ia]);
	}
}

void ReadPrisonAssumps()
{
	int ia, ir;
	ifstream file;

	file.open("PrisonAssumps.txt");
	if (file.fail()) {
		cerr << "Could not open PrisonAssumps.txt\n";
		exit(1);
	}
	file.ignore(255, '\n');
	for (ia = 0; ia < 16; ia++){
		file >> PrisonEntryRate[ia]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ir = 0; ir < 3; ir++){
		file >> PrisonEntryRace[ir];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> PrisonEntryEdu;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> PrisonReentryAdj;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> MeanUnsentencedDur;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> UnsentencedRelease;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> SentenceLength[0] >> SentenceLength[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> SentenceLength85[0] >> SentenceLength85[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ProbParole;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> PrisonEntryRisk[0] >> PrisonEntryRisk[1];
	file.close();
}

void ReadEmployment()
{
	int ia, ir, ih, is;
	ifstream file;

	file.open("Employment.txt");
	if (file.fail()) {
		cerr << "Could not open Employment.txt\n";
		exit(1);
	}
	file.ignore(255, '\n');
	file >> BaseEmployed;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ih = 0; ih < 14; ih++){
		file >> EduEffectEmployed[ih];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ir = 0; ir < 3; ir++){
		file >> RaceEffectEmployed[ir];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> UrbanEffectEmployed[0] >> UrbanEffectEmployed[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> SexEffectEmployed[0] >> SexEffectEmployed[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia < 11; ia++){
		file >> AgeEffectEmployed[ia];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ih = 0; ih < 5; ih++){
		file >> EduEffectEmployed2[ih];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ir = 0; ir < 3; ir++){
		file >> RaceEffectEmployed2[ir];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ir = 0; ir < 3; ir++){
		file >> RaceEffectEmployed3[ir];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> UrbanEffectEmployed2[0] >> UrbanEffectEmployed2[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> SexEffectEmployed2[0] >> SexEffectEmployed2[1] >> SexEffectEmployed2[2];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia < 10; ia++){
		file >> AgeEffectEmployed2[ia];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (is = 0; is < 6; is++){
		file >> HIVeffectEmployed[is];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ChildEffectEmployed;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> PrevEmployedEffect;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ORunemployedTraining;
	file.close();
}

void ReadIncome()
{
	int ia, ir, ih, iy;
	ifstream file;

	file.open("IncomeAssumps.txt");
	if (file.fail()) {
		cerr << "Could not open IncomeAssumps.txt\n";
		exit(1);
	}
	file.ignore(255, '\n');
	for (iy = 0; iy < 5; iy++){
		file >> BaseLogIncome[iy];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy < 5; iy++){
		file >> FemaleEffectLogIncome[iy];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ih = 0; ih < 5; ih++){
		for (iy = 0; iy < 5; iy++){
			file >> EduEffectLogIncome[ih][iy];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ir = 0; ir < 3; ir++){
		for (iy = 0; iy < 5; iy++){
			file >> RaceEffectLogIncome[ir][iy];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia < 9; ia++){
		for (iy = 0; iy < 5; iy++){
			file >> AgeEffectLogIncome[ia][iy];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ConscientiousEffectLogIncome;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> StdDevLogIncome;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> AnnRealEarningsGrowth;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> BaseOddsPension;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ih = 0; ih < 5; ih++){
		file >> EffectEduPension[ih];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ir = 0; ir < 3; ir++){
		file >> EffectRacePension[ir];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> PrivatePension2019;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> AnnRealPensionGrowth;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ChildWeightEquivScale;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> HHsizeAdjEquivScale;
	file.close();
}

void ReadMSM()
{
	int ir, ic, is, ia, ib;
	ifstream file;

	file.open("MSMassumps.txt");
	if (file.fail()) {
		cerr << "Could not open MSMassumps.txt\n";
		exit(1);
	}
	file.ignore(255, '\n');
	file >> MSMfraction;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> BiFraction;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> InitMalePrefBeta[0] >> InitMalePrefBeta[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> AnnChangeMalePref[0] >> AnnChangeMalePref[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> SameSexMarried;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ir = 0; ir < 16; ir++){
		for (ic = 0; ic < 16; ic++){
			file >> AgePrefMSM[ir][ic];}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> CasualEntry;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> CasualAgeAdj;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> CasualLowAdj;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> CasualHighAdj;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> CasualExit;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> CasualSexFreq;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> FreqSexST_MSM;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ORcondomCasual;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (is = 0; is<3; is++){
		file >> InitHIVtransmMSM[is][0];
	}
	for (is = 0; is<3; is++){
		file >> InitHIVtransmMSM[is][1];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (is = 0; is<3; is++){
		file >> RolePrefBi[is];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (is = 0; is<3; is++){
		file >> RolePrefHomo[is];
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> MeanDurST_MSM;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RatioInitPrevMSM;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> IntraEventVersatility;
	file.close();


	// Now allow for uncertainty in the MSM parameters
	if (MSMcalib == 1){ SimulateMSMparameters(); }
	if (MSMcalibHIV == 1){ SimulateMSM_HIV(); }

	// Set the HIV transmission probs in current year to their initial values
	HIVtransitionM.TransmProbMSM[0][0] = InitHIVtransmMSM[0][0];
	HIVtransitionM.TransmProbMSM[0][1] = InitHIVtransmMSM[0][1];
	for (is = 0; is<3; is++){
		HIVtransitionM.TransmProbMSM[is + 1][0] = InitHIVtransmMSM[1][0];
		HIVtransitionM.TransmProbMSM[is + 4][0] = InitHIVtransmMSM[2][0];
		HIVtransitionM.TransmProbMSM[is + 1][1] = InitHIVtransmMSM[1][1];
		HIVtransitionM.TransmProbMSM[is + 4][1] = InitHIVtransmMSM[2][1];
	}

	// Calculate CumAgePrefMSM
	for (ia = 0; ia<16; ia++){
		CumAgePrefMSM[ia][0] = AgePrefMSM[ia][0];
		for (ib = 1; ib<15; ib++){
			CumAgePrefMSM[ia][ib] = CumAgePrefMSM[ia][ib - 1] + AgePrefMSM[ia][ib];
		}
		CumAgePrefMSM[ia][15] = 1.0;
	}
}

void ReadMSMcalibData()
{
	int ir, ic;
	ifstream file;

	file.open("MSMbehavData.txt");
	if (file.fail()) {
		cerr << "Could not open MSMbehavData.txt\n";
		exit(1);
	}
	file.ignore(255, '\n');
	for (ir = 0; ir<MSMmarried.studies; ir++){
		for (ic = 0; ic < 4; ic++){ file >> MSMmarried.CalibData[ir][ic]; }
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ir = 0; ir < MSMrecentBi.studies; ir++){
		for (ic = 0; ic < 4; ic++){ file >> MSMrecentBi.CalibData[ir][ic]; }
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ir = 0; ir<MSMeverBi.studies; ir++){
		for (ic = 0; ic < 4; ic++){ file >> MSMeverBi.CalibData[ir][ic]; }
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ir = 0; ir<MSMaged25plus.studies; ir++){
		for (ic = 0; ic < 4; ic++){ file >> MSMaged25plus.CalibData[ir][ic]; }
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ir = 0; ir<MSMcurrentReg.studies; ir++){
		for (ic = 0; ic < 4; ic++){ file >> MSMcurrentReg.CalibData[ir][ic]; }
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ir = 0; ir<MSMmultPartners.studies; ir++){
		for (ic = 0; ic < 4; ic++){ file >> MSMmultPartners.CalibData[ir][ic]; }
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ir = 0; ir<EverMSMcurrentReg.studies; ir++){
		for (ic = 0; ic < 4; ic++){ file >> EverMSMcurrentReg.CalibData[ir][ic]; }
	}
	file.close();
}

void ReadMSM_HIVdata()
{
	int ir, ic;
	ifstream file;

	file.open("MSM_HIVdata.txt");
	if (file.fail()) {
		cerr << "Could not open MSM_HIVdata.txt\n";
		exit(1);
	}
	file.ignore(255, '\n');
	for (ir = 0; ir<MSMpositive.studies; ir++){
		for (ic = 0; ic < 4; ic++){ file >> MSMpositive.CalibData[ir][ic]; }
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ir = 0; ir<EverMSMpositive.studies; ir++){
		for (ic = 0; ic < 4; ic++){ file >> EverMSMpositive.CalibData[ir][ic]; }
	}
	file.close();
}

void ReadRatesByYear()
{
	int iy, ig, is, ia;
	ifstream file;

	file.open("RatesByYear.txt");
	if (file.fail()) {
		cerr << "Could not open RatesByYear.txt\n";
		exit(1);
	}
	file.ignore(255,'\n');
	for(iy=0; iy<81; iy++){
		file>>PropnPrivateUsingSM[iy];}
	for(iy=0; iy<81; iy++){
		file>>PropnPublicUsingSM[iy];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(iy=0; iy<81; iy++){
		file>>DrugShortage[iy];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for (iy = 0; iy<81; iy++){
		file >> HCT1stTime[iy][0];}
	for (iy = 0; iy<81; iy++){
		file >> HCT1stTime[iy][1];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> OIsDiagnosed[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for(iy=0; iy<81; iy++){
		file>>HAARTaccess[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (is = 0; is < 2; is++){
		for (iy = 0; iy<81; iy++){
			file >> ARTeligiblePreg[iy][is];}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (is = 0; is < 4; is++){
		for (iy = 0; iy<81; iy++){
			file >> ARTeligibleOI[iy][is];}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (is = 0; is < 4; is++){
		for (iy = 0; iy<81; iy++){
			file >> ARTeligibleGen[iy][is];}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> ARTuptakePreg[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> ARTuptakeOI[iy];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for (ig = 0; ig < 2; ig++){
		for (iy = 0; iy<81; iy++){
			file >> ARTinitiationU200[iy][ig];}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for(iy=0; iy<81; iy++){
		file>>PMTCTaccess[iy];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(iy=0; iy<81; iy++){
		file>>PropnCiproResistant[iy];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(iy=0; iy<81; iy++){
		file>>PropnCiproTreated[iy];}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(iy=0; iy<81; iy++){
		file>>RxPhaseIn[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> MMCtoM_ratio[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> AnnPrEPuptake[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> PrEPuptakeMSM[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> AnnVaccineUptake[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> HIVdiagNoBF[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia < 3; ia++){
		for (iy = 0; iy < 81; iy++){
			file >> ORhormonalByYear[iy][ia];}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> RRsterilizationByYear[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> PropnSTreferral[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> PropnAssistedNotif[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> HBCTfreq[iy][0];}
	for (iy = 0; iy<81; iy++){
		file >> HBCTfreq[iy][1];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> PropnST_HBCT[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> MobileTestCoverage[iy][0];}
	for (iy = 0; iy<81; iy++){
		file >> MobileTestCoverage[iy][1];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> CommunityMobilization[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> FSWtestUptake[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> MSMtestUptake[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> SchoolTestFreq[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> ANCpartnersInvited[iy][0];}
	for (iy = 0; iy<81; iy++){
		file >> ANCpartnersInvited[iy][1];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> TestingPrisons[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> TestingHIVinSTIs[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> TestingWorkplace[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> TestingFPC[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> ANCretestFreq[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> BaseEmployed2[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> ConsPriceIndex[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> CSGamount[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> CSGageLimit[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> CSGincomeLimit[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> OAPamount[iy];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> DurFSW[iy];}
	for (iy = 0; iy<81; iy++){
		file >> FSWageMean[iy];}
	for (iy = 0; iy<81; iy++){
		file >> FSWageSD[iy];}
	file.close();

	/*for (iy = 0; iy < 81; iy++){
		OIsDiagnosed[iy] = 0.0;
		TestingHIVinSTIs[iy] = 0.0;
		TestingPrisons[iy] = 0.0;
	}*/
	/*for (iy = 34; iy < 81; iy++){
		ANCpartnersInvited[iy][0] = 1.0;
		ANCpartnersInvited[iy][1] = 1.0;
		TestingFPC[iy] = 1.0;
		TestingWorkplace[iy] = 1.0;
		FSWtestUptake[iy] = 0.15;
		MSMtestUptake[iy] = 0.50;
	}*/

	//for (iy = 0; iy<41; iy++){HAARTaccess[iy] = 0.0;}
}

void ReadMortTables()
{
	int ia, iy;
	ifstream file;

	file.open("MortTables.txt");
	if (file.fail()) {
		cerr << "Could not open MortTables.txt\n";
		exit(1);
	}
	// Black mortality rates
	file.ignore(255,'\n');
	for (iy = 0; iy<81; iy++){
		file >> InfantMort1st6mM[iy][0];}
	for(ia=0; ia<91; ia++){
		for(iy=0; iy<81; iy++){
			file>>NonAIDSmortM[ia][iy][0];}
	}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for (iy = 0; iy<81; iy++){
		file >> InfantMort1st6mF[iy][0];}
	for(ia=0; ia<91; ia++){
		for(iy=0; iy<81; iy++){
			file>>NonAIDSmortF[ia][iy][0];}
	}
	// Coloured mortality rates
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for (iy = 0; iy<81; iy++){
		file >> InfantMort1st6mM[iy][1];}
	for (ia = 0; ia<91; ia++){
		for (iy = 0; iy<81; iy++){
			file >> NonAIDSmortM[ia][iy][1];}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> InfantMort1st6mF[iy][1];}
	for (ia = 0; ia<91; ia++){
		for (iy = 0; iy<81; iy++){
			file >> NonAIDSmortF[ia][iy][1];}
	}
	// White mortality rates
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> InfantMort1st6mM[iy][2];}
	for (ia = 0; ia<91; ia++){
		for (iy = 0; iy<81; iy++){
			file >> NonAIDSmortM[ia][iy][2];}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<81; iy++){
		file >> InfantMort1st6mF[iy][2];}
	for (ia = 0; ia<91; ia++){
		for (iy = 0; iy<81; iy++){
			file >> NonAIDSmortF[ia][iy][2];}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (iy = 0; iy<91; iy++){
		file >> LE_West26[iy][0] >> LE_West26[iy][1];
	}
	file.close();
}

void ReadFertTables()
{
	int ia, iy;
	ifstream file;

	file.open("FertTables.txt");
	if (file.fail()) {
		cerr << "Could not open FertTables.txt\n";
		exit(1);
	}
	file.ignore(255,'\n');
	for(ia=0; ia<35; ia++){
		for(iy=0; iy<81; iy++){
			file>>FertilityTable[ia][iy][0];}
	}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for (ia = 0; ia<35; ia++){
		for (iy = 0; iy<81; iy++){
			file >> FertilityTable[ia][iy][1];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<35; ia++){
		for (iy = 0; iy<81; iy++){
			file >> FertilityTable[ia][iy][2];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> MaleBirthPropn[0] >> MaleBirthPropn[1] >> MaleBirthPropn[2];
	file.close();
}

void ReadPopWeights()
{
	int ia, iy;
	ifstream file;

	file.open("PopWeights.txt");
	if (file.fail()) {
		cerr << "Could not open PopWeights.txt\n";
		exit(1);
	}
	file.ignore(255, '\n');
	for (ia = 0; ia<18; ia++){
		for (iy = 0; iy<81; iy++){
			file >> ThembisaTot[ia][0][iy];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<18; ia++){
		for (iy = 0; iy<81; iy++){
			file >> ThembisaTot[ia][1][iy];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<18; ia++){
		for (iy = 0; iy<41; iy++){
			file >> RaceWeights[ia][0][iy][0];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<18; ia++){
		for (iy = 0; iy<41; iy++){
			file >> RaceWeights[ia][1][iy][0];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<18; ia++){
		for (iy = 0; iy<41; iy++){
			file >> RaceWeights[ia][0][iy][1];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<18; ia++){
		for (iy = 0; iy<41; iy++){
			file >> RaceWeights[ia][1][iy][1];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<18; ia++){
		for (iy = 0; iy<41; iy++){
			file >> RaceWeights[ia][0][iy][2];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<18; ia++){
		for (iy = 0; iy<41; iy++){
			file >> RaceWeights[ia][1][iy][2];
		}
	}
	file.close();
}

void ReadOneStartProfileM(ifstream* file, int group)
{
	int ia, iz, offset; 
	double dummy;

	offset = 16 * group;
	for(ia=0; ia<16; ia++){
		*file>>dummy;
		for(iz=0; iz<5; iz++){
			*file>>HSVtransitionM.PropnByStage[ia+offset][iz];}
		for(iz=0; iz<7; iz++){
			*file>>TPtransitionM.PropnByStage[ia+offset][iz];}
		for(iz=0; iz<4; iz++){
			*file>>HDtransitionM.PropnByStage[ia+offset][iz];}
		for(iz=0; iz<4; iz++){
			*file>>NGtransitionM.PropnByStage[ia+offset][iz];}
		for(iz=0; iz<4; iz++){
			*file>>CTtransitionM.PropnByStage[ia+offset][iz];}
		for(iz=0; iz<4; iz++){
			*file>>TVtransitionM.PropnByStage[ia+offset][iz];}
	}
}

void ReadOneStartProfileF(ifstream* file, int group)
{
	int ia, iz, offset; 
	double dummy;

	offset = 16 * group;
	for(ia=0; ia<16; ia++){
		*file>>dummy;
		for(iz=0; iz<5; iz++){
			*file>>HSVtransitionF.PropnByStage[ia+offset][iz];}
		for(iz=0; iz<7; iz++){
			*file>>TPtransitionF.PropnByStage[ia+offset][iz];}
		for(iz=0; iz<4; iz++){
			*file>>HDtransitionF.PropnByStage[ia+offset][iz];}
		for(iz=0; iz<4; iz++){
			*file>>NGtransitionF.PropnByStage[ia+offset][iz];}
		for(iz=0; iz<4; iz++){
			*file>>CTtransitionF.PropnByStage[ia+offset][iz];}
		for(iz=0; iz<4; iz++){
			*file>>TVtransitionF.PropnByStage[ia+offset][iz];}
		for(iz=0; iz<3; iz++){
			*file>>VCtransitionF.PropnByStage[ia+offset][iz];}
		for(iz=0; iz<4; iz++){
			*file>>BVtransitionF.PropnByStage[ia+offset][iz];}
	}
}

void ReadStartProfile(const char* input)
{
	int ia, iz;
	double dummy;
	ifstream file;

	file.open(input);
	if (file.fail()) {
		cerr << "File open error\n";
		exit(1);
	}
	file.ignore(255,'\n');
	for(ia=0; ia<4; ia++){
		file>>dummy;}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(ia=0; ia<4; ia++){
		file>>dummy;}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(ia=0; ia<4; ia++){
		file>>dummy;
		for(iz=0; iz<3; iz++){
			file>>VCtransitionF.PropnByStage[ia][iz];}
		for(iz=0; iz<4; iz++){
			file>>BVtransitionF.PropnByStage[ia][iz];}
	}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	for(ia=0; ia<4; ia++){
		file>>dummy;
		for(iz=0; iz<3; iz++){
			file>>VCtransitionF.PropnByStage[16+ia][iz];}
		for(iz=0; iz<4; iz++){
			file>>BVtransitionF.PropnByStage[16+ia][iz];}
	}
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileM(&file, 2); // MH, 0
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileM(&file, 3); // MH, 1, S
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileM(&file, 4); // MH, 2, S
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileM(&file, 5); // MH, 1, L
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileM(&file, 6); // MH, 2, L
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileM(&file, 7); // MH, 11, S
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileM(&file, 8); // MH, 12, S
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileM(&file, 9); // MH, 22, S
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileM(&file, 10); // MH, 11, L
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileM(&file, 11); // MH, 12, L
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileM(&file, 12); // MH, 21, L
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileM(&file, 13); // MH, 22, L
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileM(&file, 14); // ML, 0
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileM(&file, 15); // ML, 1, S
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileM(&file, 16); // ML, 2, S
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileM(&file, 17); // ML, 1, L
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileM(&file, 18); // ML, 2, L
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileF(&file, 19); // FSW
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileF(&file, 2); // FH, 0
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileF(&file, 3); // FH, 1, S
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileF(&file, 4); // FH, 2, S
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileF(&file, 5); // FH, 1, L
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileF(&file, 6); // FH, 2, L
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileF(&file, 7); // FH, 11, S
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileF(&file, 8); // FH, 12, S
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileF(&file, 9); // FH, 22, S
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileF(&file, 10); // FH, 11, L
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileF(&file, 11); // FH, 12, L
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileF(&file, 12); // FH, 21, L
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileF(&file, 13); // FH, 22, L
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileF(&file, 14); // FL, 0
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileF(&file, 15); // FL, 1, S
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileF(&file, 16); // FL, 2, S
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileF(&file, 17); // FL, 1, L
	file.ignore(255,'\n');
	file.ignore(255,'\n');
	ReadOneStartProfileF(&file, 18); // FL, 2, L
	// Note that we aren't reading the rest of the file, as the initial child numbers
	// are now in StartPop and the male birth propn is in FertTables.
	file.close();
	
	// Set the proportions of uninfected virgins to 100%.
	for(ia=0; ia<4; ia++){
		if(HSVind==1){
			HSVtransitionM.PropnByStage[ia][0] = 1.0;
			HSVtransitionM.PropnByStage[ia+16][0] = 1.0;
			HSVtransitionF.PropnByStage[ia][0] = 1.0;
			HSVtransitionF.PropnByStage[ia+16][0] = 1.0;
		}
		if(TPind==1){
			TPtransitionM.PropnByStage[ia][0] = 1.0;
			TPtransitionM.PropnByStage[ia+16][0] = 1.0;
			TPtransitionF.PropnByStage[ia][0] = 1.0;
			TPtransitionF.PropnByStage[ia+16][0] = 1.0;
		}
		if(HDind==1){
			HDtransitionM.PropnByStage[ia][0] = 1.0;
			HDtransitionM.PropnByStage[ia+16][0] = 1.0;
			HDtransitionF.PropnByStage[ia][0] = 1.0;
			HDtransitionF.PropnByStage[ia+16][0] = 1.0;
		}
		if(NGind==1){
			NGtransitionM.PropnByStage[ia][0] = 1.0;
			NGtransitionM.PropnByStage[ia+16][0] = 1.0;
			NGtransitionF.PropnByStage[ia][0] = 1.0;
			NGtransitionF.PropnByStage[ia+16][0] = 1.0;
		}
		if(CTind==1){
			CTtransitionM.PropnByStage[ia][0] = 1.0;
			CTtransitionM.PropnByStage[ia+16][0] = 1.0;
			CTtransitionF.PropnByStage[ia][0] = 1.0;
			CTtransitionF.PropnByStage[ia+16][0] = 1.0;
		}
		if(TVind==1){
			TVtransitionM.PropnByStage[ia][0] = 1.0;
			TVtransitionM.PropnByStage[ia+16][0] = 1.0;
			TVtransitionF.PropnByStage[ia][0] = 1.0;
			TVtransitionF.PropnByStage[ia+16][0] = 1.0;
		}
	}
}

void ReadInitHIV()
{
	int ia;
	ifstream file;

	file.open("InitHIV.txt");
	if (file.fail()) {
		cerr << "Could not open InitHIV.txt\n";
		exit(1);
	}
	file.ignore(255, '\n');
	for (ia = 0; ia<7; ia++){
		file >> HlabisaRatio[ia][0] >> HlabisaRatio[ia][1];
	}
	file.close();
}

void ReadStartPop()
{
	int ia, ir;
	ifstream file;

	file.open("StartPop.txt");
	if (file.fail()) {
		cerr << "Could not open StartPop.txt\n";
		exit(1);
	}
	for(ia=0; ia<91; ia++){
		for (ir = 0; ir < 6; ir++){
			file >> StartPop[ia][ir];}
	}
	file.close();
}

void ReadSTDprev()
{
	if (HIVcalib == 1){
		HIVtransitionM.ReadPrevData("HIVdataM.txt");
		HIVtransitionF.ReadPrevData("HIVdataF.txt");
	}
	if (HSVcalib == 1){
		HSVtransitionM.ReadPrevData("HSVdataM.txt");
		HSVtransitionF.ReadPrevData("HSVdataF.txt");
	}
	if (TPcalib == 1){
		TPtransitionM.ReadPrevData("TPdataM.txt");
		TPtransitionF.ReadPrevData("TPdataF.txt");
	}
	if (HDcalib == 1){
		HDtransitionM.ReadPrevData("HDdataM.txt");
		HDtransitionF.ReadPrevData("HDdataF.txt");
	}
	if (NGcalib == 1){
		NGtransitionM.ReadPrevData("NGdataM.txt");
		NGtransitionF.ReadPrevData("NGdataF.txt");
	}
	if (CTcalib == 1){
		CTtransitionM.ReadPrevData("CTdataM.txt");
		CTtransitionF.ReadPrevData("CTdataF.txt");
	}
	if (TVcalib == 1){
		TVtransitionM.ReadPrevData("TVdataM.txt");
		TVtransitionF.ReadPrevData("TVdataF.txt");
	}
	if (BVcalib == 1){ BVtransitionF.ReadPrevData("BVdataF.txt"); }
	if (VCcalib == 1){ VCtransitionF.ReadPrevData("VCdataF.txt"); }
}

void ReadSTDparameters()
{
	int ic, iz, dummy;
	ifstream file, file2, file3, file4, file5, file6, file7, file8, file9, file10, file11;

	if (HIVcalib == 1){
		file.open("RandomUniformHIV.txt");
		if (file.fail()) {
			cerr << "Could not open RandomUniformHIV.txt\n";
			exit(1);
		}
		for (ic = 0; ic<ParamCombs; ic++){
			file >> SeedRecord[ic][0] >> SeedRecord[ic][1];
			for (iz = 0; iz < RandomUniformHIV.columns; iz++){
				file >> RandomUniformHIV.out[ic][iz];
			}
		}
		file.close();
	}

	if (HSVcalib == 1){
		file2.open("RandomUniformHSV.txt");
		if (file2.fail()) {
			cerr << "Could not open RandomUniformHSV.txt\n";
			exit(1);
		}
		for (ic = 0; ic<ParamCombs; ic++){
			file2 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
			for (iz = 0; iz < RandomUniformHSV.columns; iz++){
				file2 >> RandomUniformHSV.out[ic][iz];
			}
		}
		file2.close();
	}

	if (NGcalib == 1){
		file4.open("RandomUniformNG.txt");
		if (file4.fail()) {
			cerr << "Could not open RandomUniformNG.txt\n";
			exit(1);
		}
		for (ic = 0; ic<ParamCombs; ic++){
			file4 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
			for (iz = 0; iz < RandomUniformNG.columns; iz++){
				file4 >> RandomUniformNG.out[ic][iz];
			}
		}
		file4.close();
	}

	if (CTcalib == 1){
		file5.open("RandomUniformCT.txt");
		if (file5.fail()) {
			cerr << "Could not open RandomUniformCT.txt\n";
			exit(1);
		}
		for (ic = 0; ic<ParamCombs; ic++){
			file5 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
			for (iz = 0; iz < RandomUniformCT.columns; iz++){
				file5 >> RandomUniformCT.out[ic][iz];
			}
		}
		file5.close();
	}

	if (TVcalib == 1){
		file6.open("RandomUniformTV.txt");
		if (file6.fail()) {
			cerr << "Could not open RandomUniformTV.txt\n";
			exit(1);
		}
		for (ic = 0; ic<ParamCombs; ic++){
			file6 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
			for (iz = 0; iz < RandomUniformTV.columns; iz++){
				file6 >> RandomUniformTV.out[ic][iz];
			}
		}
		file6.close();
	}

	if (TPcalib == 1){
		file3.open("RandomUniformTP.txt");
		if (file3.fail()) {
			cerr << "Could not open RandomUniformTP.txt\n";
			exit(1);
		}
		for (ic = 0; ic<ParamCombs; ic++){
			file3 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
			for (iz = 0; iz < RandomUniformTP.columns; iz++){
				file3 >> RandomUniformTP.out[ic][iz];
			}
		}
		file3.close();
	}

	if (BVcalib == 1){
		file7.open("RandomUniformBV.txt");
		if (file7.fail()) {
			cerr << "Could not open RandomUniformBV.txt\n";
			exit(1);
		}
		for (ic = 0; ic<ParamCombs; ic++){
			file7 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
			for (iz = 0; iz < RandomUniformBV.columns; iz++){
				file7 >> RandomUniformBV.out[ic][iz];
			}
		}
		file7.close();
	}

	if (VCcalib == 1){
		file8.open("RandomUniformVC.txt");
		if (file8.fail()) {
			cerr << "Could not open RandomUniformVC.txt\n";
			exit(1);
		}
		for (ic = 0; ic<ParamCombs; ic++){
			file8 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
			for (iz = 0; iz < RandomUniformVC.columns; iz++){
				file8 >> RandomUniformVC.out[ic][iz];
			}
		}
		file8.close();
	}

	if (MSMcalib == 1){
		file9.open("RandomUniformMSM.txt");
		if (file9.fail()) {
			cerr << "Could not open RandomUniformMSM.txt\n";
			exit(1);
		}
		for (ic = 0; ic<ParamCombs; ic++){
			file9 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
			for (iz = 0; iz < RandomUniformMSM.columns; iz++){
				file9 >> RandomUniformMSM.out[ic][iz];
			}
		}
		file9.close();
	}

	if (MSMcalibHIV == 1){
		file10.open("RandomUniformMSMT.txt");
		if (file10.fail()) {
			cerr << "Could not open RandomUniformMSMT.txt\n";
			exit(1);
		}
		for (ic = 0; ic<ParamCombs; ic++){
			file10 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
			for (iz = 0; iz < RandomUniformMSMT.columns; iz++){
				file10 >> RandomUniformMSMT.out[ic][iz];
			}
		}
		file10.close();
	}

	if (StructuralRCTcalib == 1){
		file11.open("RandomUniformStruct.txt");
		if (file11.fail()) {
			cerr << "Could not open RandomUniformStruct.txt\n";
			exit(1);
		}
		for (ic = 0; ic<ParamCombs; ic++){
			file11 >> SeedRecord[ic][0] >> SeedRecord[ic][1];
			for (iz = 0; iz < RandomUniformStruct.columns; iz++){
				file11 >> RandomUniformStruct.out[ic][iz];
			}
		}
		file11.close();
	}
}

void ReadMCprev()
{
	int ia, ir;
	ifstream file;

	file.open("MCbaseline.txt");
	if (file.fail()) {
		cerr << "Could not open MCbaseline.txt\n";
		exit(1);
	}
	for (ia = 0; ia<91; ia++){
		for (ir = 0; ir < 3; ir++){
			file >> MCprevBaseline[ia][ir];
		}
	}
	file.close();
}

void ReadEduAssumps()
{
	int ii, ir, ig, ia, ij;
	ifstream file;

	file.open("EduAssumps.txt");
	if (file.fail()) {
		cerr << "Could not open EduAssumps.txt\n";
		exit(1);
	}

	file.ignore(255, '\n');
	for (ia = 0; ia<5; ia++){
		for (ir = 0; ir < 3; ir++){
			file >> Grade1entry[ia][ir];}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii<13; ii++){
		file >> GradeRepetition[ii][0][0] >> GradeRepetition[ii][0][1];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii<13; ii++){
		file >> GradeRepetition[ii][1][0] >> GradeRepetition[ii][1][1];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii<13; ii++){
		file >> GradeRepetition[ii][2][0] >> GradeRepetition[ii][2][1];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ConscientiousEffectGradeRep;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ParentEduEffectGradeRep;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii<13; ii++){
		file >> Dropout[ii][0][0] >> Dropout[ii][0][1];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii<13; ii++){
		file >> Dropout[ii][1][0] >> Dropout[ii][1][1];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii<13; ii++){
		file >> Dropout[ii][2][0] >> Dropout[ii][2][1];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ParentEduEffectDropout;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ir = 0; ir<3; ir++){
		file >> PregDropout[ir];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ir = 0; ir<3; ir++){
		file >> TertiaryEnrol[ir];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ConscientiousEffectTertiaryEd;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ir = 0; ir<3; ir++){
		file >> DropoutRedn[ir];}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> EduAssort;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RRdropoutSupport;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> RRdropoutIncSupport;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> ProbReturnSupport;
	file.close();

	for (ii = 0; ii < 14; ii++){
		for (ij = 0; ij < 14; ij++){
			EduMixing[ii][ij] = exp(-EduAssort * pow(1.0 * (ii - ij)/(13.5 - 0.5 * (ii + ij)), 2.0));
		}
	}

}

void ReadStartEdu()
{
	int ia, ii;
	ifstream file;

	file.open("StartEduProfile.txt");
	if (file.fail()) {
		cerr << "Could not open StartEduProfile.txt\n";
		exit(1);
	}

	file.ignore(255, '\n');
	for (ia = 0; ia<91; ia++){
		for (ii = 0; ii < 13; ii++){
			file >> CumGrade1985[ia][ii][0];}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<91; ia++){
		for (ii = 0; ii < 13; ii++){
			file >> CumGrade1985[ia][ii][1];}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<91; ia++){
		for (ii = 0; ii < 13; ii++){
			file >> CumGrade1985[ia][ii][2];}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<31; ia++){
		for (ii = 0; ii < 13; ii++){
			file >> InSchool1985[ia][ii][0];}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<31; ia++){
		for (ii = 0; ii < 13; ii++){
			file >> InSchool1985[ia][ii][1];}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ia = 0; ia<31; ia++){
		for (ii = 0; ii < 13; ii++){
			file >> InSchool1985[ia][ii][2];}
	}
	file.close();

}

void ReadCostAssumps()
{
	int im, ii;
	ifstream file;

	file.open("CostAssumps.txt");
	if (file.fail()) {
		cerr << "Could not open CostAssumps.txt\n";
		exit(1);
	}

	file.ignore(255, '\n');
	for (im = 0; im<8; im++){
		for (ii = 0; ii < 2; ii++){
			file >> BaseModalityCost[im][ii];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (im = 0; im<12; im++){
		for (ii = 0; ii < 2; ii++){
			file >> NewModalityCost[im][ii];
		}
	}
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> AnnARTcost;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> AnnCTXcost;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> CondomCost;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> MMCcost;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> PMTCTcost[0] >> PMTCTcost[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> AnnPrEPcost;
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> InfantTestCost[0] >> InfantTestCost[1];
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	file >> PalliativeCareCost;
	file.close();
}

void ReadHCTparameters()
{
	int idum, ii, ic;
	ifstream file;

	file.open("RandomUniformHCT.txt");
	if (file.fail()) {
		cerr << "Could not open RandomUniformHCT.txt\n";
		exit(1);
	}

	for (ii = 0; ii<ParamCombs; ii++){
		file >> idum >> idum;
		for (ic = 0; ic < RandomUniformHCT.columns; ic++){
			file >> RandomUniformHCT.out[ii][ic];
		}
	}
	file.close();
}

void ReadStructuralDriverData()
{
	int ii, ic;
	ifstream file;

	file.open("StructuralDriverData.txt");
	if (file.fail()) {
		cerr << "Could not open StructuralDriverData.txt\n";
		exit(1);
	}

	file.ignore(255, '\n');
	for (ii = 0; ii<4; ii++){ file >> IneqGenderBingeSvy[ii][0]; }
	for (ii = 0; ii<4; ii++){ file >> IneqGenderBingeSvy[ii][1]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii<4; ii++){ file >> IneqGenderMultSvy[ii][0]; }
	for (ii = 0; ii<4; ii++){ file >> IneqGenderMultSvy[ii][1]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii<4; ii++){ file >> BingeMultSvyM[ii][0]; }
	for (ii = 0; ii<4; ii++){ file >> BingeMultSvyM[ii][1]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii<4; ii++){ file >> BingeMultSvyF[ii][0]; }
	for (ii = 0; ii<4; ii++){ file >> BingeMultSvyF[ii][1]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii<4; ii++){ file >> BingeCondomSvyM[ii][0]; }
	for (ii = 0; ii<4; ii++){ file >> BingeCondomSvyM[ii][1]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii<4; ii++){ file >> BingeCondomSvyF[ii][0]; }
	for (ii = 0; ii<4; ii++){ file >> BingeCondomSvyF[ii][1]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii<4; ii++){ file >> EmployedTransSvy[ii][0]; }
	for (ii = 0; ii<4; ii++){ file >> EmployedTransSvy[ii][1]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii<4; ii++){ file >> EmployedMultSvyM[ii][0]; }
	for (ii = 0; ii<4; ii++){ file >> EmployedMultSvyM[ii][1]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii<4; ii++){ file >> EmployedMultSvyF[ii][0]; }
	for (ii = 0; ii<4; ii++){ file >> EmployedMultSvyF[ii][1]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii<4; ii++){ file >> EmployedHIVsvyM[ii][0]; }
	for (ii = 0; ii<4; ii++){ file >> EmployedHIVsvyM[ii][1]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii<4; ii++){ file >> EmployedHIVsvyF[ii][0]; }
	for (ii = 0; ii<4; ii++){ file >> EmployedHIVsvyF[ii][1]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii<4; ii++){ file >> EduCondomSvyM[ii][0]; }
	for (ii = 0; ii<4; ii++){ file >> EduCondomSvyM[ii][1]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii<4; ii++){ file >> EduCondomSvyF[ii][0]; }
	for (ii = 0; ii<4; ii++){ file >> EduCondomSvyF[ii][1]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii<4; ii++){ file >> EduHIVsvyM[ii][0]; }
	for (ii = 0; ii<4; ii++){ file >> EduHIVsvyM[ii][1]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii<4; ii++){ file >> EduHIVsvyF[ii][0]; }
	for (ii = 0; ii<4; ii++){ file >> EduHIVsvyF[ii][1]; }
	file.ignore(255, '\n');
	file.ignore(255, '\n');
	for (ii = 0; ii<4; ii++){ file >> SchoolMarriageSvy[ii][0]; }
	for (ii = 0; ii<4; ii++){ file >> SchoolMarriageSvy[ii][1]; }
	
	file.close();
}

void ReadStructuralRCTdata()
{
	int ii, ic, n;
	ifstream file;
	ifstream file2;
	ifstream file3;
	ifstream file4;
	ifstream file5;
	ifstream file6;
	ifstream file7;

	file.open("SingleSessionAlcoholRCTs.txt");
	if (file.fail()) {
		cerr << "Could not open SingleSessionAlcoholRCTs.txt\n";
		exit(1);
	}
	
	for (ii = 0; ii<SingleSessionAlcoholCounselling.DataPoints; ii++){
		for (ic = 0; ic < 10; ic++){
			file >> SingleSessionAlcoholCounselling.StudyDetails[ii][ic];
		}
	}
	file.close();
	SingleSessionAlcoholCounselling.GetMaxTerm();

	file2.open("MultiSessionAlcoholRCTs.txt");
	if (file2.fail()) {
		cerr << "Could not open MultiSessionAlcoholRCTs.txt\n";
		exit(1);
	}

	for (ii = 0; ii<MultiSessionAlcoholCounselling.DataPoints; ii++){
		for (ic = 0; ic < 10; ic++){
			file2 >> MultiSessionAlcoholCounselling.StudyDetails[ii][ic];
		}
	}
	file2.close();
	MultiSessionAlcoholCounselling.GetMaxTerm();

	file3.open("CashTransferRCTs.txt");
	if (file3.fail()) {
		cerr << "Could not open CashTransferRCTs.txt\n";
		exit(1);
	}

	for (ii = 0; ii<CashTransfers.DataPoints; ii++){
		for (ic = 0; ic < 10; ic++){
			file3 >> CashTransfers.StudyDetails[ii][ic];
		}
	}
	file3.close();
	CashTransfers.GetMaxTerm();

	file4.open("SchoolSupportRCTs.txt");
	if (file4.fail()) {
		cerr << "Could not open SchoolSupportRCTs.txt\n";
		exit(1);
	}

	for (ii = 0; ii<SchoolSupport.DataPoints; ii++){
		for (ic = 0; ic < 10; ic++){
			file4 >> SchoolSupport.StudyDetails[ii][ic];
		}
	}
	file4.close();
	SchoolSupport.GetMaxTerm();

	file5.open("VocationalTrainingRCTs.txt");
	if (file5.fail()) {
		cerr << "Could not open VocationalTrainingRCTs.txt\n";
		exit(1);
	}

	for (ii = 0; ii<VocationalTraining.DataPoints; ii++){
		for (ic = 0; ic < 10; ic++){
			file5 >> VocationalTraining.StudyDetails[ii][ic];
		}
	}
	file5.close();
	VocationalTraining.GetMaxTerm();

	file6.open("GenderTransformCommunRCTs.txt");
	if (file6.fail()) {
		cerr << "Could not open GenderTransformCommunRCTs.txt\n";
		exit(1);
	}

	for (ii = 0; ii<GenderTransformCommun.DataPoints; ii++){
		for (ic = 0; ic < 10; ic++){
			file6 >> GenderTransformCommun.StudyDetails[ii][ic];
		}
	}
	file6.close();
	GenderTransformCommun.GetMaxTerm();

	file7.open("GenderTransformIndivRCTs.txt");
	if (file7.fail()) {
		cerr << "Could not open GenderTransformIndivRCTs.txt\n";
		exit(1);
	}

	for (ii = 0; ii<GenderTransformIndiv.DataPoints; ii++){
		for (ic = 0; ic < 10; ic++){
			file7 >> GenderTransformIndiv.StudyDetails[ii][ic];
		}
	}
	file7.close();
	GenderTransformIndiv.GetMaxTerm();
}

void ReadAllInputFiles()
{
	ReadSexAssumps("SexAssumps.txt");
	ReadSTDepi("STDepidemiology.txt");
	if (InclMSM == 1){ ReadMSM(); }
	if (MSMcalib == 1){ ReadMSMcalibData(); }
	if (MSMcalibHIV == 1){ ReadMSM_HIVdata(); }
	ReadRatesByYear();
	ReadMortTables();
	ReadFertTables();
	if (IncludePopWeights == 1){ ReadPopWeights(); }
	ReadStartProfile("StartProfile.txt");
	ReadInitHIV();
	ReadStartPop();
	ReadSTDprev();
	ReadMCprev();
	ReadEduAssumps();
	ReadStartEdu();
	ReadConception();
	ReadBreastfeeding();
	ReadMigration();
	ReadDisclosure();
	ReadGenderNormAssumps();
	ReadAlcoholAssumps();
	ReadHouseholdAssumps();
	ReadPrisonAssumps();
	ReadEmployment();
	ReadIncome();
	if (HCTuncertainty == 1){ SimulateHCTparameters(); }
	if (FixedUncertainty == 1){ ReadCostAssumps(); }
	if (StructuralDriverCalib == 1){ ReadStructuralDriverData(); }
	if (StructuralRCTcalib == 1){ 
		ReadStructuralRCTdata(); 
		SimulateStructParams();
	}
}

/*double GetQstatistic(double Cmatrix[3][3], int MatDim)
{
	int ir, ic;
	double RowTotal, trace, Qstatistic;

	trace = 0.0;
	for(ir=0; ir<MatDim; ir++){
		RowTotal = 0.0;
		for(ic=0; ic<MatDim; ic++){
			RowTotal += Cmatrix[ir][ic];}
		trace += Cmatrix[ir][ir]/RowTotal;
	}

	Qstatistic = (trace - 1.0)/(MatDim - 1.0);

	return Qstatistic;
}*/

void SetCalibParameters()
{
	int ia;
	double AveOR;

	AveOR = 0.68; // Posterior mean from the ASSA2002 uncertainty analysis

	// Values in the 5 lines below are from the simulateparameters function in the C++
	// version of the ASSA2002 model (see derivation in section 3.3.1.3 of the uncertainty
	// analysis report)
	/*HIVtransitionF.AntenatalNlogL.BiasMult[1] = 0.885450 / (AveOR*0.044021 + 0.841431);
	HIVtransitionF.AntenatalNlogL.BiasMult[2] = 0.902648 / (AveOR*0.080730 + 0.810868);
	HIVtransitionF.AntenatalNlogL.BiasMult[3] = 0.870821 / (AveOR*0.108702 + 0.707487);
	HIVtransitionF.AntenatalNlogL.BiasMult[4] = 0.921259 / (AveOR*0.142296 + 0.698347);
	HIVtransitionF.AntenatalNlogL.BiasMult[5] = 0.937727 / (AveOR*0.117837 + 0.784181);*/

	// Setting the default values for the variance of the study effects (only relevant to
	// the sentinel surveillance data)
	if (HIVcalib == 1 && HIVtransitionF.ANClogL.VarStudyEffect == 0){
		HIVtransitionM.SetVarStudyEffect(0.09);
		HIVtransitionF.SetVarStudyEffect(0.09);
	}
	if (HSVcalib == 1 && HSVtransitionF.ANClogL.VarStudyEffect == 0){
		HSVtransitionM.SetVarStudyEffect(0.48);
		HSVtransitionF.SetVarStudyEffect(0.48);
	}
	if (TPcalib == 1 && TPtransitionF.ANClogL.VarStudyEffect == 0){
		TPtransitionM.SetVarStudyEffect(0.14);
		TPtransitionF.SetVarStudyEffect(0.14);
	}
	if (HDcalib == 1 && HDtransitionF.ANClogL.VarStudyEffect == 0){
		HDtransitionM.SetVarStudyEffect(0.30);
		HDtransitionF.SetVarStudyEffect(0.30);
	}
	if (NGcalib == 1 && NGtransitionF.ANClogL.VarStudyEffect == 0){
		NGtransitionM.SetVarStudyEffect(0.09);
		NGtransitionF.SetVarStudyEffect(0.09);
	}
	if (CTcalib == 1 && CTtransitionF.ANClogL.VarStudyEffect == 0){
		CTtransitionM.SetVarStudyEffect(0.11);
		CTtransitionF.SetVarStudyEffect(0.11);
	}
	if (TVcalib == 1 && TVtransitionF.ANClogL.VarStudyEffect == 0){
		TVtransitionM.SetVarStudyEffect(0.28);
		TVtransitionF.SetVarStudyEffect(0.28);
	}
	if (BVcalib == 1 && BVtransitionF.ANClogL.VarStudyEffect == 0){
		BVtransitionF.SetVarStudyEffect(0.0529);
	}
	if (VCcalib == 1 && VCtransitionF.ANClogL.VarStudyEffect == 0){
		VCtransitionF.SetVarStudyEffect(0.0529);
	}
}

void CalcStructuralLogL()
{
	int year;

	StructLogL = 0.0;

	for(year=0; year<4; year++){
		IneqGenderBingeAssn.CalcSvyAssnLogL(IneqGenderBingeSvy, year);
		IneqGenderMultAssn.CalcSvyAssnLogL(IneqGenderMultSvy, year);
		BingeMultAssnM.CalcSvyAssnLogL(BingeMultSvyM, year);
		BingeMultAssnF.CalcSvyAssnLogL(BingeMultSvyF, year);
		BingeCondomAssnM.CalcSvyAssnLogL(BingeCondomSvyM, year);
		BingeCondomAssnF.CalcSvyAssnLogL(BingeCondomSvyF, year);
		EmployedTransAssn.CalcSvyAssnLogL(EmployedTransSvy, year);
		EmployedMultAssnM.CalcSvyAssnLogL(EmployedMultSvyM, year);
		EmployedMultAssnF.CalcSvyAssnLogL(EmployedMultSvyF, year);
		EmployedHIVassnM.CalcSvyAssnLogL(EmployedHIVsvyM, year);
		EmployedHIVassnF.CalcSvyAssnLogL(EmployedHIVsvyF, year);
		EduCondomAssnM.CalcSvyAssnLogL(EduCondomSvyM, year);
		EduCondomAssnF.CalcSvyAssnLogL(EduCondomSvyF, year);
		EduHIVassnM.CalcSvyAssnLogL(EduHIVsvyM, year);
		EduHIVassnF.CalcSvyAssnLogL(EduHIVsvyF, year);
		SchoolMarriageAssn.CalcSvyAssnLogL(SchoolMarriageSvy, year);
	}
}


void CalcTotalLogL()
{
	int ii;
	double ModelErrorVar, TempL[7];

	TotalLogL = 0.0;
	if (HIVcalib == 1){
		HIVtransitionF.CSWlogL.CalcLogL();
		if (GetSDfromData == 1){
			HIVtransitionF.AntenatalNlogL.ModelVarEst = 0.005; // Estimate of model bias in Thembisa v12 (PLoS Med paper)
			//HIVtransitionF.AntenatalNlogL.CalcModelVar();
			/*HIVtransitionF.HouseholdNlogL.CalcModelVar();
			HIVtransitionM.HouseholdNlogL.CalcModelVar();
			ModelErrorVar = (HIVtransitionF.HouseholdNlogL.ModelVarEst * HIVtransitionF.HouseholdNlogL.Observations +
				HIVtransitionM.HouseholdNlogL.ModelVarEst * HIVtransitionM.HouseholdNlogL.Observations +
				HIVtransitionF.AntenatalNlogL.ModelVarEst * HIVtransitionF.AntenatalNlogL.Observations) /
				(HIVtransitionF.HouseholdNlogL.Observations + HIVtransitionM.HouseholdNlogL.Observations +
				HIVtransitionF.AntenatalNlogL.Observations);
			if (ModelErrorVar < 0.0){ ModelErrorVar = 0.0; }
			HIVtransitionF.HouseholdNlogL.ModelVarEst = ModelErrorVar;
			HIVtransitionM.HouseholdNlogL.ModelVarEst = ModelErrorVar;
			HIVtransitionF.AntenatalNlogL.ModelVarEst = ModelErrorVar;*/
		}
		HIVtransitionF.AntenatalNlogL.CalcLogL();
		HIVtransitionF.HouseholdNlogL.CalcLogL();
		HIVtransitionM.HouseholdNlogL.CalcLogL();
		TotalLogL += HIVtransitionF.CSWlogL.LogL + HIVtransitionF.AntenatalNlogL.LogL +
			HIVtransitionF.HouseholdNlogL.LogL + HIVtransitionM.HouseholdNlogL.LogL;
	}
	if (HSVcalib == 1){
		HSVtransitionF.ANClogL.CalcLogL();
		HSVtransitionF.FPClogL.CalcLogL();
		HSVtransitionF.CSWlogL.CalcLogL();
		HSVtransitionF.HouseholdLogL.CalcLogL();
		HSVtransitionM.HouseholdLogL.CalcLogL();
		HSVparamsLogL.out[CurrSim-1][11] = HSVtransitionF.ANClogL.LogL + HSVtransitionF.FPClogL.LogL +
			HSVtransitionF.CSWlogL.LogL + HSVtransitionF.HouseholdLogL.LogL + 
			HSVtransitionM.HouseholdLogL.LogL;
	}
	if (TPcalib == 1){
		TPtransitionF.AntenatalNlogL.CalcLogL();
		TPtransitionF.ANClogL.CalcLogL();
		TPtransitionF.FPClogL.CalcLogL();
		TPtransitionF.CSWlogL.CalcLogL();
		TPtransitionF.HouseholdLogL.CalcLogL();
		TPtransitionM.HouseholdLogL.CalcLogL();
		TPparamsLogL.out[CurrSim-1][10] = TPtransitionF.AntenatalNlogL.LogL + TPtransitionF.ANClogL.LogL +
			TPtransitionF.FPClogL.LogL + TPtransitionF.CSWlogL.LogL +
			TPtransitionF.HouseholdLogL.LogL + TPtransitionM.HouseholdLogL.LogL;
	}
	/*if (HDcalib == 1){
		HDtransitionF.GUDlogL.CalcLogL();
		HDtransitionM.GUDlogL.CalcLogL();
		TotalLogL += HDtransitionF.GUDlogL.LogL + HDtransitionM.GUDlogL.LogL;
	}*/
	if (NGcalib == 1){
		NGtransitionF.ANClogL.CalcLogL();
		NGtransitionF.FPClogL.CalcLogL();
		NGtransitionF.CSWlogL.CalcLogL();
		NGtransitionF.HouseholdLogL.CalcLogL();
		NGtransitionM.HouseholdLogL.CalcLogL();
		NGparamsLogL.out[CurrSim - 1][10] = NGtransitionF.ANClogL.LogL + NGtransitionF.FPClogL.LogL +
			NGtransitionF.CSWlogL.LogL + NGtransitionF.HouseholdLogL.LogL +
			NGtransitionM.HouseholdLogL.LogL;
	}
	if (CTcalib == 1){
		CTtransitionF.ANClogL.CalcLogL();
		CTtransitionF.FPClogL.CalcLogL();
		CTtransitionF.CSWlogL.CalcLogL();
		CTtransitionF.HouseholdLogL.CalcLogL();
		CTtransitionM.HouseholdLogL.CalcLogL();
		CTparamsLogL.out[CurrSim - 1][10] = CTtransitionF.ANClogL.LogL + CTtransitionF.FPClogL.LogL +
			CTtransitionF.CSWlogL.LogL + CTtransitionF.HouseholdLogL.LogL +
			CTtransitionM.HouseholdLogL.LogL;
	}
	if (TVcalib == 1){
		TVtransitionF.ANClogL.CalcLogL();
		TVtransitionF.FPClogL.CalcLogL();
		TVtransitionF.CSWlogL.CalcLogL();
		TVtransitionF.HouseholdLogL.CalcLogL();
		TVtransitionM.HouseholdLogL.CalcLogL();
		TVparamsLogL.out[CurrSim - 1][12] = TVtransitionF.ANClogL.LogL + TVtransitionF.FPClogL.LogL +
			TVtransitionF.CSWlogL.LogL + TVtransitionF.HouseholdLogL.LogL +
			TVtransitionM.HouseholdLogL.LogL;
	}
	/*if (BVcalib == 1){
		BVtransitionF.ANClogL.CalcLogL();
		BVtransitionF.FPClogL.CalcLogL();
		BVtransitionF.CSWlogL.CalcLogL();
		TotalLogL += BVtransitionF.ANClogL.LogL + BVtransitionF.FPClogL.LogL +
			BVtransitionF.CSWlogL.LogL;
	}
	if (VCcalib == 1){
		VCtransitionF.ANClogL.CalcLogL();
		VCtransitionF.FPClogL.CalcLogL();
		VCtransitionF.CSWlogL.CalcLogL();
		TotalLogL += VCtransitionF.ANClogL.LogL + VCtransitionF.FPClogL.LogL +
			VCtransitionF.CSWlogL.LogL;
	}*/
	if (MSMcalib == 1){
		// Calculate random effect variance
		MSMmarried.CalcRandomEffectVar();
		MSMrecentBi.CalcRandomEffectVar();
		MSMeverBi.CalcRandomEffectVar();
		ModelErrorVar = (MSMrecentBi.RandomEffectVar * MSMrecentBi.studies +
			MSMeverBi.RandomEffectVar * MSMeverBi.studies) /
			(MSMrecentBi.studies + MSMeverBi.studies);
		MSMrecentBi.RandomEffectVar = ModelErrorVar;
		MSMeverBi.RandomEffectVar = ModelErrorVar;
		//MSMaged25plus.CalcRandomEffectVar();
		MSMcurrentReg.CalcRandomEffectVar();
		EverMSMcurrentReg.CalcRandomEffectVar();
		ModelErrorVar = (MSMcurrentReg.RandomEffectVar * MSMcurrentReg.studies +
			EverMSMcurrentReg.RandomEffectVar * EverMSMcurrentReg.studies) /
			(MSMcurrentReg.studies + EverMSMcurrentReg.studies);
		MSMcurrentReg.RandomEffectVar = ModelErrorVar;
		EverMSMcurrentReg.RandomEffectVar = ModelErrorVar;
		MSMmultPartners.CalcRandomEffectVar();
		// Calculate likelihood
		MSMmarried.CalcLogL();
		MSMrecentBi.CalcLogL();
		MSMeverBi.CalcLogL();
		//MSMaged25plus.CalcLogL();
		MSMcurrentReg.CalcLogL();
		EverMSMcurrentReg.CalcLogL();
		MSMmultPartners.CalcLogL();
		MSMparamsLogL.out[CurrSim - 1][11] = MSMmarried.LogL + MSMrecentBi.LogL +
			MSMeverBi.LogL + MSMcurrentReg.LogL +
			EverMSMcurrentReg.LogL + MSMmultPartners.LogL;
		TotalLogL += MSMmarried.LogL + MSMrecentBi.LogL +
			MSMeverBi.LogL + MSMcurrentReg.LogL +
			EverMSMcurrentReg.LogL + MSMmultPartners.LogL;
	}
	if (MSMcalibHIV == 1){
		// Calculate random effect variance
		MSMpositive.CalcRandomEffectVar();
		EverMSMpositive.CalcRandomEffectVar();
		ModelErrorVar = (MSMpositive.RandomEffectVar * MSMpositive.studies +
			EverMSMpositive.RandomEffectVar * EverMSMpositive.studies) /
			(MSMpositive.studies + EverMSMpositive.studies);
		MSMpositive.RandomEffectVar = ModelErrorVar;
		EverMSMpositive.RandomEffectVar = ModelErrorVar;
		// Calculate likelihood
		MSMpositive.CalcLogL();
		EverMSMpositive.CalcLogL();
		MSMTparamsLogL.out[CurrSim - 1][8] = MSMpositive.LogL + EverMSMpositive.LogL;
		TotalLogL += MSMpositive.LogL + EverMSMpositive.LogL;
	}
	if (StructuralDriverCalib == 1){ 
		CalcStructuralLogL();
		TotalLogL += StructLogL; 
	}
	if (StructuralRCTcalib == 1){
		for (ii = 0; ii < 7; ii++){ TempL[ii] = 0.0; }
		//TempL[0] = SingleSessionAlcoholCounselling.GetLikelihood();
		//TempL[1] = MultiSessionAlcoholCounselling.GetLikelihood();
		//TempL[2] = SingleSessionAlcoholCounselling.AltLogL;
		TempL[2] = CashTransfers.GetLikelihood();
		TempL[3] = SchoolSupport.GetLikelihood();
		TempL[4] = VocationalTraining.GetLikelihood();
		//TempL[5] = GenderTransformCommun.GetLikelihood();
		//TempL[6] = GenderTransformIndiv.GetLikelihood();
		for (ii = 0; ii < 7; ii++){ 
			TotalLogL += TempL[ii]; 
			StructParamsLogL.out[CurrSim - 1][23 + ii] = TempL[ii];
		}
	}
}

void OneSimulation()
{
	int seed, SimCount2;

	for (int i = 0; i < Register.size(); ++i){
		memset(&Register[i], 0, sizeof(Register[0]));}
	ReadAllInputFiles();

	// Reset the seed
	if (FixedUncertainty == 0){
		SimCount2 = (CurrSim - 1) / IterationsPerPC;
		seed = SimCount2 * 91 + process_num * 7927 + (CurrSim - SimCount2 * IterationsPerPC);
	}
	else{
		SimCount2 = SeedRecord[CurrSim - 1][0] / IterationsPerPC;
		seed = SimCount2 * 91 + SeedRecord[CurrSim - 1][1] * 7927 + 1;
	}
	//seed = CurrSim * 91; // When FixedUncertainty = 1 but IterationsPerPC > 1 
	//seed = SeedRecord[CurrSim-1][0] * 91;
	rg.RandomInit(seed);

	RSApop.AssignAgeSex();
	RSApop.AssignPersonality();
	RSApop.GetAllPartnerRates();
	RSApop.AssignVisitFreq();
	RSApop.AssignEdu();
	RSApop.AssignBaseIncome();
	RSApop.AssignUrban();
	RSApop.SetEmployment();
	RSApop.AssignBehav();
	RSApop.AssignParents();
	RSApop.AssignBirth();
	RSApop.AssignContraception();
	RSApop.AssignPrison();
	RSApop.AssignIneqGender();
	RSApop.AssignAlcohol();
	RSApop.AssignInitHousehold();
	SetCalibParameters();
	if (SetInitPrev1990 == 0){
		RSApop.AssignHIV();}
	if(HSVind==1 || TPind==1 || HDind==1 || NGind==1 || CTind==1 || 
		TVind==1 || BVind==1 || VCind==1){
			RSApop.AssignSTIs();}
	CurrYear = StartYear;
	for(int ii=0; ii<ProjectionTerm; ii++){
		if (CurrYear == 1990 && SetInitPrev1990==1){
			RSApop.AssignHIV1990(0);
			RSApop.AssignHIV1990(1);
			RSApop.AssignHIV1990(2);
		}
		RSApop.OneYear();
	}
	//if (StructuralRCTcalib == 1){ RunStructRCTs(); }
	RSApop.UpdateAgeGroup();
}

void RunSimulations(int sims)
{
	int i, ic, ir, SimCount2, SimCount3;

	if (VaryParameters == 1 && FixedUncertainty == 1){
		ReadSTDparameters();
		if (HCTuncertainty == 1){ ReadHCTparameters(); }
	}

	if (sims > 0){ RestartFromStore(sims); }
	CurrSim = sims;
	for (i = sims; i<samplesize; i++){
		CurrSim += 1;
		if (i>0){
			Register.clear();
			Register.resize(InitPop);
			Register.reserve(MaxPop);
			HHregister.clear();
			if (StructuralRCTcalib == 1){
				TempRegister.clear();
				TempRegister.resize(InitPop);
				TempRegister.reserve(MaxPop);
				TempHHregister.clear();
			}
		}
		for (ir = 0; ir < 3; ir++){
			TotCurrFSW[ir] = 0;
			for (ic = 0; ic < MaxCSWs; ic++){
				RSApop.CSWregister[ic][ir] = 0;
			}
		}
		for (ic = 0; ic < MaxCasual; ic++){
			RSApop.CasualRegister[ic] = 0; }
		OneSimulation();
		cout << "Completed simulation " << CurrSim << endl;
		SimCount2 = CurrSim / IterationsPerPC;
		if (CurrSim == (SimCount2 * IterationsPerPC)){
			AggregateSims();
		}
		SimCount3 = CurrSim / 50;
		if (CurrSim == (SimCount3 * 50) && FixedUncertainty==0){
			StoreTempOutputs();
		}
		//cout << "TotalLogL: " << TotalLogL << endl;
	}
	StoreOutputs();
}

void StoreOutputs()
{
	if (FixedUncertainty == 0){
		HIVparamsLogL.RecordSample("ParamCombsLogL.txt");
		MSMTparamsLogL.RecordSample("MSMTparamsLogL.txt");
		if (StructuralRCTcalib == 1){ StructParamsLogL.RecordSample("StructParamsLogL.txt"); }
	}
	else{
		if (HIVind == 1){
			HIVprev15to49M.RecordSample("HIVprev15to49M.txt");
			HIVprev15to49F.RecordSample("HIVprev15to49F.txt");
			HIVprev15to49.RecordSample("HIVprev15to49.txt");
			HIVprev15to49B.RecordSample("HIVprev15to49B.txt");
			HIVprev15to49C.RecordSample("HIVprev15to49C.txt");
			HIVprev15to49W.RecordSample("HIVprev15to49W.txt");
			UrbanHIV.RecordSample("UrbanHIV.txt");
			RuralHIV.RecordSample("RuralHIV.txt");
			//HIVconcordance.RecordSample("HIVconcordance.txt");
			CD4500plus.RecordSample("CD4500plus.txt");
			CD4350to499.RecordSample("CD4350to499.txt");
			CD4200to349.RecordSample("CD4200to349.txt");
			CD4under200.RecordSample("CD4under200.txt");
			UndiagnosedAdults.RecordSample("UndiagnosedAdults.txt");
			DiagnosedUntreated.RecordSample("DiagnosedUntreated.txt");
			TreatedAdults.RecordSample("TreatedAdults.txt");
			AveHIVdurPreg.RecordSample("AveHIVdurPreg.txt");
			AveHIVdurPreg15.RecordSample("AveHIVdurPreg15.txt");
			AveHIVdurPreg20.RecordSample("AveHIVdurPreg20.txt");
			AveHIVdurPreg25.RecordSample("AveHIVdurPreg25.txt");
			AveHIVdurPreg30.RecordSample("AveHIVdurPreg30.txt");
			AveHIVdurPreg35.RecordSample("AveHIVdurPreg35.txt");
			PrevPregNoEdu.RecordSample("PrevPregNoEdu.txt");
			PrevPregPrimary.RecordSample("PrevPregPrimary.txt");
			PrevPreg2ndary.RecordSample("PrevPreg2ndary.txt");
			PrevPregTertiary.RecordSample("PrevPregTertiary.txt");
			PrevPregAllEdu.RecordSample("PrevPregAllEdu.txt");
			PrevPregRural.RecordSample("PrevPregRural.txt");
			PrevPregUrban.RecordSample("PrevPregUrban.txt");
			PrevPregB.RecordSample("PrevPregB.txt");
			PrevPregC.RecordSample("PrevPregC.txt");
			PrevPregW.RecordSample("PrevPregW.txt");
			PrevByEdu2002.RecordSample("PrevByEdu2002.txt");
			//DiscordantPropnM.RecordSample("DiscordantPropnM.txt");
			//DiscordantPropnF.RecordSample("DiscordantPropnF.txt");
			//DiscordantPropn.RecordSample("DiscordantPropn.txt");
			//DiscordantFtoM.RecordSample("DiscordantFtoM.txt");
			DiscordantMposLT.RecordSample("DiscordantMposLT.txt");
			DiscordantMposST.RecordSample("DiscordantMposST.txt");
			DiscordantFposLT.RecordSample("DiscordantFposLT.txt");
			DiscordantFposST.RecordSample("DiscordantFposST.txt");
			ConcordantPosLT.RecordSample("ConcordantPosLT.txt");
			ConcordantPosST.RecordSample("ConcordantPosST.txt");
			ConcordantNegLT.RecordSample("ConcordantNegLT.txt");
			ConcordantNegST.RecordSample("ConcordantNegST.txt");
			HIVinc15to49.RecordSample("HIVinc15to49.txt");
			HIVinc15to49_U.RecordSample("HIVinc15to49_U.txt");
			HIVinc15to49_R.RecordSample("HIVinc15to49_R.txt");
			if (HIVcalib == 1){
				PrevPreg15.RecordSample("PrevPreg15.txt");
				PrevPreg20.RecordSample("PrevPreg20.txt");
				PrevPreg25.RecordSample("PrevPreg25.txt");
				PrevPreg30.RecordSample("PrevPreg30.txt");
				PrevPreg35.RecordSample("PrevPreg35.txt");
				PrevPreg40.RecordSample("PrevPreg40.txt");
				PrevPregTotal.RecordSample("PrevPregTotal.txt");
				PrevFertile15.RecordSample("PrevFertile15.txt");
				PrevFertile20.RecordSample("PrevFertile20.txt");
				PrevFertile25.RecordSample("PrevFertile25.txt");
				PrevFertile30.RecordSample("PrevFertile30.txt");
				PrevFertile35.RecordSample("PrevFertile35.txt");
				PrevHH2005.RecordSample("PrevHH2005.txt");
				PrevHH2008.RecordSample("PrevHH2008.txt");
				PrevHH2012.RecordSample("PrevHH2012.txt");
				HIVparamsLogL.RecordSample("ParamCombsLogL.txt");
				NewHIV.RecordSample("NewHIV.txt");
				NewHIVexp.RecordSample("NewHIVexp.txt");
				HIVprevCSW.RecordSample("HIVprevCSW.txt");
			}
		}
		if (HSVind == 1){
			HSVprev15to49M.RecordSample("HSVprev15to49M.txt");
			HSVprev15to49F.RecordSample("HSVprev15to49F.txt");
			HSVprev15to49B.RecordSample("HSVprev15to49B.txt");
			HSVprev15to49C.RecordSample("HSVprev15to49C.txt");
			HSVprev15to49W.RecordSample("HSVprev15to49W.txt");
			//HSVconcordance.RecordSample("HSVconcordance.txt");
			if (HSVcalib == 1){
				HSVparamsLogL.RecordSample("HSVparamsLogL.txt");
				NewHSV.RecordSample("NewHSV.txt");
				HSVprevCSW.RecordSample("HSVprevCSW.txt");
				HSVprevANC.RecordSample("HSVprevANC.txt");
			}
		}
		if (HDind == 1){
			HDprev15to49M.RecordSample("HDprev15to49M.txt");
			HDprev15to49F.RecordSample("HDprev15to49F.txt");
			//HDconcordance.RecordSample("HDconcordance.txt");
		}
		if (TPind == 1){
			TPprev15to49M.RecordSample("TPprev15to49M.txt");
			TPprev15to49F.RecordSample("TPprev15to49F.txt");
			//TPconcordance.RecordSample("TPconcordance.txt");
			if (TPcalib == 1){
				TPparamsLogL.RecordSample("TPparamsLogL.txt");
				NewTP.RecordSample("NewTP.txt");
				TPprevCSW.RecordSample("TPprevCSW.txt");
				TPprevANC.RecordSample("TPprevANC.txt");
			}
		}
		if (NGind == 1){
			NGprev15to49M.RecordSample("NGprev15to49M.txt");
			NGprev15to49F.RecordSample("NGprev15to49F.txt");
			//NGconcordance.RecordSample("NGconcordance.txt");
			if (NGcalib == 1){
				NGparamsLogL.RecordSample("NGparamsLogL.txt");
				NewNG.RecordSample("NewNG.txt");
				NGprevCSW.RecordSample("NGprevCSW.txt");
			}
		}
		if (CTind == 1){
			CTprev15to49M.RecordSample("CTprev15to49M.txt");
			CTprev15to49F.RecordSample("CTprev15to49F.txt");
			//CTconcordance.RecordSample("CTconcordance.txt");
			if (CTcalib == 1){
				CTparamsLogL.RecordSample("CTparamsLogL.txt");
				NewCT.RecordSample("NewCT.txt");
				CTprevCSW.RecordSample("CTprevCSW.txt");
			}
		}
		if (TVind == 1){
			TVprev15to49M.RecordSample("TVprev15to49M.txt");
			TVprev15to49F.RecordSample("TVprev15to49F.txt");
			//TVconcordance.RecordSample("TVconcordance.txt");
			if (TVcalib == 1){
				TVparamsLogL.RecordSample("TVparamsLogL.txt");
				NewTV.RecordSample("NewTV.txt");
				TVprevCSW.RecordSample("TVprevCSW.txt");
			}
		}
		if (CTind == 1 && NGind == 1 && TVind == 1){
			NewCTNGTV.RecordSample("NewCTNGTV.txt");
			NewCTNGTV_M.RecordSample("NewCTNGTV_M.txt");
			NewCTNGTV_F.RecordSample("NewCTNGTV_F.txt");
		}
		if (BVind == 1){
			BVprev15to49F.RecordSample("BVprev15to49F.txt");
		}
		if (VCind == 1){
			VCprev15to49F.RecordSample("VCprev15to49F.txt");
		}
		if (MSMcalib == 1 || MSMcalibHIV == 1){
			if (MSMcalib == 1){ MSMparamsLogL.RecordSample("MSMparamsLogL.txt"); }
			if (MSMcalibHIV == 1){ MSMTparamsLogL.RecordSample("MSMTparamsLogL.txt"); }
			MSMmarried.RecordSample("MSMmarried.txt");
			MSMrecentBi.RecordSample("MSMrecentBi.txt");
			MSMeverBi.RecordSample("MSMeverBi.txt");
			MSMaged25plus.RecordSample("MSMaged25plus.txt");
			MSMpositive.RecordSample("MSMpositive.txt");
			MSMcurrentReg.RecordSample("MSMcurrentReg.txt");
			MSMmultPartners.RecordSample("MSMmultPartners.txt");
			EverMSMpropn.RecordSample("EverMSMpropn.txt");
			EverMSMcurrentReg.RecordSample("EverMSMcurrentReg.txt");
			EverMSMeverBi.RecordSample("EverMSMeverBi.txt");
			EverMSMpositive.RecordSample("EverMSMpositive.txt");
			NeverMSMpositive.RecordSample("NeverMSMpositive.txt");
			HIVageProfileMSM.RecordSample("HIVageProfileMSM.txt");
			CasualAgeProfile.RecordSample("CasualAgeProfile.txt");
			BiAgeProfile.RecordSample("BiAgeProfile.txt");
			HighRiskMSM.RecordSample("HighRiskMSM.txt");
			NewHIV_MSM.RecordSample("NewHIV_MSM.txt");
			NewHIV.RecordSample("NewHIV.txt");
			MSMpositiveBi.RecordSample("MSMpositiveBi.txt");
			MSMpositiveGay.RecordSample("MSMpositiveGay.txt");
			MSMposInsertive.RecordSample("MSMposInsertive.txt");
			MSMposReceptive.RecordSample("MSMposReceptive.txt");
			MSMposVersatile.RecordSample("MSMposVersatile.txt");
			HomoPrefPos.RecordSample("HomoPrefPos.txt");
			BiPrefPos.RecordSample("BiPrefPos.txt");
			HeteroPrefPos.RecordSample("HeteroPrefPos.txt");
		}

		// Educational outputs
		/*CompletedMatricAM.RecordSample("CompletedMatricAM.txt");
		CompletedMatricAF.RecordSample("CompletedMatricAF.txt");
		CompletedMatricCM.RecordSample("CompletedMatricCM.txt");
		CompletedMatricCF.RecordSample("CompletedMatricCF.txt");
		CompletedMatricWM.RecordSample("CompletedMatricWM.txt");
		CompletedMatricWF.RecordSample("CompletedMatricWF.txt");
		TertiaryEnrolRatioA.RecordSample("TertiaryEnrolRatioA.txt");
		TertiaryEnrolRatioC.RecordSample("TertiaryEnrolRatioC.txt");
		TertiaryEnrolRatioW.RecordSample("TertiaryEnrolRatioW.txt");
		AttainmentByAge2007.RecordSample("AttainmentByAge2007.txt");
		AttainmentByGrade2007.RecordSample("AttainmentByGrade2007.txt");
		Enrolled2001A.RecordSample("Enrolled2001A.txt");
		Enrolled2001C.RecordSample("Enrolled2001C.txt");
		Enrolled2001W.RecordSample("Enrolled2001W.txt");
		EduByParentEdu.RecordSample("EduByParentEdu.txt");*/
		DropoutDuePregnancy.RecordSample("DropoutDuePregnancy.txt");

		// Sexual behaviour outputs
		SexuallyExperienced.RecordSample("SexuallyExperienced.txt");
		SexuallyExpAdol.RecordSample("SexuallyExpAdol.txt");
		EduConcordance.RecordSample("EduConcordance.txt");
		RaceConcordance.RecordSample("RaceConcordance.txt");
		ConcurrencyPrev.RecordSample("ConcurrencyPrev.txt");
		ConcurrencyTot.RecordSample("ConcurrencyTot.txt");
		Married2001M.RecordSample("Married2001M.txt");
		Married2001F.RecordSample("Married2001F.txt");
		Married2011M.RecordSample("Married2011M.txt");
		Married2011F.RecordSample("Married2011F.txt");
		MarriedEdu2011.RecordSample("MarriedEdu2011.txt");
		MarriedYouth.RecordSample("MarriedYouth.txt");
		CondomUseByEdu.RecordSample("CondomUseByEdu.txt");
		CondomUseByRace.RecordSample("CondomUseByRace.txt");
		CondomUse15to24.RecordSample("CondomUse15to24.txt");
		CondomUse25to49.RecordSample("CondomUse25to49.txt");
		CondomUse15to49.RecordSample("CondomUse15to49.txt");
		SWsexActs.RecordSample("SWsexActs.txt");
		SWsexActsProt.RecordSample("SWsexActsProt.txt");
		CohabitTemp.RecordSample("CohabitTemp.txt");
		CohabitTempTrend.RecordSample("CohabitTempTrend.txt");
		AgeDisparate2000.RecordSample("AgeDisparate2000.txt");
		PartnerAgeDifM.RecordSample("PartnerAgeDifM.txt");
		PartnerAgeDifF.RecordSample("PartnerAgeDifF.txt");
		CasualCalib2005.RecordSample("CasualCalib2005.txt");
		CasualPrevM.RecordSample("CasualPrevM.txt");
		CasualPrevF.RecordSample("CasualPrevF.txt");

		// Intervention outputs
		TakingPrEP.RecordSample("TakingPrEP.txt");
		VaccineDoses.RecordSample("VaccineDoses.txt");
		NewHIV.RecordSample("NewHIV.txt");
		NewHIVexp.RecordSample("NewHIVexp.txt");
		NewHIVexpM.RecordSample("NewHIVexpM.txt");
		NewHIVexpF.RecordSample("NewHIVexpF.txt");

		// Demographic outputs
		TotBirths.RecordSample("TotBirths.txt");
		TeenBirths.RecordSample("TeenBirths.txt");
		TotPop.RecordSample("TotPop.txt");
		ContrPrev.RecordSample("ContrPrev.txt");
		EverUseContr98.RecordSample("EverUseContr98.txt");
		EverUseContr03.RecordSample("EverUseContr03.txt");
		CurrUseContr98.RecordSample("CurrUseContr98.txt");
		CurrUseContr03.RecordSample("CurrUseContr03.txt");
		CurrUseContr16.RecordSample("CurrUseContr16.txt");
		CurrContrRace98.RecordSample("CurrContrRace98.txt");
		CurrContrRace03.RecordSample("CurrContrRace03.txt");
		CurrContrEdu98.RecordSample("CurrContrEdu98.txt");
		CurrContrEdu03.RecordSample("CurrContrEdu03.txt");
		CurrInjectable.RecordSample("CurrInjectable.txt");
		CurrPill.RecordSample("CurrPill.txt");
		YouthHormonal.RecordSample("YouthHormonal.txt");
		BirthIntervals.RecordSample("BirthIntervals.txt");
		SexuallyExpFertB.RecordSample("SexuallyExpFertB.txt");
		SexuallyExpFertC.RecordSample("SexuallyExpFertC.txt");
		SexuallyExpFertW.RecordSample("SexuallyExpFertW.txt");
		FertByMarriage.RecordSample("FertByMarriage.txt");
		UrbanTrend.RecordSample("UrbanTrend.txt");
		UrbanAge.RecordSample("UrbanAge.txt");
		UrbanAge1996.RecordSample("UrbanAge1996.txt");
		UrbanRace.RecordSample("UrbanRace.txt");
		UrbanRace1996.RecordSample("UrbanRace1996.txt");
		UrbanEdu.RecordSample("UrbanEdu.txt");
		HHmembersByAge.RecordSample("HHmembersByAge.txt");

		// Prison outputs
		FractionInPrison.RecordSample("FractionInPrison.txt");
		PrisonYouth.RecordSample("PrisonYouth.txt");
		PrisonCompletedSchool.RecordSample("PrisonCompletedSchool.txt");
		PrisonHIV.RecordSample("PrisonHIV.txt");
		PrisonPrevious.RecordSample("PrisonPrevious.txt");
		PrisonSentenced.RecordSample("PrisonSentenced.txt");
		PrisonRace.RecordSample("PrisonRace.txt");

		// HCT outputs
		TotalTestsOI.RecordSample("TotalTestsOI.txt");
		PosTestsOI.RecordSample("PosTestsOI.txt");
		NewDiagOI.RecordSample("NewDiagOI.txt");
		TotalTestsPrEP.RecordSample("TotalTestsPrEP.txt");
		PosTestsPrEP.RecordSample("PosTestsPrEP.txt");
		TotalTestsANC.RecordSample("TotalTestsANC.txt");
		PosTestsANC.RecordSample("PosTestsANC.txt");
		NewDiagANC.RecordSample("NewDiagANC.txt");
		TotalTestsGen.RecordSample("TotalTestsGen.txt");
		PosTestsGen.RecordSample("PosTestsGen.txt");
		NewDiagGen.RecordSample("NewDiagGen.txt");
		TotalTestsMMC.RecordSample("TotalTestsMMC.txt");
		PosTestsMMC.RecordSample("PosTestsMMC.txt");
		TotalTestsPartner.RecordSample("TotalTestsPartner.txt");
		PosTestsPartner.RecordSample("PosTestsPartner.txt");
		NewDiagPartner.RecordSample("NewDiagPartner.txt");
		TotalTestsHH_U.RecordSample("TotalTestsHH_U.txt");
		PosTestsHH_U.RecordSample("PosTestsHH_U.txt");
		NewDiagHH_U.RecordSample("NewDiagHH_U.txt");
		TotalTestsHH_R.RecordSample("TotalTestsHH_R.txt");
		PosTestsHH_R.RecordSample("PosTestsHH_R.txt");
		NewDiagHH_R.RecordSample("NewDiagHH_R.txt");
		TotalTestsMobile_U.RecordSample("TotalTestsMobileU.txt");
		PosTestsMobile_U.RecordSample("PosTestsMobileU.txt");
		NewDiagMobile_U.RecordSample("NewDiagMobileU.txt");
		TotalTestsMobile_R.RecordSample("TotalTestsMobileR.txt");
		PosTestsMobile_R.RecordSample("PosTestsMobileR.txt");
		NewDiagMobile_R.RecordSample("NewDiagMobileR.txt");
		TotalTestsFSW.RecordSample("TotalTestsFSW.txt");
		PosTestsFSW.RecordSample("PosTestsFSW.txt");
		NewDiagFSW.RecordSample("NewDiagFSW.txt");
		TotalTestsMSM.RecordSample("TotalTestsMSM.txt");
		PosTestsMSM.RecordSample("PosTestsMSM.txt");
		NewDiagMSM.RecordSample("NewDiagMSM.txt");
		TotalTestsSchool.RecordSample("TotalTestsSchool.txt");
		PosTestsSchool.RecordSample("PosTestsSchool.txt");
		NewDiagSchool.RecordSample("NewDiagSchool.txt");
		TotalTestsANCpartner0.RecordSample("TotalTestsANCpartner0.txt");
		PosTestsANCpartner0.RecordSample("PosTestsANCpartner0.txt");
		NewDiagANCpartner0.RecordSample("NewDiagANCpartner0.txt");
		TotalTestsANCpartner1.RecordSample("TotalTestsANCpartner1.txt");
		PosTestsANCpartner1.RecordSample("PosTestsANCpartner1.txt");
		NewDiagANCpartner1.RecordSample("NewDiagANCpartner1.txt");
		TotalTestsPrison.RecordSample("TotalTestsPrison.txt");
		PosTestsPrison.RecordSample("PosTestsPrison.txt");
		NewDiagPrison.RecordSample("NewDiagPrison.txt");
		TotalTestsSTI.RecordSample("TotalTestsSTI.txt");
		PosTestsSTI.RecordSample("PosTestsSTI.txt");
		NewDiagSTI.RecordSample("NewDiagSTI.txt");
		TotalTestsWork.RecordSample("TotalTestsWork.txt");
		PosTestsWork.RecordSample("PosTestsWork.txt");
		NewDiagWork.RecordSample("NewDiagWork.txt");
		TotalTestsFPC.RecordSample("TotalTestsFPC.txt");
		PosTestsFPC.RecordSample("PosTestsFPC.txt");
		NewDiagFPC.RecordSample("NewDiagFPC.txt");
		AcuteTests.RecordSample("AcuteTests.txt");
		MaleDiagnosed.RecordSample("MaleDiagnosed.txt");
		FemDiagnosed.RecordSample("FemDiagnosed.txt");
		Diagnosed15to24.RecordSample("Diagnosed15to24.txt");
		Diagnosed25to49.RecordSample("Diagnosed25to49.txt");
		Diagnosed50plus.RecordSample("Diagnosed50plus.txt");
		FSWdiagnosed.RecordSample("FSWdiagnosed.txt");
		MSMdiagnosed.RecordSample("MSMdiagnosed.txt");
		LYsLost.RecordSample("LYsLost.txt");
		LYsLostExp.RecordSample("LYsLostExp.txt");
		NewART200.RecordSample("NewART200.txt");
		NewART350.RecordSample("NewART350.txt");
		NewART500.RecordSample("NewART500.txt");
		NewART500plus.RecordSample("NewART500plus.txt");
		MMCoperations.RecordSample("MMCoperations.txt");
		DiagNoART200.RecordSample("DiagNoART200.txt");
		DiagNoART350.RecordSample("DiagNoART350.txt");
		DiagNoART500.RecordSample("DiagNoART500.txt");
		DiagNoART500plus.RecordSample("DiagNoART500plus.txt");
		ProtSexActs.RecordSample("ProtSexActs.txt");
		TreatedAdults200.RecordSample("TreatedAdults200.txt");
		TreatedAdults350.RecordSample("TreatedAdults350.txt");
		TreatedAdults500.RecordSample("TreatedAdults500.txt");
		TreatedAdults500plus.RecordSample("TreatedAdults500plus.txt");
		TotBirthsHIV.RecordSample("TotBirthsHIV.txt");
		TotBirthsART.RecordSample("TotBirthsART.txt");
		AIDSdeaths.RecordSample("AIDSdeaths.txt");
		NonAIDSdeaths.RecordSample("NonAIDSdeaths.txt");
		TotalCosts.RecordSample("TotalCosts.txt");
		TotalTestCosts.RecordSample("TotalTestCosts.txt");
		//StoreHCToutputs("HCToutputs.txt");
		StoreGenOutputs("GenOutputs.txt"); 
		if (StructuralDriverCalib == 1){ StoreStructuralOutputs("StructuralOutputs.txt"); }
	}
	if (StructuralRCTcalib == 1){
		//SingleSessionAlcohol.RecordSample("SingleSessionAlcohol.txt");
		//MultiSessionAlcohol.RecordSample("MultiSessionAlcohol.txt");
		CashTransferOut.RecordSample("CashTransferOut.txt");
		SchoolSupportOut.RecordSample("SchoolSupportOut.txt");
		VocationalTrainingOut.RecordSample("VocationalTrainingOut.txt");
		if (FixedUncertainty == 1){ MeanHHincomePC.RecordSample("MeanHHincomePC.txt"); }
		//GenderCommun.RecordSample("GenderCommun.txt");
		//GenderIndiv.RecordSample("GenderIndiv.txt");
	}
}

void StoreTempOutputs()
{
	if (HIVcalib == 1){ HIVparamsLogL.RecordSample("ParamCombsLogL.txt"); }
	if (MSMcalibHIV == 1){ MSMTparamsLogL.RecordSample("MSMTparamsLogL.txt"); }
	if (StructuralRCTcalib == 1){ StructParamsLogL.RecordSample("StructParamsLogL.txt"); }
}

void RestartFromStore(int sims)
{
	int cols, idum, ir, ic;
	ifstream file1;
	ifstream file2;
	ifstream file3;

	if (HIVcalib == 1){ 
		cols = HIVparamsLogL.columns;
		file1.open("ParamCombsLogL.txt");
		if (file1.fail()) {
			cerr << "Could not open ParamCombsLogL.txt\n";
			exit(1);
		}
		for (ir = 0; ir<sims; ir++){
			file1 >> idum >> idum;
			for (ic = 0; ic <= cols; ic++){
				file1 >> HIVparamsLogL.out[ir][ic];
			}
		}
		file1.close();
	}
	if (MSMcalibHIV == 1){
		cols = MSMTparamsLogL.columns;
		file2.open("MSMTparamsLogL.txt");
		if (file2.fail()) {
			cerr << "Could not open MSMTparamsLogL.txt\n";
			exit(1);
		}
		for (ir = 0; ir<sims; ir++){
			file2 >> idum >> idum;
			for (ic = 0; ic <= cols; ic++){
				file2 >> MSMTparamsLogL.out[ir][ic];
			}
		}
		file2.close();
	}
	if (StructuralRCTcalib == 1){ 
		cols = StructParamsLogL.columns;
		file3.open("StructParamsLogL.txt");
		if (file3.fail()) {
			cerr << "Could not open StructParamsLogL.txt\n";
			exit(1);
		}
		for (ir = 0; ir<sims; ir++){
			file3 >> idum >> idum >> idum;
			for (ic = 0; ic <= cols; ic++){
				file3 >> StructParamsLogL.out[ir][ic];
			}
		}
		file3.close();
	}
}

void StoreHCToutputs(const char* filout)
{
	int ic, iy;
	double temp1, temp2;

	for (ic = 0; ic < samplesize; ic++){
		for (iy = 15; iy < CurrYear - StartYear; iy++){
			TotalTests.out[ic][iy] = TotalTestsGen.out[ic][iy] + TotalTestsANC.out[ic][iy] +
				TotalTestsOI.out[ic][iy] + TotalTestsPrison.out[ic][iy] + TotalTestsPartner.out[ic][iy] +
				TotalTestsMMC.out[ic][iy] + TotalTestsPrEP.out[ic][iy] + TotalTestsSTI.out[ic][iy] + 
				TotalTestsHH_U.out[ic][iy] + TotalTestsHH_R.out[ic][iy] + TotalTestsMobile_U.out[ic][iy] + 
				TotalTestsMobile_R.out[ic][iy] + TotalTestsMSM.out[ic][iy] + TotalTestsFSW.out[ic][iy] +
				TotalTestsANCpartner0.out[ic][iy] + TotalTestsANCpartner1.out[ic][iy]  +
				TotalTestsSchool.out[ic][iy] + TotalTestsWork.out[ic][iy] + TotalTestsFPC.out[ic][iy];
			NewDiagnoses.out[ic][iy] = NewDiagGen.out[ic][iy] + NewDiagANC.out[ic][iy] +
				NewDiagOI.out[ic][iy] + NewDiagPrison.out[ic][iy] + NewDiagPartner.out[ic][iy] +
				PosTestsMMC.out[ic][iy] + PosTestsPrEP.out[ic][iy] + NewDiagSTI.out[ic][iy] +
				NewDiagHH_U.out[ic][iy] + NewDiagHH_R.out[ic][iy] + NewDiagMobile_U.out[ic][iy] +
				NewDiagMobile_R.out[ic][iy] + NewDiagMSM.out[ic][iy] + NewDiagFSW.out[ic][iy] +
				NewDiagANCpartner0.out[ic][iy] + NewDiagANCpartner1.out[ic][iy] +
				NewDiagSchool.out[ic][iy] + NewDiagWork.out[ic][iy] + NewDiagFPC.out[ic][iy];
		}
	}

	YieldGen.GetSummHCTout2(1, &NewDiagGen, &TotalTestsGen);
	YieldANC.GetSummHCTout2(89, &NewDiagANC, &TotalTestsANC);
	YieldOI.GetSummHCTout2(3, &NewDiagOI, &TotalTestsOI);
	YieldPrison.GetSummHCTout2(4, &NewDiagPrison, &TotalTestsPrison);
	YieldPartner.GetSummHCTout2(5, &NewDiagPartner, &TotalTestsPartner);
	YieldMMC.GetSummHCTout2(6, &PosTestsMMC, &TotalTestsMMC);
	YieldPrEP.GetSummHCTout2(7, &PosTestsPrEP, &TotalTestsPrEP);
	YieldSTI.GetSummHCTout2(8, &NewDiagSTI, &TotalTestsSTI);
	YieldHH_U.GetSummHCTout2(9, &NewDiagHH_U, &TotalTestsHH_U);
	YieldHH_R.GetSummHCTout2(10, &NewDiagHH_R, &TotalTestsHH_R);
	YieldMobile_U.GetSummHCTout2(11, &NewDiagMobile_U, &TotalTestsMobile_U);
	YieldMobile_R.GetSummHCTout2(12, &NewDiagMobile_R, &TotalTestsMobile_R);
	YieldMSM.GetSummHCTout2(13, &NewDiagMSM, &TotalTestsMSM);
	YieldFSW.GetSummHCTout2(14, &NewDiagFSW, &TotalTestsFSW);
	YieldSchool.GetSummHCTout2(15, &NewDiagSchool, &TotalTestsSchool);
	YieldANCpartner0.GetSummHCTout2(16, &NewDiagANCpartner0, &TotalTestsANCpartner0);
	YieldANCpartner1.GetSummHCTout2(17, &NewDiagANCpartner1, &TotalTestsANCpartner1);
	YieldWork.GetSummHCTout2(18, &NewDiagWork, &TotalTestsWork);
	YieldFPC.GetSummHCTout2(19, &NewDiagFPC, &TotalTestsFPC);

	TotalTestsGen.GetSummHCTout1(22);
	TotalTestsANC.GetSummHCTout1(90);
	TotalTestsOI.GetSummHCTout1(24);
	TotalTestsPrison.GetSummHCTout1(25);
	TotalTestsPartner.GetSummHCTout1(26);
	TotalTestsMMC.GetSummHCTout1(27);
	TotalTestsPrEP.GetSummHCTout1(28);
	TotalTestsSTI.GetSummHCTout1(29);
	TotalTestsHH_U.GetSummHCTout1(30);
	TotalTestsHH_R.GetSummHCTout1(31);
	TotalTestsMobile_U.GetSummHCTout1(32);
	TotalTestsMobile_R.GetSummHCTout1(33);
	TotalTestsMSM.GetSummHCTout1(34);
	TotalTestsFSW.GetSummHCTout1(35);
	TotalTestsSchool.GetSummHCTout1(36);
	TotalTestsANCpartner0.GetSummHCTout1(37);
	TotalTestsANCpartner1.GetSummHCTout1(38);
	TotalTestsWork.GetSummHCTout1(39);
	TotalTestsFPC.GetSummHCTout1(40);
	TotalTests.GetSummHCTout1(41);

	PosTestsGen.GetSummHCTout1(46);
	PosTestsANC.GetSummHCTout1(47);
	PosTestsOI.GetSummHCTout1(48);
	PosTestsPrison.GetSummHCTout1(49);
	PosTestsPartner.GetSummHCTout1(50);
	PosTestsMMC.GetSummHCTout1(51);
	PosTestsPrEP.GetSummHCTout1(52);
	PosTestsSTI.GetSummHCTout1(53);
	PosTestsHH_U.GetSummHCTout1(54);
	PosTestsHH_R.GetSummHCTout1(55);
	PosTestsMobile_U.GetSummHCTout1(56);
	PosTestsMobile_R.GetSummHCTout1(57);
	PosTestsMSM.GetSummHCTout1(58);
	PosTestsFSW.GetSummHCTout1(59);
	PosTestsSchool.GetSummHCTout1(60);
	PosTestsANCpartner0.GetSummHCTout1(61);
	PosTestsANCpartner1.GetSummHCTout1(62);
	PosTestsWork.GetSummHCTout1(63);
	PosTestsFPC.GetSummHCTout1(64);

	NewDiagGen.GetSummHCTout1(67);
	NewDiagANC.GetSummHCTout1(68);
	NewDiagOI.GetSummHCTout1(69);
	NewDiagPrison.GetSummHCTout1(70);
	NewDiagPartner.GetSummHCTout1(71);
	for (iy = 0; iy < 56; iy++){
		SummOut[71][iy] = SummOut[50][iy];
		SummOut[72][iy] = SummOut[51][iy];
	}
	NewDiagSTI.GetSummHCTout1(74);
	NewDiagHH_U.GetSummHCTout1(75);
	NewDiagHH_R.GetSummHCTout1(76);
	NewDiagMobile_U.GetSummHCTout1(77);
	NewDiagMobile_R.GetSummHCTout1(78);
	NewDiagMSM.GetSummHCTout1(79);
	NewDiagFSW.GetSummHCTout1(80);
	NewDiagSchool.GetSummHCTout1(81);
	NewDiagANCpartner0.GetSummHCTout1(82);
	NewDiagANCpartner1.GetSummHCTout1(83);
	NewDiagWork.GetSummHCTout1(84);
	NewDiagFPC.GetSummHCTout1(85);
	NewDiagnoses.GetSummHCTout1(86);

	// Calculate additional outputs for the ANC retest scenarios
	for (iy = 0; iy < 40; iy++){
		SummOut[22][iy] = SummOut[89][iy] + (SummOut[89][iy] - SummOut[46][iy]) * ANCretestFreq[iy + 15];
		SummOut[90][iy] = SummOut[46][iy]; // Positives
		SummOut[91][iy] = SummOut[67][iy]; // New diagnoses
		SummOut[1][iy] = SummOut[91][iy] / SummOut[22][iy];
	}
	SummOut[90][41] = SummOut[46][41];
	SummOut[90][42] = SummOut[46][42];
	SummOut[91][41] = SummOut[67][41];
	SummOut[91][42] = SummOut[67][42];
	temp1 = 0.0;
	temp2 = 0.0;
	for (iy = 19; iy < 39; iy++){
		temp1 += SummOut[91][iy];
		temp2 += SummOut[22][iy];
	}
	SummOut[22][41] = temp2;
	SummOut[1][41] = 1.0 * temp1 / temp2;
	SummOut[22][42] = SummOut[89][42] * SummOut[22][41] / SummOut[89][41];
	SummOut[1][42] = SummOut[88][42] * SummOut[1][41] / SummOut[88][41];

	NewHIV.GetSummHCTout1(94);
	LYsLost.GetSummHCTout1(95);
	TreatedAdults.GetSummHCTout1(96);
	TreatedAdults200.GetSummHCTout1(97);
	TreatedAdults350.GetSummHCTout1(98);
	TreatedAdults500.GetSummHCTout1(99);
	TreatedAdults500plus.GetSummHCTout1(100);
	NewART200.GetSummHCTout1(101);
	NewART350.GetSummHCTout1(102);
	NewART500.GetSummHCTout1(103);
	NewART500plus.GetSummHCTout1(104);
	MMCoperations.GetSummHCTout1(105);
	DiagNoART200.GetSummHCTout1(106);
	DiagNoART350.GetSummHCTout1(107);
	DiagNoART500.GetSummHCTout1(108);
	DiagNoART500plus.GetSummHCTout1(109);
	TakingPrEP.GetSummHCTout1(110);
	ProtSexActs.GetSummHCTout1(111);
	TotBirthsHIV.GetSummHCTout1(112);
	TotBirthsART.GetSummHCTout1(113);
	AIDSdeaths.GetSummHCTout1(114);
	NonAIDSdeaths.GetSummHCTout1(115);
	TotalCosts.GetSummHCTout1(117);

	// Calculate increases relative to baseline
	NewHIV_B.ReadBaseline("NewHIV_B.txt");
	LYsLostB.ReadBaseline("LYsLostB.txt");
	NewDiagnosesB.ReadBaseline("NewDiagnosesB.txt");
	TotalTestsB.ReadBaseline("TotalTestsB.txt");
	TotalCostsB.ReadBaseline("TotalCostsB.txt");
	NewHIV.CalcPeriodIncrease(2019, 2038, &NewHIV_B);
	LYsLost.CalcPeriodIncrease(2019, 2038, &LYsLostB);
	NewDiagnoses.CalcPeriodIncrease(2019, 2038, &NewDiagnosesB);
	TotalTests.CalcPeriodIncrease(2019, 2038, &TotalTestsB);
	TotalCosts.CalcPeriodIncrease(2019, 2038, &TotalCostsB);
	SummOut[93][43] = NewHIV.PeriodIncrease;
	SummOut[93][44] = NewHIV.PeriodIncreaseSE;
	SummOut[94][43] = LYsLost.PeriodIncrease;
	SummOut[94][44] = LYsLost.PeriodIncreaseSE;
	SummOut[85][43] = NewDiagnoses.PeriodIncrease;
	SummOut[85][44] = NewDiagnoses.PeriodIncreaseSE;
	SummOut[40][43] = TotalTests.PeriodIncrease;
	SummOut[40][44] = TotalTests.PeriodIncreaseSE;
	SummOut[116][43] = TotalCosts.PeriodIncrease;
	SummOut[116][44] = TotalCosts.PeriodIncreaseSE;

	int i, c;
	ostringstream s;

	if (process_num > 0){
		s << process_num << filout;
	}
	else{
		s << filout;
	}

	//ofstream file(filout);
	ofstream file(s.str().c_str()); // Converts s to a C string

	for (i = 0; i<120; i++){
		file << setw(6) << right;
		for (c = 0; c<45; c++){
			file << "	" << setw(10) << right << SummOut[i][c];
		}
		file << endl;
	}
	file.close();
}

void StoreGenOutputs(const char* filout)
{
	int ic, iy;
	double temp1, temp2;

	// HIV incidence outputs
	SummOutRow = 0;
	NewHIV.GetSummGenOut(0);
	NewHIV_U.GetSummGenOut(0);
	NewHIV_R.GetSummGenOut(0);
	NewHIV_MSM.GetSummGenOut(0);
	HIVinc15to49.GetSummGenOut(0);
	HIVinc15to49_U.GetSummGenOut(0);
	HIVinc15to49_R.GetSummGenOut(0);
	
	// HIV prevalence outputs
	SummOutRow += 3;
	HIVprev15to49M.GetSummGenOut(0);
	HIVprev15to49F.GetSummGenOut(0);
	HIVprevCSW.GetSummGenOut(0);
	MSMprev.GetSummGenOut(0);
	PrevPreg15.GetSummGenOut(5);
	PrevPreg20.GetSummGenOut(5);
	PrevPreg25.GetSummGenOut(5);
	PrevPreg30.GetSummGenOut(5);
	PrevPreg35.GetSummGenOut(5);
	PrevPregTotal.GetSummGenOut(5);
	HIVprev15to49B.GetSummGenOut(0);
	HIVprev15to49C.GetSummGenOut(0);
	HIVprev15to49W.GetSummGenOut(0);
	UrbanHIV.GetSummGenOut(0);
	RuralHIV.GetSummGenOut(0);
	PrisonHIV.GetSummGenOut(0);

	// AIDS mortality and other mortality outputs
	SummOutRow += 3;
	AIDSdeaths.GetSummGenOut(0);
	ARTdeathsExp.GetSummGenOut(0);
	NonAIDSdeaths.GetSummGenOut(0);
	LYsLost.GetSummGenOut(0);

	// HIV testing outputs
	SummOutRow += 3;
	for (ic = 0; ic < samplesize; ic++){
		for (iy = 0; iy < CurrYear - StartYear; iy++){
			TotalTests.out[ic][iy] = TotalTestsGen.out[ic][iy] + TotalTestsANC.out[ic][iy] +
				TotalTestsOI.out[ic][iy] + TotalTestsPrison.out[ic][iy] + TotalTestsPartner.out[ic][iy] +
				TotalTestsMMC.out[ic][iy] + TotalTestsPrEP.out[ic][iy] + TotalTestsSTI.out[ic][iy] +
				TotalTestsHH_U.out[ic][iy] + TotalTestsHH_R.out[ic][iy] + TotalTestsMobile_U.out[ic][iy] +
				TotalTestsMobile_R.out[ic][iy] + TotalTestsMSM.out[ic][iy] + TotalTestsFSW.out[ic][iy] +
				TotalTestsANCpartner0.out[ic][iy] + TotalTestsANCpartner1.out[ic][iy] +
				TotalTestsSchool.out[ic][iy] + TotalTestsWork.out[ic][iy] + TotalTestsFPC.out[ic][iy];
			TotPosTests.out[ic][iy] = PosTestsGen.out[ic][iy] + PosTestsANC.out[ic][iy] +
				PosTestsOI.out[ic][iy] + PosTestsPrison.out[ic][iy] + PosTestsPartner.out[ic][iy] +
				PosTestsMMC.out[ic][iy] + PosTestsPrEP.out[ic][iy] + PosTestsSTI.out[ic][iy] +
				PosTestsHH_U.out[ic][iy] + PosTestsHH_R.out[ic][iy] + PosTestsMobile_U.out[ic][iy] +
				PosTestsMobile_R.out[ic][iy] + PosTestsMSM.out[ic][iy] + PosTestsFSW.out[ic][iy] +
				PosTestsANCpartner0.out[ic][iy] + PosTestsANCpartner1.out[ic][iy] +
				PosTestsSchool.out[ic][iy] + PosTestsWork.out[ic][iy] + PosTestsFPC.out[ic][iy];
			NewDiagnoses.out[ic][iy] = NewDiagGen.out[ic][iy] + NewDiagANC.out[ic][iy] +
				NewDiagOI.out[ic][iy] + NewDiagPrison.out[ic][iy] + NewDiagPartner.out[ic][iy] +
				PosTestsMMC.out[ic][iy] + PosTestsPrEP.out[ic][iy] + NewDiagSTI.out[ic][iy] +
				NewDiagHH_U.out[ic][iy] + NewDiagHH_R.out[ic][iy] + NewDiagMobile_U.out[ic][iy] +
				NewDiagMobile_R.out[ic][iy] + NewDiagMSM.out[ic][iy] + NewDiagFSW.out[ic][iy] +
				NewDiagANCpartner0.out[ic][iy] + NewDiagANCpartner1.out[ic][iy] +
				NewDiagSchool.out[ic][iy] + NewDiagWork.out[ic][iy] + NewDiagFPC.out[ic][iy];
		}
	}
	TotalTests.GetSummGenOut(0);
	TotPosTests.GetSummGenOut(0);
	NewDiagnoses.GetSummGenOut(0);
	DiagNoART200.GetSummGenOut(0);
	DiagNoART350.GetSummGenOut(0);
	DiagNoART500.GetSummGenOut(0);
	DiagNoART500plus.GetSummGenOut(0);
	DiagnosedUntreated.GetSummGenOut(0);
	UndiagnosedAdults.GetSummGenOut(0);

	// ART outputs
	SummOutRow += 3;
	TreatedAdults.GetSummGenOut(0);
	AdultARTcoverage.GetSummGenOut(0);
	TreatedAdults200.GetSummGenOut(0);
	TreatedAdults350.GetSummGenOut(0);
	TreatedAdults500.GetSummGenOut(0);
	TreatedAdults500plus.GetSummGenOut(0);
	NewART200.GetSummGenOut(0);
	NewART350.GetSummGenOut(0);
	NewART500.GetSummGenOut(0);
	NewART500plus.GetSummGenOut(0);

	// Other HIV programme outputs
	SummOutRow += 3;
	MMCoperations.GetSummGenOut(0);
	TakingPrEP.GetSummGenOut(0);
	ProtSexActs.GetSummGenOut(0);
	TotalCosts.GetSummGenOut(0);

	// STI incidence and prevalence
	SummOutRow += 3;
	NewHSV.GetSummGenOut(0);
	NewTP.GetSummGenOut(0);
	NewNG.GetSummGenOut(0);
	NewCT.GetSummGenOut(0);
	NewTV.GetSummGenOut(0);
	HSVprev15to49M.GetSummGenOut(0);
	HSVprev15to49F.GetSummGenOut(0);
	TPprev15to49M.GetSummGenOut(0);
	TPprev15to49F.GetSummGenOut(0);
	NGprev15to49M.GetSummGenOut(0);
	NGprev15to49F.GetSummGenOut(0);
	CTprev15to49M.GetSummGenOut(0);
	CTprev15to49F.GetSummGenOut(0);
	TVprev15to49M.GetSummGenOut(0);
	TVprev15to49F.GetSummGenOut(0);
	BVprev15to49F.GetSummGenOut(0);
	VCprev15to49F.GetSummGenOut(0);

	// Fertility and contraception
	SummOutRow += 3;
	TotBirths.GetSummGenOut(0);
	TeenBirths.GetSummGenOut(0);
	TotBirthsHIV.GetSummGenOut(0);
	TotBirthsART.GetSummGenOut(0);
	CurrInjectable.GetSummGenOut(0);
	CurrPill.GetSummGenOut(0);

	// Other demographic outputs
	SummOutRow += 3;
	TotPop.GetSummGenOut(0);
	TotPop15to49_U.GetSummGenOut(0);
	TotPop15to49_R.GetSummGenOut(0);
	UrbanTrend.GetSummGenOut(0);
	
	// Sexual behaviour
	SummOutRow += 3;
	CondomUse15to24.GetSummGenOut(0);
	CondomUse25to49.GetSummGenOut(0);
	MultPartnersM.GetSummGenOut(0);
	MultPartnersF.GetSummGenOut(0);
	MarriedTrendM.GetSummGenOut(0);
	MarriedTrendF.GetSummGenOut(0);
	EverMSMpropn.GetSummGenOut(0);

	// Education
	SummOutRow += 3;
	CompletedMatric.GetSummGenOut(0);
	CompletedMatricAM.GetSummGenOut(0);
	CompletedMatricAF.GetSummGenOut(0);
	CompletedMatricCM.GetSummGenOut(0);
	CompletedMatricCF.GetSummGenOut(0);
	CompletedMatricWM.GetSummGenOut(0);
	CompletedMatricWF.GetSummGenOut(0);
	TertiaryEnrolRatioA.GetSummGenOut(0);
	TertiaryEnrolRatioC.GetSummGenOut(0);
	TertiaryEnrolRatioW.GetSummGenOut(0);
	DropoutDuePregnancy.GetSummGenOut(0);

	// Employment
	SummOutRow += 3;
	EmployedPropn.GetSummGenOut(0);
	EmployedPropnAM.GetSummGenOut(0);
	EmployedPropnAF.GetSummGenOut(0);
	EmployedPropnCM.GetSummGenOut(0);
	EmployedPropnCF.GetSummGenOut(0);
	EmployedPropnWM.GetSummGenOut(0);
	EmployedPropnWF.GetSummGenOut(0);

	// Income and income inequality
	SummOutRow += 3;
	MedianIncome.GetSummGenOut(0);
	ZeroIncomeHH.GetSummGenOut(0);
	PalmaRatio.GetSummGenOut(0);
	GiniIncome.GetSummGenOut(0);
	GiniIncomeB.GetSummGenOut(0);
	GiniIncomeC.GetSummGenOut(0);
	GiniIncomeW.GetSummGenOut(0);

	// Incarceration
	SummOutRow += 3;
	FractionInPrison.GetSummGenOut(0);
	PrisonYouth.GetSummGenOut(0);
	PrisonCompletedSchool.GetSummGenOut(0);
	PrisonPrevious.GetSummGenOut(0);
	PrisonSentenced.GetSummGenOut(0);

	// Alcohol use
	SummOutRow += 3;
	BingeDrinking15to24M.GetSummGenOut(0);
	BingeDrinking15to24F.GetSummGenOut(0);
	BingeDrinking15to49M.GetSummGenOut(0);
	BingeDrinking15to49F.GetSummGenOut(0);
	BingeDrinkingHSchoolM.GetSummGenOut(0);
	BingeDrinkingHSchoolF.GetSummGenOut(0);
	AlcoholPast12moM.GetSummGenOut(0);
	AlcoholPast12moF.GetSummGenOut(0);
	AveDrinksPerDay.GetSummGenOut(0);

	// Household structures
	SummOutRow += 3;
	IndivsInHHsized1or2.GetSummGenOut(0);
	IndivsInHHsized3or4.GetSummGenOut(0);
	IndivsInHHsized5or6.GetSummGenOut(0);
	IndivsInHHsized7to9.GetSummGenOut(0);
	IndivsInHHsized10plus.GetSummGenOut(0);
	HeadOfHH.GetSummGenOut(0);
	PartnerOfHead.GetSummGenOut(0);
	ChildOfHead.GetSummGenOut(0);
	GrandchildOfHead.GetSummGenOut(0);

	// Homelessness outputs
	SummOutRow += 3;
	Homeless.GetSummGenOut(0);
	HomelessAdult.GetSummGenOut(0);
	HomelessChild.GetSummGenOut(0);
	HomelessChildPropn.GetSummGenOut(0);
	HomelessFem.GetSummGenOut(0);
	HomelessEmployed.GetSummGenOut(0);
	HomelessAveDur.GetSummGenOut(0);
	HomelessDrinkGT2perWeek.GetSummGenOut(0);

	// Inequitable gender norms
	SummOutRow += 3;
	IneqGenderAllM.GetSummGenOut(0);
	IneqGender15to24M.GetSummGenOut(0);
	IneqGender25to34M.GetSummGenOut(0);
	IneqGender35plusM.GetSummGenOut(0);

	int i, c;
	ostringstream s;

	if (process_num > 0){
		s << process_num << filout;
	}
	else{
		s << filout;
	}

	//ofstream file(filout);
	ofstream file(s.str().c_str()); // Converts s to a C string

	for (i = 0; i<1000; i++){
		file << setw(6) << right;
		for (c = 0; c<56; c++){
			//file << "	" << setw(10) << right << SummOut[i][c];
			if (SummOut[i][c]<0.0 || SummOut[i][c]>0.0){
				file << setw(10) << right << SummOut[i][c] << "	";
			}
			else{
				file << setw(10) << right << "	";
			}
		}
		file << endl;
	}
	file.close();
}

void StoreStructuralOutputs(const char* filout)
{
	int ic, iy;
	double temp1, temp2;

	SummOutRow = 0;
	IneqGenderBingeAssn.GetSummGenOut(0);
	IneqGenderMultAssn.GetSummGenOut(0);
	BingeMultAssnM.GetSummGenOut(0);
	BingeMultAssnF.GetSummGenOut(0);
	BingeCondomAssnM.GetSummGenOut(0);
	BingeCondomAssnF.GetSummGenOut(0);
	EmployedTransAssn.GetSummGenOut(0);
	EmployedMultAssnM.GetSummGenOut(0);
	EmployedMultAssnF.GetSummGenOut(0);
	EmployedHIVassnM.GetSummGenOut(0);
	EmployedHIVassnF.GetSummGenOut(0);
	EduCondomAssnM.GetSummGenOut(0);
	EduCondomAssnF.GetSummGenOut(0);
	EduHIVassnM.GetSummGenOut(0);
	EduHIVassnF.GetSummGenOut(0);
	SchoolMarriageAssn.GetSummGenOut(0);

	int i, c;
	ostringstream s;

	if (process_num > 0){
		s << process_num << filout;
	}
	else{
		s << filout;
	}

	//ofstream file(filout);
	ofstream file(s.str().c_str()); // Converts s to a C string

	for (i = 0; i<32; i++){
		// Increase upper limit on i if you add more associations.
		file << setw(6) << right;
		for (c = 0; c<56; c++){
			//file << "	" << setw(10) << right << SummOut[i][c];
			if (SummOut[i][c]<0.0 || SummOut[i][c]>0.0){
				file << setw(10) << right << SummOut[i][c] << "	";
			}
			else{
				file << setw(10) << right << "	";
			}
		}
		file << endl;
	}
	file.close();
}

void AggregateSims()
{
	int SimCount2, iy;

	SimCount2 = (CurrSim - 1) / IterationsPerPC;
	if (HIVcalib == 1){
		HIVtransitionM.GetPrev();
		HIVtransitionF.GetPrev();
		HIVtransitionF.GetCSWprev();
	}
	if (HSVcalib == 1){ 
		HSVtransitionF.GetCSWprev(); 
		HSVtransitionF.GetANCprevSmooth();
	}
	if (TPcalib == 1){ 
		TPtransitionF.GetCSWprev(); 
		TPtransitionF.GetANCprevSmooth();
	}
	if (NGcalib == 1){ 
		NGtransitionF.GetCSWprev(); 
		NGtransitionF.GetANCprevSmooth();
	}
	if (CTcalib == 1){ 
		CTtransitionF.GetCSWprev(); 
		CTtransitionF.GetANCprevSmooth();
	}
	if (TVcalib == 1){ 
		TVtransitionF.GetCSWprev(); 
		TVtransitionF.GetANCprevSmooth();
	}
	if (FixedUncertainty == 1){
		for (iy = 3; iy < ProjectionTerm-3; iy++){
			if (HIVind == 1){
				HIVprevCSW.out[CurrSim - 1][iy] = 0.05 * HIVtransitionF.CSWprevUnsmoothed[iy - 3] +
					0.12 * HIVtransitionF.CSWprevUnsmoothed[iy - 2] + 0.20 * HIVtransitionF.CSWprevUnsmoothed[iy - 1] +
					0.26 * HIVtransitionF.CSWprevUnsmoothed[iy] + 0.20 * HIVtransitionF.CSWprevUnsmoothed[iy + 1] +
					0.12 * HIVtransitionF.CSWprevUnsmoothed[iy + 2] + 0.05 * HIVtransitionF.CSWprevUnsmoothed[iy + 3];
			}
			if (HSVcalib == 1){
				//HSVprevCSW.out[CurrSim - 1][iy] = HSVtransitionF.CSWprevUnsmoothed[iy];
				HSVprevCSW.out[CurrSim - 1][iy] = 0.05 * HSVtransitionF.CSWprevUnsmoothed[iy - 3] +
					0.12 * HSVtransitionF.CSWprevUnsmoothed[iy - 2] + 0.20 * HSVtransitionF.CSWprevUnsmoothed[iy - 1] +
					0.26 * HSVtransitionF.CSWprevUnsmoothed[iy] + 0.20 * HSVtransitionF.CSWprevUnsmoothed[iy + 1] +
					0.12 * HSVtransitionF.CSWprevUnsmoothed[iy + 2] + 0.05 * HSVtransitionF.CSWprevUnsmoothed[iy + 3];
			}
			if (TPcalib == 1){
				TPprevCSW.out[CurrSim - 1][iy] = 0.05 * TPtransitionF.CSWprevUnsmoothed[iy - 3] +
					0.12 * TPtransitionF.CSWprevUnsmoothed[iy - 2] + 0.20 * TPtransitionF.CSWprevUnsmoothed[iy - 1] +
					0.26 * TPtransitionF.CSWprevUnsmoothed[iy] + 0.20 * TPtransitionF.CSWprevUnsmoothed[iy + 1] +
					0.12 * TPtransitionF.CSWprevUnsmoothed[iy + 2] + 0.05 * TPtransitionF.CSWprevUnsmoothed[iy + 3];
			}
			if (NGcalib == 1){
				//NGprevCSW.out[CurrSim - 1][iy] = NGtransitionF.CSWprevUnsmoothed[iy];
				NGprevCSW.out[CurrSim - 1][iy] = 0.05 * NGtransitionF.CSWprevUnsmoothed[iy - 3] +
					0.12 * NGtransitionF.CSWprevUnsmoothed[iy - 2] + 0.20 * NGtransitionF.CSWprevUnsmoothed[iy - 1] +
					0.26 * NGtransitionF.CSWprevUnsmoothed[iy] + 0.20 * NGtransitionF.CSWprevUnsmoothed[iy + 1] +
					0.12 * NGtransitionF.CSWprevUnsmoothed[iy + 2] + 0.05 * NGtransitionF.CSWprevUnsmoothed[iy + 3]; 
			}
			if (CTcalib == 1){
				//CTprevCSW.out[CurrSim - 1][iy] = CTtransitionF.CSWprevUnsmoothed[iy];
				CTprevCSW.out[CurrSim - 1][iy] = 0.05 * CTtransitionF.CSWprevUnsmoothed[iy - 3] +
					0.12 * CTtransitionF.CSWprevUnsmoothed[iy - 2] + 0.20 * CTtransitionF.CSWprevUnsmoothed[iy - 1] +
					0.26 * CTtransitionF.CSWprevUnsmoothed[iy] + 0.20 * CTtransitionF.CSWprevUnsmoothed[iy + 1] +
					0.12 * CTtransitionF.CSWprevUnsmoothed[iy + 2] + 0.05 * CTtransitionF.CSWprevUnsmoothed[iy + 3]; 
			}
			if (TVcalib == 1){
				//TVprevCSW.out[CurrSim - 1][iy] = TVtransitionF.CSWprevUnsmoothed[iy];
				TVprevCSW.out[CurrSim - 1][iy] = 0.05 * TVtransitionF.CSWprevUnsmoothed[iy - 3] +
					0.12 * TVtransitionF.CSWprevUnsmoothed[iy - 2] + 0.20 * TVtransitionF.CSWprevUnsmoothed[iy - 1] +
					0.26 * TVtransitionF.CSWprevUnsmoothed[iy] + 0.20 * TVtransitionF.CSWprevUnsmoothed[iy + 1] +
					0.12 * TVtransitionF.CSWprevUnsmoothed[iy + 2] + 0.05 * TVtransitionF.CSWprevUnsmoothed[iy + 3]; 
			}
		}
	}
	CalcTotalLogL();
	if (HIVcalib == 1){
		HIVparamsLogL.out[SimCount2][9] = TotalLogL;
	}
}

void SimulateParameters()
{
	// This function gets called after the heterosexual behaviour and STI parameters have
	// been read in, so it should not be used to simulate other sources of uncertainty.
	if (HIVcalib == 1){ SimulateHIVparams(); }
	if (TPcalib == 1){ SimulateTPparameters(); }
	if (HSVcalib == 1){ SimulateHSVparameters(); }
	if (NGcalib == 1){ SimulateNGparameters(); }
	if (CTcalib == 1){ SimulateCTparameters(); }
	if (TVcalib == 1){ SimulateTVparameters(); }
	if (BVcalib == 1){ SimulateBVparameters(); }
	if (VCcalib == 1){ SimulateVCparameters(); }
	// Note that we don't call SimulateMSMparameters and SimulateMSM_HIV here because 
	// they both get called in the ReadMSM function, which gets called later.
}

void SimulateHIVparams()
{
	int SimCount2, seed, ii, ind, is;
	double x, y, a, b, p, q, bound;
	double r[9], MeanSurvival;

	SimCount2 = (CurrSim - 1) / IterationsPerPC;
	if (FixedUncertainty == 0){
		seed = SimCount2 * 91 + process_num * 7927;
		CRandomMersenne rg(seed);
		for (ii = 0; ii < 9; ii++){
			r[ii] = rg.Random();
		}
	}
	else{
		SimCount2 = (CurrSim - 1) / IterationsPerPC;
		for (ii = 0; ii < 9; ii++){
			r[ii] = RandomUniformHIV.out[SimCount2][ii];
			//r[ii] = RandomUniformHIV.out[CurrSim-1][ii];
		}
	}

	// Simulate relative infectiousness during acute HIV
	a = 7.1111;
	b = 0.4444;
	p = r[0];
	q = 1 - r[0];
	bound = 0.0;
	ind = 2;
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, &bound);
	HIVtransitionM.HIVinfecIncrease[0] = x - 1.0;
	HIVtransitionF.HIVinfecIncrease[0] = x - 1.0;

	// Simulate M->F transmission prob per sex act in SW-client contacts
	a = 11.099;
	b = 11087.9;
	p = r[1];
	q = 1 - r[1];
	bound = 0.0;
	ind = 2;
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
	InitHIVtransm[0][0] = x;

	// Simulate relative infectiousness during advanced HIV
	a = 9.00;
	b = 5.00;
	p = r[8];
	q = 1 - r[8];
	bound = 0.0;
	ind = 2;
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, &bound);
	HIVtransitionM.HIVinfecIncrease[3] = x;
	HIVtransitionF.HIVinfecIncrease[3] = x;

	// Simulate M->F transmission prob per sex act in short-term partnerships
	a = 24.94;
	b = 9949.1;
	//a = 5.6789;
	//b = 467.56;
	p = r[2];
	q = 1 - r[2];
	bound = 0.0;
	ind = 2;
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
	InitHIVtransm[1][0] = x;

	// Simulate F->M transmission prob per sex act in short-term partnerships
	a = 17.34;
	b = 13853.2;
	//a = 8.9400;
	//b = 1481.06;
	p = r[3];
	q = 1 - r[3];
	bound = 0.0;
	ind = 2;
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
	InitHIVtransm[1][1] = x;
	InitHIVtransm[0][1] = x;

	// Simulate M->F transmission prob per sex act in spousal partnerships (low risk spouse)
	a = 15.986;
	b = 19967.0;
	//a = 7.0949;
	//b = 3540.3;
	p = r[4];
	q = 1 - r[4];
	bound = 0.0;
	ind = 2;
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
	//InitHIVtransm[2][0] = x;
	InitHIVtransm[2][0] = r[4] * InitHIVtransm[1][0];

	// Simulate F->M transmission prob per sex act in spousal partnerships (low risk spouse)
	a = 6.2464;
	b = 12486.5;
	//a = 7.0949;
	//b = 3540.3;
	p = r[5];
	q = 1 - r[5];
	bound = 0.0;
	ind = 2;
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
	//InitHIVtransm[2][1] = x;
	InitHIVtransm[2][1] = r[5] * InitHIVtransm[1][1];

	// Simulate mean survival in absence of HIV
	/*a = 144.0;
	b = 12.0;
	p = r[6];
	q = 1 - r[6];
	bound = 0.0;
	ind = 2;
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, &bound);
	MeanSurvival = x;
	HIVtransitionM.AveDuration[2] = MeanSurvival * 52.0 * 0.358;
	HIVtransitionF.AveDuration[2] = MeanSurvival * 52.0 * 0.358;
	HIVtransitionM.AveDuration[3] = MeanSurvival * 52.0 * 0.174;
	HIVtransitionF.AveDuration[3] = MeanSurvival * 52.0 * 0.174;
	HIVtransitionM.AveDuration[1] = MeanSurvival * 52.0 * (1.0 - 0.358 - 0.174) -
		HIVtransitionM.AveDuration[0];
	HIVtransitionF.AveDuration[1] = MeanSurvival * 52.0 * (1.0 - 0.358 - 0.174) -
		HIVtransitionF.AveDuration[0];*/

	// Simulate bias in reporting of condom use
	//CondomScaling = r[8];
	CondomScaling = 0.20; // Same value as in Thembisa

	// Simulate initial HIV prevalence
	InitHIVprevHigh = (r[6] * 0.02) + 0.01;

	// Simulate sexual mixing parameter
	AssortativeM = r[7];
	AssortativeF = r[7];

	// Simulate the effect of VL on HIV transmission
	/*a = 12.09;
	b = 15.90;
	p = r[8];
	q = 1 - r[8];
	bound = 0.0;
	ind = 2;
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
	VLeffectTransm[0] = x;*/

	// Simulate the standard deviation of the variation in desired partner acquisition rates
	/*a = 4.0;
	b = 8.0;
	p = r[8];
	q = 1 - r[8];
	bound = 0.0;
	ind = 2;
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, &bound);
	SDpartnerRateAdj = x;*/

	// Simulate the relative rate of partner acquisition in married women
	//PartnerEffectNew[1][1] = PartnerEffectNew[0][1] * r[8];

	// Store the parameter combination
	if (CurrSim == (SimCount2 * IterationsPerPC + 1)){
		for (ii = 0; ii < 9; ii++){
			HIVparamsLogL.out[SimCount2][ii] = r[ii];
		}
	}
}

void SimulateTPparameters()
{
	int ind, i, is, seed;
	double x, y, a, b, p, q;
	double r[10]; // Random variables from U(0, 1)

	ind = 2;

	if(FixedUncertainty==0){
		seed = (CurrSim - 1) * 91 + process_num * 7927;
		CRandomMersenne rg(seed);
		for(i=0; i<10; i++){
			r[i] = rg.Random();
			RandomUniformTP.out[CurrSim-1][i] = r[i];
		}
	}
	else{
		for(i=0; i<10; i++){
			r[i] = RandomUniformTP.out[CurrSim-1][i];}
	}

	// Simulate M->F transmission prob
	a = 4.44;
	b = 13.31;
	p = r[0];
	q = 1 - r[0];
	cdfbet(&ind,&p,&q,&x,&y,&a,&b,0,0);
	x = 0.228; // Median of 100-best fitting parameters
	TPtransitionM.TransmProb = x;

	//TPtransitionM.RelTransmCSW = TPtransitionM.TransmProbSW/TPtransitionM.TransmProb;
	TPtransitionM.RelTransmCSW = 1.0;
	TPtransitionM.TransmProbSW = x;

	// Simulate F->M transmission prob
	a = 7.5;
	b = 42.5;
	p = r[1];
	q = 1 - r[1];
	cdfbet(&ind,&p,&q,&x,&y,&a,&b,0,0);
	x = 0.155; // Median of 100-best fitting parameters
	TPtransitionF.TransmProb = x;
	
	// Simulate the average duration of primary syphilis (mean 6.6, std dev 2)
	a = 10.89;
	b = 1.65;
	p = r[2];
	q = 1 - r[2];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	x = 6.59; // Median of 100-best fitting parameters
	TPtransitionM.AveDuration[1] = x;
	TPtransitionF.AveDuration[1] = x;

	// Simulate the average duration of 2ndary syphilis (mean 15.6, std dev 4)
	a = 15.21;
	b = 0.975;
	p = r[3];
	q = 1 - r[3];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	x = 15.1; // Median of 100-best fitting parameters
	TPtransitionM.AveDuration[2] = x;
	TPtransitionF.AveDuration[2] = x;

	// Simulate the average duration of latent syphilis (mean 520, std dev 150)
	a = 12.02;
	b = 0.02311;
	p = r[4];
	q = 1 - r[4];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	x = 595.0; // Median of 100-best fitting parameters
	TPtransitionM.AveDuration[3] = x;
	TPtransitionF.AveDuration[3] = x;

	// Simulate the average duration of early immunity (mean 26, std dev 8)
	a = 10.56;
	b = 0.4063;
	p = r[5];
	q = 1 - r[5];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	x = 23.0; // Median of 100-best fitting parameters
	TPtransitionM.AveDuration[4] = x;
	TPtransitionF.AveDuration[4] = x;

	// Simulate the average duration of late immunity (mean 52, std dev 16)
	a = 10.56;
	b = 0.2031;
	p = r[6];
	q = 1 - r[6];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	x = 53.8; // Median of 100-best fitting parameters
	TPtransitionM.AveDuration[5] = x;
	TPtransitionF.AveDuration[5] = x;

	// Simulate the std deviation of the study effect
	/*a = 4.0;
	b = 13.333;
	p = r[7];
	q = 1 - r[7];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	TPtransitionM.HouseholdLogL.VarStudyEffect = pow(x, 2.0);
	TPtransitionF.HouseholdLogL.VarStudyEffect = pow(x, 2.0);
	TPtransitionF.ANClogL.VarStudyEffect = pow(x, 2.0);
	TPtransitionF.FPClogL.VarStudyEffect = pow(x, 2.0);
	TPtransitionF.CSWlogL.VarStudyEffect = pow(x, 2.0);
	TPparamsLogL.out[CurrSim - 1][7] = pow(x, 2.0);*/
	//SecondaryRxMult = r[7];
	//TPparamsLogL.out[CurrSim - 1][7] = r[7];
	SecondaryRxMult = 0.682; // Median of 100-best fitting parameters

	// Simulate the propn of cases correctly treated prior to introduction of SM
	a = 14.0;
	b = 6.0;
	p = r[8];
	q = 1 - r[8];
	cdfbet(&ind,&p,&q,&x,&y,&a,&b,0,0);
	x = 0.711; // Median of 100-best fitting parameters
	TPtransitionF.CorrectRxPreSM = x;
	TPtransitionM.CorrectRxPreSM = x;

	// Simulate the propn of primary syphilis cases that are RPR-neg after Rx
	a = 9.2;
	b = 13.8;
	p = r[9];
	q = 1 - r[9];
	cdfbet(&ind,&p,&q,&x,&y,&a,&b,0,0);
	x = 0.398; // Median of 100-best fitting parameters
	TPtransitionF.PropnSuscepAfterRx = x;
	TPtransitionM.PropnSuscepAfterRx = x;

	// Store parameters
	for (i = 0; i < 10; i++){
		// Note that for i=7 the simulated random number doesn't actually get used
		TPparamsLogL.out[CurrSim - 1][i] = r[i];
	}
}

void SimulateHSVparameters()
{
	int ind, i, is, seed;
	double x, y, a, b, p, q;
	double r[11]; // Random variables from U(0, 1)

	ind = 2;

	if (FixedUncertainty == 0){
		seed = (CurrSim - 1) * 91 + process_num * 7927;
		CRandomMersenne rg(seed);
		for (i = 0; i<11; i++){
			r[i] = rg.Random();
			RandomUniformHSV.out[CurrSim - 1][i] = r[i];
		}
	}
	else{
		for (i = 0; i<11; i++){
			r[i] = RandomUniformHSV.out[CurrSim - 1][i];
		}
	}

	// Simulate M->F transmission prob in ST relationships (mean 0.0095, SD 0.0038)
	a = 6.1811;
	b = 644.46;
	p = r[0];
	q = 1 - r[0];
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
	x = 0.016; // Median of 100-best fitting parameters
	HSVtransitionM.TransmProb = x;
	
	// Simulate F->M transmission prob in ST relationships (mean 0.0065, SD 0.0026)
	a = 6.2029;
	b = 948.09;
	p = r[1];
	q = 1 - r[1];
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
	x = 0.0035; // Median of 100-best fitting parameters
	HSVtransitionF.TransmProb = x;
	
	// Simulate M->F transmission prob in LT relationships (mean 0.0009, SD 0.00036)
	a = 6.2435;
	b = 6931.0;
	p = r[2];
	q = 1 - r[2];
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
	x = 0.0009; // Median of 100-best fitting parameters
	HSVtransitionM.RelTransmLT = x / HSVtransitionM.TransmProb;
	
	// Simulate F->M transmission prob in LT relationships (mean 0.00015, SD 0.00006)
	a = 6.2489;
	b = 41653.2;
	p = r[3];
	q = 1 - r[3];
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
	x = 0.0001; // Median of 100-best fitting parameters
	HSVtransitionM.RelTransmLT = x / HSVtransitionM.TransmProb;
	
	// Simulate M->F transmission prob in CSW-client relationships (mean 0.002, SD 0.0005)
	a = 15.966;
	b = 7967.0;
	p = r[4];
	q = 1 - r[4];
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
	x = 0.0016; // Median of 100-best fitting parameters
	HSVtransitionM.TransmProbSW = x;
	HSVtransitionM.RelTransmCSW = x / HSVtransitionM.TransmProb;
	
	// Simulate the symptomatic proportion (mean 0.15, SD 0.05)
	a = 7.5;
	b = 42.5;
	p = r[5];
	q = 1 - r[5];
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
	x = 0.14; // Median of 100-best fitting parameters
	HSVtransitionM.SymptomaticPropn = x;
	HSVtransitionF.SymptomaticPropn = x;
	
	// Simulate the average male frequency of reactivation pa (mean 6, std dev 1)
	a = 36.0;
	b = 6.0;
	p = r[6];
	q = 1 - r[6];
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
	x = 5.95; // Median of 100-best fitting parameters
	HSVtransitionM.RecurrenceRate = x / 52.0;
	//InitRecurrenceRateM = x/52.0;

	// Simulate the average female frequency of reactivation pa (mean 3, std dev 0.5)
	a = 36.0;
	b = 12.0;
	p = r[7];
	q = 1 - r[7];
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
	x = 2.88; // Median of 100-best fitting parameters
	HSVtransitionF.RecurrenceRate = x / 52.0;
	//InitRecurrenceRateF = x/52.0;

	// Simulate the increase in infectiousness during symptomatic episodes (mean 15, SD 5)
	a = 9.0;
	b = 0.6;
	p = r[8];
	q = 1 - r[8];
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
	x = 13.9; // Median of 100-best fitting parameters
	HSVsymptomInfecIncrease = x;

	// Simulate the rate of transition out of early HSV stage (mean 0.1, std dev 0.02)
	a = 25.0;
	b = 250.0;
	p = r[9];
	q = 1 - r[9];
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
	x = 0.098; // Median of 100-best fitting parameters
	HSVtransitionM.AveDuration[1] = 52.0 / x;
	HSVtransitionF.AveDuration[1] = 52.0 / x;

	// Simulate the std deviation of the study effect
	/*a = 4.0;
	b = 13.333;
	p = r[10];
	q = 1 - r[10];
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
	HSVtransitionM.HouseholdLogL.VarStudyEffect = pow(x, 2.0);
	HSVtransitionF.HouseholdLogL.VarStudyEffect = pow(x, 2.0);
	HSVtransitionF.ANClogL.VarStudyEffect = pow(x, 2.0);
	HSVtransitionF.FPClogL.VarStudyEffect = pow(x, 2.0);
	HSVtransitionF.CSWlogL.VarStudyEffect = pow(x, 2.0);
	HSVparamsLogL.out[CurrSim - 1][10] = pow(x, 2.0);*/

	// Store parameters
	for (i = 0; i < 11; i++){
		// Note that for i=10 the simulated random number doesn't actually get used
		HSVparamsLogL.out[CurrSim - 1][i] = r[i];
	}
}

void SimulateNGparameters()
{
	int ind, i, is, seed;
	double x, y, a, b, p, q;
	double r[10]; // Random variables from U(0, 1)

	ind = 2;

	if(FixedUncertainty==0){
		seed = (CurrSim - 1) * 91 + process_num * 7927;
		CRandomMersenne rg(seed);
		for(i=0; i<10; i++){
			r[i] = rg.Random();
			RandomUniformNG.out[CurrSim-1][i] = r[i];
		}
	}
	else{
		for(i=0; i<10; i++){
			r[i] = RandomUniformNG.out[CurrSim-1][i];}
	}

	// Simulate M->F transmission prob
	a = 9.2;
	b = 13.8;
	p = r[0];
	q = 1 - r[0];
	cdfbet(&ind,&p,&q,&x,&y,&a,&b,0,0);
	x = 0.459; // Median of 100-best fitting parameters
	NGtransitionM.TransmProb = x;
	
	//NGtransitionM.RelTransmCSW = NGtransitionM.TransmProbSW/NGtransitionM.TransmProb;
	NGtransitionM.RelTransmCSW = 1.0;
	NGtransitionM.TransmProbSW = x;

	// Simulate F->M transmission prob
	a = 12.6;
	b = 50.4;
	p = r[1];
	q = 1 - r[1];
	cdfbet(&ind,&p,&q,&x,&y,&a,&b,0,0);
	x = 0.237; // Median of 100-best fitting parameters
	NGtransitionF.TransmProb = x;
	
	// Simulate % of male NG cases that become symptomatic
	a = 31.5;
	b = 3.5;
	p = r[2];
	q = 1 - r[2];
	cdfbet(&ind,&p,&q,&x,&y,&a,&b,0,0);
	x = 0.866; // Median of 100-best fitting parameters
	NGtransitionM.SymptomaticPropn = x;
	
	// Simulate % of female NG cases that become symptomatic
	a = 3.867;
	b = 5.8;
	p = r[3];
	q = 1 - r[3];
	cdfbet(&ind,&p,&q,&x,&y,&a,&b,0,0);
	x = 0.299; // Median of 100-best fitting parameters
	NGtransitionF.SymptomaticPropn = x;
	
	// Simulate the average duration of untreated NG in males
	a = 4.0;
	b = 0.2;
	p = r[4];
	q = 1 - r[4];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	x = 34.0; // Median of 100-best fitting parameters
	NGtransitionM.AveDuration[0] = x;
	NGtransitionM.AveDuration[1] = x;
	
	// Simulate the average duration of untreated NG in females
	a = 4.0;
	b = 0.2;
	p = r[5];
	q = 1 - r[5];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	x = 33.6; // Median of 100-best fitting parameters
	NGtransitionF.AveDuration[0] = x;
	NGtransitionF.AveDuration[1] = x;
	
	// Simulate the std deviation of the study effect
	/*a = 4.0;
	b = 13.333;
	p = r[6];
	q = 1 - r[6];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	NGtransitionM.HouseholdLogL.VarStudyEffect = pow(x, 2.0);
	NGtransitionF.HouseholdLogL.VarStudyEffect = pow(x, 2.0);
	NGtransitionF.ANClogL.VarStudyEffect = pow(x, 2.0);
	NGtransitionF.FPClogL.VarStudyEffect = pow(x, 2.0);
	NGtransitionF.CSWlogL.VarStudyEffect = pow(x, 2.0);
	NGparamsLogL.out[CurrSim - 1][6] = pow(x, 2.0);*/
	
	// Simulate the proportion of treated NG cases that are immune to reinfection
	r[7] = 0.401; // Median of 100-best fitting parameters
	NGtransitionM.PropnImmuneAfterRx = r[7];
	NGtransitionF.PropnImmuneAfterRx = r[7];
	
	// Simulate the average duration of immunity
	a = 4.0;
	b = 4.0;
	p = r[8];
	q = 1 - r[8];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	//NGtransitionM.AveDuration[2] = 52.0 * x;
	//NGtransitionF.AveDuration[2] = 52.0 * x;
	//NGparamsLogL.out[CurrSim - 1][8] = 52.0 * x;
	NGtransitionM.AveDuration[2] = 48.8;
	NGtransitionF.AveDuration[2] = 48.8;
	
	// Simulate the proportion of cases correctly treated prior to syndromic mngt
	a = 14.0;
	b = 6.0;
	p = r[9];
	q = 1 - r[9];
	cdfbet(&ind,&p,&q,&x,&y,&a,&b,0,0);
	x = 0.704; // Median of 100-best fitting parameters
	NGtransitionM.CorrectRxPreSM = x;
	NGtransitionF.CorrectRxPreSM = x;
	
	// Store parameters
	for (i = 0; i < 10; i++){
		NGparamsLogL.out[CurrSim - 1][i] = r[i];
	}
}

void SimulateCTparameters()
{
	int ind, i, is, seed;
	double x, y, a, b, p, q;
	double r[10]; // Random variables from U(0, 1)

	ind = 2;

	if(FixedUncertainty==0){
		seed = (CurrSim - 1) * 91 + process_num * 7927;
		CRandomMersenne rg(seed);
		for(i=0; i<10; i++){
			r[i] = rg.Random();
			RandomUniformCT.out[CurrSim-1][i] = r[i];
		}
	}
	else{
		for(i=0; i<10; i++){
			r[i] = RandomUniformCT.out[CurrSim-1][i];}
	}

	// Simulate M->F transmission prob
	a = 3.4;
	b = 24.933;
	p = r[0];
	q = 1 - r[0];
	cdfbet(&ind,&p,&q,&x,&y,&a,&b,0,0);
	x = 0.162; // Median of 100-best fitting parameters
	CTtransitionM.TransmProb = x;
	
	//CTtransitionM.RelTransmCSW = CTtransitionM.TransmProbSW/CTtransitionM.TransmProb;
	CTtransitionM.RelTransmCSW = 1.0;
	CTtransitionM.TransmProbSW = x;

	// Simulate F->M transmission prob
	a = 1.99;
	b = 10.45;
	p = r[1];
	q = 1 - r[1];
	cdfbet(&ind,&p,&q,&x,&y,&a,&b,0,0);
	x = 0.0975; // Median of 100-best fitting parameters
	CTtransitionF.TransmProb = x;
	
	// Simulate % of male CT cases that become symptomatic
	a = 2.5;
	b = 5.8333;
	p = r[2];
	q = 1 - r[2];
	cdfbet(&ind,&p,&q,&x,&y,&a,&b,0,0);
	x = 0.367; // Median of 100-best fitting parameters
	CTtransitionM.SymptomaticPropn = x;
	
	// Simulate % of female CT cases that become symptomatic
	a = 2.838;
	b = 16.084;
	p = r[3];
	q = 1 - r[3];
	cdfbet(&ind,&p,&q,&x,&y,&a,&b,0,0);
	x = 0.116; // Median of 100-best fitting parameters
	CTtransitionF.SymptomaticPropn = x;
	
	// Simulate the average duration of untreated symptomatic CT
	a = 4.0;
	b = 0.25;
	p = r[4];
	q = 1 - r[4];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	x = 15.0; // Median of 100-best fitting parameters
	CTtransitionM.AveDuration[0] = x;
	CTtransitionF.AveDuration[0] = x;
	
	// Simulate the average duration of untreated asymptomatic CT
	a = 36.0;
	b = 0.4;
	p = r[5];
	q = 1 - r[5];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	x = 106.6; // Median of 100-best fitting parameters
	CTtransitionM.AveDuration[1] = x;
	CTtransitionF.AveDuration[1] = x;
	
	// Simulate the std deviation of the study effect
	/*a = 4.0;
	b = 13.333;
	p = r[6];
	q = 1 - r[6];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	CTtransitionM.HouseholdLogL.VarStudyEffect = pow(x, 2.0);
	CTtransitionF.HouseholdLogL.VarStudyEffect = pow(x, 2.0);
	CTtransitionF.ANClogL.VarStudyEffect = pow(x, 2.0);
	CTtransitionF.FPClogL.VarStudyEffect = pow(x, 2.0);
	CTtransitionF.CSWlogL.VarStudyEffect = pow(x, 2.0);
	CTparamsLogL.out[CurrSim - 1][6] = pow(x, 2.0);*/
	
	// Simulate the proportion of treated CT cases that are immune to reinfection
	r[7] = 0.732; // Median of 100-best fitting parameters
	CTtransitionM.PropnImmuneAfterRx = r[7];
	CTtransitionF.PropnImmuneAfterRx = r[7];
	
	// Simulate the average duration of immunity
	a = 6.76;
	b = 0.013;
	p = r[8];
	q = 1 - r[8];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	x = 295.0; // Median of 100-best fitting parameters
	CTtransitionM.AveDuration[2] = x;
	CTtransitionF.AveDuration[2] = x;
	
	// Simulate the proportion of cases correctly treated prior to syndromic mngt
	a = 14.0;
	b = 6.0;
	p = r[9];
	q = 1 - r[9];
	cdfbet(&ind,&p,&q,&x,&y,&a,&b,0,0);
	x = 0.711; // Median of 100-best fitting parameters
	CTtransitionM.CorrectRxPreSM = x;
	CTtransitionF.CorrectRxPreSM = x;
	
	// Store parameters
	for (i = 0; i < 10; i++){
		// Note that for i=6 the simulated value doesn't get used
		CTparamsLogL.out[CurrSim - 1][i] = r[i];
	}
}

void SimulateTVparameters()
{
	int ind, i, is, seed;
	double x, y, a, b, p, q;
	double r[12]; // Random variables from U(0, 1)

	ind = 2;

	if(FixedUncertainty==0){
		seed = (CurrSim - 1) * 91 + process_num * 7927;
		CRandomMersenne rg(seed);
		for(i=0; i<12; i++){
			r[i] = rg.Random();
			RandomUniformTV.out[CurrSim-1][i] = r[i];
		}
	}
	else{
		for(i=0; i<12; i++){
			r[i] = RandomUniformTV.out[CurrSim-1][i];}
	}

	// Simulate M->F transmission prob
	a = 2.8383;
	b = 16.084;
	p = r[0];
	q = 1 - r[0];
	cdfbet(&ind,&p,&q,&x,&y,&a,&b,0,0);
	x = 0.192; // Median of 100-best fitting parameters
	TVtransitionM.TransmProb = x;
	
	//TVtransitionM.RelTransmCSW = TVtransitionM.TransmProbSW/TVtransitionM.TransmProb;
	TVtransitionM.RelTransmCSW = 1.0;
	TVtransitionM.TransmProbSW = x;

	// Simulate F->M transmission prob
	a = 3.8;
	b = 91.2;
	p = r[1];
	q = 1 - r[1];
	cdfbet(&ind,&p,&q,&x,&y,&a,&b,0,0);
	x = 0.0386; // Median of 100-best fitting parameters
	TVtransitionF.TransmProb = x;
	
	// Simulate % of male TV cases that become symptomatic
	a = 9.2;
	b = 13.8;
	p = r[2];
	q = 1 - r[2];
	cdfbet(&ind,&p,&q,&x,&y,&a,&b,0,0);
	x = 0.391; // Median of 100-best fitting parameters
	TVtransitionM.SymptomaticPropn = x;
	
	// Simulate % of female TV cases that become symptomatic
	a = 6.0;
	b = 14.0;
	p = r[3];
	q = 1 - r[3];
	cdfbet(&ind,&p,&q,&x,&y,&a,&b,0,0);
	x = 0.275; // Median of 100-best fitting parameters
	TVtransitionF.SymptomaticPropn = x;
	
	// Simulate the average duration of untreated symptomatic TV in men
	a = 8.1633;
	b = 4.0816;
	p = r[4];
	q = 1 - r[4];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	x = 2.08; // Median of 100-best fitting parameters
	TVtransitionM.AveDuration[0] = x;
	
	// Simulate the average duration of untreated symptomatic TV in women
	a = 9.0;
	b = 0.6;
	p = r[5];
	q = 1 - r[5];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	x = 14.6; // Median of 100-best fitting parameters
	TVtransitionF.AveDuration[0] = x;
	
	// Simulate the average duration of untreated asymptomatic TV in men
	a = 8.1633;
	b = 0.4082;
	p = r[6];
	q = 1 - r[6];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	x = 20.1; // Median of 100-best fitting parameters
	TVtransitionM.AveDuration[1] = x;
	
	// Simulate the average duration of untreated asymptomatic TV in women
	a = 9.0;
	b = 0.06;
	p = r[7];
	q = 1 - r[7];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	x = 249.0; // Median of 100-best fitting parameters
	TVtransitionF.AveDuration[1] = x;
	
	// Simulate the std deviation of the study effect
	/*a = 4.0;
	b = 13.333;
	p = r[8];
	q = 1 - r[8];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	TVtransitionM.HouseholdLogL.VarStudyEffect = pow(x, 2.0);
	TVtransitionF.HouseholdLogL.VarStudyEffect = pow(x, 2.0);
	TVtransitionF.ANClogL.VarStudyEffect = pow(x, 2.0);
	TVtransitionF.FPClogL.VarStudyEffect = pow(x, 2.0);
	TVtransitionF.CSWlogL.VarStudyEffect = pow(x, 2.0);
	TVparamsLogL.out[CurrSim - 1][8] = pow(x, 2.0);*/
	
	// Simulate the average duration of immunity
	a = 4.0;
	b = 4.0;
	p = r[9];
	q = 1 - r[9];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	//TVtransitionM.AveDuration[2] = 52.0 * x;
	//TVtransitionF.AveDuration[2] = 52.0 * x;
	//TVparamsLogL.out[CurrSim - 1][9] = 52.0 * x;
	TVtransitionM.AveDuration[2] = 30.0; // Median of 100-best fitting parameters
	TVtransitionF.AveDuration[2] = 30.0;
	
	// Simulate propn of cases correctly treated pre-SM
	a = 3.8667;
	b = 5.8;
	p = r[10];
	q = 1 - r[10];
	cdfbet(&ind,&p,&q,&x,&y,&a,&b,0,0);
	x = 0.447; // Median of 100-best fitting parameters
	TVtransitionM.CorrectRxPreSM = x;
	TVtransitionF.CorrectRxPreSM = x;
	
	// Simulate the proportion of treated TV cases that are immune to reinfection
	r[11] = 0.741; // Median of 100-best fitting parameters
	TVtransitionM.PropnImmuneAfterRx = r[11];
	TVtransitionF.PropnImmuneAfterRx = r[11];

	// Store parameters
	for (i = 0; i < 12; i++){
		// Note that for i=8 the simulated random number doesn't actually get used
		TVparamsLogL.out[CurrSim - 1][i] = r[i];
	}
}

void SimulateBVparameters()
{
	int ind, i, is;
	double x, y, a, b, p, q;
	double r[8]; // Random variables from U(0, 1)

	ind = 2;

	if (FixedUncertainty == 0){
		int seed = (CurrSim - 1) * 91 + process_num * 7927;
		CRandomMersenne rg(seed);
		for (i = 0; i<8; i++){
			r[i] = rg.Random();
			RandomUniformBV.out[CurrSim - 1][i] = r[i];
		}
	}
	else{
		for (i = 0; i<8; i++){
			r[i] = RandomUniformBV.out[CurrSim - 1][i];
		}
	}

	// Simulate % of female BV cases that become symptomatic
	a = 4.4375;
	b = 13.313;
	p = r[0];
	q = 1 - r[0];
	cdfbet(&ind,&p,&q,&x,&y,&a,&b,0,0);
	BVtransitionF.SymptomaticPropn = x;
	BVparamsLogL.out[CurrSim-1][0] = x;

	// Simulate the weekly incidence of BV in women with intermediate flora (mean 0.1, SD 0.03)
	a = 11.111;
	b = 111.11;
	p = r[1];
	q = 1 - r[1];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	BVtransitionF.Incidence1 = x;
	BVparamsLogL.out[CurrSim-1][1] = x;

	// Simulate the transition rate from BV to normal flora (mean 0.008, SD 0.003)
	a = 7.1111;
	b = 888.89;
	p = r[2];
	q = 1 - r[2];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	BVtransitionF.CtsTransition[2][0] = x;
	BVtransitionF.CtsTransition[3][0] = x;
	BVparamsLogL.out[CurrSim-1][2] = x;

	// Simulate the transition rate from BV to intermediate flora (mean 0.051, SD 0.015)
	a = 11.56;
	b = 226.67;
	p = r[3];
	q = 1 - r[3];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	BVtransitionF.CtsTransition[2][1] = x;
	BVtransitionF.CtsTransition[3][1] = x;
	BVparamsLogL.out[CurrSim-1][3] = x;

	// Simulate the transition rate from normal to intermediate flora (mean 0.03, SD 0.01)
	a = 9.0;
	b = 300.0;
	p = r[4];
	q = 1 - r[4];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	BVtransitionF.CtsTransition[0][1] = x;
	BVparamsLogL.out[CurrSim-1][4] = x;

	// Simulate the transition rate from intermediate to normal flora (mean 0.069, SD 0.02)
	a = 11.903;
	b = 172.5;
	p = r[5];
	q = 1 - r[5];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	BVtransitionF.CtsTransition[1][0] = x;
	BVparamsLogL.out[CurrSim-1][5] = x;

	// Simulate the std deviation of the study effect
	/*a = 4.0;
	b = 13.333;
	p = r[6];
	q = 1 - r[6];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	BVtransitionF.HouseholdLogL.VarStudyEffect = pow(x, 2.0);
	BVtransitionF.ANClogL.VarStudyEffect = pow(x, 2.0);
	BVtransitionF.FPClogL.VarStudyEffect = pow(x, 2.0);
	BVtransitionF.CSWlogL.VarStudyEffect = pow(x, 2.0);
	BVparamsLogL.out[CurrSim-1][6] = pow(x, 2.0);*/

	// Simulate propn of cases correctly treated pre-SM
	a = 3.8667;
	b = 5.8;
	p = r[7];
	q = 1 - r[7];
	cdfbet(&ind,&p,&q,&x,&y,&a,&b,0,0);
	BVtransitionF.CorrectRxPreSM = x;
	BVparamsLogL.out[CurrSim-1][7] = x;

	// Calculate remaining elements of the CtsTransition matrix and AveDuration for BV
	/*BVtransitionF.CtsTransition[1][2] = BVtransitionF.Incidence1 *
		BVtransitionF.SymptomaticPropn;
	BVtransitionF.CtsTransition[1][3] = BVtransitionF.Incidence1 *
		(1.0 - BVtransitionF.SymptomaticPropn);
	BVtransitionF.AveDuration[0] = 1.0/BVtransitionF.CtsTransition[0][1];
	BVtransitionF.AveDuration[1] = 1.0/(BVtransitionF.CtsTransition[1][0] +
		BVtransitionF.CtsTransition[1][2] + BVtransitionF.CtsTransition[1][3]);
	BVtransitionF.AveDuration[2] = 1.0/(BVtransitionF.CtsTransition[2][0] +
		BVtransitionF.CtsTransition[2][1]);
	BVtransitionF.AveDuration[3] = 1.0/(BVtransitionF.CtsTransition[3][0] +
		BVtransitionF.CtsTransition[3][1]);*/
}

void SimulateVCparameters()
{
	int ind, i, is;
	double x, y, a, b, p, q;
	double r[6]; // Random variables from U(0, 1)

	ind = 2;

	if (FixedUncertainty == 0){
		int seed = (CurrSim - 1) * 91 + process_num * 7927;
		CRandomMersenne rg(seed);
		for (i = 0; i<6; i++){
			r[i] = rg.Random();
			RandomUniformVC.out[CurrSim - 1][i] = r[i];
		}
	}
	else{
		for (i = 0; i<6; i++){
			r[i] = RandomUniformVC.out[CurrSim - 1][i];
		}
	}

	// Simulate the ave dur of symptomatic VC (mean 12, SD 3)
	a = 16.0;
	b = 1.3333;
	p = r[0];
	q = 1 - r[0];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	VCtransitionF.AveDuration[1] = x;
	VCparamsLogL.out[CurrSim-1][0] = x;

	// Simulate the ave dur of asymptomatic VC (mean 26, SD 6)
	a = 18.778;
	b = 0.7222;
	p = r[1];
	q = 1 - r[1];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	VCtransitionF.AveDuration[0] = x;
	VCparamsLogL.out[CurrSim-1][1] = x;

	// Simulate the incidence of Candida colonization (mean 0.8, SD 0.2)
	a = 16.0;
	b = 20.0;
	p = r[2];
	q = 1 - r[2];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	VCtransitionF.Incidence = x;
	VCparamsLogL.out[CurrSim-1][2] = x;

	// Simulate the incidence of symptomatic VC in all asymp women (mean 0.15, SD 0.05)
	a = 9.0;
	b = 60.0;
	p = r[3];
	q = 1 - r[3];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	VCtransitionF.RecurrenceRate = x * (VCtransitionF.Incidence + 1.0/
		VCtransitionF.AveDuration[0])/VCtransitionF.Incidence;
	VCparamsLogL.out[CurrSim-1][3] = x;

	// Simulate the std deviation of the study effect
	/*a = 4.0;
	b = 13.333;
	p = r[4];
	q = 1 - r[4];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	VCtransitionF.HouseholdLogL.VarStudyEffect = pow(x, 2.0);
	VCtransitionF.ANClogL.VarStudyEffect = pow(x, 2.0);
	VCtransitionF.FPClogL.VarStudyEffect = pow(x, 2.0);
	VCtransitionF.CSWlogL.VarStudyEffect = pow(x, 2.0);
	VCparamsLogL.out[CurrSim-1][4] = pow(x, 2.0);*/

	// Simulate the proportion of VC cases correctly treated (mean 0.5, SD 0.2)
	a = 2.625;
	b = 2.625;
	p = r[5];
	q = 1 - r[5];
	cdfbet(&ind,&p,&q,&x,&y,&a,&b,0,0);
	VCtransitionF.CorrectRxPreSM = x;
	VCtransitionF.CorrectRxWithSM = x;
	VCparamsLogL.out[CurrSim-1][5] = x;
}

void SimulateHCTparameters()
{
	int ind, i, iy;
	double x, y, a, b, p, q;
	double r[15]; // Random variables from U(0, 1)

	ind = 2;

	for (i = 0; i<15; i++){ r[i] = RandomUniformHCT.out[CurrSim - 1][i]; }

	// Simulate the RR retesting in diagnosed (for new testing modalities)
	RetestAdjDiagnosed[3] = r[0]; // HBCT
	RetestAdjDiagnosed[4] = r[0]; // Mobile
	for (i = 6; i <= 11; i++){ RetestAdjDiagnosed[i] = r[0]; }

	// Simulate the RR retesting on ART (for new testing modalities)
	a = 2.2;
	b = 3.911;
	p = r[1];
	q = 1 - r[1];
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
	RetestAdjART[3] = x; // HBCT
	RetestAdjART[4] = x; // Mobile
	for (i = 6; i <= 11; i++){ RetestAdjART[i] = x; }

	// Simulate the prob of referral for testing in context of assisted partner notification
	a = 18.37;
	b = 6.122;
	p = r[2];
	q = 1 - r[2];
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
	ProbReferral[1] = (1.0 / (1.0 + ((1.0 - 0.3) / 0.3) / x)) / 0.71;
	if (ProbReferral[1]>1.0){ ProbReferral[1] = 1.0; }

	// Simulate the average uptake of testing through home-based testing
	a = 29.3;
	b = 12.56;
	p = r[3];
	q = 1 - r[3];
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
	HBCTuptake[0] = x * 0.80;
	HBCTuptake[1] = x * 1.20;

	// Simulate the annual uptake of testing in family planning clinics
	a = 20.25;
	b = 45.00;
	p = r[4];
	q = 1 - r[4];
	cdfgam(&ind,&p,&q,&x,&a,&b,0,0);
	FPCtestUptake = x;

	// Simulate the annual uptake of testing in mobile clinics
	a = 7.563;
	b = 137.5;
	p = r[5];
	q = 1 - r[5];
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
	MobileTestUptake[0] = x;
	MobileTestUptake[1] = x * 2.0;

	// Simulate the annual uptake of testing in MSM
	a = 4.0;
	b = 10.0;
	p = r[6];
	q = 1 - r[6];
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
	if (MSMtestUptake[35] > 0.0){
		for (iy = 35; iy < 81; iy++){
			MSMtestUptake[iy] = x; }
	}

	// Simulate the annual uptake of testing in FSW
	a = 1.823;
	b = 1.35;
	p = r[7];
	q = 1 - r[7];
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
	if (FSWtestUptake[35] > 0.0){
		for (iy = 35; iy < 81; iy++){
			FSWtestUptake[iy] = x; }
	}

	// Simulate the RR testing in virgins in schools
	a = 2.625;
	b = 2.625;
	p = r[8];
	q = 1 - r[8];
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
	SchoolTestVirginRR = x;
	
	// Simulate the fraction of the employed population that can be reached by workplace testing
	a = 1.644;
	b = 3.837;
	p = r[9];
	q = 1 - r[9];
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
	EmployedReachable = x;

	// Simulate ANC partner uptake after invitation letter (assuming OR of 4 for married)
	a = 6.966;
	b = 14.14;
	p = r[10];
	q = 1 - r[10];
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
	ANCpartnerTested[0][1] = x;
	ANCpartnerTested[1][1] = x;
	ANCpartnerTested[0][0] = 1.0 / (1.0 + 4.0 * (1.0 - x) / x);
	ANCpartnerTested[1][0] = 1.0 / (1.0 + 4.0 * (1.0 - x) / x);

	// Simulate effect of self-testing on odds of testing
	a = 66.75;
	b = 8.17;
	p = r[11];
	q = 1 - r[11];
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, 0);
	OR_STuptake = x;
	iy = CurrYear - StartYear;
	if (PropnSTreferral[iy] > 0.0){
		ANCpartnerTested[0][1] = PropnSTreferral[iy] / (1.0 + (1.0 - ANCpartnerTested[0][1]) / (OR_STuptake * ANCpartnerTested[0][1])) +
			(1.0 - PropnSTreferral[iy]) * ANCpartnerTested[0][1];
		ANCpartnerTested[1][1] = PropnSTreferral[iy] / (1.0 + (1.0 - ANCpartnerTested[1][1]) / (OR_STuptake * ANCpartnerTested[1][1])) +
			(1.0 - PropnSTreferral[iy]) * ANCpartnerTested[1][1];
		ANCpartnerTested[0][0] = PropnSTreferral[iy] / (1.0 + (1.0 - ANCpartnerTested[0][0]) / (OR_STuptake * ANCpartnerTested[0][0])) +
			(1.0 - PropnSTreferral[iy]) * ANCpartnerTested[0][0];
		ANCpartnerTested[1][0] = PropnSTreferral[iy] / (1.0 + (1.0 - ANCpartnerTested[1][0]) / (OR_STuptake * ANCpartnerTested[1][0])) +
			(1.0 - PropnSTreferral[iy]) * ANCpartnerTested[1][0];
	}

	// Simulate rate of linkage associated with community-based testing
	a = 5.896;
	b = 2.775;
	p = r[12];
	q = 1 - r[12];
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, 0);
	RR_ARTstartCommunity = x;
}

void SimulateStructParams()
{
	int SimCount2, seed, ii, ind, is;
	double x, y, a, b, p, q, bound;
	double r[24];

	SimCount2 = (CurrSim - 1) / IterationsPerPC;
	if (FixedUncertainty == 0){
		seed = SimCount2 * 91 + process_num * 7927;
		CRandomMersenne rg(seed);
		for (ii = 0; ii < 24; ii++){
			r[ii] = rg.Random();
		}
	}
	else{
		SimCount2 = (CurrSim - 1) / IterationsPerPC;
		for (ii = 0; ii < 24; ii++){
			r[ii] = RandomUniformStruct.out[SimCount2][ii];
		}
	}

	// Store the parameter combination
	if (CurrSim == (SimCount2 * IterationsPerPC + 1)){
		for (ii = 0; ii < 23; ii++){
			StructParamsLogL.out[SimCount2][ii] = r[ii];
		}
	}

	// Simulate OR of condom use per day of binge drinking per week (beta hurdle)
	/*if (r[0] < 0.5){
		r[0] *= 1.0 / 0.5;
		a = 3.05;
		b = 1.64;
		p = r[0];
		q = 1 - r[0];
		bound = 0.0;
		ind = 2;
		cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
		ORcondomBingePW = x;
	}
	else{ ORcondomBingePW = 1.0; }

	// Simulate RR for entry into casual sex for women who binge drink (gamma hurdle)
	if (r[1] < 0.5){ RRcasualBinge[1] = 1.0; }
	else{ 
		r[1] = (r[1] - 0.5) / 0.5;
		a = 5.29;
		b = 2.30;
		p = r[1];
		q = 1 - r[1];
		bound = 0.0;
		ind = 2;
		cdfgam(&ind, &p, &q, &x, &a, &b, 0, &bound);
		RRcasualBinge[1] = x + 1.0;
	}

	// Simulate RR for entry into casual sex for men who binge drink (gamma hurdle)
	if (r[2] < 0.5){ RRcasualBinge[0] = 1.0; }
	else{
		r[2] = (r[2] - 0.5) / 0.5;
		a = 6.25;
		b = 8.33;
		p = r[2];
		q = 1 - r[2];
		bound = 0.0;
		ind = 2;
		cdfgam(&ind, &p, &q, &x, &a, &b, 0, &bound);
		RRcasualBinge[0] = x + 1.0;
	}

	// Simulate prob of confounding between drinks per DD and risk group
	ConfoundingAlcSex = r[3];

	// RR binge drinking after single-session alcohol counselling intervention
	a = 12.0;
	b = 3.00;
	p = r[4];
	q = 1 - r[4];
	bound = 0.0;
	ind = 2;
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
	RRalcoholSingle = x;

	// RR binge drinking after multi-session alcohol counselling intervention
	a = 2.85;
	b = 2.33;
	p = r[5];
	q = 1 - r[5];
	bound = 0.0;
	ind = 2;
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
	RRalcoholMultiple = x;

	// Annual prob of reverting to alcohol consumption before intervention
	ProbAlcoholReversion = r[6];*/

	// ln(OR) for consistent condom use per year of education (gamma hurdle)
	if (r[7] < 0.5){ ORcondomPerYrSchool = 1.0; }
	else{
		r[7] = (r[7] - 0.5) / 0.5;
		a = 4.00;
		b = 40.0;
		p = r[7];
		q = 1 - r[7];
		bound = 0.0;
		ind = 2;
		cdfgam(&ind, &p, &q, &x, &a, &b, 0, &bound);
		ORcondomPerYrSchool = exp(x);
	}
	CondomEduMedian = pow(ORcondomPerYrSchool / 1.165, 1 / 3.5);

	// Factor by which time to behaviour change reduces per year of schooling (beta hurdle)
	/*if (r[8] < 0.25){
		r[8] *= 1.0 / 0.25;
		a = 91.2;
		b = 3.80;
		p = r[8];
		q = 1 - r[8];
		bound = 0.0;
		ind = 2;
		cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
		CondomEduMedian = x;
	}
	else{ CondomEduMedian = 1.0; }*/

	// RR of sexual debut if in school (boys) - beta hurdle
	if (r[9] < 0.50){
		r[9] *= 1.0 / 0.50;
		a = 2.12;
		b = 1.41;
		p = r[9];
		q = 1 - r[9];
		bound = 0.0;
		ind = 2;
		cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
		RRdebutInSchool[0] = x;
	}
	else{ RRdebutInSchool[0] = 1.0; }

	// RR of entry into casual sex for men who are employed (gamma hurdle)
	if (r[10] < 0.50){ RRcasualEmployedM = 1.0; }
	else{
		r[10] = (r[10] - 0.50) / 0.50;
		a = 2.25;
		b = 3.57;
		p = r[10];
		q = 1 - r[10];
		bound = 0.0;
		ind = 2;
		cdfgam(&ind, &p, &q, &x, &a, &b, 0, &bound);
		RRcasualEmployedM = x + 1.0;
	}

	// RR of FSW contact for men who are employed (gamma hurdle)
	if (r[11] < 0.50){ RR_FSWcontactEmployedM = 1.0; }
	else{
		r[11] = (r[11] - 0.50) / 0.50;
		a = 4.00;
		b = 8.00;
		p = r[11];
		q = 1 - r[11];
		bound = 0.0;
		ind = 2;
		cdfgam(&ind, &p, &q, &x, &a, &b, 0, &bound);
		RR_FSWcontactEmployedM = x + 1.0;
	}

	// RR of entry into casual sex for women per log reduction in household per capita income (gamma hurdle)
	if (r[12] < 0.50){ RRcasualLogDropIncomeF = 1.00; }
	else{
		r[12] = (r[12] - 0.50) / 0.50;
		a = 3.27;
		b = 3.85;
		p = r[12];
		q = 1 - r[12];
		bound = 0.0;
		ind = 2;
		cdfgam(&ind, &p, &q, &x, &a, &b, 0, &bound);
		RRcasualLogDropIncomeF = x + 1.0;
	}

	// RR of ST partnership formation in men who are employed (gamma hurdle)
	if (r[13] < 0.25){ RR_STemployedM = 1.0; }
	else{ 
		r[13] = (r[13] - 0.25) / 0.75;
		a = 2.96;
		b = 6.88;
		p = r[13];
		q = 1 - r[13];
		bound = 0.0;
		ind = 2;
		cdfgam(&ind, &p, &q, &x, &a, &b, 0, &bound);
		RR_STemployedM = x + 1.0;
	}

	// RR of entry into marriage while in school (beta hurdle)
	if (r[14] < 0.90){
		r[14] *= 1.0 / 0.90;
		a = 1.28;
		b = 2.98;
		p = r[14];
		q = 1 - r[14];
		bound = 0.0;
		ind = 2;
		cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
		ORmarriageInSchool = x;
	}
	else{ ORmarriageInSchool = 1.0; }

	// RR of sexual debut for females per log reduction in household per capita income (gamma hurdle)
	if (r[15] < 0.10){ RRdebutLogDropIncomeF = 1.0; }
	else{
		r[15] = (r[15] - 0.10) / 0.90;
		a = 2.39;
		b = 7.02;
		p = r[15];
		q = 1 - r[15];
		bound = 0.0;
		ind = 2;
		cdfgam(&ind, &p, &q, &x, &a, &b, 0, &bound);
		RRdebutLogDropIncomeF = x + 1.0;
	}

	// RR school dropout if receiving school support
	a = 3.02;
	b = 1.42;
	p = r[16];
	q = 1 - r[16];
	bound = 0.0;
	ind = 2;
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
	RRdropoutSupport = x;

	// OR unemployment if receiving vocational training/microfinance
	a = 8.97;
	b = 1.34;
	p = r[17];
	q = 1 - r[17];
	bound = 0.0;
	ind = 2;
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
	ORunemployedTraining = x;

	// RR of VMMC per log increase in household income
	// NB: the parameter we simulate is the proportional increase in the rate of VMMC.
	if (r[18] < 0.25){ RR_VMMClogIncome = 1.0; }
	else{
		r[18] = (r[18] - 0.25) / 0.75;
		a = 4.00;
		b = 20.0;
		p = r[18];
		q = 1 - r[18];
		bound = 0.0;
		ind = 2;
		cdfgam(&ind, &p, &q, &x, &a, &b, 0, &bound);
		RR_VMMClogIncome = exp(x);
	}

	// RR of sexual debut if in school (girls)
	if (r[19] < 0.90){
		r[19] *= 1.0 / 0.90;
		a = 1.41;
		b = 2.12;
		p = r[19];
		q = 1 - r[19];
		bound = 0.0;
		ind = 2;
		cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
		RRdebutInSchool[1] = x;
	}
	else{ RRdebutInSchool[1] = 1.0; }

	// RR school dropout if for R800 income support
	a = 9.49;
	b = 1.55;
	p = r[20];
	q = 1 - r[20];
	bound = 0.0;
	ind = 2;
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
	RRdropoutIncSupport = x;

	// Effect of inequitable gender norms on condom use
	/*if (r[15] < 0.50){
		r[15] *= 1.0 / 0.50;
		a = 1.83;
		b = 5.50;
		p = r[15];
		q = 1 - r[15];
		bound = 0.0;
		ind = 2;
		cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
		ORcondomGenderIneq = x;
	}
	else{ ORcondomGenderIneq = 1.0; }

	// Effect of inequitable gender norms on entry into casual sex
	if (r[16] < 0.50){
		r[16] *= 1.0 / 0.50;
		a = 2.77;
		b = 0.92;
		p = r[16];
		q = 1 - r[16];
		bound = 0.0;
		ind = 2;
		cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
		EffectIneqGenderCasual = x;
	}
	else{ EffectIneqGenderCasual = 1.0; }

	// Effect of high-risk group on inequitable gender norms
	a = 2.63;
	b = 2.63;
	p = r[17];
	q = 1 - r[17];
	bound = 0.0;
	ind = 2;
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
	HighRiskEffectIneqGender = 1.0 / x;

	// Increase in concurrency rate if endorse ineq gender norms (gamma hurdle)
	if (r[18] < 0.50){ EffectIneqGenderConcurrency = 0.0; }
	else{
		r[18] = (r[18] - 0.50) / 0.50;
		a = 4.00;
		b = 0.80;
		p = r[18];
		q = 1 - r[18];
		bound = 0.0;
		ind = 2;
		cdfgam(&ind, &p, &q, &x, &a, &b, 0, &bound);
		EffectIneqGenderConcurrency = x;
	}

	// Effect of endorsing inequitable gender norms on # drinks per drinking day (gamma hurdle)
	if (r[19] < 0.25){ GenderIneqEffectDrinksPerDD = 0.0; }
	else{
		r[19] = (r[19] - 0.25) / 0.75;
		a = 4.34;
		b = 0.52;
		p = r[19];
		q = 1 - r[19];
		bound = 0.0;
		ind = 2;
		cdfgam(&ind, &p, &q, &x, &a, &b, 0, &bound);
		GenderIneqEffectDrinksPerDD = x;
	}

	// RR of endoring inequitable gender norms after participating in gender-transformative interventions, at individual level
	RRgenderIneqIndiv = r[20];

	// RR of endoring inequitable gender norms after participating in gender-transformative interventions, at community level
	RRgenderIneqComm = r[21];

	// Annual probability of reverting to pre-intervention gender norms
	ProbGenderIneqReversion = r[22];*/
}

void SaveBaselinePop()
{
	// There is a fair amount of white space in this function, which is deliberate. This is so that if you paste this function
	// and the corresponding section of the header file into a spreadsheet (in different columns), you can easily see where
	// variables have been omitted.

	int ic, ih, ii, ir, ig, CurrPop, CurrHHs;
	HouseholdGroup HouseA(1); // Argument is arbitrary

	memset(&HouseA, 0, sizeof(HouseA));

	CurrPop = Register.size();
	TempRegister.resize(CurrPop);
	for (ic = 0; ic < CurrPop; ic++){
		TempRegister[ic].AliveInd = Register[ic].AliveInd; 
		TempRegister[ic].SexInd = Register[ic].SexInd;
		TempRegister[ic].RiskGroup = Register[ic].RiskGroup;
		TempRegister[ic].PopWeight = Register[ic].PopWeight;
		TempRegister[ic].DOB = Register[ic].DOB;
		TempRegister[ic].DOD = Register[ic].DOD;
		TempRegister[ic].DOLC = Register[ic].DOLC;
		TempRegister[ic].DOLB = Register[ic].DOLB;
		TempRegister[ic].DOLW = Register[ic].DOLW;
		TempRegister[ic].DOAI = Register[ic].DOAI;
		TempRegister[ic].DOUD = Register[ic].DOUD;
		TempRegister[ic].Fecundability = Register[ic].Fecundability;
		TempRegister[ic].Conscientiousness = Register[ic].Conscientiousness;
		TempRegister[ic].AgeGroup = Register[ic].AgeGroup;

		TempRegister[ic].CurrAge = Register[ic].CurrAge;		
		TempRegister[ic].Race = Register[ic].Race;
		TempRegister[ic].HouseholdID = Register[ic].HouseholdID;
		TempRegister[ic].DateHomeless = Register[ic].DateHomeless;
		for (ii = 0; ii < 2; ii++){ TempRegister[ic].ParentID[ii] = Register[ic].ParentID[ii]; }
		for (ii = 0; ii < 20; ii++){ TempRegister[ic].ChildIDs[ii] = Register[ic].ChildIDs[ii]; }
		TempRegister[ic].FatherBaby = Register[ic].FatherBaby;
		TempRegister[ic].CircInd = Register[ic].CircInd;
		TempRegister[ic].HighestGrade = Register[ic].HighestGrade;
		TempRegister[ic].InSchool = Register[ic].InSchool;
		TempRegister[ic].CurrContr = Register[ic].CurrContr;
		TempRegister[ic].PrevContr = Register[ic].PrevContr;
		TempRegister[ic].CurrUrban = Register[ic].CurrUrban;
		TempRegister[ic].Imprisoned = Register[ic].Imprisoned;
		TempRegister[ic].ReleaseDate = Register[ic].ReleaseDate;
		TempRegister[ic].Employed = Register[ic].Employed;
		TempRegister[ic].LogIncomeDif = Register[ic].LogIncomeDif;
		TempRegister[ic].PrivatePension = Register[ic].PrivatePension;
		TempRegister[ic].MarriedInd = Register[ic].MarriedInd;
		TempRegister[ic].FSWind = Register[ic].FSWind;
		TempRegister[ic].CasualInd = Register[ic].CasualInd;
		TempRegister[ic].HetCasualInd = Register[ic].HetCasualInd;
		TempRegister[ic].VirginInd = Register[ic].VirginInd;
		TempRegister[ic].IDprimary = Register[ic].IDprimary;
		TempRegister[ic].ID2ndary = Register[ic].ID2ndary;
		
		
		TempRegister[ic].CurrPartners = Register[ic].CurrPartners;
		TempRegister[ic].LifetimePartners = Register[ic].LifetimePartners;
		TempRegister[ic].AnnPartners = Register[ic].AnnPartners;
		TempRegister[ic].IneqGender = Register[ic].IneqGender;
		TempRegister[ic].NonHIVmortProb = Register[ic].NonHIVmortProb;
		TempRegister[ic].NonHIVfertRate = Register[ic].NonHIVfertRate;
		TempRegister[ic].CumSelectionProb = Register[ic].CumSelectionProb;
		TempRegister[ic].MalePref = Register[ic].MalePref;
		TempRegister[ic].ChangeMalePref = Register[ic].ChangeMalePref;
		TempRegister[ic].InsertivePref = Register[ic].InsertivePref;
		TempRegister[ic].CD4 = Register[ic].CD4;
		TempRegister[ic].BaselineCD4 = Register[ic].BaselineCD4;
		TempRegister[ic].logVL = Register[ic].logVL;
		TempRegister[ic].SPVL = Register[ic].SPVL;
		TempRegister[ic].HIVstage = Register[ic].HIVstage;
		TempRegister[ic].HIVstageE = Register[ic].HIVstageE;
		TempRegister[ic].CTstage = Register[ic].CTstage;
		TempRegister[ic].CTstageE = Register[ic].CTstageE;
		TempRegister[ic].HDstage = Register[ic].HDstage;
		TempRegister[ic].HDstageE = Register[ic].HDstageE;
		TempRegister[ic].NGstage = Register[ic].NGstage;
		TempRegister[ic].NGstageE = Register[ic].NGstageE;
		TempRegister[ic].TVstage = Register[ic].TVstage;
		TempRegister[ic].TVstageE = Register[ic].TVstageE;
		TempRegister[ic].TPstage = Register[ic].TPstage;
		TempRegister[ic].TPstageE = Register[ic].TPstageE;
		TempRegister[ic].HSVstage = Register[ic].HSVstage;
		
		TempRegister[ic].HSVstageE = Register[ic].HSVstageE;
		TempRegister[ic].BVstage = Register[ic].BVstage;
		TempRegister[ic].BVstageE = Register[ic].BVstageE;
		TempRegister[ic].VCstage = Register[ic].VCstage;
		TempRegister[ic].VCstageE = Register[ic].VCstageE;
		TempRegister[ic].SuscepHIVadj = Register[ic].SuscepHIVadj;
		TempRegister[ic].DateInfect = Register[ic].DateInfect;
		TempRegister[ic].PartnerRateAdj = Register[ic].PartnerRateAdj;
		
		TempRegister[ic].DesiredNewPartners = Register[ic].DesiredNewPartners;
		TempRegister[ic].CondomPref = Register[ic].CondomPref;
		TempRegister[ic].TempCondomAdj = Register[ic].TempCondomAdj;
		TempRegister[ic].CondomPrimary = Register[ic].CondomPrimary;
		TempRegister[ic].Condom2ndary = Register[ic].Condom2ndary;
		TempRegister[ic].DisclosedPrimary = Register[ic].DisclosedPrimary;
		TempRegister[ic].Disclosed2ndary = Register[ic].Disclosed2ndary;
		TempRegister[ic].NewStatus = Register[ic].NewStatus;
		TempRegister[ic].UVIprimary = Register[ic].UVIprimary;
		TempRegister[ic].PVIprimary = Register[ic].PVIprimary;
		TempRegister[ic].UVI2ndary = Register[ic].UVI2ndary;
		TempRegister[ic].PVI2ndary = Register[ic].PVI2ndary;
		TempRegister[ic].UVICSW = Register[ic].UVICSW;
		TempRegister[ic].PVICSW = Register[ic].PVICSW;
		TempRegister[ic].UAIcasual = Register[ic].UAIcasual;
		TempRegister[ic].PAIcasual = Register[ic].PAIcasual;
		TempRegister[ic].IDofCSW = Register[ic].IDofCSW;
		TempRegister[ic].IDofCasual = Register[ic].IDofCasual;
		TempRegister[ic].VCThistory = Register[ic].VCThistory;
		TempRegister[ic].Visiting = Register[ic].Visiting;
		TempRegister[ic].VisitFreq = Register[ic].VisitFreq;
		TempRegister[ic].OnPrEP = Register[ic].OnPrEP;
		TempRegister[ic].VaccineEff = Register[ic].VaccineEff;
		TempRegister[ic].Date1stVacc = Register[ic].Date1stVacc;
		TempRegister[ic].DailyDrinkProb = Register[ic].DailyDrinkProb;
		TempRegister[ic].DrinksPerDD = Register[ic].DrinksPerDD;
		TempRegister[ic].DrinkProbConstant = Register[ic].DrinkProbConstant;
		TempRegister[ic].DrinksPerDDconstant = Register[ic].DrinksPerDDconstant;
		
		
		TempRegister[ic].EverMSM = Register[ic].EverMSM;
		TempRegister[ic].EverBi = Register[ic].EverBi;
		TempRegister[ic].RecentMSM = Register[ic].RecentMSM;
		TempRegister[ic].RecentBi = Register[ic].RecentBi;
		
		TempRegister[ic].EverInjectable = Register[ic].EverInjectable;
		TempRegister[ic].EverPill = Register[ic].EverPill;
		TempRegister[ic].PrevImprisoned = Register[ic].PrevImprisoned;
		TempRegister[ic].BirthInterval = Register[ic].BirthInterval;
		TempRegister[ic].DateMig = Register[ic].DateMig;
		TempRegister[ic].EverCasual = Register[ic].EverCasual;
		TempRegister[ic].CasualLast3mo = Register[ic].CasualLast3mo;
		TempRegister[ic].PartnersLast3mo = Register[ic].PartnersLast3mo;
	}

	CurrHHs = HHregister.size();
	for (ih = 0; ih < CurrHHs; ih++){
		HouseA.IDhead = HHregister[ih].IDhead;
		for (ii = 0; ii < MaxHHsize; ii++){ HouseA.Members[ii] = HHregister[ih].Members[ii]; }
		HouseA.Size = HHregister[ih].Size;
		HouseA.Active = HHregister[ih].Active;
		HouseA.Urban = HHregister[ih].Urban;
		HouseA.PerCapitaIncome = HHregister[ih].PerCapitaIncome;
		HouseA.PerCapitaIncomeAdj = HHregister[ih].PerCapitaIncomeAdj;
		TempHHregister.push_back(HouseA);
	}

	for (ii = 0; ii < 3; ii++){
		for (ic = 0; ic < MaxCSWs; ic++){
			TempRSApop.CSWregister[ic][ii] = RSApop.CSWregister[ic][ii];
		}
		for (ic = 0; ic < MaxCasual; ic++){
			TempRSApop.CasualRegisterM[ic][ii] = RSApop.CasualRegisterM[ic][ii];
			TempRSApop.CasualRegisterF[ic][ii] = RSApop.CasualRegisterF[ic][ii];
		}
	}
	for (ic = 0; ic < MaxCasual; ic++){
		TempRSApop.CasualRegister[ic] = RSApop.CasualRegister[ic];
	}

	TempTotCurrCasual = TotCurrCasual;
	for (ir = 0; ir < 3; ir++){
		TempTotCurrFSW[ir] = TotCurrFSW[ir];
		for (ig = 0; ig < 2; ig++){
			TempTotCurrCasualHet[ir][ig] = TotCurrCasualHet[ir][ig];
		}
	}
}

void RestartFromBaseline()
{
	// There is a fair amount of white space in this function, which is deliberate. This is so that if you paste this function
	// and the corresponding section of the header file into a spreadsheet (in different columns), you can easily see where
	// variables have been omitted.

	int ic, ih, ii, ir, ig, CurrPop, CurrHHs;
	HouseholdGroup HouseA(1); // Argument is arbitrary

	memset(&HouseA, 0, sizeof(HouseA));
	
	CurrPop = TempRegister.size();
	Register.resize(CurrPop);
	for (ic = 0; ic < CurrPop; ic++){
		Register[ic].AliveInd = TempRegister[ic].AliveInd;
		Register[ic].SexInd = TempRegister[ic].SexInd;
		Register[ic].RiskGroup = TempRegister[ic].RiskGroup;
		Register[ic].PopWeight = TempRegister[ic].PopWeight;
		Register[ic].DOB = TempRegister[ic].DOB;
		Register[ic].DOD = TempRegister[ic].DOD;
		Register[ic].DOLC = TempRegister[ic].DOLC;
		Register[ic].DOLB = TempRegister[ic].DOLB;
		Register[ic].DOLW = TempRegister[ic].DOLW;
		Register[ic].DOAI = TempRegister[ic].DOAI;
		Register[ic].DOUD = TempRegister[ic].DOUD;
		Register[ic].Fecundability = TempRegister[ic].Fecundability;
		Register[ic].Conscientiousness = TempRegister[ic].Conscientiousness;
		Register[ic].AgeGroup = TempRegister[ic].AgeGroup;

		Register[ic].CurrAge = TempRegister[ic].CurrAge;
		Register[ic].Race = TempRegister[ic].Race;
		Register[ic].HouseholdID = TempRegister[ic].HouseholdID;
		Register[ic].DateHomeless = TempRegister[ic].DateHomeless;
		for (ii = 0; ii < 2; ii++){Register[ic].ParentID[ii] = TempRegister[ic].ParentID[ii]; }
		for (ii = 0; ii < 20; ii++){Register[ic].ChildIDs[ii] = TempRegister[ic].ChildIDs[ii]; }
		Register[ic].FatherBaby = TempRegister[ic].FatherBaby;
		Register[ic].CircInd = TempRegister[ic].CircInd;
		Register[ic].HighestGrade = TempRegister[ic].HighestGrade;
		Register[ic].InSchool = TempRegister[ic].InSchool;
		Register[ic].CurrContr = TempRegister[ic].CurrContr;
		Register[ic].PrevContr = TempRegister[ic].PrevContr;
		Register[ic].CurrUrban = TempRegister[ic].CurrUrban;
		Register[ic].Imprisoned = TempRegister[ic].Imprisoned;
		Register[ic].ReleaseDate = TempRegister[ic].ReleaseDate;
		Register[ic].Employed = TempRegister[ic].Employed;
		Register[ic].LogIncomeDif = TempRegister[ic].LogIncomeDif;
		Register[ic].PrivatePension = TempRegister[ic].PrivatePension;
		Register[ic].MarriedInd = TempRegister[ic].MarriedInd;
		Register[ic].FSWind = TempRegister[ic].FSWind;
		Register[ic].CasualInd = TempRegister[ic].CasualInd;
		Register[ic].HetCasualInd = TempRegister[ic].HetCasualInd;
		Register[ic].VirginInd = TempRegister[ic].VirginInd;
		Register[ic].IDprimary = TempRegister[ic].IDprimary;
		Register[ic].ID2ndary = TempRegister[ic].ID2ndary;


		Register[ic].CurrPartners = TempRegister[ic].CurrPartners;
		Register[ic].LifetimePartners = TempRegister[ic].LifetimePartners;
		Register[ic].AnnPartners = TempRegister[ic].AnnPartners;
		Register[ic].IneqGender = TempRegister[ic].IneqGender;
		Register[ic].NonHIVmortProb = TempRegister[ic].NonHIVmortProb;
		Register[ic].NonHIVfertRate = TempRegister[ic].NonHIVfertRate;
		Register[ic].CumSelectionProb = TempRegister[ic].CumSelectionProb;
		Register[ic].MalePref = TempRegister[ic].MalePref;
		Register[ic].ChangeMalePref = TempRegister[ic].ChangeMalePref;
		Register[ic].InsertivePref = TempRegister[ic].InsertivePref;
		Register[ic].CD4 = TempRegister[ic].CD4;
		Register[ic].BaselineCD4 = TempRegister[ic].BaselineCD4;
		Register[ic].logVL = TempRegister[ic].logVL;
		Register[ic].SPVL = TempRegister[ic].SPVL;
		Register[ic].HIVstage = TempRegister[ic].HIVstage;
		Register[ic].HIVstageE = TempRegister[ic].HIVstageE;
		Register[ic].CTstage = TempRegister[ic].CTstage;
		Register[ic].CTstageE = TempRegister[ic].CTstageE;
		Register[ic].HDstage = TempRegister[ic].HDstage;
		Register[ic].HDstageE = TempRegister[ic].HDstageE;
		Register[ic].NGstage = TempRegister[ic].NGstage;
		Register[ic].NGstageE = TempRegister[ic].NGstageE;
		Register[ic].TVstage = TempRegister[ic].TVstage;
		Register[ic].TVstageE = TempRegister[ic].TVstageE;
		Register[ic].TPstage = TempRegister[ic].TPstage;
		Register[ic].TPstageE = TempRegister[ic].TPstageE;
		Register[ic].HSVstage = TempRegister[ic].HSVstage;

		Register[ic].HSVstageE = TempRegister[ic].HSVstageE;
		Register[ic].BVstage = TempRegister[ic].BVstage;
		Register[ic].BVstageE = TempRegister[ic].BVstageE;
		Register[ic].VCstage = TempRegister[ic].VCstage;
		Register[ic].VCstageE = TempRegister[ic].VCstageE;
		Register[ic].SuscepHIVadj = TempRegister[ic].SuscepHIVadj;
		Register[ic].DateInfect = TempRegister[ic].DateInfect;
		Register[ic].PartnerRateAdj = TempRegister[ic].PartnerRateAdj;

		Register[ic].DesiredNewPartners = TempRegister[ic].DesiredNewPartners;
		Register[ic].CondomPref = TempRegister[ic].CondomPref;
		Register[ic].TempCondomAdj = TempRegister[ic].TempCondomAdj;
		Register[ic].CondomPrimary = TempRegister[ic].CondomPrimary;
		Register[ic].Condom2ndary = TempRegister[ic].Condom2ndary;
		Register[ic].DisclosedPrimary = TempRegister[ic].DisclosedPrimary;
		Register[ic].Disclosed2ndary = TempRegister[ic].Disclosed2ndary;
		Register[ic].NewStatus = TempRegister[ic].NewStatus;
		Register[ic].UVIprimary = TempRegister[ic].UVIprimary;
		Register[ic].PVIprimary = TempRegister[ic].PVIprimary;
		Register[ic].UVI2ndary = TempRegister[ic].UVI2ndary;
		Register[ic].PVI2ndary = TempRegister[ic].PVI2ndary;
		Register[ic].UVICSW = TempRegister[ic].UVICSW;
		Register[ic].PVICSW = TempRegister[ic].PVICSW;
		Register[ic].UAIcasual = TempRegister[ic].UAIcasual;
		Register[ic].PAIcasual = TempRegister[ic].PAIcasual;
		Register[ic].IDofCSW = TempRegister[ic].IDofCSW;
		Register[ic].IDofCasual = TempRegister[ic].IDofCasual;
		Register[ic].VCThistory = TempRegister[ic].VCThistory;
		Register[ic].Visiting = TempRegister[ic].Visiting;
		Register[ic].VisitFreq = TempRegister[ic].VisitFreq;
		Register[ic].OnPrEP = TempRegister[ic].OnPrEP;
		Register[ic].VaccineEff = TempRegister[ic].VaccineEff;
		Register[ic].Date1stVacc = TempRegister[ic].Date1stVacc;
		Register[ic].DailyDrinkProb = TempRegister[ic].DailyDrinkProb;
		Register[ic].DrinksPerDD = TempRegister[ic].DrinksPerDD;
		Register[ic].DrinkProbConstant = TempRegister[ic].DrinkProbConstant;
		Register[ic].DrinksPerDDconstant = TempRegister[ic].DrinksPerDDconstant;


		Register[ic].EverMSM = TempRegister[ic].EverMSM;
		Register[ic].EverBi = TempRegister[ic].EverBi;
		Register[ic].RecentMSM = TempRegister[ic].RecentMSM;
		Register[ic].RecentBi = TempRegister[ic].RecentBi;

		Register[ic].EverInjectable = TempRegister[ic].EverInjectable;
		Register[ic].EverPill = TempRegister[ic].EverPill;
		Register[ic].PrevImprisoned = TempRegister[ic].PrevImprisoned;
		Register[ic].BirthInterval = TempRegister[ic].BirthInterval;
		Register[ic].DateMig = TempRegister[ic].DateMig;
		Register[ic].EverCasual = TempRegister[ic].EverCasual;
		Register[ic].CasualLast3mo = TempRegister[ic].CasualLast3mo;
		Register[ic].PartnersLast3mo = TempRegister[ic].PartnersLast3mo;
	}

	CurrHHs = TempHHregister.size();
	HHregister.clear();
	for (ih = 0; ih < CurrHHs; ih++){
		HouseA.IDhead = TempHHregister[ih].IDhead;
		for (ii = 0; ii < MaxHHsize; ii++){ HouseA.Members[ii] = TempHHregister[ih].Members[ii]; }
		HouseA.Size = TempHHregister[ih].Size;
		HouseA.Active = TempHHregister[ih].Active;
		HouseA.Urban = TempHHregister[ih].Urban;
		HouseA.PerCapitaIncome = TempHHregister[ih].PerCapitaIncome;
		HouseA.PerCapitaIncomeAdj = TempHHregister[ih].PerCapitaIncomeAdj;
		HHregister.push_back(HouseA);
	}

	for (ii = 0; ii < 3; ii++){
		for (ic = 0; ic < MaxCSWs; ic++){
			RSApop.CSWregister[ic][ii] = TempRSApop.CSWregister[ic][ii];
		}
		for (ic = 0; ic < MaxCasual; ic++){
			RSApop.CasualRegisterM[ic][ii] = TempRSApop.CasualRegisterM[ic][ii];
			RSApop.CasualRegisterF[ic][ii] = TempRSApop.CasualRegisterF[ic][ii];
		}
	}
	for (ic = 0; ic < MaxCasual; ic++){
		RSApop.CasualRegister[ic] = TempRSApop.CasualRegister[ic];
	}

	TotCurrCasual = TempTotCurrCasual;
	for (ir = 0; ir < 3; ir++){
		TotCurrFSW[ir] = TempTotCurrFSW[ir];
		for (ig = 0; ig < 2; ig++){
			TotCurrCasualHet[ir][ig] = TempTotCurrCasualHet[ir][ig];
		}
	}
}

void UpdateMalePref()
{
	int ic, ExactAge;

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 0 &&
			Register[ic].MalePref>0.0 && Register[ic].MalePref<1.0){
				ExactAge = CurrYear + 0.5 - Register[ic].DOB;
				if (ExactAge>20){
					Register[ic].MalePref += Register[ic].ChangeMalePref;
					if (Register[ic].MalePref <= 0.0){ Register[ic].MalePref = 0.001; }
					if (Register[ic].MalePref >= 1.0){ Register[ic].MalePref = 0.999; }
				}
		}
	}
}

void ResetRecentMSM()
{
	int ic, pID1, pID2;

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 0 && Register[ic].MalePref>0.0){
			Register[ic].RecentMSM = 0;
			Register[ic].RecentBi = 0;
			if (Register[ic].CurrPartners > 0){
				pID1 = Register[ic].IDprimary;
				if (Register[pID1-1].SexInd == 0){ Register[ic].RecentMSM = 1; }
				else{ Register[ic].RecentBi = 1; }
			}
			if (Register[ic].CurrPartners == 2){
				pID2 = Register[ic].ID2ndary;
				if (Register[pID2 - 1].SexInd == 0){ Register[ic].RecentMSM += 1; }
				else{ Register[ic].RecentBi = 1; }
			}
		}
	}
}

void GetMSMcalib()
{
	int ic, ia;
	double denom, married, recentBi, everBi, aged25plus, positive, currentReg, multPartners;
	double denomIns, denomRecep, posBi, posIns, posRecep, everMSM, positiveMSW;
	int pID1, pID2;
	double ExactAge;

	// Calculate outputs for men who've recently been active with other men (last 6 months)

	denom = 0.0;
	married = 0.0;
	recentBi = 0.0;
	everBi = 0.0;
	aged25plus = 0.0;
	positive = 0.0;
	currentReg = 0.0;
	multPartners = 0.0;
	denomIns = 0.0;
	denomRecep = 0.0;
	posBi = 0.0;
	posIns = 0.0;
	posRecep = 0.0;

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 0 &&
			Register[ic].MalePref>0.0 && Register[ic].VirginInd == 0){
				ExactAge = 0.5 + CurrYear - Register[ic].DOB;
				if (ExactAge >= 18.0 && Register[ic].RecentMSM>0){
					denom += Register[ic].PopWeight;
					if (Register[ic].MarriedInd == 1){ married += Register[ic].PopWeight; }
					if (Register[ic].RecentBi > 0){ recentBi += Register[ic].PopWeight; }
					if (Register[ic].EverBi == 1){ everBi += Register[ic].PopWeight; }
					if (ExactAge >= 25.0){ aged25plus += Register[ic].PopWeight; }
					if (Register[ic].HIVstage > 0){ positive += Register[ic].PopWeight; }
					if (Register[ic].CurrPartners > 0){
						pID1 = Register[ic].IDprimary;
						if (Register[pID1 - 1].SexInd == 0){ currentReg += Register[ic].PopWeight; }
						else{
							if (Register[ic].CurrPartners == 2){
								pID2 = Register[ic].ID2ndary;
								if (Register[pID2 - 1].SexInd == 0){ currentReg += Register[ic].PopWeight; }
							}
						}
					}
					if (Register[ic].RecentMSM > 1){ multPartners += Register[ic].PopWeight; }
					if (Register[ic].RecentBi > 0){ 
						if (Register[ic].HIVstage > 0){ posBi += Register[ic].PopWeight; }
					}
					if (Register[ic].InsertivePref == 1.0){
						denomIns += Register[ic].PopWeight;
						if (Register[ic].HIVstage > 0){ posIns += Register[ic].PopWeight; }
					}
					if (Register[ic].InsertivePref == 0.0){
						denomRecep += Register[ic].PopWeight;
						if (Register[ic].HIVstage > 0){ posRecep += Register[ic].PopWeight; }
					}
				}
		}
	}
	MSMmarried.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * married / denom;
	MSMrecentBi.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * recentBi / denom;
	MSMeverBi.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * everBi / denom;
	MSMaged25plus.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * aged25plus / denom;
	MSMpositive.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * positive / denom;
	MSMprev.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * positive / denom; // Same as MSMpositive
	MSMcurrentReg.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * currentReg / denom;
	MSMmultPartners.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * multPartners / denom;
	MSMpositiveBi.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * posBi / recentBi;
	MSMpositiveGay.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * (positive - posBi) / (denom - recentBi);
	MSMposInsertive.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * posIns / denomIns;
	MSMposReceptive.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * posRecep / denomRecep;
	MSMposVersatile.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * (positive - posIns - posRecep) / 
		(denom - denomIns - denomRecep);

	// Calculate outputs for men who've EVER been active with other men
	// Outputs defined to correspond to sampling strategy used by Dunkle et al (2013).

	denom = 0.0;
	everMSM = 0.0;
	everBi = 0.0;
	positive = 0.0;
	positiveMSW = 0.0;
	currentReg = 0.0;

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 0){
			ExactAge = 0.5 + CurrYear - Register[ic].DOB;
			if (ExactAge >= 18.0 && ExactAge<50){
				denom += Register[ic].PopWeight;
				if (Register[ic].EverMSM == 1){ 
					everMSM += Register[ic].PopWeight;
					if (Register[ic].EverBi == 1){ everBi += Register[ic].PopWeight; }
					if (Register[ic].CurrPartners > 0){
						pID1 = Register[ic].IDprimary;
						if (Register[pID1 - 1].SexInd == 0){ currentReg += Register[ic].PopWeight; }
						else{
							if (Register[ic].CurrPartners == 2){
								pID2 = Register[ic].ID2ndary;
								if (Register[pID2 - 1].SexInd == 0){ currentReg += Register[ic].PopWeight; }
							}
						}
					}
					if (Register[ic].HIVstage > 0){ positive += Register[ic].PopWeight; }
				}
				else{
					if (Register[ic].HIVstage > 0){ positiveMSW += Register[ic].PopWeight; }
				}
			}
		}
	}

	EverMSMpropn.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * everMSM / denom;
	EverMSMcurrentReg.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * currentReg / everMSM;
	EverMSMeverBi.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * everBi / everMSM;
	EverMSMpositive.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * positive / everMSM;
	NeverMSMpositive.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * positiveMSW / (denom - everMSM);

	// Calculate age profiles in 2010

	if (FixedUncertainty == 1 && CurrYear == 2010){
		for (ia = 0; ia < 36; ia++){HIVageProfileMSM.out[CurrSim - 1][ia] = 0;}
		for (ia = 0; ia < 9; ia++){ CasualAgeProfile.out[CurrSim - 1][ia] = 0; }
		for (ia = 0; ia < 9; ia++){ BiAgeProfile.out[CurrSim - 1][ia] = 0; }
		HighRiskMSM.out[CurrSim - 1][0] = 0;
		HighRiskMSM.out[CurrSim - 1][1] = 0;
		for (ic = 0; ic < Register.size(); ic++){
			if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 0 && 
				Register[ic].AgeGroup>2 && Register[ic].AgeGroup<12 && Register[ic].VirginInd==0){
				ia = Register[ic].AgeGroup - 3;
				if (Register[ic].RecentMSM >0){
					HIVageProfileMSM.out[CurrSim - 1][ia] += Register[ic].PopWeight;
					if (Register[ic].HIVstage > 0){
						HIVageProfileMSM.out[CurrSim - 1][ia + 9] += Register[ic].PopWeight;
					}
					if (Register[ic].CasualInd == 1){
						CasualAgeProfile.out[CurrSim - 1][ia] += Register[ic].PopWeight;
					}
					if (Register[ic].RecentBi > 0){
						BiAgeProfile.out[CurrSim - 1][ia] += Register[ic].PopWeight;
					}
					if (Register[ic].RiskGroup == 1){
						if (Register[ic].RecentBi > 0){ HighRiskMSM.out[CurrSim - 1][0] += Register[ic].PopWeight; }
						else{ HighRiskMSM.out[CurrSim - 1][1] += Register[ic].PopWeight; }
					}
				}
				else{
					HIVageProfileMSM.out[CurrSim - 1][ia + 18] += Register[ic].PopWeight;
					if (Register[ic].HIVstage > 0){
						HIVageProfileMSM.out[CurrSim - 1][ia + 27] += Register[ic].PopWeight;
					}
				}
			}
		}
	}

	// Calculate HIV prevalence according to sexual preference, in men aged 15-49
	if (FixedUncertainty == 1){
		positive = 0.0;
		posBi = 0.0;
		positiveMSW = 0.0;
		everMSM = 0.0;
		everBi = 0.0;
		denom = 0.0;
		for (ic = 0; ic < Register.size(); ic++){
			if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 0 &&
				Register[ic].AgeGroup>2 && Register[ic].AgeGroup < 10){
				denom += Register[ic].PopWeight;
				if (Register[ic].HIVstage > 0){ positive += Register[ic].PopWeight; }
				if (Register[ic].MalePref > 0.0 && Register[ic].MalePref < 1.0){
					everBi += Register[ic].PopWeight;
					if (Register[ic].HIVstage > 0){ posBi += Register[ic].PopWeight; }
				}
				if (Register[ic].MalePref > 0.0){ everMSM += Register[ic].PopWeight; }
				if (Register[ic].MalePref == 0.0 && Register[ic].HIVstage > 0){ 
					positiveMSW += Register[ic].PopWeight;
				}
			}
		}
		HomoPrefPos.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * (positive - posBi - positiveMSW) / 
			(everMSM - everBi);
		BiPrefPos.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * posBi / everBi;
		HeteroPrefPos.out[CurrSim - 1][CurrYear - StartYear] = 1.0 * positiveMSW / (denom - everMSM);
	}
}

void GetMSMcalib2()
{
	// Same as GemMSMcalib, but allowing for age standardization

	int ic, ia, ii;
	double denom[2], married[2], recentBi[2], everBi[2];
	double positive[2], currentReg[2], multPartners[2], everMSM[2];
	int pID1, pID2;
	double ExactAge;

	// Calculate outputs for men who've recently been active with other men (last 6 months)

	for (ii = 0; ii < 2; ii++){
		denom[ii] = 0.0;
		married[ii] = 0.0;
		recentBi[ii] = 0.0;
		everBi[ii] = 0.0;
		positive[ii] = 0.0;
		currentReg[ii] = 0.0;
		multPartners[ii] = 0.0;
		everMSM[ii] = 0.0;
	}

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 0 &&
			Register[ic].MalePref>0.0 && Register[ic].VirginInd == 0){
			ExactAge = 0.5 + CurrYear - Register[ic].DOB;
			if (ExactAge >= 18.0 && Register[ic].RecentMSM>0){
				if (ExactAge >= 25.0){ii = 1;}
				else{ii = 0;}
				denom[ii] += Register[ic].PopWeight;
				if (Register[ic].MarriedInd == 1){ married[ii] += Register[ic].PopWeight; }
				if (Register[ic].RecentBi > 0){ recentBi[ii] += Register[ic].PopWeight; }
				if (Register[ic].EverBi == 1){ everBi[ii] += Register[ic].PopWeight; }
				if (Register[ic].HIVstage > 0){ positive[ii] += Register[ic].PopWeight; }
				if (Register[ic].CurrPartners > 0){
					pID1 = Register[ic].IDprimary;
					if (Register[pID1 - 1].SexInd == 0){ currentReg[ii] += 1; }
					else{
						if (Register[ic].CurrPartners == 2){
							pID2 = Register[ic].ID2ndary;
							if (Register[pID2 - 1].SexInd == 0){ currentReg[ii] += Register[ic].PopWeight; }
						}
					}
				}
				if (Register[ic].RecentMSM > 1){ multPartners[ii] += Register[ic].PopWeight; }
			}
		}
	}
	for (ii = 0; ii < 2; ii++){
		MSMmarried.out2[CurrSim - 1][CurrYear - StartYear][ii] = 1.0 * married[ii] / denom[ii];
		MSMrecentBi.out2[CurrSim - 1][CurrYear - StartYear][ii] = 1.0 * recentBi[ii] / denom[ii];
		MSMeverBi.out2[CurrSim - 1][CurrYear - StartYear][ii] = 1.0 * everBi[ii] / denom[ii];
		MSMpositive.out2[CurrSim - 1][CurrYear - StartYear][ii] = 1.0 * positive[ii] / denom[ii];
		MSMcurrentReg.out2[CurrSim - 1][CurrYear - StartYear][ii] = 1.0 * currentReg[ii] / denom[ii];
		MSMmultPartners.out2[CurrSim - 1][CurrYear - StartYear][ii] = 1.0 * multPartners[ii] / denom[ii];
	}

	// Calculate outputs for men who've EVER been active with other men
	// Outputs defined to correspond to sampling strategy used by Dunkle et al (2013).

	for (ii = 0; ii < 2; ii++){
		everMSM[ii] = 0.0;
		positive[ii] = 0.0;
		currentReg[ii] = 0.0;
	}

	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].SexInd == 0){
			ExactAge = 0.5 + CurrYear - Register[ic].DOB;
			if (ExactAge >= 18.0 && ExactAge<50){
				if (ExactAge >= 25.0){ii = 1;}
				else{ii = 0;}
				if (Register[ic].EverMSM == 1){
					everMSM[ii] += Register[ic].PopWeight;
					if (Register[ic].CurrPartners > 0){
						pID1 = Register[ic].IDprimary;
						if (Register[pID1 - 1].SexInd == 0){ currentReg[ii] += Register[ic].PopWeight; }
						else{
							if (Register[ic].CurrPartners == 2){
								pID2 = Register[ic].ID2ndary;
								if (Register[pID2 - 1].SexInd == 0){ currentReg[ii] += Register[ic].PopWeight; }
							}
						}
					}
					if (Register[ic].HIVstage > 0){ positive[ii] += Register[ic].PopWeight; }
				}
			}
		}
	}

	for (ii = 0; ii < 2; ii++){
		EverMSMcurrentReg.out2[CurrSim - 1][CurrYear - StartYear][ii] = 1.0 * currentReg[ii] / everMSM[ii];
		EverMSMpositive.out2[CurrSim - 1][CurrYear - StartYear][ii] = 1.0 * positive[ii] / everMSM[ii];
	}
}

void SimulateMSMparameters()
{
	int SimCount2, seed, ii, ind, is, s;
	double x, y, a, b, p, q, bound;
	double r[12], betaVar, CasualAgeAdj10, PrevCasual;

	SimCount2 = (CurrSim - 1) / IterationsPerPC;
	if (FixedUncertainty == 0){
		seed = SimCount2 * 84 + process_num * 8523;
		CRandomMersenne rg(seed);
		for (ii = 0; ii < 12; ii++){
			r[ii] = rg.Random();
		}
	}
	else{
		for (ii = 0; ii < 11; ii++){
			r[ii] = RandomUniformMSM.out[CurrSim - 1][ii];
		}
	}

	// Simulate BiFraction
	a = 12.0;
	b = 3.0;
	p = r[0];
	q = 1 - r[0];
	bound = 0.0;
	ind = 2;
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
	BiFraction = x;

	// Simulate InitMalePrefBeta (2 parameters)
	// Applying the method of moments to calculate the beta parameters from the mean (r[1]) 
	// and the std deviation (0.4*r[1] for r[1] <0.5, 0.4*(1-r[1]) otherwise).
	if (r[1] < 0.5){ betaVar = pow(0.4*r[1], 2.0); }
	else{ betaVar = pow(0.4*(1.0 - r[1]), 2.0); }
	InitMalePrefBeta[0] = r[1]*(r[1] * (1.0 - r[1]) / betaVar - 1.0);
	InitMalePrefBeta[1] = (1.0 - r[1])*(r[1] * (1.0 - r[1]) / betaVar - 1.0);

	// Simulate the AVERAGE annual change in male preference
	a = 0.0;
	b = 0.05;
	p = r[2];
	q = 1 - r[2];
	bound = 0.0;
	ind = 2;
	s = 0;
	cdfnor(&ind, &p, &q, &x, &a, &b, &s, &bound);
	AnnChangeMalePref[0] = x;

	// Simulate the std deviation of the annual change in male preference
	a = 7.50;
	b = 42.5;
	p = r[3];
	q = 1 - r[3];
	bound = 0.0;
	ind = 2;
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
	AnnChangeMalePref[1] = x;

	// Simulate the casual sex adjustment factor for low-risk MSM
	CasualLowAdj = r[4];

	// Simulate the casual sex adjustment factor for high-risk MSM in regular partnerships
	CasualHighAdj = r[5];

	// Simulate the casual sex age adjustment factor
	a = 16.0;
	b = 20.0;
	p = r[6];
	q = 1 - r[6];
	bound = 0.0;
	ind = 2;
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, &bound);
	CasualAgeAdj10 = x;
	CasualAgeAdj = pow(CasualAgeAdj10, 0.1);

	// Simulate the annual rate of exit from casual sex
	a = 4.0;
	b = 4.0;
	p = r[7];
	q = 1 - r[7];
	bound = 0.0;
	ind = 2;
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, &bound);
	CasualExit = 1.0 / x;

	// Simulate the annual rate of entry into casual sex for single high risk gay men aged 20
	a = 12.0;
	b = 3.0;
	p = r[8];
	q = 1 - r[8];
	bound = 0.0;
	ind = 2;
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
	PrevCasual = x;
	CasualEntry = CasualExit * PrevCasual / (1.0 - PrevCasual);

	// Simulate the average duration of ST MSM relationships
	a = 7.716;
	b = 15.4321;
	p = r[9];
	q = 1 - r[9];
	bound = 0.0;
	ind = 2;
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, &bound);
	MeanDurST_MSM = x;

	// Simulate the relative rate of entry into marriage for same-sex relationships
	SameSexMarried = r[10];

	// Store parameters
	if (CurrSim == (SimCount2 * IterationsPerPC + 1)){
		for (ii = 0; ii < 11; ii++){
			MSMparamsLogL.out[SimCount2][ii] = r[ii];
		}
	}
}

void SimulateMSM_HIV()
{
	int SimCount2, seed, ii, ind, is, s;
	double x, y, a, b, p, q, bound;
	double r[12];

	SimCount2 = (CurrSim - 1) / IterationsPerPC;
	if (FixedUncertainty == 0){
		seed = SimCount2 * 103 + process_num * 4077;
		CRandomMersenne rg(seed);
		for (ii = 0; ii < 12; ii++){
			r[ii] = rg.Random();
		}
	}
	else{
		for (ii = 0; ii < 8; ii++){
			r[ii] = RandomUniformMSMT.out[CurrSim - 1][ii];
		}
	}

	// Simulate InitHIVtransmMSM receptive ST/casual
	if (CofactorType == 0){
		a = 3.945;
		b = 354.7;
	}
	else{ // Mean 0.005 and std deviation 0.001
		a = 24.87;
		b = 4949.13;
	}
	p = r[0];
	q = 1 - r[0];
	bound = 0.0;
	ind = 2;
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
	InitHIVtransmMSM[0][0] = x;
	InitHIVtransmMSM[1][0] = x;

	// Simulate InitHIVtransmMSM insertive ST/casual
	/*if (CofactorType == 0){
		a = 3.984;
		b = 1241.0;
	}
	else{
		a = 8.127;
		b = 2023.5;
	}
	p = r[1];
	q = 1 - r[1];
	bound = 0.0;
	ind = 2;
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
	InitHIVtransmMSM[0][1] = x;
	InitHIVtransmMSM[1][1] = x;*/
	InitHIVtransmMSM[0][1] = 0.25 * InitHIVtransmMSM[0][0];
	InitHIVtransmMSM[1][1] = 0.25 * InitHIVtransmMSM[1][0];

	// Simulate InitHIVtransmMSM receptive LT
	/*if (CofactorType == 0){
		a = 3.990;
		b = 1991.0;
	}
	else{
		a = 11.09;
		b = 5532.4;
	}
	p = r[2];
	q = 1 - r[2];
	bound = 0.0;
	ind = 2;
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
	InitHIVtransmMSM[2][0] = x;*/
	InitHIVtransmMSM[2][0] = r[2] * InitHIVtransmMSM[1][0];
	
	// Simulate InitHIVtransmMSM insertive LT
	/*a = 3.997;
	b = 5705.3;
	p = r[3];
	q = 1 - r[3];
	bound = 0.0;
	ind = 2;
	cdfbet(&ind, &p, &q, &x, &y, &a, &b, 0, &bound);
	InitHIVtransmMSM[2][1] = x;*/
	InitHIVtransmMSM[2][1] = 0.25 * InitHIVtransmMSM[2][0];

	// Simulate the freq of sex in ST relationships
	/*a = 13.44;
	b = 2.444;
	p = r[4];
	q = 1 - r[4];
	bound = 0.0;
	ind = 2;
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, &bound);
	FreqSexST_MSM = x;

	// Simulate the freq of casual sex while in the casual sex state
	a = 16.00;
	b = 4.00;
	p = r[5];
	q = 1 - r[5];
	bound = 0.0;
	ind = 2;
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, &bound);
	CasualSexFreq = x;

	// Simulate the OR of condom use in casual sex (relative to ST relationships)
	a = 16.0;
	b = 8.8889;
	p = r[6];
	q = 1 - r[6];
	bound = 0.0;
	ind = 2;
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, &bound);
	ORcondomCasual = x;

	// Simulate the initial HIV prevalence ratio in MSM
	a = 5.325;
	b = 1.775;
	p = r[7];
	q = 1 - r[7];
	bound = 0.0;
	ind = 2;
	cdfgam(&ind, &p, &q, &x, &a, &b, 0, &bound);
	RatioInitPrevMSM = x;*/

	// Store parameters
	if (CurrSim == (SimCount2 * IterationsPerPC + 1)){
		// Note that for ii = 1 and ii = 3 to 7 we don't use the simulated values
		for (ii = 0; ii < 8; ii++){
			MSMTparamsLogL.out[SimCount2][ii] = r[ii];
		}
	}
}

void CalcStructuralAssns()
{
	int ic, ih, ig, im, ib, ii, ip, it, ie, iw, is, ir, iy, year;
	double TempIB[2][2], TempIM[2][2], TempBM_M[2][2], TempBP_M[2][2];
	double TempWT[2][2], TempWM_M[2][2], TempWH_M[2][2], TempEP_M[2][2];
	double TempEH_M[2][2], TempBM_F[2][2], TempBP_F[2][2], TempWM_F[2][2];
	double TempWH_F[2][2], TempEP_F[2][2], TempEH_F[2][2], TempSR[2][2];

	iy = CurrYear - StartYear;

	for (ic = 0; ic < 2; ic++){
		for (ii = 0; ii < 2; ii++){
			TempIB[ic][ii] = 0.0; 
			TempIM[ic][ii] = 0.0; 
			TempBM_M[ic][ii] = 0.0; 
			TempBP_M[ic][ii] = 0.0;
			TempWT[ic][ii] = 0.0; 
			TempWM_M[ic][ii] = 0.0; 
			TempWH_M[ic][ii] = 0.0; 
			TempEP_M[ic][ii] = 0.0;
			TempEH_M[ic][ii] = 0.0; 
			TempBM_F[ic][ii] = 0.0; 
			TempBP_F[ic][ii] = 0.0; 
			TempWM_F[ic][ii] = 0.0;
			TempWH_F[ic][ii] = 0.0; 
			TempEP_F[ic][ii] = 0.0; 
			TempEH_F[ic][ii] = 0.0; 
			TempSR[ic][ii] = 0.0;
		}
	}
	for (ic = 0; ic < Register.size(); ic++){
		if (Register[ic].AliveInd == 1 && Register[ic].Race == 0 && Register[ic].CurrAge >= 15 &&
			Register[ic].CurrAge < 50 && Register[ic].VirginInd == 0){
			ih = 0; // HIV status
			if (Register[ic].HIVstage > 0){ ih = 1; }
			im = 0; // Multiple partners indicator
			if (Register[ic].AnnPartners > 1 || Register[ic].HetCasualInd == 1 || Register[ic].CasualInd == 1){ im = 1; }
			ib = 0; // Binge drinking indicator
			if (Register[ic].DailyDrinkProb > 1.0 / 30.0 && Register[ic].DrinksPerDD >= 5.0){ ib = 1; }
			ii = 0; // Inequitable gender norm indicator
			if (Register[ic].SexInd == 0 && Register[ic].IneqGender >= 0.4){ ii = 1; }
			// The 0.4 threshold is arbitrary - it's 2 times the mean inequitable norm score.
			ip = 0; // Condom use at last sex
			if (Register[ic].CurrPartners > 0 && Register[ic].CondomPrimary > 0.5){ ip = 1; }
			// Not strictly correct - should ideally also take into account other partner types.
			it = 0; // Transactional sex indicator
			if (Register[ic].CasualInd == 1 || Register[ic].HetCasualInd == 1){ it = 1; }
			ie = 0; // Completed 2ndary edu indicator
			if (Register[ic].HighestGrade >= 12){ ie = 1; }
			iw = Register[ic].Employed; // Employment status
			is = Register[ic].InSchool; // Current schooling
			ir = Register[ic].MarriedInd; 
			if (Register[ic].SexInd == 0){
				TempIB[ii][ib] += Register[ic].PopWeight;
				TempIM[ii][im] += Register[ic].PopWeight;
				TempBM_M[ib][im] += Register[ic].PopWeight;
				if (Register[ic].CurrPartners > 0){ TempBP_M[ib][ip] += Register[ic].PopWeight; }
				TempWT[iw][it] += Register[ic].PopWeight;
				TempWM_M[iw][im] += Register[ic].PopWeight;
				TempWH_M[iw][ih] += Register[ic].PopWeight;
				if (Register[ic].CurrPartners > 0){ TempEP_M[ie][ip] += Register[ic].PopWeight; }
				TempEH_M[ie][ih] += Register[ic].PopWeight;
			}
			else{
				TempBM_F[ib][im] += Register[ic].PopWeight;
				if (Register[ic].CurrPartners > 0){ TempBP_F[ib][ip] += Register[ic].PopWeight; }
				TempWM_F[iw][im] += Register[ic].PopWeight;
				TempWH_F[iw][ih] += Register[ic].PopWeight;
				if (Register[ic].CurrPartners > 0){ TempEP_F[ie][ip] += Register[ic].PopWeight; }
				TempEH_F[ie][ih] += Register[ic].PopWeight;
				if (Register[ic].CurrAge < 25){ TempSR[is][ir] += Register[ic].PopWeight; }
			}
		}
	}

	if (TempIB[0][0] == 0 || TempIB[0][1] == 0 || TempIB[1][0] == 0 || TempIB[1][1] == 0){
		IneqGenderBingeAssn.out[CurrSim - 1][iy] = 0.0; // Temporary fix
	}
	else{
		IneqGenderBingeAssn.out[CurrSim - 1][iy] = log(TempIB[0][0] * TempIB[1][1] /
			(TempIB[0][1] * TempIB[1][0]));
	}
	if (TempIM[0][0] == 0 || TempIM[0][1] == 0 || TempIM[1][0] == 0 || TempIM[1][1] == 0){
		IneqGenderMultAssn.out[CurrSim - 1][iy] = 1.0; // Temporary fix
	}
	else{
		IneqGenderMultAssn.out[CurrSim - 1][iy] = log(TempIM[0][0] * TempIM[1][1] /
			(TempIM[0][1] * TempIM[1][0]));
	}
	BingeMultAssnM.out[CurrSim - 1][iy] = log(TempBM_M[0][0] * TempBM_M[1][1] /
		(TempBM_M[0][1] * TempBM_M[1][0]));
	BingeMultAssnF.out[CurrSim - 1][iy] = log(TempBM_F[0][0] * TempBM_F[1][1] /
		(TempBM_F[0][1] * TempBM_F[1][0]));
	BingeCondomAssnM.out[CurrSim - 1][iy] = log(TempBP_M[0][0] * TempBP_M[1][1] /
		(TempBP_M[0][1] * TempBP_M[1][0]));
	BingeCondomAssnF.out[CurrSim - 1][iy] = log(TempBP_F[0][0] * TempBP_F[1][1] /
		(TempBP_F[0][1] * TempBP_F[1][0]));
	EmployedTransAssn.out[CurrSim - 1][iy] = log(TempWT[0][0] * TempWT[1][1] /
		(TempWT[0][1] * TempWT[1][0]));
	EmployedMultAssnM.out[CurrSim - 1][iy] = log(TempWM_M[0][0] * TempWM_M[1][1] /
		(TempWM_M[0][1] * TempWM_M[1][0]));
	EmployedMultAssnF.out[CurrSim - 1][iy] = log(TempWM_F[0][0] * TempWM_F[1][1] /
		(TempWM_F[0][1] * TempWM_F[1][0]));
	EmployedHIVassnM.out[CurrSim - 1][iy] = log(TempWH_M[0][0] * TempWH_M[1][1] /
		(TempWH_M[0][1] * TempWH_M[1][0]));
	EmployedHIVassnF.out[CurrSim - 1][iy] = log(TempWH_F[0][0] * TempWH_F[1][1] /
		(TempWH_F[0][1] * TempWH_F[1][0]));
	EduCondomAssnM.out[CurrSim - 1][iy] = log(TempEP_M[0][0] * TempEP_M[1][1] /
		(TempEP_M[0][1] * TempEP_M[1][0]));
	EduCondomAssnF.out[CurrSim - 1][iy] = log(TempEP_F[0][0] * TempEP_F[1][1] /
		(TempEP_F[0][1] * TempEP_F[1][0]));
	EduHIVassnM.out[CurrSim - 1][iy] = log(TempEH_M[0][0] * TempEH_M[1][1] /
		(TempEH_M[0][1] * TempEH_M[1][0]));
	EduHIVassnF.out[CurrSim - 1][iy] = log(TempEH_F[0][0] * TempEH_F[1][1] /
		(TempEH_F[0][1] * TempEH_F[1][0]));
	if (TempSR[0][0] == 0 || TempSR[0][1] == 0 || TempSR[1][0] == 0 || TempSR[1][1] == 0){
		SchoolMarriageAssn.out[CurrSim - 1][iy] = 0.0;
	}
	else{
		SchoolMarriageAssn.out[CurrSim - 1][iy] = log(TempSR[0][0] * TempSR[1][1] /
			(TempSR[0][1] * TempSR[1][0]));
	}
}

void RunStructRCTs()
{
	int ii, iy, ic, Term;

	for (ii = 1; ii <= 7; ii++){
		if (ii > 2 && ii < 6){
			RestartFromBaseline();
			CurrYear = BaselineStart;
			StructIntScenario = ii;
			Term = BaselineStart - StartYear;
			if (ii == 1){ Term += SingleSessionAlcoholCounselling.MaxTerm; }
			if (ii == 2){ Term += MultiSessionAlcoholCounselling.MaxTerm; }
			if (ii == 3){ Term += CashTransfers.MaxTerm; }
			if (ii == 4){ Term += SchoolSupport.MaxTerm; }
			if (ii == 5){ Term += VocationalTraining.MaxTerm; }
			if (ii == 6){ Term += GenderTransformCommun.MaxTerm; }
			if (ii == 7){ Term += GenderTransformIndiv.MaxTerm; }
			for (iy = BaselineStart - StartYear; iy < Term; iy++){
				RSApop.OneYear();
			}
		}
	}
	//SingleSessionAlcoholCounselling.CalcModelEffects();
	//MultiSessionAlcoholCounselling.CalcModelEffects();
	CashTransfers.CalcModelEffects();
	SchoolSupport.CalcModelEffects();
	VocationalTraining.CalcModelEffects();
	//GenderTransformCommun.CalcModelEffects();
	//GenderTransformIndiv.CalcModelEffects();
	//if (FixedUncertainty == 1){
		/*for (ii = 0; ii < SingleSessionAlcohol.columns; ii++){
			SingleSessionAlcohol.out[CurrSim - 1][ii] = SingleSessionAlcoholCounselling.ModelEffects[ii];}
		for (ii = 0; ii < MultiSessionAlcohol.columns; ii++){
			MultiSessionAlcohol.out[CurrSim - 1][ii] = MultiSessionAlcoholCounselling.ModelEffects[ii];}*/
		for (ii = 0; ii < CashTransferOut.columns; ii++){
			CashTransferOut.out[CurrSim - 1][ii] = CashTransfers.ModelEffects[ii];}
		for (ii = 0; ii < SchoolSupportOut.columns; ii++){
			SchoolSupportOut.out[CurrSim - 1][ii] = SchoolSupport.ModelEffects[ii];}
		for (ii = 0; ii < VocationalTrainingOut.columns; ii++){
			VocationalTrainingOut.out[CurrSim - 1][ii] = VocationalTraining.ModelEffects[ii];}
		/*for (ii = 0; ii < GenderCommun.columns; ii++){
			GenderCommun.out[CurrSim - 1][ii] = GenderTransformCommun.ModelEffects[ii];}
		for (ii = 0; ii < GenderIndiv.columns; ii++){
			GenderIndiv.out[CurrSim - 1][ii] = GenderTransformIndiv.ModelEffects[ii];}*/
	//}

	// Reset to zero for the next simulation.
	StructIntScenario = 0;
}