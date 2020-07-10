
functions{


  real tform(real parin, int transform, data real multiplier, data real meanscale, data real offset, data real inneroffset){
    real param=parin;
    if(meanscale!=1.0) param *= meanscale; 
if(inneroffset != 0.0) param += inneroffset; 
if(transform==0) param = param;
if(transform==1) param = (log1p_exp(param));
if(transform==2) param = (exp(param));
if(transform==3) param = (1/(1+exp(-param)));
if(transform==4) param = ((param)^3);
if(transform==5) param = log1p(param);
if(transform==50) param = meanscale;
if(transform==51) param = 1/(1+exp(-param));
if(transform==52) param = exp(param);
if(transform==53) param = 1/(1+exp(-param))-(exp(param)^2)/(1+exp(param))^2;
if(transform==54) param = 3*param^2;
if(transform==55) param = 1/(1+param);

if(multiplier != 1.0) param *=multiplier;
if(transform < 49 && offset != 0.0) param+=offset;
    return param;
  }
  
}
data {
  int<lower=0> ndatapoints;
  int<lower=1> nmanifest;
  int<lower=1> nlatent;
  int nlatentpop;
  int nsubjects;
  int<lower=0> ntipred; // number of time independent covariates
  int<lower=0> ntdpred; // number of time dependent covariates
  matrix[ntipred ? nsubjects : 0, ntipred ? ntipred : 0] tipredsdata;
  int nmissingtipreds;
  int ntipredeffects;
  real<lower=0> tipredsimputedscale;
  real<lower=0> tipredeffectscale;

  vector[nmanifest] Y[ndatapoints];
  int nopriors;
  int nldynamics;
  vector[ntdpred] tdpreds[ndatapoints];
  
  real maxtimestep;
  real time[ndatapoints];
  int subject[ndatapoints];
  int<lower=0> nparams;
  int continuoustime; // logical indicating whether to incorporate timing information
  int nindvarying; // number of subject level parameters that are varying across subjects
  int nindvaryingoffdiagonals; //number of off diagonal parameters needed for popcov matrix
  vector[nindvarying] sdscale;
  int indvaryingindex[nindvarying];
  int notindvaryingindex[nparams-nindvarying];

  int nt0varstationary;
  int nt0meansstationary;
  int t0varstationary [nt0varstationary, 2];
  int t0meansstationary [nt0meansstationary, 2];

  int nobs_y[ndatapoints];  // number of observed variables per observation
  int whichobs_y[ndatapoints, nmanifest]; // index of which variables are observed per observation
  int ndiffusion; //number of latents involved in system noise calcs
  int derrind[ndiffusion]; //index of which latent variables are involved in system noise calculations
  int drcintoffdiag[nlatent+1];

  int manifesttype[nmanifest];
  int nbinary_y[ndatapoints];  // number of observed binary variables per observation
  int whichbinary_y[ndatapoints, nmanifest]; // index of which variables are observed and binary per observation
  int ncont_y[ndatapoints];  // number of observed continuous variables per observation
  int whichcont_y[ndatapoints, nmanifest]; // index of which variables are observed and continuous per observation
  
  int intoverpop;
  int statedep[10];
  int choleskymats;
  int intoverstates;
  int verbose; //level of printing during model fit
  int subindices[10];
  int TIPREDEFFECTsetup[nparams, ntipred];
  int nrowmatsetup;
  int matsetup[nrowmatsetup,9];
  real matvalues[nrowmatsetup,6];
  int matrixdims[10,2];
  int savescores;
  int savesubjectmatrices;
  int fixedsubpars;
  vector[fixedsubpars ? nindvarying : 0] fixedindparams[fixedsubpars ? nsubjects : 0];
  int dokalman;
  int dokalmanrowsdata[ndatapoints];
  real Jstep;
  real dokalmanpriormodifier;
  int intoverpopindvaryingindex[intoverpop ? nindvarying : 0];
  int nsJAxfinite;
  int sJAxfinite[nsJAxfinite];
  int nJyfinite;
  int sJyfinite[nJyfinite];
  int taylorheun;
  int difftype;
  int jacoffdiag[nlatentpop];
  int njacoffdiagindex;
  int jacoffdiagindex[njacoffdiagindex];
  int popcovn;
  int llsinglerow;
  int doonesubject;
}
      
transformed data{
  matrix[nlatent,nlatent] IIlatent= diag_matrix(rep_vector(1,nlatent));
  matrix[nlatent*nlatent,nlatent*nlatent] IIlatent2 = diag_matrix(rep_vector(1,nlatent*nlatent));
  matrix[nlatentpop,nlatentpop] IIlatentpop = diag_matrix(rep_vector(1,nlatentpop));
  matrix[nindvarying,nindvarying] IIindvar = diag_matrix(rep_vector(1,nindvarying));
  vector[nlatentpop-nlatent] nlpzerovec = rep_vector(0,nlatentpop-nlatent);
  vector[nlatent+1] nlplusonezerovec = rep_vector(0,nlatent+1);
  int nsubjects2 = doonesubject ? 1 : nsubjects;
  int tieffectindices[nparams]=rep_array(0,nparams);
  int ntieffects = 0;
  
  if(ntipred >0){
    for(pi in 1:nparams){
      if(sum(TIPREDEFFECTsetup[pi,]) > .5){
      ntieffects+=1;
      tieffectindices[ntieffects] = pi;
      }
    }
  }
  
}
      
parameters{
  vector[nparams] rawpopmeans; // population level means 

  vector[nindvarying] rawpopsdbase; //population level std dev
  vector[nindvaryingoffdiagonals] sqrtpcov; // unconstrained basis of correlation parameters
  vector[fixedsubpars ? 0 : (intoverpop ? 0 : nindvarying)] baseindparams[fixedsubpars ? 0 : (intoverpop ? 0 :  (doonesubject ? 1 : nsubjects2) )]; //vector of subject level deviations, on the raw scale
  
  vector[ntipredeffects] tipredeffectparams; // effects of time independent covariates
  vector[nmissingtipreds] tipredsimputed;
  //vector[ (( (ntipredeffects-1) * (1-nopriors) ) > 0) ? 1 : 0] tipredglobalscalepar;
  
  vector[intoverstates ? 0 : nlatentpop*ndatapoints] etaupdbasestates; //sampled latent states posterior
  real onesubject[doonesubject ? doonesubject : 0]; //allows multiple specific

  // transformed parameters
  vector[nindvarying] rawpopsd; //population level std dev
  matrix[nindvarying, nindvarying] rawpopcovsqrt; 
  matrix[nindvarying, nindvarying] rawpopcovchol; 
  matrix[nindvarying,nindvarying] rawpopcorr;
  matrix[nindvarying,nindvarying] rawpopcov;

  real ll;
  vector[ndatapoints] llrow;
  matrix[nlatentpop,nlatentpop] etapriorcov[savescores ? ndatapoints : 0];
  matrix[nlatentpop,nlatentpop] etaupdcov[savescores ? ndatapoints : 0];
  matrix[nlatentpop,nlatentpop] etasmoothcov[savescores ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] ypriorcov[savescores ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] yupdcov[savescores ? ndatapoints : 0];
  matrix[nmanifest,nmanifest] ysmoothcov[savescores ? ndatapoints : 0];
  vector[nlatentpop] etaprior[savescores ? ndatapoints : 0];
  vector[nlatentpop] etaupd[savescores ? ndatapoints : 0];
  vector[nlatentpop] etasmooth[savescores ? ndatapoints : 0];
  vector[nmanifest] yprior[savescores ? ndatapoints : 0];
  vector[nmanifest] yupd[savescores ? ndatapoints : 0];
  vector[nmanifest] ysmooth[savescores ? ndatapoints : 0];
  matrix[matrixdims[1, 1], matrixdims[1, 2] ] T0MEANS[subindices[1]  ? (savesubjectmatrices ? nsubjects2 : 1) : 1]; 
  matrix[matrixdims[2, 1], matrixdims[2, 2] ] LAMBDA[subindices[2]  ? (savesubjectmatrices ? nsubjects2 : 1) : 1]; 
  matrix[matrixdims[3, 1], matrixdims[3, 2] ] DRIFT[subindices[3]  ? (savesubjectmatrices ? nsubjects2 : 1) : 1]; 
  matrix[matrixdims[4, 1], matrixdims[4, 2] ] DIFFUSION[subindices[4]  ? (savesubjectmatrices ? nsubjects2 : 1) : 1]; 
  matrix[matrixdims[5, 1], matrixdims[5, 2] ] MANIFESTVAR[subindices[5]  ? (savesubjectmatrices ? nsubjects2 : 1) : 1]; 
  matrix[matrixdims[6, 1], matrixdims[6, 2] ] MANIFESTMEANS[subindices[6]  ? (savesubjectmatrices ? nsubjects2 : 1) : 1]; 
  matrix[matrixdims[7, 1], matrixdims[7, 2] ] CINT[subindices[7]  ? (savesubjectmatrices ? nsubjects2 : 1) : 1]; 
  matrix[matrixdims[8, 1], matrixdims[8, 2] ] T0VAR[subindices[8]  ? (savesubjectmatrices ? nsubjects2 : 1) : 1]; 
  matrix[matrixdims[9, 1], matrixdims[9, 2] ] TDPREDEFFECT[subindices[9]  ? (savesubjectmatrices ? nsubjects2 : 1) : 1]; 
  matrix[matrixdims[10, 1], matrixdims[10, 2] ] PARS[subindices[10]  ? (savesubjectmatrices ? nsubjects2 : 1) : 1];

  matrix[nlatent,nlatent] asymDIFFUSION[(subindices[3] + subindices[4])  ? (savesubjectmatrices ? nsubjects2 : 1) : 1]; //stationary latent process variance
  vector[nlatent] asymCINT[(subindices[3] + subindices[7])  ? (savesubjectmatrices ? nsubjects2 : 1) : 1]; // latent process asymptotic level
  matrix[nlatent, nlatent] DIFFUSIONcov[subindices[4] ? (savesubjectmatrices ? nsubjects2 : 1) : 1];
  matrix[matrixdims[1, 1], matrixdims[1, 2] ] pop_T0MEANS; 
  matrix[matrixdims[2, 1], matrixdims[2, 2] ] pop_LAMBDA; 
  matrix[matrixdims[3, 1], matrixdims[3, 2] ] pop_DRIFT; 
  matrix[matrixdims[4, 1], matrixdims[4, 2] ] pop_DIFFUSION; 
  matrix[matrixdims[5, 1], matrixdims[5, 2] ] pop_MANIFESTVAR; 
  matrix[matrixdims[6, 1], matrixdims[6, 2] ] pop_MANIFESTMEANS; 
  matrix[matrixdims[7, 1], matrixdims[7, 2] ] pop_CINT; 
  matrix[matrixdims[8, 1], matrixdims[8, 2] ] pop_T0VAR; 
  matrix[matrixdims[9, 1], matrixdims[9, 2] ] pop_TDPREDEFFECT; 
  matrix[matrixdims[10, 1], matrixdims[10, 2] ] pop_PARS;

  matrix[nlatent,nlatent] pop_asymDIFFUSION; //stationary latent process variance
  vector[nlatent] pop_asymCINT; // latent process asymptotic level
  matrix[nlatent, nlatent] pop_DIFFUSIONcov;

  matrix[ntipred ? (nmissingtipreds ? nsubjects : 0) : 0, ntipred ? (nmissingtipreds ? ntipred : 0) : 0] tipreds; //tipred values to fill from data and, when needed, imputation vector
  matrix[nparams, ntipred] TIPREDEFFECT; //design matrix of individual time independent predictor effects

}
generated quantities{
  vector[nparams] popmeans;
  vector[nindvarying] popsd; // = rep_vector(0,nparams);
  matrix[nindvarying,nindvarying] popcov;
  matrix[nparams,ntipred] linearTIPREDEFFECT;


  {
    matrix[popcovn, nindvarying] x;
    if(nindvarying){
      for(ri in 1:rows(x)){
        x[ri,] = (rawpopcovchol * 
          to_vector(normal_rng(rep_vector(0,nindvarying),rep_vector(1,nindvarying))) + 
          rawpopmeans[indvaryingindex])';
      }
    }
    
    for(pi in 1:nparams){
      int found=0;
      int pr1;
      int pr2;
      real rawpoppar = rawpopmeans[pi];
      while(!found){
        for(ri in 1:size(matsetup)){
          if(matsetup[ri,9] <=0 && matsetup[ri,3]==pi && matsetup[ri,8]==0) { //if a free parameter 
            pr1 = ri; 
            pr2=ri;// unless intoverpop, pop matrix row reference is simply current row
            found=1;
            if(intoverpop && matsetup[ri,5]) { //check if shifted
              for(ri2 in 1:size(matsetup)){ //check when state reference param of matsetup corresponds to row of t0means in current matsetup row
                if(matsetup[ri2,8]  && matsetup[ri2,3] == matsetup[ri,1] && 
                matsetup[ri2,3] > nlatent && matsetup[ri2,7] < 20 &&
                matsetup[ri,9] <=0) pr2 = ri2; //if param is dynamic and matches row (state ref) and is not in jacobian
                //print("ri = ",ri, " pr2 = ",pr2, " ri2 = ",ri2);
              }
            }
          }
        }
      }
        
      popmeans[pi] = tform(rawpoppar, matsetup[pr2,4], matvalues[pr2,2], matvalues[pr2,3], matvalues[pr2,4], matvalues[pr2,6] ); 
      if(matsetup[pr1,5]){ //if indvarying, transform random sample
        for(ri in 1:rows(x)){
          x[ri,matsetup[pr1,5]] = tform(x[ri,matsetup[pr1,5]],matsetup[pr2,4],matvalues[pr2,2],matvalues[pr2,3],matvalues[pr2,4],matvalues[pr2,6]);
        }
        x[,matsetup[pr1,5]] += rep_vector(-mean(x[,matsetup[pr1,5]]),rows(x));
      }
      if(ntipred > 0){
      for(tij in 1:ntipred){
        if(TIPREDEFFECTsetup[matsetup[pr1,3],tij] ==0){
          linearTIPREDEFFECT[matsetup[pr1,3],tij] = 0;
        } else {
        linearTIPREDEFFECT[matsetup[pr1,3],tij] = ( //tipred reference is from row pr1, tform reference from row pr2 in case of intoverpop
          tform(rawpoppar + TIPREDEFFECT[matsetup[pr1,3],tij] * .01, matsetup[pr2,4], matvalues[pr2,2], matvalues[pr2,3], matvalues[pr2,4], matvalues[pr2,6] ) -
          tform(rawpoppar - TIPREDEFFECT[matsetup[pr1,3],tij] * .01, matsetup[pr2,4], matvalues[pr2,2], matvalues[pr2,3], matvalues[pr2,4], matvalues[pr2,6] )
          ) /2 * 100;
        }
      }
    }
    } //end nparams loop
  
  if(nindvarying){
    popcov = crossprod(x) /(rows(x)-1);
    popsd = sqrt(diagonal(popcov));
  }
  }
}
