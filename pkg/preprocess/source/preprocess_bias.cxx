
//#include "bichip.hxx"
#include "imagestacksnifs.hxx"
#include "catorfile.hxx"
#include "kombinator.hxx"
#include "preprocessor.hxx"

int main(int argc, char **argv) {

  char **argval, **arglabel;
  
  set_arglist("-in none -out none -stackout null -sigcut 3.0 -nlines 5000 -all");
  init_session(argv,argc,&arglabel,&argval);

  CatOrFile catIn(argval[0]);
  CatOrFile catOut(argval[1]);
  int nlines;
  char inName[lg_name+1],outName[lg_name+1];
  Preprocessor P;
  if (is_true(argval[5]))
    P.SetAllImage(1);

  get_argval(4,"%d",&nlines);
  
  if (is_set(argval[2])) {
    /*
    BiChipStackSnifs * outStack = new BiChipStackSnifs(nlines);
  
    while (catIn.NextFile(inName) && catOut.NextFile(outName)) {
      print_msg("Opening %s",inName);
      P.SetIoMethod(kIoSlice);
      BiChipSnifs * tmpOut = P.PreprocessBias(inName,outName);
      outStack->AddBiChip(tmpOut);
    }
    
    double sigma;
    get_argval(3,"%lf", &sigma);
    KGauss k(sigma);
    BiChipSnifs* out = outStack->Kombine(argval[2],&k);

    delete outStack;
    delete out;
  */
  
    
    ImageStackSnifs * outStack = new ImageStackSnifs(nlines);
  
    while (catIn.NextFile(inName) && catOut.NextFile(outName)) {
      print_msg("Opening %s",inName);
      P.SetIoMethod(kIoSlice);
      ImageSnifs * tmpOut = P.PreprocessAssemble(inName,outName);
      outStack->AddImage(tmpOut);
    }
    
    double sigma;
    get_argval(3,"%lf", &sigma);
    KGauss k(sigma);
    ImageSnifs* out = outStack->Kombine(argval[2],&k);

    delete outStack;
    delete out;
  
  } else { // no need to store temporary information
  
    while (catIn.NextFile(inName) && catOut.NextFile(outName)) {
      print_msg("Opening %s",inName);
      ImageSnifs * out = P.PreprocessAssemble(inName,outName);
      delete out;  
    }
  }

  exit_session(0);
}
